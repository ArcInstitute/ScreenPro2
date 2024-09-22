## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

## This part of the software is conceptualized and developed by Abolfazl Arab (@abearab)
## with support from the Nick Youngblut (@nick-youngblut).

'''Scripts to work with NGS data

This module provides functions to process FASTQ files from screens with single or dual guide
libraries. In general, the algorithm is fairly simple:

1. Read the FASTQ file and extract the proper sequences
2. Count the exact number of occurrences for each unique sequence
3. Map the counted sequences to the reference sequence library
4. Return the counted mapped or unmapped events as a dataframe(s)

For single-guide screens, the sequences are counted as single protospacer
from a single-end read file (R1). Then, these sequences are mapped to the reference
library of protospacer sequences.

For dual-guide screens, the sequences are counted as pairs of protospacer A and B
from paired-end read files (R1 and R2). Then, sequences are mapped to the reference
library of protospacer A and B pairs.

Theoretically, the algorithm is able to detect any observed sequence since it is counting first
and then mapping. Therefore, the recombination events can be detected. In dual-guide design
protospacer A and B are not the same pairs as in the reference library. These events include:

- Protospacer A and B pairs are present in the reference library but paired differently
- Only one of the protospacer A and B is present in the reference library
- None of the protospacer A and B is present in the reference library
'''

import pandas as pd
import polars as pl
import anndata as ad
import os

from . import cas9
from . import cas12
from ..load import load_cas9_sgRNA_library
from simple_colors import green


class GuideCounter:
    '''Class to count sequences from FASTQ files
    '''

    def __init__(self, cas_type, library_type):
        self.cas_type = cas_type
        self.library_type = library_type
        self.counts_dict = None
        self.counts_mat = None
        self.recombinants = None

    def load_library(self, library_path, sep='\t', index_col=0, protospacer_length=19, verbose=False, **args):
        '''Load library file
        '''
        if self.cas_type == 'cas9':

            library = load_cas9_sgRNA_library(library_path, library_type=self.library_type, sep=sep, index_col=index_col, protospacer_length=protospacer_length, verbose=verbose, **args)

            # Check if the library has duplicate sequences and remove them
            if library.duplicated('sequence').any():
                shape_before_dedup = library.shape[0]
                library = library.drop_duplicates(subset='sequence', keep='first')
                shape_after_dedup = library.shape[0]
                if verbose:
                    print(f"Warning: {shape_before_dedup - shape_after_dedup} duplicate sgRNA sequences found and removed.")

            if self.library_type == "single_guide_design":
                sgRNA_table = library[['target','sgID','protospacer']].set_index('sgID')

            elif self.library_type == "dual_guide_design":
                sgRNA_table = pd.concat([
                    library[['target','sgID_A','protospacer_A']].rename(columns={'sgID_A':'sgID','protospacer_A':'protospacer'}),
                    library[['target','sgID_B', 'protospacer_B']].rename(columns={'sgID_B':'sgID','protospacer_B':'protospacer'})
                ])
                # drop duplicates and set index
                sgRNA_table = sgRNA_table.drop_duplicates(keep='first')

            if verbose: print('total # of cas9 sgRNAs:', sgRNA_table.shape[0])

        elif self.cas_type == 'cas12':
            raise NotImplementedError("Cas12 library is not yet implemented.")
        
        # covert to polar DataFrame
        library = pl.from_pandas(library)
        sgRNA_table = pl.from_pandas(sgRNA_table)

        self.library = library
        self.sgRNA_table = sgRNA_table

    def _process_cas9_single_guide_sample(self, fastq_dir, sample_id, trim_first_g, protospacer_length, write, verbose=False):
        if verbose: print(green(sample_id, ['bold']))
        get_counts = True

        # check if df_count is already available
        if os.path.exists(f'{fastq_dir}/{sample_id}_count.arrow'):
            if verbose: print('count file exists ...')
            if write != "force":
                df_count = pl.read_ipc_stream(f'{fastq_dir}/{sample_id}_count.arrow')
                get_counts = False
            else:
                if verbose: print('skip loading count file, force write is set ...')

        if get_counts:
            if trim_first_g:
                trim5p_start = 2
            else:
                trim5p_start = 1
            df_count = cas9.fastq_to_count_single_guide(
                fastq_file_path=f'{fastq_dir}/{sample_id}.fastq.gz',
                trim5p_start=trim5p_start,
                trim5p_length=protospacer_length,
                verbose=verbose
            )
            if write == "force" or write == True:
                # write df_count to file
                df_count.write_ipc_stream(f'{fastq_dir}/{sample_id}_count.arrow', compression='lz4')
                if verbose: print('count file written ...')

        out = cas9.map_to_library_single_guide(
            df_count=df_count,
            library=self.library,
            return_type='all',
            verbose=verbose
        )
        
        return out
    
    def _process_cas9_dual_guide_sample(self, fastq_dir, sample_id, get_recombinant, trim_first_g, protospacer_A_length, protospacer_B_length, write, verbose=False):
        if verbose: print(green(sample_id, ['bold']))
        get_counts = True

        # check if df_count is already available
        if os.path.exists(f'{fastq_dir}/{sample_id}_count.arrow'):
            if verbose: print('count file exists ...')
            if write != "force":
                df_count = pl.read_ipc_stream(f'{fastq_dir}/{sample_id}_count.arrow')
                get_counts = False
            else:
                if verbose: print('skip loading count file, force write is set ...')
        
        if get_counts:
            if get_counts:
                if trim_first_g == True or trim_first_g == {'A':True, 'B':True}:
                    trim5p_pos1_start = 2
                    trim5p_pos2_start = 2
                elif trim_first_g == False or trim_first_g == {'A':False, 'B':False}:
                    trim5p_pos1_start = 1
                    trim5p_pos2_start = 1
                elif trim_first_g == {'A':True, 'B':False}:
                    trim5p_pos1_start = 2
                    trim5p_pos2_start = 1
                elif trim_first_g == {'A':False, 'B':True}:
                    trim5p_pos1_start = 1
                    trim5p_pos2_start = 2
                else:
                    raise ValueError("Invalid trim_first_g argument. Please provide a boolean or a dictionary with 'A' and 'B' keys.")

            df_count = cas9.fastq_to_count_dual_guide(
                R1_fastq_file_path=f'{fastq_dir}/{sample_id}_R1.fastq.gz',
                R2_fastq_file_path=f'{fastq_dir}/{sample_id}_R2.fastq.gz',
                trim5p_pos1_start=trim5p_pos1_start,
                trim5p_pos1_length=protospacer_A_length,
                trim5p_pos2_start=trim5p_pos2_start,
                trim5p_pos2_length=protospacer_B_length,
                verbose=verbose
            )
            if write == "force" or write == True:
                # write df_count to file
                df_count.write_ipc_stream(f'{fastq_dir}/{sample_id}_count.arrow', compression='lz4')
                if verbose: print('count file written ...')

        out = cas9.map_to_library_dual_guide(
            df_count=df_count,
            library=self.library,
            get_recombinant=get_recombinant,
            return_type='all',
            verbose=verbose
        )
        
        return out        

    def get_counts_matrix(self, fastq_dir, samples, get_recombinant=False, cas_type='cas9', protospacer_length='auto', trim_first_g=False, write=True, verbose=False):
        '''Get count matrix for given samples
        '''
        if self.cas_type == 'cas9':
            counts = {}

            if self.library_type == "single_guide_design":
                if get_recombinant:
                    raise ValueError("Recombinants are not applicable for single guide design!")
                if protospacer_length == 'auto':
                    protospacer_length = self.library['protospacer'].str.len_bytes().unique().to_list()[0]

                for sample_id in samples:
                    cnt = self._process_cas9_single_guide_sample(
                        fastq_dir=fastq_dir, 
                        sample_id=sample_id, 
                        trim_first_g=trim_first_g,
                        protospacer_length=protospacer_length,
                        write=write,
                        verbose=verbose
                    )
                    
                    counts[sample_id] = cnt['mapped']
        
                counts_mat = pd.concat([
                    counts[sample_id].to_pandas().set_index('sgID')['count'].rename(sample_id)
                    for sample_id in counts.keys()
                ],axis=1).fillna(0)

            elif self.library_type == "dual_guide_design":
                if get_recombinant: recombinants = {}

                if protospacer_length == 'auto':
                    protospacer_A_length = self.library['protospacer_A'].str.len_bytes().unique().to_list()[0]
                    protospacer_B_length = self.library['protospacer_B'].str.len_bytes().unique().to_list()[0]
                elif isinstance(protospacer_length, dict):
                    protospacer_A_length = protospacer_length['protospacer_A']
                    protospacer_B_length = protospacer_length['protospacer_B']
                elif isinstance(protospacer_length, int):
                    protospacer_A_length = protospacer_length
                    protospacer_B_length = protospacer_length
                else:
                    raise ValueError("Invalid protospacer_length argument. If not 'auto', please provide an integer or a dictionary with 'protospacer_A' and 'protospacer_B' keys.")

                for sample_id in samples:
                    cnt = self._process_cas9_dual_guide_sample(
                        fastq_dir=fastq_dir, 
                        sample_id=sample_id, 
                        get_recombinant=get_recombinant, 
                        trim_first_g=trim_first_g,
                        protospacer_A_length=protospacer_A_length,
                        protospacer_B_length=protospacer_B_length,
                        write=write, 
                        verbose=verbose
                    )
                    counts[sample_id] = cnt['mapped']
                    if get_recombinant:
                        recombinants[sample_id] = cnt['recombinant']
                
                counts_mat = pd.concat([
                    counts[sample_id].to_pandas().set_index('sgID_AB')['count'].rename(sample_id) 
                    for sample_id in counts.keys()
                ],axis=1).fillna(0)
            
            else:
                raise ValueError("Invalid library type. Please choose from 'single_guide_design' or 'dual_guide_design'.")

        if cas_type == 'cas12':
            # TODO: Implement codes to build count matrix for given samples
            raise NotImplementedError("Cas12 count matrix is not yet implemented.")
        
        self.counts_dict = counts
        self.counts_mat = counts_mat
        if get_recombinant:
            self.recombinants = recombinants
    
    def load_counts_matrix(self, counts_mat_path, **kwargs):
        '''Load count matrix from file
        '''
        self.counts_mat = pd.read_csv(counts_mat_path, **kwargs)
    
    def _build_cas9_dual_guide_var_table(self, counts_table, source, ctrl_label='negative_control'):
        '''Build variant table for dual guide design

        Args:
            counts_table (pd.DataFrame): count table for dual guide design (e.g. main library mapped counts or recombinant counts)
        '''
        if source=='library':
            var_table = pd.DataFrame(
                counts_table.index.str.split('|').to_list(),
                index = counts_table.index.to_list(),
                columns=['sgID_A','sgID_B']
            )
        
        elif source=='recombinant':
            var_table = pd.DataFrame(
                counts_table.index.to_list(),
                index = ['|'.join(i) for i in counts_table.index.to_list()],
                columns=['sgID_A','sgID_B']
            )
        var_table.index.name = 'sgID_AB'

        sgRNA_table = self.sgRNA_table.to_pandas().set_index('sgID')
        
        #TODO: extract "target" values from protospacer IDs
        var_table = pd.concat([
            var_table.reset_index().reset_index(drop=True),
            sgRNA_table.loc[var_table['sgID_A']].rename(columns={'target':'target_A', 'protospacer':'protospacer_A'}).reset_index(drop=True),
            sgRNA_table.loc[var_table['sgID_B']].rename(columns={'target':'target_B', 'protospacer':'protospacer_B'}).reset_index(drop=True),
        ], axis=1).set_index('sgID_AB')

        var_table['targetType'] = ''
        var_table['target'] = ''

        ### assign target types: negative_control
        control_targets = (var_table.target_A.eq(ctrl_label)) & (var_table.target_B.eq(ctrl_label))
        var_table.loc[control_targets,'targetType']  = 'negative_control'
        var_table.loc[control_targets,'target']  = ctrl_label

        ### assign target types: gene
        same_gene_targets = (var_table.target_A == var_table.target_B) & ~(var_table.target_A.eq(ctrl_label)) & ~(var_table.target_B.eq(ctrl_label))
        var_table.loc[same_gene_targets,'targetType']  = 'gene'
        var_table.loc[same_gene_targets,'target']  = var_table.target_A # or target_B

        ### assign target types: gene-negative_control
        gene_control_targets = ~(var_table.target_A.eq(ctrl_label)) & (var_table.target_B.eq(ctrl_label))
        var_table.loc[gene_control_targets,'targetType']  = 'gene--negative_control'
        var_table.loc[gene_control_targets,'target']  = var_table.target_A + '|' + var_table.target_B

        ### assign target types: negative_control-gene
        control_gene_targets = (var_table.target_A.eq(ctrl_label)) & ~(var_table.target_B.eq(ctrl_label))
        var_table.loc[control_gene_targets,'targetType']  = 'negative_control--negative_control'
        var_table.loc[control_gene_targets,'target']  = var_table.target_A + '|' + var_table.target_B

        ### assign target types: gene-gene
        gene_gene_targets = (var_table.target_A != var_table.target_B) & ~(var_table.target_A.eq(ctrl_label)) & ~(var_table.target_B.eq(ctrl_label))
        var_table.loc[gene_gene_targets,'targetType']  = 'gene--gene'
        var_table.loc[gene_gene_targets,'target']  = var_table.target_A + '|' + var_table.target_B

        var_table.index.name = None
        var_table.targetType = pd.Categorical(
            var_table.targetType, categories=[
                'gene','gene--gene',
                'gene--negative_control','negative_control--gene',
                'negative_control'
            ]
        ).remove_unused_categories()

        var_table['sequence'] = var_table['protospacer_A'] + ';' + var_table['protospacer_B']

        return var_table

    def build_counts_anndata(self, source='library', verbose=False):
        '''Build AnnData object from count matrix
        '''
        if source == 'recombinant' and self.library_type == "single_guide_design":
            raise ValueError("Recombinants are not applicable for single guide design!")
        if source == 'recombinant' and self.recombinants is None:
            raise ValueError("Recombinants are not available. If applicable, please set get_recombinant=True in get_counts_matrix method.")

        if self.library_type == "single_guide_design":
            adata = ad.AnnData(
                X = self.counts_mat.T, 
                var = self.library.to_pandas().set_index('sgID')
            )
        
        elif self.library_type == "dual_guide_design":
            
            adata = ad.AnnData(
                X = self.counts_mat.T, 
                var = self._build_cas9_dual_guide_var_table(self.counts_mat, source='library')
            )
            
            if source == 'recombinant':
                counts_recombinants = {}

                for sample in self.recombinants.keys():
                    if verbose: print(green(sample, ['bold']))
                    d = self.recombinants[sample].drop_nulls()
                    d = d.to_pandas()
                    counts_recombinants[sample] = d.set_index(['sgID_A','sgID_B'])['count']
                    if verbose: print('recombinant count added ...')

                counts_recombinants = pd.concat(counts_recombinants,axis=1).fillna(0)

                if verbose: print('recombinant count matrix built ...')

                var_table = self._build_cas9_dual_guide_var_table(counts_recombinants, source='recombinant')

                rdata = ad.AnnData(
                    X = counts_recombinants.T.to_numpy(),
                    var = var_table,
                    obs = adata.obs
                )

                if verbose: print('recombinant AnnData created.')

        if source == 'mapped' or source == 'library':
            return adata
        elif source == 'recombinant':
            return rdata
        else:
            raise ValueError("Invalid source argument. Please choose from 'mapped', 'recombinant' or 'library'. Note: 'mapped' and 'library' act the same way.")
