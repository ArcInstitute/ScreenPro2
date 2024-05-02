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

from . import cas9
from . import cas12
from ..load import load_cas9_sgRNA_library
from simple_colors import green


class Counter:
    '''Class to count sequences from FASTQ files
    '''

    def __init__(self, cas_type, library_type):
        self.cas_type = cas_type
        self.library_type = library_type

    def load_library(self, library_path, sep='\t', verbose=False):
        '''Load library file
        '''
        if self.cas_type == 'cas9':
            library = load_cas9_sgRNA_library(library_path, library_type=self.library_type, sep=sep, verbose=verbose)

        elif self.cas_type == 'cas12':
            raise NotImplementedError("Cas12 library is not yet implemented.")
        
        # covert to polar DataFrame
        library = pl.from_pandas(library)

        self.library = library
        
    def get_counts_matrix(self, fastq_dir, samples,get_recombinant=False, cas_type='cas9',verbose=False):
        '''Get count matrix for given samples
        '''
        if self.cas_type == 'cas9':
            counts = {}

            if self.library_type == "single_guide_design":
                # TODO: Implement codes to build count matrix for given samples
                pass
            elif self.library_type == "dual_guide_design":
                if get_recombinant: recombinants = {}
                for sample_id in samples:
                    if verbose: print(green(sample_id, ['bold']))
                    df_count = cas9.fastq_to_count_dual_guide(
                        R1_fastq_file_path=f'{fastq_dir}/{sample_id}_R1.fastq.gz',
                        R2_fastq_file_path=f'{fastq_dir}/{sample_id}_R2.fastq.gz',
                        trim5p_pos1_start=1,
                        trim5p_pos1_length=19,
                        trim5p_pos2_start=1,
                        trim5p_pos2_length=19,
                        verbose=verbose
                    )
                    if not get_recombinant:
                        counts[sample_id] = cas9.map_to_library_dual_guide(
                            df_count=df_count,
                            library=self.library,
                            get_recombinant=False,
                            return_type='mapped',
                            verbose=verbose
                        )
                    elif get_recombinant:
                        cnt = cas9.map_to_library_dual_guide(
                            df_count=df_count,
                            library=self.library,
                            get_recombinant=False,
                            return_type='all',
                            verbose=verbose
                        )
                        counts[sample_id] = cnt['mapped']
                        recombinants[sample_id] = cnt['recombinant']
            
            counts_mat = pd.concat([
                counts[sample_id].to_pandas()['count'].rename(sample_id)
                for sample_id in counts.keys()
            ],axis=1).fillna(0)
        if cas_type == 'cas12':
            # TODO: Implement codes to build count matrix for given samples
            raise NotImplementedError("Cas12 count matrix is not yet implemented.")
        
        self.counts_dict = counts
        self.counts_mat = counts_mat
        if get_recombinant:
            self.recombinants = recombinants
