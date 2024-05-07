import pandas as pd
import polars as pl
import anndata as ad
# import multiprocess as mp

from . import cas9
from . import cas12
from ..load import load_cas9_sgRNA_library
from ..utils import find_low_counts
from simple_colors import green


class Counter:
    '''Class to count sequences from FASTQ files
    '''

    def __init__(self, cas_type, library_type):
        self.cas_type = cas_type
        self.library_type = library_type
        self.counts_dict = None
        self.counts_mat = None
        self.recombinants = None

    def load_library(self, library_path, sep='\t', index_col=0, verbose=False):
        '''Load library file
        '''
        if self.cas_type == 'cas9':
            library = load_cas9_sgRNA_library(library_path, library_type=self.library_type, sep=sep, index_col=index_col, verbose=verbose)

        elif self.cas_type == 'cas12':
            raise NotImplementedError("Cas12 library is not yet implemented.")
        
        # covert to polar DataFrame
        library = pl.from_pandas(library)

        self.library = library
        
    def _get_sgRNA_table(self):

        if self.cas_type == 'cas9':
            if self.library_type == "single_guide_design":
                sgRNA_table = self.library.to_pandas()[['target','sgID','protospacer']].set_index('sgID')

            elif self.library_type == "dual_guide_design":
                sgRNA_table = pd.concat([
                    self.library.to_pandas()[['target','sgID_A','protospacer_A']].rename(columns={'sgID_A':'sgID','protospacer_A':'protospacer'}),
                    self.library.to_pandas()[['target','sgID_B', 'protospacer_B']].rename(columns={'sgID_B':'sgID','protospacer_B':'protospacer'})
                ]).set_index('sgID')

        return sgRNA_table

    def get_counts_matrix(self, fastq_dir, samples,get_recombinant=False, cas_type='cas9', parallel=False, verbose=False):
        '''Get count matrix for given samples
        '''
        if self.cas_type == 'cas9':
            counts = {}

            if self.library_type == "single_guide_design":
                # TODO: Implement codes to build count matrix for given samples
                pass
            elif self.library_type == "dual_guide_design":
                if get_recombinant: recombinants = {}
                
                # TODO: Fix the function for parallel processing
                def process_sample(sample_id):
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
                            get_recombinant=True,
                            return_type='all',
                            verbose=verbose
                        )
                        counts[sample_id] = cnt['mapped']
                        recombinants[sample_id] = cnt['recombinant']
                
                if parallel:
                    raise NotImplementedError("Parallel processing is not yet implemented.")
                    # pool = mp.Pool(len(samples))
                    # pool.map(process_sample, samples)
                    # pool.close()
                    # pool.join()
                
                else:
                    for sample_id in samples:
                        process_sample(sample_id)
                
            counts_mat = pd.concat([
                counts[sample_id].to_pandas().set_index('sgID_AB')['count'].rename(sample_id) 
                for sample_id in counts.keys()
            ],axis=1).fillna(0)

        if cas_type == 'cas12':
            # TODO: Implement codes to build count matrix for given samples
            raise NotImplementedError("Cas12 count matrix is not yet implemented.")
        
        self.counts_dict = counts
        self.counts_mat = counts_mat
        if get_recombinant:
            self.recombinants = recombinants
    
    def load_counts_matrix(self, counts_mat_path, **kwargs):
        '''Load count matrix
        '''
        self.counts_mat = pd.read_csv(counts_mat_path, **kwargs)
    
    def build_counts_anndata(self, source='library'):
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
                var = self.library.to_pandas().set_index('sgID_AB')
            )
            
            if source == 'recombinant':
                counts_recombinants = {}

                for sample in self.recombinants.keys():
                    d = self.recombinants[sample].drop_nulls()
                    d = d.to_pandas()
                    counts_recombinants[sample] = d.set_index(['sgID_A','sgID_B'])['count']

                counts_recombinants = pd.concat(counts_recombinants,axis=1).fillna(0)

                counts_recombinants = pd.concat([
                    counts_recombinants,
                    pd.concat([
                        # add non-targeting counts from the main count matrix
                        self.counts_mat[self.counts_mat.index.str.contains('non-targeting')],
                        # add non-targeting counts from the main count matrix
                        self.library.to_pandas().set_index('sgID_AB')[['sgID_A','sgID_B']][self.counts_mat.index.str.contains('non-targeting')]
                    ],axis=1).set_index(['sgID_A','sgID_B'])
                ])

                var_table = pd.DataFrame(
                    counts_recombinants.index.to_list(),
                    index = ['|'.join(i) for i in counts_recombinants.index.to_list()],
                    columns=['sgID_A','sgID_B'])

                var_table.index.name = 'sgID_AB'

                sgRNA_table = self._get_sgRNA_table()
                
                var_table = pd.concat([
                    var_table.reset_index().reset_index(drop=True),
                    sgRNA_table.loc[var_table['sgID_A']].rename(columns={'target':'target_A', 'protospacer':'protospacer_A'}).reset_index(drop=True),
                    sgRNA_table.loc[var_table['sgID_B']].rename(columns={'target':'target_B', 'protospacer':'protospacer_B'}).reset_index(drop=True),
                ], axis=1).set_index('sgID_AB')

                var_table['targetType'] = ''

                var_table.loc[
                    (var_table.target_A.eq('negative_control')) & 
                    (var_table.target_B.eq('negative_control')),'targetType']  = 'negCtrl'

                var_table.loc[
                    ~(var_table.target_A.eq('negative_control')) & 
                    (var_table.target_B.eq('negative_control')),'targetType']  = 'gene-ctrl'

                var_table.loc[
                    (var_table.target_A.eq('negative_control')) & 
                    ~(var_table.target_B.eq('negative_control')),'targetType']  = 'ctrl-gene'

                var_table.loc[
                    ~(var_table.target_A.eq('negative_control')) & 
                    ~(var_table.target_B.eq('negative_control')),'targetType']  = 'gene-gene'

                var_table.index.name = None

                var_table['target'] = var_table['target_A'] + '|' + var_table['target_B']
                var_table['sequence'] = var_table['protospacer_A'] + ';' + var_table['protospacer_B']

                rdata = ad.AnnData(
                    X = counts_recombinants.T.to_numpy(),
                    var = var_table,
                    obs = adata.obs
                )

                rdata = rdata[:,~rdata.var.low_count]

        if source == 'mapped' or source == 'library':
            return adata
        elif source == 'recombinant':
            return rdata
        else:
            raise ValueError("Invalid source argument. Please choose from 'mapped', 'recombinant' or 'library'. Note: 'mapped' and 'library' act the same way.")