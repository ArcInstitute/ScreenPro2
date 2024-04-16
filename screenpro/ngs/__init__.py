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
from simple_colors import green

from . import cas9
from . import cas12

class Counter:
    '''Class to count sequences from FASTQ files
    '''

    def __init__(self, cas_type, library_type):
        self.cas_type = cas_type
        self.library_type = library_type        

    def load_library(self, library_path):
        '''Load library file
        '''
        if self.cas_type == 'cas9':
            library = pd.read_csv(
                library_path,
                sep='\t',
                index_col=0,
            )

            # I would like to name the target column 'target' if it is named 'gene'!
            if 'gene' in library.columns:
                # rename gene column to target
                library = library.rename(columns={'gene': 'target'})

            # Evaluate library table
            if self.library_type == "single_guide_design":
                pass

            elif self.library_type == "dual_guide_design":
                # reformat columns for downstream analysis
                if 'sgID_AB' in library.columns:
                    library = library.set_index('sgID_AB')
                library = library.rename(
                    columns={'protospacer_A':'protospacer_a','protospacer_B':'protospacer_b'}
                )

                if 'sequence' not in library.columns:
                    # TODO: Enable trimming of protospacer sequences through command line arguments.
                    library.protospacer_a = library.protospacer_a.str.upper()
                    library.protospacer_b = library.protospacer_b.str[1:].str.upper()
                    library['sequence'] = library.protospacer_a + ';' + library.protospacer_b

                # check required columns are present:
                eval_columns = ['target', 'sequence', 'protospacer_A', 'protospacer_B', 'sequence']
                for col in eval_columns:
                    if col not in library.columns:
                        raise ValueError(f"Column '{col}' not found in library table.")

                library = library[eval_columns]

        elif self.cas_type == 'cas12':
            raise NotImplementedError("Cas12 library is not yet implemented.")
        
        self.library = library
        
    def get_matrix(self, fastq_dir, samples,library,get_recombinant=False, cas_type='cas9',verbose=False):
        '''Get count matrix for given samples
        '''
        if self.cas_type == 'cas9':
            if self.library_type == "single_guide_design":
                # TODO: Implement codes to build count matrix for given samples
                pass
            elif self.library_type == "dual_guide_design":
                counts = {}
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
                            library=library,
                            get_recombinant=False,
                            return_type='mapped',
                            verbose=verbose
                        )
                    elif get_recombinant:
                        raise NotImplementedError("Recombinant count matrix are not yet implemented.")
            
            counts_mat = pd.concat([
                counts[sample_id].to_pandas()['count'].rename(sample_id)
                for sample_id in counts.keys()
            ],axis=1).fillna(0)
        if cas_type == 'cas12':
            # TODO: Implement codes to build count matrix for given samples
            raise NotImplementedError("Cas12 count matrix is not yet implemented.")
        
        return counts_mat
