"""
Module for loading screen datasets
"""

import pickle
import pandas as pd

from .utils import check_protospacer_length, trim_protospacer


def load_cas9_sgRNA_library(library_path, library_type, sep='\t', index_col=0, protospacer_length=19, verbose=True, **args):
    '''Load Cas9 sgRNA library table for single or dual guide design.
    '''
    library = pd.read_csv(
        library_path,
        sep=sep,
        index_col=index_col,
        **args
    )

    ## Evaluate library table and reformat columns for downstream analysis
    # I would like to name the target column 'target' if it is named 'gene'!
    
    if library_type == "single_guide_design":
        eval_columns = ['target', 'sgID', 'protospacer', 'sequence']

        # reformating columns as needed
        if 'gene' in library.columns:
            # rename gene column to target
            library = library.rename(columns={'gene': 'target'})
        if 'sequence' in library.columns and 'protospacer' not in library.columns:
            library.rename(columns={'sequence': 'protospacer'}, inplace=True)
        if 'sgId' in library.columns:
            library.rename(columns={'sgId': 'sgID'}, inplace=True)

        # Upper case protospacer sequences
        library['protospacer'] = library['protospacer'].str.upper()

        protospacer_col = 'protospacer'
        in_length = check_protospacer_length(library, 'protospacer')
        if in_length == protospacer_length:
            pass
        elif in_length > protospacer_length:
            if verbose: print(f"Trimming protospacer sequences in '{protospacer_col}' column.")
            library = trim_protospacer(
                library, protospacer_col, 
                '5prime', 
                in_length - protospacer_length
            )

        elif in_length < protospacer_length:
            raise ValueError(
                f"Input protospacer length for '{protospacer_col}' is less than {protospacer_length}"
            )
        
        # write `sequence` column as `protospacer` (after trimming)
        library['sequence'] = library['protospacer']

        for col in eval_columns:
            if col not in library.columns:
                raise ValueError(f"Column '{col}' not found in library table.")
        
        library = library[eval_columns]

    elif library_type == "dual_guide_design":
        eval_columns = [
            'target', 'sgID_AB',
            'sgID_A', 'protospacer_A', 
            'sgID_B', 'protospacer_B', 
            'sequence'
        ]
        
        # reformating columns as needed
        if 'gene' in library.columns:
            # rename gene column to target
            library = library.rename(columns={'gene': 'target'})

        # Upper case protospacer sequences
        library['protospacer_A'] = library['protospacer_A'].str.upper()
        library['protospacer_B'] = library['protospacer_B'].str.upper()

        # # TODO: Enable trimming of protospacer sequences through command line arguments.
        for protospacer_col in ['protospacer_A', 'protospacer_B']:
            in_length = check_protospacer_length(library, protospacer_col) 
            if in_length == protospacer_length:
                pass
            elif in_length > protospacer_length:
                if verbose: print(f"Trimming protospacer sequences in '{protospacer_col}' column.")
                library = trim_protospacer(
                    library, protospacer_col, 
                    '5prime', 
                    in_length - protospacer_length
                )

            elif in_length < protospacer_length:
                raise ValueError(
                    f"Input protospacer length for '{protospacer_col}' is less than {protospacer_length}"
                )
    
        # write `sequence` column as `protospacer_A;protospacer_B` (after trimming)
        library['sequence'] = library['protospacer_A'] + ';' + library['protospacer_B']

        for col in eval_columns:
            if col not in library.columns:
                raise ValueError(f"Column '{col}' not found in library table.")

        library = library[eval_columns]

    else:
        raise ValueError(f"Invalid library type: {library_type}. Please choose 'single_guide_design' or 'dual_guide_design'.")
    
    if verbose: print("Library table successfully loaded.")

    return library


def loadScreenProcessingData(experimentName, collapsedToTranscripts=True, premergedCounts=False):
    """
    Load ScreenProcessing outputs
    (see original code `here <https://github.com/mhorlbeck/ScreenProcessing/blob/master/screen_analysis.py#L70>`__)
    Input files:
        * `*_librarytable.txt` => library table
        * `*_mergedcountstable.txt` => merged counts table
        * `*_phenotypetable.txt` => phenotype table

    Parameters:        
        experimentName (str): name of the experiment
        collapsedToTranscripts (bool): whether the gene scores are collapsed to transcripts
        premergedCounts (bool): whether the counts are premerged
    
    Returns:
        dict: dictionary of dataframes
    """
    # dict of dataframes
    dataDict = {
        'library': pd.read_csv(
            experimentName + '_librarytable.txt', 
            sep='\t', 
            header=0, 
            index_col=0
        ),
        'counts': pd.read_csv(
            experimentName + '_mergedcountstable.txt', 
            sep='\t', 
            header=list(range(2)),
            index_col=list(range(1))
        ),
        'phenotypes': pd.read_csv(
            experimentName + '_phenotypetable.txt', 
            sep='\t', 
            header=list(range(2)),
            index_col=list(range(1))
        )
    }

    if premergedCounts:
        # add premerged counts
        dataDict['premerged counts'] = pd.read_csv(
            experimentName + '_rawcountstable.txt', 
            sep='\t',
            header=list(range(3)), 
            index_col=list(range(1))
        )

    if collapsedToTranscripts:
        # add transcript scores
        dataDict['transcript scores'] = pd.read_csv(
            experimentName + '_genetable.txt', 
            sep='\t', 
            header=list(range(3)),
            index_col=list(range(2))
        )
        dataDict['gene scores'] = pd.read_csv(
            experimentName + '_genetable_collapsed.txt', 
            sep='\t',
            header=list(range(3)), 
            index_col=list(range(1))
        )
    else:
        # add gene scores
        dataDict['gene scores'] = pd.read_csv(
            experimentName + '_genetable.txt', 
            sep='\t', 
            header=list(range(3)),
            index_col=list(range(1))
        )

    return dataDict


def write_screen_pkl(screen, name):
    """
    Write AnnData object to a pickle file
    
    Parameters:
        screen (object): ScreenPro object to save
        name (str): name of the output file (.pkl extension will be added)
    """
    file_name = f'{name}.pkl'
    with open(file_name, 'wb') as file:
        pickle.dump(screen, file)
        print(f'Object successfully saved to "{file_name}"')


def read_screen_pkl(name):
    """
    Read ScreenPro object from a pickle file
    
    Parameters:
        name (str): name of the input file (.pkl extension will be added)
    """
    file_name = f'{name}.pkl'
    with open(file_name, 'rb') as f:
        screen = pickle.load(f)
    return screen
