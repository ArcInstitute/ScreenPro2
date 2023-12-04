"""
Module for loading screen datasets
"""

# import
import pickle
import pandas as pd


# functions
def loadScreenProcessingData(experimentName, collapsedToTranscripts=True, premergedCounts=False):
    """
    Load ScreenProcessing outputs
    (see original code `here <https://github.com/mhorlbeck/ScreenProcessing/blob/master/screen_analysis.py#L70>`__)
    Input files:
        * `*_librarytable.txt` => library table
        * `*_mergedcountstable.txt` => merged counts table
        * `*_phenotypetable.txt` => phenotype table
    Args:
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
    Args:
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
    Args:
        name (str): name of the input file (.pkl extension will be added)
    """
    file_name = f'{name}.pkl'
    with open(file_name, 'rb') as f:
        screen = pickle.load(f)
    return screen
