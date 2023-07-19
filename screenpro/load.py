"""load screen datasets
"""
import pickle
import pandas as pd


def loadScreenProcessingData(experimentName, collapsedToTranscripts=True, premergedCounts=False):
    """load ScreenProcessing outputs
    (see original code `here <https://github.com/mhorlbeck/ScreenProcessing/blob/master/screen_analysis.py#L70>`__)
    """
    dataDict = {'library': pd.read_csv(experimentName + '_librarytable.txt', sep='\t', header=0, index_col=0),
                'counts': pd.read_csv(experimentName + '_mergedcountstable.txt', sep='\t', header=list(range(2)),
                                      index_col=list(range(1))),
                'phenotypes': pd.read_csv(experimentName + '_phenotypetable.txt', sep='\t', header=list(range(2)),
                                          index_col=list(range(1)))}

    if premergedCounts:
        dataDict['premerged counts'] = pd.read_csv(experimentName + '_rawcountstable.txt', sep='\t',
                                                   header=list(range(3)), index_col=list(range(1)))

    if collapsedToTranscripts:
        dataDict['transcript scores'] = pd.read_csv(experimentName + '_genetable.txt', sep='\t', header=list(range(3)),
                                                    index_col=list(range(2)))
        dataDict['gene scores'] = pd.read_csv(experimentName + '_genetable_collapsed.txt', sep='\t',
                                              header=list(range(3)), index_col=list(range(1)))
    else:
        dataDict['gene scores'] = pd.read_csv(experimentName + '_genetable.txt', sep='\t', header=list(range(3)),
                                              index_col=list(range(1)))

    return dataDict


def write_adata_pkl(adata, name):
    file_name = f'{name}.pkl'
    with open(file_name, 'wb') as file:
        pickle.dump(adata, file)
        print(f'Object successfully saved to "{file_name}"')


def read_adata_pkl(name):
    file_name = f'{name}.pkl'
    with open(file_name, 'rb') as f:
        adata = pickle.load(f)
    return adata
