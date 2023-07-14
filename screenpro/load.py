"""load screen datasets
"""
import pandas as pd


def loadScreenProcessingData(experimentName, collapsedToTranscripts=True, premergedCounts=False):
    """load ScreenProcessing outputs (original code)
    https://github.com/mhorlbeck/ScreenProcessing/blob/ddf27b1d2e968984f2073bc9f77969524ac18d23/screen_analysis.py#L70
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

