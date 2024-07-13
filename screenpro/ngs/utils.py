import pandas as pd
import numpy as np


def find_low_counts(adata, filter_type='either', minimum_reads=50):
    """
    Label variables with low counts in either or all samples.

    Parameters:
        adata (AnnData): AnnData object
        filter_type (str): either or all
        minimum_reads (int): minimum number of reads

    Returns:
        None
    """
    count_bin = adata.X >= minimum_reads
    if filter_type == 'either':
        out = adata[:, ~(~count_bin.all(axis=0))].copy()
    elif filter_type == 'all':
        out = adata[:, count_bin.all(axis=0)].copy()
    elif filter_type == 'sum':
        out = adata[:, adata.to_df().sum(axis=0) >= minimum_reads].copy()
    else:
        raise ValueError(f'filter_type "{filter_type}" not recognized. Use "either", "all", or "sum".')
    
    # print the number of removed variables
    n_removed = adata.shape[1] - out.shape[1]
    print(
        f"{n_removed} variables with less than {minimum_reads} reads (filter_type: '{filter_type}')"
    )

    adata.var['low_count'] = ~adata.var.index.isin(out.var.index.to_list())


def addPseudoCount(counts, pseudocountBehavior, pseudocountValue):
    pass

    ## possible pseudocount behaviors
    # 1. remove 0
    # 2. add pseudocount
    # 3. impute 0 (it's hard)

    # # pseudocount
    # if pseudocountBehavior == 'default' or pseudocountBehavior == 'zeros only':
    #     def defaultBehavior(row):
    #         return row if min(
    #             row) != 0 else row + pseudocountValue
    
    #     combinedCountsPseudo = combinedCounts.apply(defaultBehavior, axis=1)
    # elif pseudocountBehavior == 'all values':
    #     combinedCountsPseudo = combinedCounts.apply(
    #         lambda row: row + pseudocountValue, axis=1)
    # elif pseudocountBehavior == 'filter out':
    #     combinedCountsPseudo = combinedCounts.copy()
    #     zeroRows = combinedCounts.apply(lambda row: min(row) <= 0, axis=1)
    #     combinedCountsPseudo.loc[zeroRows, :] = np.nan
    # else:
    #     raise ValueError(
    #         'Pseudocount behavior not recognized or not implemented')