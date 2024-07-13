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


def add_pseudo_count(counts, behavior, value):
    """
    Add pseudocounts to the given counts based on the specified behavior.

    Args:
        counts (DataFrame): The counts to which pseudocounts will be added.
        behavior (str): The behavior for adding pseudocounts. Possible values are:
            - 'default' or 'zeros_only': Add pseudocounts only to rows with at least one zero value.
            - 'all_values': Add pseudocounts to all rows.
            - 'filter_out': Set rows with at least one zero value to NaN.
        value (float): The value of the pseudocount to be added.

    Returns:
        DataFrame: The counts DataFrame with pseudocounts added based on the specified behavior.

    Raises:
        ValueError: If the pseudocount behavior is not recognized or not implemented.
    """
    ## possible pseudocount behaviors
    # 1. remove 0
    # 2. add pseudocount
    # 3. impute 0 (it's hard)
    
    # Source:
    # https://github.com/mhorlbeck/ScreenProcessing/blob/0ee5192ecc17348665bd1387ddfa9037efb7964f/process_experiments.py#L485
    
    # pseudocount
    if behavior == 'default' or behavior == 'zeros_only':
        counts_pseudo = counts.apply(
            lambda row: row if min(row) != 0 else row + value,
            axis=1
        )
    elif behavior == 'all_values':
        counts_pseudo = counts.apply(
            lambda row: row + value,
            axis=1
        )
    elif behavior == 'filter_out':
        counts_pseudo = counts.copy()
        zero_rows = counts.apply(lambda row: min(row) <= 0, axis=1)
        counts_pseudo.loc[zero_rows, :] = np.nan
    else:
        raise ValueError(
            'Pseudocount behavior not recognized or not implemented')

    return counts_pseudo
