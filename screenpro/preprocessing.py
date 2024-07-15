import pandas as pd
import numpy as np
from pydeseq2 import preprocessing


def findLowCounts(adata, filter_type, minimum_reads, verbose=True):
    """
    Label variables with low counts in either or all samples.

    Parameters:
        adata (AnnData): AnnData object containing the counts to be filtered.
        filter_type (str): specify the filter type. Possible values are: 'all', or 'sum'.
        minimum_reads (int): minimum number of reads.
        verbose (bool): print the number of removed variables. Default is True.

    Returns:
        None
    """
    count_bin = adata.X >= minimum_reads
 
    if filter_type == 'all':
        out = adata[:, count_bin.all(axis=0)].copy()
    elif filter_type == 'sum':
        out = adata[:, adata.to_df().sum(axis=0) >= minimum_reads].copy()
    else:
        raise ValueError(f'filter_type "{filter_type}" not recognized. Use "all", or "sum".')
    
    if verbose:
        n_removed = adata.shape[1] - out.shape[1]
        # print the number of removed variables
        print(
            f"{n_removed} variables with less than {minimum_reads} reads (filter_type: '{filter_type}')"
        )

    adata.var['low_count'] = ~adata.var.index.isin(out.var.index.to_list())


def addPseudoCount(adata, behavior, value, inplace=True):
    """
    Add pseudocounts to the given counts based on the specified behavior.

    Args:
        adata (AnnData): AnnData object containing the counts to which pseudocounts will be added.
        behavior (str): The behavior for adding pseudocounts. Possible values are:
            - 'default' or 'zeros_only': Add pseudocounts only to rows with at least one zero value.
            - 'all_values': Add pseudocounts to all rows.
            - 'filter_out': Set rows with at least one zero value to NaN.
        value (float): The value of the pseudocount to be added.
        inplace (bool): If True, the pseudocounts will replace the original counts in the AnnData object.

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
    
    counts = adata.to_df()
    
    # pseudocount
    if behavior == 'default' or behavior == 'zeros_only':
        counts_pseudo = counts.replace(0, value)

    elif behavior == 'all_values':
        counts_pseudo = counts + value
    elif behavior == 'filter_out':
        counts_pseudo = counts.replace(0, np.nan)
    else:
        raise ValueError(
            'Pseudocount behavior not recognized or not implemented')

    if inplace:
        adata.X = counts_pseudo.to_numpy()
    else:
        return counts_pseudo


def normalizeSeqDepth(adata):
    """
    Normalize counts by sequencing depth and update the adata object.
    This function uses the PyDESeq2 normalization method.

    Args:
        adata (AnnData): AnnData object containing the counts to be normalized.
    """
    # normalize counts by sequencing depth
    norm_counts, size_factors = preprocessing.deseq2_norm(adata.X)
    # update adata object
    adata.obs['size_factors'] = size_factors
    adata.layers['seq_depth_norm'] = norm_counts
    adata.X = adata.layers['seq_depth_norm']
