"""
phenostat module
"""

from scipy.stats import ttest_rel
import numpy as np
from statsmodels.stats.multitest import multipletests


def matrixStat(x, y, test, level):
    """
    Get p-values comparing `y` vs `x` matrices.

    Parameters:
        x (np.array): array of values
        y (np.array): array of values
        test (str): test to use for calculating p-value
        level (str): level at which to calculate p-value
    
    Returns:
        np.array: array of p-values
    """
    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test
        raise ValueError('Mann-Whitney U rank test not implemented')
    
    elif test == 'KS':
        # run Kolmorogov-Smirnov test
        raise ValueError('Kolmorogov-Smirnov test not implemented')
    
    elif test == 'ttest':
        # run ttest
        if level == 'col':
            p_value = ttest_rel(y, x, axis=1)[1]
        elif level == 'row':
            p_value = ttest_rel(y, x, axis=0)[1]
        elif level == 'all':
            # average across all values
            p_value = ttest_rel(y, x)[1]
        else:
            raise ValueError(f'Level "{level}" not recognized')
        return p_value
    else:
        raise ValueError(f'Test "{test}" not recognized')


def getFDR(p_values, method='fdr_bh'):
    """
    Calculate FDR.

    Parameters:
        p_values (np.array): array of p-values
        method (str): method to use for calculating FDR
    
    Returns:
        np.array: array of adjusted p-values
    """
    # fill na with 1
    p_values[np.isnan(p_values)] = 1
    # Calculate the adjusted p-values using the Benjamini-Hochberg method
    if p_values is None:
        raise ValueError('p_values is None')
    _, adj_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

    return adj_p_values
