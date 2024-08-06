"""
phenostat module: internal module for statistical analysis of phenoscore data.
"""

from scipy.stats import ttest_rel
import numpy as np
from statsmodels.stats.multitest import multipletests


def matrixStat(x, y, test, level, transform='log10'):
    """
    Get p-values comparing `y` vs `x` matrices.

    Parameters:
        x (np.array): array of values
        y (np.array): array of values
        test (str): test to use for calculating p-value
        level (str): level at which to calculate p-value
        transform (str): transformation to apply to values before running test
    
    Returns:
        np.array: array of p-values
    """
    # log-transform values
    if transform == None:
        pass
    elif transform == 'log10':
        x = np.log10(x)
        y = np.log10(y)
    else:
        raise ValueError(f'Transform "{transform}" not recognized')
    
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
            p_value = ttest_rel(y, x, axis=None)[1]
        else:
            raise ValueError(f'Level "{level}" not recognized')
        return p_value
    else:
        raise ValueError(f'Test "{test}" not recognized')


def multipleTestsCorrection(p_values, method='fdr_bh'):
    """
    Calculate adjusted p-values using multiple testing correction.

    Parameters:
        p_values (np.array): array of p-values
        method (str): method to use for multiple testing correction
    
    Returns:
        np.array: array of adjusted p-values
    """
    if method == 'fdr_bh':
        # fill na with 1
        p_values[np.isnan(p_values)] = 1
        # Calculate the adjusted p-values using the Benjamini-Hochberg method
        if p_values is None:
            raise ValueError('p_values is None')
        _, adj_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    else:
        raise ValueError(f'Method "{method}" not recognized')

    return adj_p_values


def empiricalFDR():
    pass
