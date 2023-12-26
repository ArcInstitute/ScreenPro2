"""
phenoStats module
"""

from scipy.stats import ttest_rel
import numpy as np
from statsmodels.stats.multitest import multipletests


def matrixStat(x, y, test, ave_reps):
    """
    Get p-values comparing `y` vs `x` matrices.
    Args:
        x (np.array): array of values
        y (np.array): array of values
        test (str): test to use for calculating p-value
        ave_reps (bool): average replicates
    Returns:
        np.array: array of p-values
    """
    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test
        pass
    elif test == 'ttest':
        # run ttest
        if ave_reps:
            # average across replicates
            p_value = ttest_rel(y, x, axis=1)[1]
        else:
            # average across all values
            p_value = ttest_rel(y, x)[1]
        return p_value
    else:
        raise ValueError(f'Test "{test}" not recognized')


def getFDR(p_values, method='fdr_bh'):
    """
    Calculate FDR.
    """
    # fill na with 1
    p_values[np.isnan(p_values)] = 1
    # Calculate the adjusted p-values using the Benjamini-Hochberg method
    if p_values is None:
        raise ValueError('p_values is None')
    _, adj_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

    return adj_p_values
