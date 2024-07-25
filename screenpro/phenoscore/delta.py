"""
delta module
"""
import numpy as np
from .phenostat import matrixStat


def calculateDelta(x, y, transformation, level):
    """Calculate log ratio of y / x.
    `level` == 'all' (i.e. averaged across all values, sgRNA library elements and replicates)
    `level` == 'col' (i.e. averaged across columns, replicates)

    Args:
        x (np.array): array of values
        y (np.array): array of values
        transformation (str): transformation to use for calculating score
        level (str): level to use for calculating score
    
    Returns:
        np.array: array of log ratio values
    """
    # check if transformation is implemented
    if transformation not in ['log2', 'log10', 'log1p']:
        raise ValueError(f'transformation "{transformation}" not recognized')
    
    if level == 'all':
        # average across all values
        if transformation == 'log2':
            return np.mean(np.log2(y) - np.log2(x))
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x))
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x))
    elif level == 'row':
        # average across rows
        if transformation == 'log2':
            return np.mean(np.log2(y) - np.log2(x), axis=0)
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=0)
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=0)
    elif level == 'col':
        # average across columns
        if transformation == 'log2':
            return np.mean(np.log2(y) - np.log2(x), axis=1)
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=1)
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=1)


def calculatePhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, transformation, level):
    """Calculate phenotype score normalized by negative control and growth rate.
    
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
        transformation (str): transformation to use for calculating score
        level (str): level to use for calculating score
    
    Returns:
        np.array: array of scores
    """
    # calculate control median and std
    ctrl_median = np.median(calculateDelta(x=x_ctrl, y=y_ctrl, transformation=transformation, level=level))

    # calculate delta
    delta = calculateDelta(x=x, y=y, transformation=transformation, level=level)

    # calculate score
    return (delta - ctrl_median) / growth_rate


def matrixTest(x, y, x_ctrl, y_ctrl, transformation, level, test = 'ttest', growth_rate = 1):
    """Calculate phenotype score and p-values comparing `y` vs `x` matrices.

    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        transformation (str): transformation to use for calculating score
        level (str): level to use for calculating score and p-value
        test (str): test to use for calculating p-value
        growth_rate (int): growth rate
    
    Returns:
        np.array: array of scores
        np.array: array of p-values
    """
    # calculate growth score
    scores = calculatePhenotypeScore(
        x = x, y = y, x_ctrl = x_ctrl, y_ctrl = y_ctrl,
        growth_rate = growth_rate, transformation = transformation,
        level = level
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, level = level)

    return scores, p_values
