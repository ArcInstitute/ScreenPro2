"""
phenoScore module
"""

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests


def runPhenoScore(adata, cond1, cond2, math, test, score_level,
                  growth_rate=1, n_reps=2, ctrl_label='negCtrl'):
    """
    Calculate phenotype score and p-values when comparing `cond2` vs `cond1`.
    Args:
        adata (AnnData): AnnData object
        cond1 (str): condition 1
        cond2 (str): condition 2
        math (str): math to use for calculating score
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        score_level (str): score level
        growth_rate (int): growth rate
        n_reps (int): number of replicates
        ctrl_label (str): control label
    Returns:
        str: result name
        pd.DataFrame: result dataframe
    """
    # format result name
    result_name = f'{cond2}_vs_{cond1}'
    print(f'\t{cond2} vs {cond1}')

    count_layer = 'seq_depth_norm'
    # check if count_layer exists
    if 'seq_depth_norm' not in adata.layers.keys():
        seqDepthNormalization(adata)

    # calc phenotype score and p-value
    if score_level == 'compare_reps':
        # prep counts for phenoScore calculation
        df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index[:n_reps],].to_df(count_layer).T
        df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index[:n_reps],].to_df(count_layer).T

        # convert to numpy arrays
        x = df_cond1.to_numpy()
        y = df_cond2.to_numpy()

        # get control values
        x_ctrl = df_cond1[adata.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond2[adata.var.targetType.eq(ctrl_label)].to_numpy()

        # calculate growth score and p_value
        scores, p_values = matrixTest(
            x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
            math=math, ave_reps=True, test=test, 
            growth_rate=growth_rate
        )

        # Calculate the adjusted p-values using the Benjamini-Hochberg method
        if p_values is None:
            raise ValueError('p_values is None')
        _, adj_pvalues, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        # get targets
        targets = adata.var.index.str.split('_[-,+]_').str[0].to_list()
        # combine results into a dataframe
        result = pd.concat([
            pd.Series(targets, index=adata.var.index, name='target'),
            pd.Series(scores, index=adata.var.index, name='score'),
            pd.Series(p_values, index=adata.var.index, name=f'{test} pvalue'),
            pd.Series(adj_pvalues, index=adata.var.index, name='BH adj_pvalue'),
        ], axis=1)

        return result_name, result
    elif score_level in ['compare_oligos', 'compare_oligos_and_reps']:
        return None
    else:
        raise ValueError(f'score_level "{score_level}" not recognized')


def seqDepthNormalization(adata):
    """
    Normalize counts by sequencing depth.
    Args:
        adata (AnnData): AnnData object
    """
    # normalize counts by sequencing depth
    norm_counts, size_factors = preprocessing.deseq2_norm(adata.X)
    # update adata object
    adata.obs['size_factors'] = size_factors
    adata.layers['seq_depth_norm'] = norm_counts


def addPseudoCount():
    pass
    # # pseudocount
    # if pseudocountBehavior == 'default' or pseudocountBehavior == 'zeros only':
    #     def defaultBehavior(row):
    #         return row if min(
    #             row) != 0 else row + pseudocountValue
    #
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


def getDelta(x, y, math, ave):
    """
    Calculate log ratio of y / x.
    `ave` == 'all' (i.e. averaged across all values, oligo and replicates)
    `ave` == 'col' (i.e. averaged across columns, replicates)
    Args:
        x (np.array): array of values
        y (np.array): array of values
        math (str): math to use for calculating score
        ave (str): average method
    Returns:
        np.array: array of log ratio values
    """
    if ave == 'all':
        # average across all values
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1))
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x))
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x))
    elif ave == 'row':
        # average across rows
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=0)
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=0)
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=0)
    elif ave == 'col':
        # average across columns
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=1)
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=1)
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=1)


def getPhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, math, ave):
    """
    Calculate phenotype score normalized by negative control and growth rate.
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
        math (str): math to use for calculating score
        ave (str): average method
    Returns:
        np.array: array of scores
    """
    # calculate control median and std
    ctrl_median = np.median(getDelta(x=x_ctrl, y=y_ctrl, math=math, ave=ave))

    # calculate delta
    delta = getDelta(x=x, y=y, math=math, ave=ave)

    # calculate score
    return (delta - ctrl_median) / growth_rate


def getZScorePhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, math, ave):
    pass
    # # calculate control median and std
    # ctrl_std = np.std(getDelta(x=x_ctrl, y=y_ctrl, math=math, ave=ave))
    # ctrl_median = np.median(getDelta(x=x_ctrl, y=y_ctrl, math=math, ave=ave))
    #
    # # calculate delta
    # delta = getDelta(x=x, y=y, math=math, ave=ave)
    #
    # # calculate score
    # return ((delta - ctrl_median) / growth_rate) / ctrl_std


def generatePseudoGeneLabels(adata, num_pseudogenes=None, ctrl_label='negCtrl'):
    """
    Generate new labels per `num_pseudogenes` randomly selected non targeting oligo in `adata.var`.
    Args:
        adata (AnnData): AnnData object
        num_pseudogenes (int): number of pseudogenes
        ctrl_label (str): control label
    """
    # check if num_pseudogenes is defined. 
    ## If not, set to half of the number of non-targeting oligos
    if num_pseudogenes is None:
        num_pseudogenes = len(adata.var[adata.var.targetType.eq(ctrl_label)]) // 2
    # Get non-targeting oligos
    ctrl_oligos = adata.var[adata.var.targetType.eq(ctrl_label)].index
    adata.var['pseudoLabel'] = ''
    # Check if there are more than 1 non-targeting oligos to label as pseudogenes
    if len(ctrl_oligos) / 2 <= num_pseudogenes:
        # raise error if `num_pseudogenes` is greater than (total number of non-targeting oligos) / 2
        raise TypeError(
            "Define `num_pseudogenes` to be less than (total number of non-targeting oligos) / 2"
        )
    else:
        # randomly select `num` non-targeting oligos
        while len(ctrl_oligos) > num_pseudogenes:
            # randomly select `num` non-targeting oligos
            pseudo_oligos = np.random.choice(ctrl_oligos, num_pseudogenes, replace=False)
            # generate new labels
            pseudo_labels = [f'pseudo_{i}' for i in range(num_pseudogenes)]
            # update adata.var
            adata.var.loc[pseudo_oligos, 'pseudoLabel'] = pseudo_labels
            # remove selected oligos from ctrl_oligos
            ctrl_oligos = ctrl_oligos.drop(pseudo_oligos)

        # label remaining non-targeting oligos as pseudogenes
        adata.var.loc[adata.var.targetType.eq('gene'), 'pseudoLabel'] = 'gene'
        adata.var.loc[adata.var.pseudoLabel.eq(''), 'pseudoLabel'] = np.nan


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


def matrixTest(x, y, x_ctrl, y_ctrl, math, ave_reps, test = 'ttest', growth_rate = 1):
    """
    Calculate phenotype score and p-values comparing `y` vs `x` matrices.
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        math (str): math to use for calculating score
        ave_reps (bool): average replicates
        test (str): test to use for calculating p-value
        growth_rate (int): growth rate
    Returns:
        np.array: array of scores
        np.array: array of p-values
    """
    # check if average across replicates
    ave = 'col' if ave_reps else 'all'

    # calculate growth score
    scores = getPhenotypeScore(
        x = x, y = y, x_ctrl = x_ctrl, y_ctrl = y_ctrl,
        growth_rate = growth_rate, math = math,
        ave = ave
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, ave_reps=ave_reps)

    return scores, p_values
