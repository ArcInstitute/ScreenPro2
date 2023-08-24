"""phenoScore module
"""

import numpy as np
import pandas as pd

from pydeseq2 import preprocessing
from scipy.stats import ttest_rel
# from scipy.stats import mannwhitneyu
# from statsmodels.stats.multitest import multipletests


def seqDepthNormalization(adata):
    """Normalize counts by sequencing depth
    """
    norm_counts, size_factors = preprocessing.deseq2_norm(adata.X)

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
    """log ratio of y / x
    `ave` == 'all' – i.e. averaged across all values, oligo and replicates
    `ave` == 'col' – i.e. averaged across columns, replicates)
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


def getScore(x, y, x_ctrl, y_ctrl, growth_rate, math, ave):
    """Calculate phenotype score normalized by negative control and growth rate
    """
    ctrl_std = np.std(getDelta(x=x_ctrl, y=y_ctrl, math=math, ave=ave))
    ctrl_median = np.median(getDelta(x=x_ctrl, y=y_ctrl, math=math, ave=ave))
    delta = getDelta(x=x, y=y, math=math, ave=ave)

    return ((delta - ctrl_median) / growth_rate) / ctrl_std


def generatePseudoGeneLabels(adata, num_pseudogenes=None, ctrl_label='negCtrl'):
    """Generate new labels per `num_pseudogenes` randomly selected non targeting oligo in `adata.var`
    """
    if num_pseudogenes is None:
        num_pseudogenes = len(adata.var[adata.var.targetType.eq(ctrl_label)]) // 2
    # get non-targeting oligos
    ctrl_oligos = adata.var[adata.var.targetType.eq(ctrl_label)].index
    adata.var['pseudoLabel'] = ''
    # check if there are more than 1 non-targeting oligos to label as pseudogenes
    if len(ctrl_oligos) / 2 <= num_pseudogenes:
        raise TypeError("Define `num_pseudogenes` to be less than total number of non-targeting oligos / 2")
    else:
        while len(ctrl_oligos) > num_pseudogenes:
            # randomly select `num` non-targeting oligos
            pseudo_oligos = np.random.choice(ctrl_oligos, num_pseudogenes, replace=False)
            # generate new labels
            pseudo_labels = [f'pseudo_{i}' for i in range(num_pseudogenes)]
            # update adata.var
            adata.var.loc[pseudo_oligos, 'pseudoLabel'] = pseudo_labels
            # ...
            ctrl_oligos = ctrl_oligos.drop(pseudo_oligos)

        adata.var.loc[adata.var.targetType.eq('gene'), 'pseudoLabel'] = 'gene'
        adata.var.loc[adata.var.pseudoLabel.eq(''), 'pseudoLabel'] = np.nan


def matrixStat(x, y, test, ave_reps):
    """Get p-values comparing `y` vs `x` matrices
    """
    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test
        pass
    elif test == 'ttest':
        # run ttest
        if ave_reps:
            p_value = ttest_rel(y, x, axis=1)[1]

        else:
            p_value = ttest_rel(y, x)[1]

        return p_value
    else:
        raise ValueError(f'Test "{test}" not recognized')


def matrixTest(x, y, x_ctrl, y_ctrl, math, ave_reps, test = 'ttest', growth_rate = 1):
    """Calculate phenotype score and p-values comparing `y` vs `x` matrices
    """
    if ave_reps:
        ave = 'col'
    else:
        ave = 'all'

    # calculate growth score
    scores = getScore(
        x = x, y = y, x_ctrl = x_ctrl, y_ctrl = y_ctrl,
        growth_rate = growth_rate, math = math,
        ave = ave
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, ave_reps=ave_reps)

    return scores, p_values


def ann_score_df(df_in, up_hit='resistance_hit', down_hit='sensitivity_hit', ctrl_label='non-targeting', threshold=10):
    """Annotate score dataframe with hit labels using given `threshold`
       (i.e. `score/pseudo_sd * -np.log10(pvalue) >= threshold`)
    """
    df = df_in.copy()

    df.columns = ['target', 'score', 'pvalue']
    df['score'] = df['score'].astype(float)
    df['pvalue'] = df['pvalue'].astype(float)

    pseudo_sd = df[df['target'].str.contains(ctrl_label)]['score'].tolist()
    pseudo_sd = np.std(pseudo_sd)

    df['label'] = '.'

    df.loc[
        (df['score'] > 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) >= threshold), 'label'
    ] = up_hit

    df.loc[
        (df['score'] < 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) <= -threshold), 'label'
    ] = down_hit

    df.loc[df['target'].str.contains(ctrl_label), 'label'] = ctrl_label

    df.loc[df['label'] == '.', 'label'] = 'target_non_hit'

    # reorder factors
    df['label'] = pd.Categorical(
        df['label'],
        categories=[down_hit, up_hit, ctrl_label, 'target_non_hit']
    )

    return df
