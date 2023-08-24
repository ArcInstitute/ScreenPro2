from .__version__ import __version__
from .phenoScore import seqDepthNormalization, getDelta, getScore
from scipy.stats import ttest_rel
# from scipy.stats import mannwhitneyu
# from statsmodels.stats.multitest import multipletests


def runPhenoScoreByReps(adata, cond1, cond2, growth_rate=1, n_reps=2, ctrl_label='negCtrl', test='ttest'):
    """Calculate phenotype score and p-values comparing `cond2` vs `cond1`
    """
    print(f'\t{cond2} vs {cond1}')

    count_layer = 'seq_depth_norm'
    # check if count_layer exists
    if 'seq_depth_norm' not in adata.layers.keys():
        seqDepthNormalization(adata)
    # prep counts for phenoScore calculation
    df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index[:n_reps], ].to_df(count_layer).T
    df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index[:n_reps], ].to_df(count_layer).T

    x = df_cond1.to_numpy()
    y = df_cond2.to_numpy()

    x_ctrl = df_cond1[adata.var.targetType.eq(ctrl_label)].to_numpy()
    y_ctrl = df_cond2[adata.var.targetType.eq(ctrl_label)].to_numpy()

    # calculate growth score
    phenotype_score = getScore(x, y, x_ctrl, y_ctrl, growth_rate)

    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test on replicates
        pass
    if test == 'ttest':
        # run ttest on replicates
        pvalues = ttest_rel(y, x, axis=1)[1]

    return phenotype_score, pvalues
