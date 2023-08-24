import pandas as pd
from .__version__ import __version__
from .phenoScore import seqDepthNormalization, getDelta, getScore
from scipy.stats import ttest_rel
# from scipy.stats import mannwhitneyu
# from statsmodels.stats.multitest import multipletests


def runPhenoScoreByReps(adata,
                        cond1, cond2, growth_rate=1, n_reps=2,
                        ctrl_label='negCtrl', test='ttest', math='log2(x+1)'):
    """Calculate phenotype score and p-values comparing `cond2` vs `cond1`
    """
    result_name = f'{cond2}_vs_{cond1}'
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
    phenotype_score = getScore(
        x = x, y = y, x_ctrl = x_ctrl, y_ctrl = y_ctrl,
        growth_rate = growth_rate, math = math
    )

    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test on replicates
        pass
    elif test == 'ttest':
        # run ttest on replicates
        pvalues = ttest_rel(y, x, axis=1)[1]

        return phenotype_score, pvalues, result_name
    else:
        raise ValueError(f'Test "{test}" not recognized')


def runPhenoScoreByGuideSet(adata,
                            cond1, cond2, growth_rate=1, n_reps=2,
                            ctrl_label='negCtrl', test='ttest', math='log2(x+1)'
                            ):
    """Calculate phenotype score and p-values comparing `cond2` vs `cond1`
    """
    print(f'\t{cond2} vs {cond1}')

    count_layer = 'seq_depth_norm'
    # check if count_layer exists
    if 'seq_depth_norm' not in adata.layers.keys():
        seqDepthNormalization(adata)
    # prep counts for phenoScore calculation
    pass


def convertResultsToDataFrame(adata, targets, scores, pvalues):
    return pd.concat([
        pd.Series(targets, index=adata.var.index, name='target'),
        pd.Series(scores, index=adata.var.index, name='score'),
        pd.Series(pvalues, index=adata.var.index, name='pvalue')
    ], axis=1)


class ScreenPro(object):
    """`ScreenPro` class for processing CRISPR screen datasets
    """

    def __init__(self, adata, math='log2(x+1)'):
        self.phenotypes = {}
        self.adata = adata
        self.math = math

    def __repr__(self):
        descriptions = ''
        for scoreLevel in self.phenotypes.keys():
            scores = "', '".join(self.phenotypes[scoreLevel].columns.get_level_values(0).unique().to_list())
            descriptions += f"Phenotypes in scoreLevel: '{scoreLevel}':\n    scores: '{scores}'\n"

        return f'obs->samples\nvar->oligos\n\n{self.adata.__repr__()}\n\n{descriptions}'

    def calculateDrugScreen(self,
                            t0, untreated, treated, growth_rate,
                            scoreLevel, method):
        """Calculate gamma, rho, and tau phenotype scores using given `method` in `scoreLevel`
        """
        if method == 'ByReps':
            gamma, gamma_pv, gamma_name= runPhenoScoreByReps(
                self.adata, cond1=t0, cond2=untreated, growth_rate=growth_rate, math=self.math)
            tau, tau_pv, tau_name = runPhenoScoreByReps(
                self.adata, cond1=t0, cond2=treated, growth_rate=growth_rate, math=self.math)
            rho, rho_pv, rho_name = runPhenoScoreByReps(
                self.adata, cond1=untreated, cond2=treated, growth_rate=growth_rate, math=self.math)
            targets = self.adata.var.index.str.split('_[-,+]_').str[0].to_list()

            self.phenotypes[scoreLevel] = pd.concat({
                f'rho:{rho_name}': convertResultsToDataFrame(self.adata, targets, rho, rho_pv),
                f'gamma:{gamma_name}': convertResultsToDataFrame(self.adata, targets, gamma, gamma_pv),
                f'tau:{tau_name}': convertResultsToDataFrame(self.adata, targets, tau, tau_pv)
            }, axis=1)

        elif method == 'ByGuideSet':
            pass

        else:
            raise ValueError('Method not recognized')

    # adata.var[f'condition_{cond2}_vs_{cond1}_pvalue'] = pvalues
    # adata.var[f'condition_{cond2}_vs_{cond1}_delta'] = phenotype_score
    ## calculate FDR
    # Calculate the adjusted p-values using the Benjamini-Hochberg method
    # _, adj_pvalues, _, _ = multipletests(adata.var[f'condition_{cond1}_vs_{cond2}_pvalue'], alpha=0.05, method='fdr_bh')
    # adata.var[f'condition_{cond1}_vs_{cond2}_adj_pvalue'] = adj_pvalues


