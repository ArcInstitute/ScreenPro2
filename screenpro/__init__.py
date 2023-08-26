import pandas as pd
from .__version__ import __version__
from .phenoScore import seqDepthNormalization, matrixStat, matrixTest
from statsmodels.stats.multitest import multipletests


def runPhenoScore(adata,
                  cond1, cond2, math, test,
                  score_level,
                  growth_rate=1, n_reps=2,
                  ctrl_label='negCtrl',
                  ):
    """Calculate phenotype score and p-values comparing `cond2` vs `cond1`
    """
    result_name = f'{cond2}_vs_{cond1}'
    print(f'\t{cond2} vs {cond1}')

    count_layer = 'seq_depth_norm'
    # check if count_layer exists
    if 'seq_depth_norm' not in adata.layers.keys():
        seqDepthNormalization(adata)

    if score_level == 'compare_reps':
        # prep counts for phenoScore calculation
        df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index[:n_reps],].to_df(count_layer).T
        df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index[:n_reps],].to_df(count_layer).T

        x = df_cond1.to_numpy()
        y = df_cond2.to_numpy()

        x_ctrl = df_cond1[adata.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond2[adata.var.targetType.eq(ctrl_label)].to_numpy()

        # calculate growth score and p_value
        scores, p_values = matrixTest(
            x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
            math=math, ave_reps=True, test=test, growth_rate=growth_rate
        )

        # Calculate the adjusted p-values using the Benjamini-Hochberg method
        _, adj_pvalues, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

        targets = adata.var.index.str.split('_[-,+]_').str[0].to_list()

        result = pd.concat([
            pd.Series(targets, index=adata.var.index, name='target'),
            pd.Series(scores, index=adata.var.index, name='score'),
            pd.Series(p_values, index=adata.var.index, name=f'{test} pvalue'),
            pd.Series(adj_pvalues, index=adata.var.index, name='BH adj_pvalue'),
        ], axis=1)

        return result_name, result

    elif score_level == 'compare_oligos':
        pass
    elif score_level == 'compare_oligos_and_reps':
        pass
    else:
        raise ValueError(f'score_level "{score_level}" not recognized')


class ScreenPro(object):
    """`ScreenPro` class for processing CRISPR screen datasets
    """

    def __init__(self, adata, math='log2(x+1)', test='ttest'):
        self.adata = adata
        self.math = math
        self.test = test
        self.phenotypes = None
        self.phenotypes_score_level = None

    def __repr__(self):
        descriptions = ''
        for scoreLevel in self.phenotypes.keys():
            scores = "', '".join(self.phenotypes[scoreLevel].columns.get_level_values(0).unique().to_list())
            descriptions += f"Phenotypes in scoreLevel = '{scoreLevel}':\n    scores: '{scores}'\n"

        return f'obs->samples\nvar->oligos\n\n{self.adata.__repr__()}\n\n{descriptions}'

    def calculateDrugScreen(self,
                            t0, untreated, treated, growth_rate,
                            score_level):
        """Calculate gamma, rho, and tau phenotype scores for a drug screen dataset in a given `score_level`
        """
        # calculate phenotype scores: gamma, tau, rho
        gamma_name, gamma = runPhenoScore(
            self.adata, cond1=t0, cond2=untreated, growth_rate=growth_rate,
            math=self.math, test=self.test, score_level=score_level
        )
        tau_name, tau = runPhenoScore(
            self.adata, cond1=t0, cond2=treated, growth_rate=growth_rate,
            math=self.math, test=self.test, score_level=score_level
        )
        rho_name, rho = runPhenoScore(
            self.adata, cond1=untreated, cond2=treated, growth_rate=growth_rate,
            math=self.math, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)

        self.phenotypes_score_level = score_level
