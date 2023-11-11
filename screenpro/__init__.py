import pandas as pd
from .__version__ import __version__
from .phenoScore import seqDepthNormalization, matrixStat, matrixTest
from .phenoScore import runPhenoScore
from copy import copy


class ScreenPro(object):
    """
    `ScreenPro` class for processing CRISPR screen datasets
    """

    def __init__(self, adata, math='log2(x+1)', test='ttest'):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            math (str): math transformation to apply to the data before calculating phenotype scores
            test (str): statistical test to use for calculating phenotype scores
        """
        self.adata = adata
        self.math = math
        self.test = test
        self.phenotypes = {}

    def __repr__(self):
        descriptions = ''
        for score_level in self.phenotypes.keys():
            scores = "', '".join(self.phenotypes[score_level].columns.get_level_values(0).unique().to_list())
            descriptions += f"Phenotypes in score_level = '{score_level}':\n    scores: '{scores}'\n"

        return f'obs->samples\nvar->oligos\n\n{self.adata.__repr__()}\n\n{descriptions}'

    def copy(self):
        return copy(self)

    def calculateDrugScreen(self, t0, untreated, treated, growth_rate, score_level, n_reps=2):
        """
        Calculate gamma, rho, and tau phenotype scores for a drug screen dataset in a given `score_level`
        see this issue for discussion https://github.com/abearab/ScreenPro2/issues/15.
        Args:
            t0 (str): name of the untreated condition
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            growth_rate (float): growth rate of the untreated condition
            score_level (str): name of the score level
        """
        # calculate phenotype scores: gamma, tau, rho
        gamma_name, gamma = runPhenoScore(
            self.adata, cond1=t0, cond2=untreated, growth_rate=growth_rate, n_reps=n_reps,
            math=self.math, test=self.test, score_level=score_level
        )
        tau_name, tau = runPhenoScore(
            self.adata, cond1=t0, cond2=treated, growth_rate=growth_rate, n_reps=n_reps,
            math=self.math, test=self.test, score_level=score_level
        )
        rho_name, rho = runPhenoScore(
            self.adata, cond1=untreated, cond2=treated, growth_rate=growth_rate, n_reps=n_reps,
            math=self.math, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes[score_level] = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)

    def calculateFlowBasedScreen(self):
        """
        Calculate phenotype scores for a flow-based screen dataset
        see this issue for discussion https://github.com/abearab/ScreenPro2/issues/17
        """
        # TODO
        pass
