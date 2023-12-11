import numpy as np
import pandas as pd
from . import plotting as pl
from . import phenoScore as ps
from . import utils

from .__version__ import __version__
from copy import copy


class ScreenPro(object):
    """
    `ScreenPro` class for processing CRISPR screen datasets
    """

    def __init__(self, adata, math='log2(x+1)', test='ttest', n_reps=3):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            math (str): math transformation to apply to the data before calculating phenotype scores
            test (str): statistical test to use for calculating phenotype scores
        """
        self.adata = adata
        self.math = math
        self.test = test
        self.n_reps = n_reps
        self.phenotypes = {}

    def __repr__(self):
        descriptions = ''
        for score_level in self.phenotypes.keys():
            scores = "', '".join(self.phenotypes[score_level].columns.get_level_values(0).unique().to_list())
            descriptions += f"Phenotypes in score_level = '{score_level}':\n    scores: '{scores}'\n"

        return f'obs->samples\nvar->oligos\n\n{self.adata.__repr__()}\n\n{descriptions}'

    def copy(self):
        return copy(self)

    def calculateDrugScreen(self, t0, untreated, treated, db_untreated, db_treated, score_level):
        """
        Calculate gamma, rho, and tau phenotype scores for a drug screen dataset in a given `score_level`
        see this issue for discussion https://github.com/abearab/ScreenPro2/issues/15.
        Args:
            t0 (str): name of the untreated condition
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            db_untreated (float): doubling rate of the untreated condition
            db_treated (float): doubling rate of the treated condition
            score_level (str): name of the score level
        """
        # calculate phenotype scores: gamma, tau, rho
        gamma_name, gamma = ps.runPhenoScore(
            self.adata, cond1=t0, cond2=untreated, growth_rate=db_untreated,
            n_reps=self.n_reps,
            math=self.math, test=self.test, score_level=score_level
        )
        tau_name, tau = ps.runPhenoScore(
            self.adata, cond1=t0, cond2=treated, growth_rate=db_treated,
            n_reps=self.n_reps,
            math=self.math, test=self.test, score_level=score_level
        )
        # TO-DO: warning / error if db_untreated and db_treated are too close, i.e. growth_rate ~= 0.
        rho_name, rho = ps.runPhenoScore(
            self.adata, cond1=untreated, cond2=treated, growth_rate=np.abs(db_untreated - db_treated),
            n_reps=self.n_reps,
            math=self.math, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes[score_level] = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)

    def calculateFlowBasedScreen(self, low_bin, high_bin, score_level):
        """
        Calculate phenotype scores for a flow-based screen dataset
        see this issue for discussion https://github.com/abearab/ScreenPro2/issues/17
        """
        # calculate phenotype scores
        phenotype_name, phenotype = ps.runPhenoScore(
            self.adata, cond1=low_bin, cond2=high_bin, n_reps=self.n_reps,
            math=self.math, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes[score_level] = pd.concat({
            f'phenotype:{phenotype_name}': phenotype
        }, axis=1)
