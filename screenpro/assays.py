"""
Assays module
"""

import numpy as np
import pandas as pd

from .phenoscore import runPhenoScore

from copy import copy


class PooledScreens(object):
    """
    pooledScreens class for processing CRISPR screen datasets
    """

    def __init__(self, adata, transformation='log2(x+1)', test='ttest', n_reps=3):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            transformation (str): transformation to apply to the data before calculating phenotype scores
            test (str): statistical test to use for calculating phenotype scores
        """
        self.adata = adata
        self.transformation = transformation
        self.test = test
        self.n_reps = n_reps
        self.phenotypes = {}

    # def __repr__(self):
    #     descriptions = ''
    #     for score_level in self.phenotypes.keys():
    #         scores = "', '".join(self.phenotypes[score_level].columns.get_level_values(0).unique().to_list())
    #         descriptions += f"Phenotypes in score_level = '{score_level}':\n    scores: '{scores}'\n"

    #     return f'obs->samples\nvar->elementss\n\n{self.__repr__()}\n\n{descriptions}'

    def copy(self):
        return copy(self)

    def calculateDrugScreen(self, t0, untreated, treated, db_untreated, db_treated, score_level):
        """
        Calculate `gamma`, `rho`, and `tau` phenotype scores for a drug screen dataset in a given `score_level`.
        To normalize by growth rate, the doubling rate of the untreated and treated conditions are required.

        Args:
            t0 (str): name of the untreated condition
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            db_untreated (float): doubling rate of the untreated condition
            db_treated (float): doubling rate of the treated condition
            score_level (str): name of the score level
        """
        # calculate phenotype scores: gamma, tau, rho
        gamma_name, gamma = runPhenoScore(
            self.adata, cond1=t0, cond2=untreated, growth_rate=db_untreated,
            n_reps=self.n_reps,
            transformation=self.transformation, test=self.test, score_level=score_level
        )
        tau_name, tau = runPhenoScore(
            self.adata, cond1=t0, cond2=treated, growth_rate=db_treated,
            n_reps=self.n_reps,
            transformation=self.transformation, test=self.test, score_level=score_level
        )
        # TO-DO: warning / error if db_untreated and db_treated are too close, i.e. growth_rate ~= 0.
        rho_name, rho = runPhenoScore(
            self.adata, cond1=untreated, cond2=treated, growth_rate=np.abs(db_untreated - db_treated),
            n_reps=self.n_reps,
            transformation=self.transformation, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes[score_level] = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)

    def calculateFlowBasedScreen(self, low_bin, high_bin, score_level):
        """
        Calculate phenotype scores for a flow-based screen dataset.

        Args:
            low_bin (str): name of the low bin condition
            high_bin (str): name of the high bin condition
            score_level (str): name of the score level
        """
        # calculate phenotype scores
        phenotype_name, phenotype = runPhenoScore(
            self.adata, cond1=low_bin, cond2=high_bin, n_reps=self.n_reps,
            transformation=self.transformation, test=self.test, score_level=score_level
        )

        # save all results into a multi-index dataframe
        self.phenotypes[score_level] = pd.concat({
            f'phenotype:{phenotype_name}': phenotype
        }, axis=1)


class GImaps(object):
    pass