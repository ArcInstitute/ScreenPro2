## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

"""Assays module

"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from pydeseq2 import preprocessing
from .phenoscore.deseq import runDESeq, extractDESeqResults
from .phenoscore import runPhenoScore, runPhenoScoreForReplicate
from .phenoscore.annotate import ann_score_df
from copy import copy


class PooledScreens(object):
    """
    pooledScreens class for processing CRISPR screen datasets
    """

    def __init__(self, adata, fc_transformation='log2', test='ttest', n_reps=3):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            fc_transformation (str): fold change transformation to apply for calculating phenotype scores
            test (str): statistical test to use for calculating phenotype scores
        """
        self.adata = adata.copy()
        self.pdata = None
        self.fc_transformation = fc_transformation
        self.test = test
        self.n_reps = n_reps
        self.phenotypes = {}
        self.phenotype_names = []

    # def __repr__(self):
    #     descriptions = ''
    #     for score_level in self.phenotypes.keys():
    #         scores = "', '".join(self.phenotypes[score_level].columns.get_level_values(0).unique().to_list())
    #         descriptions += f"Phenotypes in score_level = '{score_level}':\n    scores: '{scores}'\n"

    #     return f'obs->samples\nvar->elementss\n\n{self.__repr__()}\n\n{descriptions}'

    def copy(self):
        return copy(self)
    
    def _add_phenotype_results(self, phenotype_name):
        if phenotype_name in self.phenotype_names:
            raise ValueError(f"Phenotype '{phenotype_name}' already exists in self.phenotype_names!")
        self.phenotype_names.append(phenotype_name)

    def _calculateGrowthFactor(self, untreated, treated, db_rate_col):
        """
        Calculate growth factor for gamma, tau, or rho score per replicates.

        Parameters:
            untreated (str): untreated condition
            treated (str): treated condition
            db_rate_col (str): column name for doubling rate
        
        Returns:
            pd.DataFrame: growth factor dataframe
        """
        adat = self.adata.copy()
        growth_factors = []
        # calculate growth factor for gamma, tau, or rho score per replicates
        for replicate in adat.obs.replicate.unique():
            db_untreated = adat.obs.query(f'condition == "{untreated}" & replicate == {str(replicate)}')[db_rate_col][0]
            db_treated = adat.obs.query(f'condition == "{treated}" & replicate == {str(replicate)}')[db_rate_col][0]

            growth_factors.append(('gamma', db_untreated, replicate, f'gamma_replicate_{replicate}'))
            growth_factors.append(('tau', db_treated, replicate, f'tau_replicate_{replicate}'))
            growth_factors.append(('rho', np.abs(db_untreated - db_treated), replicate, f'rho_replicate_{replicate}'))

        out = pd.DataFrame(growth_factors, columns=['score', 'growth_factor', 'replicate', 'index']).set_index('index')
        out.index.name = None

        return out

    def countNormalization(self):
        """
        Normalize the counts data in adata.X
        """
        self.adata.layers['raw_counts'] = self.adata.X.copy()
        
        # normalize counts by sequencing depth
        norm_counts, size_factors = preprocessing.deseq2_norm(self.adata.X)
        # update adata object
        self.adata.obs['size_factors'] = size_factors
        self.adata.layers['seq_depth_norm'] = norm_counts
        self.adata.X = self.adata.layers['seq_depth_norm']
    
    def calculateDrugScreenDESeq(self, t0, untreated, treated, run_name=None, **kwargs):
        """
        Calculate DESeq2 results for a given drug screen dataset.

        Args:
            design (str): design matrix for DESeq2
            run_name (str): name for the DESeq2 calculation run
            **kwargs: additional arguments to pass to runDESeq
        """
        dds = runDESeq(self.adata, 'condition', **kwargs)
        
        # Calculate `gamma`, `rho`, and `tau` phenotype scores
        gamma_name, gamma = extractDESeqResults(
            dds, 'condition', t0, untreated, **kwargs
        )

        tau_name, tau = extractDESeqResults(
            dds, 'condition', t0, treated, **kwargs
        )

        rho_name, rho = extractDESeqResults(
            dds, 'condition', untreated, treated, **kwargs
        )

        if not run_name: run_name = 'pyDESeq2'

        self.phenotypes[run_name] = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)        

    def calculateDrugScreen(self, t0, untreated, treated, db_untreated, db_treated, score_level, db_rate_col='pop_doublings', run_name=None, **kwargs):
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
            db_rate_col (str): column name for the doubling rate, default is 'pop_doublings'
            run_name (str): name for the phenotype calculation run
            **kwargs: additional arguments to pass to runPhenoScore
        """
        # calculate phenotype scores: gamma, tau, rho
        gamma_name, gamma = runPhenoScore(
            self.adata, cond1=t0, cond2=untreated, growth_rate=db_untreated,
            n_reps=self.n_reps,
            transformation=self.fc_transformation, test=self.test, score_level=score_level,
            **kwargs
        )
        tau_name, tau = runPhenoScore(
            self.adata, cond1=t0, cond2=treated, growth_rate=db_treated,
            n_reps=self.n_reps,
            transformation=self.fc_transformation, test=self.test, score_level=score_level,
            **kwargs
        )
        # TO-DO: warning / error if db_untreated and db_treated are too close, i.e. growth_rate ~= 0.
        rho_name, rho = runPhenoScore(
            self.adata, cond1=untreated, cond2=treated, growth_rate=np.abs(db_untreated - db_treated),
            n_reps=self.n_reps,
            transformation=self.fc_transformation, test=self.test, score_level=score_level,
            **kwargs
        )

        if not run_name: run_name = score_level
        # save all results into a multi-index dataframe
        self.phenotypes[run_name] = pd.concat({
            f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        }, axis=1)

        # save phenotype name for reference
        self._add_phenotype_results(f'gamma:{gamma_name}')
        self._add_phenotype_results(f'tau:{tau_name}')
        self._add_phenotype_results(f'rho:{rho_name}')

        growth_factor_table = self._calculateGrowthFactor(
            untreated = untreated, treated = treated, db_rate_col = db_rate_col
        )

        # get replicate level phenotype scores
        pdata_df = pd.concat([
            runPhenoScoreForReplicate(
                self.adata, x_label = x_label, y_label = y_label, score = score_label,
                transformation=self.fc_transformation, 
                growth_factor_table=growth_factor_table,
                **kwargs
            ).add_prefix(f'{score_label}_')

            for x_label, y_label, score_label in [
                ('T0', untreated, 'gamma'),
                ('T0', treated, 'tau'),
                (untreated, treated, 'rho')
            ]
        ],axis=1).T
        # add .pdata
        self.pdata = ad.AnnData(
            X = pdata_df,
            obs = growth_factor_table.loc[pdata_df.index,:],
            var=self.adata.var
        )
        
    def calculateFlowBasedScreen(self, low_bin, high_bin, score_level, run_name=None, **kwargs):
        """
        Calculate phenotype scores for a flow-based screen dataset.

        Args:
            low_bin (str): name of the low bin condition
            high_bin (str): name of the high bin condition
            score_level (str): name of the score level
            run_name (str): name for the phenotype calculation run
            **kwargs: additional arguments to pass to runPhenoScore
        """
        # calculate phenotype scores
        delta_name, delta = runPhenoScore(
            self.adata, cond1=low_bin, cond2=high_bin, n_reps=self.n_reps,
            transformation=self.fc_transformation, test=self.test, score_level=score_level,
            **kwargs
        )

        if not run_name: run_name = score_level
        # save all results into a multi-index dataframe
        self.phenotypes[run_name] = pd.concat({
            f'delta:{delta_name}': delta
        }, axis=1)

        # save phenotype name for reference
        self._add_phenotype_results(f'delta:{delta_name}')

    def getPhenotypeScores(self, score_name, run_name='auto', threshold=5, ctrl_label='control', target_col='target',pvalue_column='ttest pvalue', score_column='score'):
        """
        Get phenotype scores for a given score level

        Args:
            run_name (str): name of the phenotype calculation run to retrieve
            score_name (str): name of the score to retrieve, e.g. 'gamma', 'tau', 'rho', 'delta'
            threshold (float): threshold for filtering significant hits, default is 5
            ctrl_label (str): label for the negative control, default is 'control'
            target_col (str): column name for the target gene, default is 'target'
            pvalue_column (str): column name for the p-value, default is 'ttest pvalue'
            score_column (str): column name for the score, default is 'score'
        """
        hit_dict = {
            'gamma':{
                'up_hit':'up_hit',
                'down_hit':'essential_hit'
            },
            'tau':{
                'up_hit':'up_hit', 
                'down_hit':'down_hit'
            },
            'rho':{
                'up_hit':'resistance_hit', 
                'down_hit':'sensitivity_hit'
            }
        }

        if run_name == 'auto':
            if len(list(self.phenotypes.keys())) == 1:
                run_name = list(self.phenotypes.keys())[0]
            else:
                raise ValueError(
                    'Multiple phenotype calculation runs found.'
                    'Please specify run_name. Available runs: '
                    '' + ', '.join(self.phenotypes.keys())
                )
        
        if score_name not in self.phenotype_names:
            raise ValueError(f"Phenotype '{score_name}' not found in self.phenotype_names")

        keep_col = [target_col, score_column, pvalue_column]
        score_tag = score_name.split(':')[0]
        out = ann_score_df(
            self.phenotypes[run_name][score_name].loc[:,keep_col],
            ctrl_label=ctrl_label, 
            up_hit=hit_dict[score_tag]['up_hit'],
            down_hit=hit_dict[score_tag]['down_hit'],
            threshold=threshold
        )

        return out

    def getAnnotatedTable(self, run_name='auto', threshold=5, ctrl_label='control', target_col='target',pvalue_column='ttest pvalue', score_column='score'):
        hit_dict = {
            'gamma':{
                'up_hit':'up_hit',
                'down_hit':'essential_hit'
            },
            'tau':{
                'up_hit':'up_hit', 
                'down_hit':'down_hit'
            },
            'rho':{
                'up_hit':'resistance_hit', 
                'down_hit':'sensitivity_hit'
            }
        }
        
        if run_name == 'auto':
            if len(list(self.phenotypes.keys())) == 1:
                run_name = list(self.phenotypes.keys())[0]
            else:
                raise ValueError(
                    'Multiple phenotype calculation runs found.'
                    'Please specify run_name. Available runs: '
                    '' + ', '.join(self.phenotypes.keys())
                )

        keep_col = [target_col, score_column, pvalue_column]

        score_names = {s for s, col in self.phenotypes[run_name].columns}
        sort_var = self.adata.var.sort_values(['targetType','target']).index.to_list()
        
        df_list = {}
        for score_name in score_names:
            score_tag = score_name.split(':')[0]
            # get label
            df_label = ann_score_df(
                self.phenotypes[run_name][score_name].loc[:,keep_col],
                up_hit=hit_dict[score_tag]['up_hit'],
                down_hit=hit_dict[score_tag]['down_hit'],
                ctrl_label=ctrl_label,
                threshold=threshold
            )['label']
            # get replicate phe
            df_phe_reps = self.pdata[self.pdata.obs.score.eq(score_tag)].to_df().T
            
            # make table
            df = pd.concat([
                self.phenotypes['compare_reps'][score_name], df_phe_reps, df_label
            ],axis=1).loc[sort_var,:]
            
            df_list.update({score_name:df})
        
        out = pd.concat(df_list,axis=1)
        
        return out


class GImaps(object):
    pass