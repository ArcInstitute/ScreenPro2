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

from ..phenoscore import runDESeq, extractDESeqResults
from ..phenoscore import runPhenoScore, runPhenoScoreForReplicate
from ..preprocessing import addPseudoCount, findLowCounts, normalizeSeqDepth
from ..phenoscore.annotate import annotateScoreTable, hit_dict

import warnings
from copy import copy


class PooledScreens(object):
    """
    pooledScreens class for processing CRISPR screen datasets
    """

    def __init__(self, adata, fc_transformation='log2', test='ttest', n_reps=3, verbose=False):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            fc_transformation (str): fold change transformation to apply for calculating phenotype scores
            test (str): statistical test to use for calculating phenotype scores
            n_reps (int): number of replicates to use for calculating phenotype scores
            verbose (bool): whether to print verbose output
        """
        self.adata = adata.copy()
        self.pdata = None
        self.fc_transformation = fc_transformation
        self.test = test
        self.n_reps = n_reps
        self.phenotypes = {}
        self.verbose = verbose

    def copy(self):
        return copy(self)
    
    def _add_phenotype_results(self, run_name, phenotype_name, phenotype_table):
        if phenotype_name in self.phenotypes[run_name]['results'].keys():
            raise ValueError(f"Phenotype '{phenotype_name}' already exists in self.phenotypes['results']!")
        self.phenotypes[run_name]['results'][phenotype_name] = phenotype_table

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

    def _getTreatmentDoublingRate(self, untreated, treated, db_rate_col):
        if 'pop_doublings' not in self.adata.obs.columns or db_rate_col == None:
            warnings.warn('No doubling rate information provided.')
            db_untreated = 1
            db_treated = 1
            db_diff = 1
            growth_factor_table = None
        
        else:
            growth_factor_table = self._calculateGrowthFactor(
                untreated = untreated, treated = treated, db_rate_col = db_rate_col
            )

            db_untreated=growth_factor_table.query(f'score=="gamma"')['growth_factor'].mean()
            db_treated=growth_factor_table.query(f'score=="tau"')['growth_factor'].mean()
            db_diff = np.abs(db_untreated - db_treated)
        
        return db_untreated, db_treated, db_diff

    def filterLowCounts(self, filter_type='all', minimum_reads=50):
        """
        Filter low counts in adata.X
        """
        findLowCounts(
            self.adata, 
            filter_type=filter_type, 
            minimum_reads=minimum_reads,
            verbose=self.verbose
        )

        self.adata = self.adata[:,~self.adata.var.low_count].copy()

    def countNormalization(self, pseudo_count_value=0.5):
        """
        Preprocess and normalize the counts data in adata.X

        Steps:
            1. Add pseudocount to counts
            2. Normalize counts by sequencing depth
            3. Log10 transformation
        
        """
        self.adata.layers['raw_counts'] = self.adata.X.copy()
        
        # add pseudocount
        addPseudoCount(self.adata, behavior='default', value=pseudo_count_value)
        
        if self.verbose: print('Pseudocount added to counts.')
        
        # normalize counts by sequencing depth
        normalizeSeqDepth(self.adata)

        if self.verbose: print('Counts normalized by sequencing depth.')

        # log scale the counts
        self.adata.X = np.log10(self.adata.X)
        
        if self.verbose: print('`log10` transformation applied to counts.')

    def calculateDrugScreenDESeq(self, untreated, treated, t0=None, run_name='pyDESeq2', **kwargs):
        """
        Calculate DESeq2 results for a given drug screen dataset.

        Args:
            design (str): design matrix for DESeq2-based analysis
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            t0 (str): name of the untreated condition
            run_name (str): name for the phenotype calculation run
            **kwargs: additional arguments to pass to runDESeq
        """
        if run_name in self.phenotypes.keys():
            raise ValueError(f"Phenotype calculation run '{run_name}' already exists in self.phenoypes!")
        else:
            self.phenotypes[run_name] = {}
            
            self.phenotypes[run_name]['config'] = {
                    'method':'pyDESeq2',
                    'untreated':untreated,
                    'treated':treated,
                    't0':t0,
                    'n_reps':self.n_reps,
            }
            self.phenotypes[run_name]['results'] = {}

        if type(treated) != list: treated = [treated]

        # run pyDESeq2 analysis
        dds = runDESeq(self.adata, 'condition', **kwargs)

        # extract comparison results
        if t0 != None and type(t0) == str:
            # Calculate `gamma`, `rho`, and `tau` phenotype scores
            gamma_name, gamma = extractDESeqResults(
                dds, 'condition', t0, untreated, **kwargs
            )
            self._add_phenotype_results(run_name, f'gamma:{gamma_name}', gamma)

            for tr in treated:
                tau_name, tau = extractDESeqResults(
                    dds, 'condition', t0, treated, **kwargs
                )
                self._add_phenotype_results(run_name, f'tau:{tau_name}', tau)

        for tr in treated:
            rho_name, rho = extractDESeqResults(
                dds, 'condition', untreated, tr, **kwargs
            )
            self._add_phenotype_results(run_name, f'rho:{rho_name}', rho)

    def calculateDrugScreen(self, score_level, untreated, treated, t0=None, db_rate_col='pop_doublings', run_name=None, **kwargs):
        """
        Calculate `gamma`, `rho`, and `tau` phenotype scores for a drug screen dataset in a given `score_level`.

        Args:
            score_level (str): name of the score level
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            t0 (str): name of the untreated condition
            db_rate_col (str): column name for the doubling rate, default is 'pop_doublings'
            run_name (str): name for the phenotype calculation run
            **kwargs: additional arguments to pass to runPhenoScore
        """
        if not run_name: run_name = score_level
        if run_name in self.phenotypes.keys():
            raise ValueError(f"Phenotype calculation run '{run_name}' already exists in self.phenoypes!")
        else:
            self.phenotypes[run_name] = {}
            self.phenotypes[run_name]['config'] = {
                'method':'ScreenPro2 - phenoscore',
                'untreated':untreated,
                'treated':treated,
                't0':t0,
                'n_reps':self.n_reps,
                'test':self.test,
                'score_level':score_level,
            }
            self.phenotypes[run_name]['results'] = {}

        if type(treated) != list: treated = [treated]

        # calculate phenotype scores: gamma, tau, rho
        if t0 != None and type(t0) == str:
            db_untreated,_,_ = self._getTreatmentDoublingRate(untreated, treated[0], db_rate_col)
            gamma_name, gamma = runPhenoScore(
                self.adata, cond_ref=t0, cond_test=untreated, growth_rate=db_untreated,
                n_reps=self.n_reps,
                transformation=self.fc_transformation, test=self.test, score_level=score_level,
                **kwargs
            )
            self._add_phenotype_results(run_name, f'gamma:{gamma_name}', gamma)

        for tr in treated:
            _, db_tr, db_diff = self._getTreatmentDoublingRate(untreated, tr, db_rate_col)

            if t0 != None and type(t0) == str:
                tau_name, tau = runPhenoScore(
                    self.adata, cond_ref=t0, cond_test=tr, growth_rate=db_tr,
                    n_reps=self.n_reps,
                    transformation=self.fc_transformation, test=self.test, score_level=score_level,
                    **kwargs
                )
                self._add_phenotype_results(run_name, f'tau:{tau_name}', tau)
            
            #TODO: warning / error if db_untreated and db_treated are too close, i.e. growth_rate ~= 0.
            rho_name, rho = runPhenoScore(
                self.adata, cond_ref=untreated, cond_test=tr, growth_rate=db_diff,
                n_reps=self.n_reps,
                transformation=self.fc_transformation, test=self.test, score_level=score_level,
                **kwargs
            )
            self._add_phenotype_results(run_name, f'rho:{rho_name}', rho)

        # gnerate replicate level phenotype scores
        pdata_dict = {}
        for score_name in self.phenotypes[score_level]['results'].keys():
            score_label, comparison = score_name.split(':')
            y_label, x_label = comparison.split('_vs_')
            
            #TODO: get growth rates for replicate level scores
            
        
            pdata_dict.update({
                score_name: runPhenoScoreForReplicate(
                    self.adata, x_label = x_label, y_label = y_label,
                    transformation=self.fc_transformation, 
                    # growth_factor_reps=
                    # **kwargs
                ).add_prefix(f'{score_label}_').T # transpose to match pdata format
            })

        pdata_df = pd.concat(pdata_dict, axis=0)

        #TODO: fix `_calculateGrowthFactor` and `_getTreatmentDoublingRate` to maintain same format
        # add .pdata
        self.pdata = ad.AnnData(
            X = pdata_df,
            # obs = growth_factor_table.loc[pdata_df.index,:],
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
            self.adata, cond_ref=low_bin, cond_test=high_bin, n_reps=self.n_reps,
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

    def listPhenotypeScores(self, run_name='auto'):
        """
        List available phenotype scores for a given run_name

        Args:
            run_name (str): name of the phenotype calculation run to retrieve
        """
        if run_name == 'auto':
            if len(list(self.phenotypes.keys())) == 1:
                run_name = list(self.phenotypes.keys())[0]
            else:
                raise ValueError(
                    'Multiple phenotype calculation runs found.'
                    'Please specify run_name. Available runs: '
                    '' + ', '.join(self.phenotypes.keys())
                )

        out = list(self.phenotypes[run_name]['results'].keys())

        return out
    
    def getPhenotypeScores(self, score_name, threshold, run_name='auto', ctrl_label='negative_control', target_col='target',pvalue_col='ttest pvalue', score_col='score'):
        """
        Get phenotype scores for a given score level

        Args:
            score_name (str): name of the score to retrieve, e.g. 'gamma', 'tau', 'rho', 'delta'
            threshold (float): threshold for filtering significant hits, default is 5
            run_name (str): name of the phenotype calculation run to retrieve
            ctrl_label (str): label for the negative control, default is 'negative_control'
            target_col (str): column name for the target gene, default is 'target'
            pvalue_column (str): column name for the p-value, default is 'ttest pvalue'
            score_column (str): column name for the score, default is 'score'
        """

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

        keep_col = [target_col, score_col, pvalue_col]
        score_tag = score_name.split(':')[0]
        out = annotateScoreTable(
            self.phenotypes[run_name][score_name].loc[:,keep_col],
            threshold=threshold,
            up_hit=hit_dict[score_tag]['up_hit'],
            down_hit=hit_dict[score_tag]['down_hit'],
            ctrl_label=ctrl_label, 
            score_col=score_col,
            pvalue_col=pvalue_col
        )

        return out

    def getAnnotatedTable(self, threshold, run_name='auto', ctrl_label='negative_control', target_col='target', pvalue_col='ttest pvalue', score_col='score'):
        """
        Returns an annotated table with scores, labels, and replicate phenotypes.

        Args:
            threshold (int, optional): The threshold value for determining hits. Defaults to 5.
            run_name (str, optional): The name of the phenotype calculation run. Defaults to 'auto'.
            ctrl_label (str, optional): The label for the control group. Defaults to 'negative_control'.
            target_col (str, optional): The column name for the target. Defaults to 'target'.
            pvalue_column (str, optional): The column name for the p-value. Defaults to 'ttest pvalue'.
            score_column (str, optional): The column name for the score. Defaults to 'score'.

        Returns:
            pandas.DataFrame: An annotated table with scores, labels, and replicate phenotypes.
        """
        if run_name == 'auto':
            if len(list(self.phenotypes.keys())) == 1:
                run_name = list(self.phenotypes.keys())[0]
            else:
                raise ValueError(
                    'Multiple phenotype calculation runs found.'
                    'Please specify run_name. Available runs: '
                    '' + ', '.join(self.phenotypes.keys())
                )

        keep_col = [target_col, score_col, pvalue_col]
        
        # self.phenotypes[run_name] = pd.concat({
        #     f'gamma:{gamma_name}': gamma, f'tau:{tau_name}': tau, f'rho:{rho_name}': rho
        # }, axis=1)

        score_names = set(self.phenotypes[run_name]['results'].keys())
        sort_var = self.adata.var.sort_values(['targetType','target']).index.to_list()

        df_list = {}
        for score_name in score_names:
            score_tag = score_name.split(':')[0]
            
            # get annotated table
            df_ann = annotateScoreTable(
                self.phenotypes[run_name][score_name].loc[:,keep_col],
                up_hit=hit_dict[score_tag]['up_hit'],
                down_hit=hit_dict[score_tag]['down_hit'],
                score_col=score_col,
                pvalue_col=pvalue_col,
                ctrl_label=ctrl_label,
                threshold=threshold
            )

            # get replicate phe
            df_phe_reps = self.pdata[self.pdata.obs.score.eq(score_tag)].to_df().T

            # make table
            df = pd.concat([
                df_ann.drop(columns=['label']),
                df_phe_reps, 
                df_ann['label']
            ],axis=1).loc[sort_var,:]

            df_list.update({score_name:df})

        out = pd.concat(df_list,axis=1)

        return out


class GImaps(object):
    pass