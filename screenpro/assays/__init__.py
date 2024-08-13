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

from ..phenoscore import (
    runPhenoScore, getPhenotypeData,
    runDESeq, extractDESeqResults,
)
from ..preprocessing import addPseudoCount, findLowCounts, normalizeSeqDepth
from ..phenoscore._annotate import annotateScoreTable, hit_dict
from ..plotting import volcano_plot, label_resistance_hit, label_sensitivity_hit

import warnings
from copy import copy


class PooledScreens(object):
    """
    pooledScreens class for processing CRISPR screen datasets
    """

    def __init__(self, adata, test='ttest', n_reps=3, verbose=False):
        """
        Args:
            adata (AnnData): AnnData object with adata.X as a matrix of sgRNA counts
            test (str): statistical test to use for calculating phenotype scores
            n_reps (int): number of replicates to use for calculating phenotype scores
            verbose (bool): whether to print verbose output
        """
        self.adata = adata.copy()
        self.pdata = None
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
        if 'pop_doubling' not in self.adata.obs.columns or db_rate_col == None:
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

    def _auto_run_name(self):
        if len(list(self.phenotypes.keys())) == 1:
            run_name = list(self.phenotypes.keys())[0]
        else:
            raise ValueError(
                'Multiple phenotype calculation runs found.'
                'Please specify run_name. Available runs: '
                '' + ', '.join(self.phenotypes.keys())
            )
        return run_name

    def filterLowCounts(self, filter_type='all', minimum_reads=1):
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
        
        """
        self.adata.layers['raw_counts'] = self.adata.X.copy()
        
        # add pseudocount
        addPseudoCount(self.adata, behavior='default', value=pseudo_count_value)
        
        if self.verbose: print('Pseudocount added to counts.')
        
        # normalize counts by sequencing depth
        normalizeSeqDepth(self.adata)

        if self.verbose: print('Counts normalized by sequencing depth.')

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
        adt = self.adata.copy()
        adt.X = adt.layers['raw_counts']
        dds = runDESeq(adt, 'condition', **kwargs)

        # extract comparison results
        if t0 != None and type(t0) == str:
            # Calculate `gamma`, `rho`, and `tau` phenotype scores
            gamma_name, gamma = extractDESeqResults(
                dds, design='condition', ref_level=t0, tested_level=untreated, **kwargs
            )
            self._add_phenotype_results(run_name, f'gamma:{gamma_name}', gamma)

            for tr in treated:
                tau_name, tau = extractDESeqResults(
                    dds, design='condition', ref_level=t0, tested_level=tr, **kwargs
                )
                self._add_phenotype_results(run_name, f'tau:{tau_name}', tau)

        for tr in treated:
            rho_name, rho = extractDESeqResults(
                dds, design='condition', ref_level=untreated, tested_level=tr, **kwargs
            )
            self._add_phenotype_results(run_name, f'rho:{rho_name}', rho)

    def calculateDrugScreen(self, score_level, untreated, treated, t0=None, db_rate_col='pop_doubling', run_name=None, **kwargs):
        """
        Calculate `gamma`, `rho`, and `tau` phenotype scores for a drug screen dataset in a given `score_level`.

        Args:
            score_level (str): name of the score level
            untreated (str): name of the untreated condition
            treated (str): name of the treated condition
            t0 (str): name of the untreated condition
            db_rate_col (str): column name for the doubling rate, default is 'pop_doubling'
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
                test=self.test, score_level=score_level,
                **kwargs
            )
            self._add_phenotype_results(run_name, f'gamma:{gamma_name}', gamma)

        for tr in treated:
            _, db_tr, db_diff = self._getTreatmentDoublingRate(untreated, tr, db_rate_col)

            if t0 != None and type(t0) == str:
                tau_name, tau = runPhenoScore(
                    self.adata, cond_ref=t0, cond_test=tr, growth_rate=db_tr,
                    n_reps=self.n_reps,
                    test=self.test, score_level=score_level,
                    **kwargs
                )
                self._add_phenotype_results(run_name, f'tau:{tau_name}', tau)
            
            #TODO: warning / error if db_untreated and db_treated are too close, i.e. growth_rate ~= 0.
            rho_name, rho = runPhenoScore(
                self.adata, cond_ref=untreated, cond_test=tr, growth_rate=db_diff,
                n_reps=self.n_reps,
                test=self.test, score_level=score_level,
                **kwargs
            )
            self._add_phenotype_results(run_name, f'rho:{rho_name}', rho)
        
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
            test=self.test, score_level=score_level,
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
        if run_name == 'auto': run_name = self._auto_run_name()

        out = list(self.phenotypes[run_name]['results'].keys())

        return out
    
    def getPhenotypeScores(self, phenotype_name, threshold, run_name='auto', **kwargs):
        """
        Get phenotype scores for a given phenotype_name

        Args:
            phenotype_name (str): name of the phenotype score
            run_name (str): name of the phenotype calculation run to retrieve
        """
        if run_name == 'auto': run_name = self._auto_run_name()
        
        score_tag, _ = phenotype_name.split(':')

        out = annotateScoreTable(
            self.phenotypes[run_name]['results'][phenotype_name],
            up_hit=hit_dict[score_tag]['up_hit'],
            down_hit=hit_dict[score_tag]['down_hit'],
            threshold=threshold,
            **kwargs
        )

        return out

    def buildPhenotypeData(self, run_name='auto',db_rate_col='pop_doubling', **kwargs):
        if run_name == 'auto': run_name = self._auto_run_name()
        if run_name=='compare_reps':
            pass
        else:
            raise ValueError('Only `compare_reps` run_name is supported for now!')
        
        untreated = self.phenotypes[run_name]['config']['untreated']
        treated = self.phenotypes[run_name]['config']['treated']

        if type(treated) != list: treated = [treated]

        if db_rate_col:
            #TODO: fix `_calculateGrowthFactor` and `_getTreatmentDoublingRate`
            growth_factor_table = self._calculateGrowthFactor(
                untreated = untreated, treated = treated, 
                db_rate_col = db_rate_col
            )
        
        pdata_list = []

        for phenotype_name in self.listPhenotypeScores(run_name=run_name):

            score_tag, comparison = phenotype_name.split(':')
            cond_test, cond_ref = comparison.split('_vs_')

            if db_rate_col:
                growth_rate_reps=growth_factor_table.query(
                    f'score=="{score_tag}"'
                ).set_index('replicate')['growth_factor'].to_dict()
            else:
                growth_rate_reps=None
            
            pdata = getPhenotypeData(
                self.adata, score_tag=score_tag, 
                cond_ref=cond_ref, cond_test=cond_test, 
                growth_rate_reps=growth_rate_reps,
                **kwargs
            )
            # obs = growth_factor_table.loc[pdata_df.index,:],

            pdata_list.append(pdata)

        self.pdata = ad.concat(pdata_list, axis=0)
        self.pdata.var = self.adata.var.copy()
    
    def drawVolcano(
            self, ax,
            phenotype_name,
            threshold,
            dot_size=1,
            run_name='auto',
            score_col='score', 
            pvalue_col='pvalue',
            xlabel='auto',
            ylabel='-log10(pvalue)',
            xlims='auto',
            ylims='auto',
            ctrl_label='negative_control',
            resistance_hits=None,
            sensitivity_hits=None,
            size_txt=None,
            t_x=0, t_y=0,
            **args
            ):
        if run_name == 'auto': run_name = self._auto_run_name()
        
        score_tag, _ = phenotype_name.split(':')

        df = self.phenotypes[run_name]['results'][phenotype_name].dropna()

        df = annotateScoreTable(
            df, 
            up_hit=hit_dict[score_tag]['up_hit'], 
            down_hit=hit_dict[score_tag]['down_hit'],
            score_col=score_col, pvalue_col=pvalue_col,
            ctrl_label=ctrl_label,
            threshold=threshold,
        )

        df['-log10(pvalue)'] = -np.log10(df[pvalue_col])

        if xlabel == 'auto':
            xlabel = phenotype_name.replace(':', ': ').replace('_', ' ')
        
        volcano_plot(ax, df, 
                      up_hit=hit_dict[score_tag]['up_hit'], 
                      down_hit=hit_dict[score_tag]['down_hit'],
                      score_col=score_col, pvalue_col=pvalue_col,
                      xlabel=xlabel, ylabel=ylabel,
                      dot_size=dot_size, xlims=xlims, ylims=ylims,
                      ctrl_label=ctrl_label,
                      **args)
        
        if resistance_hits != None:
            if type(resistance_hits) != list: resistance_hits = [resistance_hits]
            for hit in resistance_hits:
                label_resistance_hit(
                    ax=ax, df_in=df, label=hit,
                    x_col=score_col,
                    y_col='-log10(pvalue)',
                    size=dot_size * 2,
                    size_txt=size_txt,
                    t_x=t_x, t_y=t_y
                )
        
        if sensitivity_hits != None:
            if type(sensitivity_hits) != list: sensitivity_hits = [sensitivity_hits]
            for hit in sensitivity_hits:
                label_sensitivity_hit(
                    ax=ax, df_in=df, label=hit,
                    x_col=score_col,
                    y_col='-log10(pvalue)',
                    size=dot_size * 2,
                    size_txt=size_txt,
                    t_x=t_x, t_y=t_y
            )


class GImaps(object):
    pass
