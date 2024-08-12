## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

"""phenoscore module

This module contains functions for calculating relative phenotypes from CRISPR screens
datasets.
"""

import numpy as np
import anndata as ad
import pandas as pd

from .delta import (
    compareByReplicates, compareByTargetGroup,
    getPhenotypeData,
    calculateDelta,
    getBestTargetByTSS,
    generatePseudoGeneAnnData
)
from .deseq import runDESeq, extractDESeqResults
from ._annotate import annotateScoreTable
from .phenostat import matrixStat, multipleTestsCorrection


def runPhenoScore(adata, cond_ref, cond_test, score_level, var_names='target', test='ttest',
                  growth_rate=1, n_reps='auto', keep_top_n = None, collapse_var=False,
                  num_pseudogenes='auto', pseudogene_size='auto',
                  count_layer=None, count_filter_type='mean', count_filter_threshold=40,
                  ctrl_label='negative_control'
                  ):
    """Calculate phenotype score and p-values when comparing `cond_test` vs `cond_ref`.

    Args:
        adata (AnnData): AnnData object
        cond_ref (str): condition reference
        cond_test (str): condition test
        score_level (str): score level
        var_names (str): variable names to use as index in the result dataframe
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        growth_rate (int): growth rate
        n_reps (int): number of replicates
        keep_top_n (int): number of top guides to keep per target
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        count_layer (str): count layer to use for calculating score, default is None (use default count layer in adata.X)
        count_filter_type (str): filter type for counts, default is 'mean'
        count_filter_threshold (int): filter threshold for counts, default is 40
        ctrl_label (str): control label, default is 'negative_control'
    
    Returns:
        str: result name
        pd.DataFrame: result dataframe
    """
    adat = adata.copy()

    # format result name
    result_name = f'{cond_test}_vs_{cond_ref}'
    print(f'\t{cond_test} vs {cond_ref}')
    
    # set n_reps if not provided
    if n_reps == 'auto':
        n_reps = min(
            adat.obs.query(f'condition=="{cond_ref}"').shape[0], 
            adat.obs.query(f'condition=="{cond_test}"').shape[0]
        )

    # check if count_layer exists
    if count_layer is None:
        pass
    elif count_layer not in adat.layers.keys():
        raise ValueError(f"Layer '{count_layer}' not found in adata.layers.keys().")
    elif count_layer in adat.layers.keys():
        adat.X = adat.layers[count_layer].copy()
    
    # calc phenotype score and p-value
    if score_level in ['compare_reps']:

        # prep counts for phenoScore calculation
        df_cond_ref = adat[adat.obs.query(f'condition=="{cond_ref}"').index[:n_reps],].to_df(count_layer).T
        df_cond_test = adat[adat.obs.query(f'condition=="{cond_test}"').index[:n_reps],].to_df(count_layer).T

        result = compareByReplicates(
            adata=adat, 
            df_cond_ref=df_cond_ref, 
            df_cond_test=df_cond_test,
            var_names=var_names,
            test=test,
            ctrl_label=ctrl_label,
            growth_rate=growth_rate,
            filter_type=count_filter_type,
            filter_threshold=count_filter_threshold
        )
    
    elif score_level in ['compare_guides']:

        # prep counts for phenoScore calculation
        df_cond_ref = adat[adat.obs.query(f'condition=="{cond_ref}"').index].to_df().T
        df_cond_test = adat[adat.obs.query(f'condition=="{cond_test}"').index].to_df().T
        del df_cond_ref, df_cond_test
        
        adat_pseudo = generatePseudoGeneAnnData(adat, num_pseudogenes=num_pseudogenes, pseudogene_size=pseudogene_size, ctrl_label=ctrl_label)
        if 'transcript' in var_names: adat_pseudo.var['transcript'] = 'na'

        adat_test = ad.concat([adat[:,~adat.var.targetType.eq(ctrl_label)], adat_pseudo], axis=1)
        adat_test.obs = adat.obs.copy()

        # prep counts for phenoScore calculation
        df_cond_ref = adat_test[adat_test.obs.query(f'condition=="{cond_ref}"').index].to_df().T
        df_cond_test = adat_test[adat_test.obs.query(f'condition=="{cond_test}"').index].to_df().T
        
        result = compareByTargetGroup(
            adata=adat_test, 
            df_cond_ref=df_cond_ref, 
            df_cond_test=df_cond_test,
            keep_top_n=keep_top_n,
            var_names=var_names,
            test=test,
            ctrl_label=ctrl_label,
            growth_rate=growth_rate,
            filter_type=count_filter_type,
            filter_threshold=count_filter_threshold
        )

        # get best best transcript as lowest p-value for each target
        if collapse_var not in [False, None]:
            if collapse_var not in result.columns:
                raise ValueError(f'collapse_var "{collapse_var}" not found in result columns.')
            else:
                result = getBestTargetByTSS(
                    score_df=result, target_col=collapse_var, pvalue_col=f'{test} pvalue'
                )
                result.index.name = None
        
        # change target name to control label if it is a pseudo gene
        result['target'] = result['target'].apply(lambda x: ctrl_label if 'pseudo' in x else x).to_list()
        
    else:
        raise ValueError(f'score_level "{score_level}" not recognized. Currently, "compare_reps" and "compare_guides" are supported.')
    
    return result_name, result
