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

from .delta import calculatePhenotypeScore, matrixTest
from .deseq import runDESeq, extractDESeqResults
from .annotate import annotateScoreTable
from .phenostat import multipleTestsCorrection


def _generatePseudoGeneAnnData(adata, num_pseudogenes='auto', pseudogene_size='auto', ctrl_label='negative_control'):
    """Generate pseudogenes from negative control elements in the library.

    Args:
        adata (AnnData): AnnData object
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        ctrl_label (str): control label, default is 'negative_control'
    
    Returns:
        AnnData: AnnData object with pseudogenes
    """
    if pseudogene_size == 'auto':
        # sgRNA elements / target in the library
        pseudogene_size = int(adata.var[~adata.var.targetType.eq(ctrl_label)].groupby('target').size().mean())

    if num_pseudogenes == 'auto':
        # approx number of target in the library
        num_pseudogenes = len(adata.var.loc[~adata.var.targetType.eq(ctrl_label),'target'].unique()) * pseudogene_size
    
    adata_ctrl = adata[:,adata.var.targetType.eq(ctrl_label)].copy()
    ctrl_elements = adata_ctrl.var.index.to_list()
    
    adata_pseudo_list = []
    pseudo_source_sgrna = []

    for pseudo_num in range(0, num_pseudogenes, pseudogene_size):
        pseudo_elements = np.random.choice(ctrl_elements, pseudogene_size, replace=False)
        pseudo_labels = [f'pseudo_{pseudo_num}_{i}' for i in range(1,pseudogene_size+1)]
        
        adata_pseudo = ad.AnnData(
            X = adata_ctrl.X[:,adata_ctrl.var.index.isin(pseudo_elements)],
            obs = adata_ctrl.obs   
        )
        adata_pseudo.var_names = pseudo_labels
        adata_pseudo_list.append(adata_pseudo)
        
        for element in pseudo_elements: 
            pseudo_source_sgrna.append(element)
        
    out = ad.concat(adata_pseudo_list, axis=1)
    out.var['target'] = out.var.index.str.split('_').str[:-1].str.join('_')
    out.var['targetType'] = ctrl_label
    out.var['source'] = pseudo_source_sgrna
    out.obs = adata_ctrl.obs.copy()
    
    return out


def runPhenoScore(adata, cond_ref, cond_test, score_level, test, transformation='log2',
                  growth_rate=1, n_reps='auto', keep_top_n = None,num_pseudogenes='auto', pseudogene_size='auto',
                  count_layer=None, ctrl_label='negative_control'):
    """Calculate phenotype score and p-values when comparing `cond_test` vs `cond_ref`.

    Args:
        adata (AnnData): AnnData object
        cond_ref (str): condition reference
        cond_test (str): condition test
        score_level (str): score level
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        transformation (str): transformation to use for calculating score
        growth_rate (int): growth rate
        n_reps (int): number of replicates
        keep_top_n (int): number of top guides to keep per target
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        count_layer (str): count layer to use for calculating score, default is None (use default count layer in adata.X)
        ctrl_label (str): control label, default is 'negative_control'
    
    Returns:
        str: result name
        pd.DataFrame: result dataframe
    """
    # format result name
    result_name = f'{cond_test}_vs_{cond_ref}'
    print(f'\t{cond_test} vs {cond_ref}')

    # set n_reps if not provided
    if n_reps == 'auto':
        n_reps = min(
            adata.obs.query(f'condition=="{cond_ref}"').shape[0], 
            adata.obs.query(f'condition=="{cond_test}"').shape[0]
        )

    # check if count_layer exists
    if count_layer is None:
        pass
    elif count_layer not in adata.layers.keys():
        raise ValueError(f"Layer '{count_layer}' not found in adata.layers.keys().")
    elif count_layer in adata.layers.keys():
        adata.X = adata.layers[count_layer].copy()
    
    # evaluate library table to get targets and riase error if not present
    required_columns = ['target'] #, 'sequence']
    missing_columns = list(set(required_columns) - set(adata.var.columns))
    if len(missing_columns) > 0:
        raise ValueError(f"Missing required columns in library table: {missing_columns}")

    # calc phenotype score and p-value
    if score_level in ['compare_reps']:
        # prep counts for phenoScore calculation
        df_cond_ref = adata[adata.obs.query(f'condition=="{cond_ref}"').index[:n_reps],].to_df(count_layer).T
        df_cond_test = adata[adata.obs.query(f'condition=="{cond_test}"').index[:n_reps],].to_df(count_layer).T

        # convert to numpy arrays
        x = df_cond_ref.to_numpy()
        y = df_cond_test.to_numpy()

        # get control values
        x_ctrl = df_cond_ref[adata.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond_test[adata.var.targetType.eq(ctrl_label)].to_numpy()

        # calculate growth score and p_value
        scores, p_values = matrixTest(
            x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
            transformation=transformation, level='col', test=test, 
            growth_rate=growth_rate
        )
        # get adjusted p-values
        adj_p_values = multipleTestsCorrection(p_values)
                
        # get targets
        targets = adata.var['target'].to_list()

        # combine results into a dataframe
        result = pd.concat([
            pd.Series(targets, index=adata.var.index, name='target'),
            pd.Series(scores, index=adata.var.index, name='score'),
        ], axis=1)
        
        # add p-values
        result[f'{test} pvalue'] = p_values
        result['BH adj_pvalue'] = adj_p_values
    
    elif score_level in ['compare_guides']:
        # keep original adata for later use
        adata0 = adata.copy()

        # prep counts for phenoScore calculation
        df_cond_ref = adata0[adata0.obs.query(f'condition=="{cond_ref}"').index].to_df().T
        df_cond_test = adata0[adata0.obs.query(f'condition=="{cond_test}"').index].to_df().T
        # get control values
        x_ctrl = df_cond_ref[adata0.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond_test[adata0.var.targetType.eq(ctrl_label)].to_numpy()
        del df_cond_ref, df_cond_test
        
        adata_pseudo = _generatePseudoGeneAnnData(adata0, num_pseudogenes=num_pseudogenes, pseudogene_size=pseudogene_size, ctrl_label=ctrl_label)
        adata = ad.concat([adata0[:,~adata0.var.targetType.eq(ctrl_label)], adata_pseudo], axis=1)
        adata.obs = adata0.obs.copy()

        # prep counts for phenoScore calculation
        df_cond_ref = adata[adata.obs.query(f'condition=="{cond_ref}"').index].to_df().T
        df_cond_test = adata[adata.obs.query(f'condition=="{cond_test}"').index].to_df().T
        
        targets = []
        scores = []
        p_values = []
        
        # group by target genes or pseudogenes to aggregate counts for score calculation
        for target_name, target_group in adata.var.groupby('target'):
            # convert to numpy arrays
            x = df_cond_ref.loc[target_group.index,:]
            y = df_cond_test.loc[target_group.index,:]
            # Sort and find top n guide per target, see #18
            if keep_top_n:
                x = x.sort_values(x.columns.to_list(), ascending=False)
                y = y.sort_values(y.columns.to_list(), ascending=False)
                x = x.iloc[:keep_top_n, :]
                y = y.iloc[:keep_top_n, :]
            # convert to numpy arrays
            x = x.to_numpy()
            y = y.to_numpy()

            # calculate growth score and p_value
            target_scores, target_p_values = matrixTest(
                x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
                transformation=transformation, 
                level='all', # test across all guides and replicates per target
                test=test,
                growth_rate=growth_rate
            )
            
            scores.append(target_scores)
            p_values.append(target_p_values)
            targets.append(target_name)

        # get adjusted p-values
        adj_p_values = multipleTestsCorrection(np.array(p_values))
        
        # combine results into a dataframe
        result = pd.concat([
            pd.Series(targets, index=targets, name='target'),
            pd.Series(scores, index=targets, name='score'),
            pd.Series(p_values, index=targets, name=f'{test} pvalue'),
            pd.Series(adj_p_values, index=targets, name='BH adj_pvalue'),
        ], axis=1)

        # rename pseudo genes in target column to `ctrl_label`
        result['target'] = result['target'].apply(lambda x: ctrl_label if 'pseudo' in x else x)

    
    else:
        raise ValueError(f'score_level "{score_level}" not recognized. Currently, "compare_reps" and "compare_guides" are supported.')
    
    return result_name, result


def runPhenoScoreForReplicate(adata, x_label, y_label, growth_rate_reps=None, transformation='log2', ctrl_label='negative_control'):
    """Calculate phenotype score for each pair of replicates.

    Args:
        adata (AnnData): AnnData object
        x_label: name of the first condition in column `condition` of `screen.adata.obs`
        y_label: name of the second condition in column `condition` of `screen.adata.obs`
        growth_rate_reps (dict): dictionary of growth rates for each replicate
        transformation (str): transformation to use for calculating score
        ctrl_label: string to identify labels of negative control elements in sgRNA library (default is 'negative_control')

    Returns:
        pd.DataFrame: dataframe of phenotype scores
    """
    adat = adata.copy()

    adat_ctrl = adat[:, adat.var.targetType.eq(ctrl_label)].copy()

    results = {}

    if growth_rate_reps is None:
        growth_rate_reps = dict([(replicate, 1) for replicate in adat.obs.replicate.unique()])

    for replicate in adat.obs.replicate.unique():

        res = calculatePhenotypeScore(
            x=adat[adat.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y=adat[adat.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,
            x_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,
            growth_rate=growth_rate_reps[replicate],
            transformation=transformation,
            level='row'  # there is only one column so `row` option here is equivalent to the value before averaging.
        )

        results.update({f'replicate_{replicate}': res.reshape(-1)})

    out = pd.DataFrame(
        results,
        index=adat.var.index
    )

    return out
