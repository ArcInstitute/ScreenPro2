"""
delta module
"""
import numpy as np
import anndata as ad
import pandas as pd
from .phenostat import matrixStat


### Key functions for calculating delta phenotype score


def calculatePhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, level):
    """Calculate phenotype score normalized by negative control and growth rate.
    
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
        level (str): level to use for calculating averaged score
    
    Returns:
        np.array: array of scores
    """
    # calculate control median and std 
    ctrl_median = np.median(
        calcLog2e(x=x_ctrl, y=y_ctrl), 
        axis=0 # for each individual sample (i.e. replicate)
    )

    # calculate log2e (i.e. log2 fold change enrichment y / x)
    log2e = calcLog2e(x=x, y=y)

    # calculate delta score normalized by control median and growth rate
    delta = (log2e - ctrl_median) / growth_rate

    # average delta score by given level (replicate or target)
    if level == 'replicate' or level == 'col':
        delta = np.mean(delta, axis=1)
    elif level == 'target' or level == 'row':
        delta = np.mean(delta, axis=0)
    elif level == 'all':
        delta = np.mean(delta)
    else:
        raise ValueError(f'Invalid level: {level}')
    
    return delta


def matrixTest(x, y, x_ctrl, y_ctrl, level, test = 'ttest', growth_rate = 1):
    """Calculate phenotype score and p-values comparing `y` vs `x` matrices.

    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        level (str): level to use for calculating score and p-value
        test (str): test to use for calculating p-value
        growth_rate (int): growth rate
    
    Returns:
        np.array: array of scores
        np.array: array of p-values
    """
    # calculate growth score
    scores = calculatePhenotypeScore(
        x = x, y = y, x_ctrl = x_ctrl, y_ctrl = y_ctrl,
        growth_rate = growth_rate,
        level = level
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, level = level)

    return scores, p_values


### Utility functions

def calcLog2e(x, y):
    return np.log2(y) - np.log2(x)


def generatePseudoGeneAnnData(adata, num_pseudogenes='auto', pseudogene_size='auto', ctrl_label='negative_control'):
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


def averageBestN(target_group, df_cond_ref, df_cond_test, keep_top_n):
    if keep_top_n and keep_top_n>0:
        df = pd.concat({
            'ref':df_cond_ref.loc[target_group.index,:],
            'test': df_cond_test.loc[target_group.index,:],
        }, axis=1)
        
        # Sort and find top n guide per target, see #18
        df = df.sort_values(df.columns.to_list(), ascending=False).iloc[:keep_top_n, :]

        df_cond_ref  = df['ref']
        df_cond_test = df['test']
    
    else:
        pass
    
    return df_cond_ref, df_cond_test
