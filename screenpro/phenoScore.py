"""phenoScore module
"""

import numpy as np
import pandas as pd
import anndata as ad

from pydeseq2 import preprocessing
from scipy.stats import ttest_rel
# from scipy.stats import mannwhitneyu
# from statsmodels.stats.multitest import multipletests


def seqDepthNormalization(adata):
    norm_counts, size_factors = preprocessing.deseq2_norm(adata.X)

    adata.obs['size_factors'] = size_factors
    adata.layers['seq_depth_norm'] = norm_counts

    
def getDelta(x, y, math='log2'):
    """log ratio of y / x, averaged across replicates 
    """
    if math == 'log2':
        return np.mean(np.log2(y) - np.log2(x), axis=1)
    elif math == 'log10':
        return np.mean(np.log10(y) - np.log10(x), axis=1)
    elif math == 'log1p':
        return np.mean(np.log1p(y) - np.log1p(x), axis=1)


def getScore(x, y, x_ctrl, y_ctrl, growth_rate):
    """
    """
    ctrl_std = np.std(getDelta(x_ctrl, y_ctrl))
    ctrl_median = np.median(getDelta(x_ctrl, y_ctrl))

    return ((getDelta(x,y) - ctrl_median) / growth_rate) / ctrl_std


def runPhenoScore(adata, cond1, cond2, growth_rate=1, n_reps=2, test='ttest', layer='seq_depth_norm'):
    # prep fqcounter
    df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index[:n_reps], ].to_df(layer).T
    df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index[:n_reps], ].to_df(layer).T

    x = df_cond1.to_numpy()
    y = df_cond2.to_numpy()

    x_ctrl = df_cond1[adata.var.targetType.eq('negCtrl')].to_numpy()
    y_ctrl = df_cond2[adata.var.targetType.eq('negCtrl')].to_numpy()
    
    # calculate growth score
    phenotype_score = getScore(x, y, x_ctrl, y_ctrl, growth_rate)
    
    adata.var[f'condition_{cond2}_vs_{cond1}_delta'] = phenotype_score
    
    # calculate p-values
    if test == 'MW':
        # run Mann-Whitney U rank test on replicates
        pass
    if test == 'ttest':
        # run ttest on replicates
        pvalues = ttest_rel(y, x, axis=1)[1]
        adata.var[f'condition_{cond2}_vs_{cond1}_pvalue'] = pvalues

    ## calculate FDR
    # Calculate the adjusted p-values using the Benjamini-Hochberg method
    # _, adj_pvalues, _, _ = multipletests(adata.var[f'condition_{cond1}_vs_{cond2}_pvalue'], alpha=0.05, method='fdr_bh')
    # adata.var[f'condition_{cond1}_vs_{cond2}_adj_pvalue'] = adj_pvalues
    
    
def ann_score_df(df_in, up_hit='resistance_hit', down_hit='sensitivity_hit', ctrl_label='non-targeting', threshold=10):
    df = df_in.copy()

    df.columns = ['target', 'score', 'pvalue']
    df['score'] = df['score'].astype(float)
    df['pvalue'] = df['pvalue'].astype(float)

    pseudo_sd = df[df['target'].str.contains(ctrl_label)]['score'].tolist()
    pseudo_sd = np.std(pseudo_sd)
    # print (pseudo_sd)
    
    df['label'] = '.'

    df.loc[
        (df['score'] > 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) >= threshold), 'label'
    ] = up_hit

    df.loc[
        (df['score'] < 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) <= -threshold), 'label'
    ] = down_hit

    df.loc[df['target'].str.contains(ctrl_label), 'label'] = ctrl_label

    df.loc[df['label'] == '.', 'label'] = 'target_non_hit'

    # reorder factors
    df['label'] = pd.Categorical(
        df['label'],
        categories=[down_hit, up_hit, ctrl_label, 'target_non_hit']
    )

    return df
