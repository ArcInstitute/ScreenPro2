"""
delta module
"""
import numpy as np
import anndata as ad
import pandas as pd

from .phenostat import (
    matrixStat, multipleTestsCorrection
)


### Key functions for calculating delta phenotype score

def compareByReplicates(adata, df_cond_ref, df_cond_test, var_names='target', test='ttest', ctrl_label='negative_control', growth_rate=1, filter_type='mean', filter_threshold=40):
    """Calculate phenotype score and p-values comparing `cond_test` vs `cond_ref`.

    In this function, the phenotype calculation is done by comparing multiple replicates of `cond_test` vs `cond_ref`.

    Args:
        adata (AnnData): AnnData object
        df_cond_ref (pd.DataFrame): dataframe of condition reference
        df_cond_test (pd.DataFrame): dataframe of condition test
        var_names (str): variable names to use as index in the result dataframe
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        ctrl_label (str): control label, default is 'negative_control'
        growth_rate (int): growth rate
        filter_type (str): filter type to apply to low counts ('mean', 'both', 'either')
        filter_threshold (int): filter threshold for low counts (default is 40)
    
    Returns:
        pd.DataFrame: result dataframe
    """
    adat = adata.copy()

    # apply NA to low counts
    df_cond_ref, df_cond_test = applyNAtoLowCounts(
        df_cond_ref=df_cond_ref, df_cond_test=df_cond_test,
        filter_type=filter_type, 
        filter_threshold=filter_threshold
    )

    # convert to numpy arrays
    x = df_cond_ref.to_numpy()
    y = df_cond_test.to_numpy()

    # get control values
    x_ctrl = df_cond_ref[adat.var.targetType.eq(ctrl_label)].dropna().to_numpy()
    y_ctrl = df_cond_test[adat.var.targetType.eq(ctrl_label)].dropna().to_numpy()

    # calculate phenotype scores
    scores = calculateDelta(
        x = x, y = y, 
        x_ctrl = x_ctrl, y_ctrl = y_ctrl, 
        growth_rate = growth_rate,
    )

    # average scores across replicates
    scores = [np.mean(s) for s in scores]
    
    # compute p-value
    p_values = matrixStat(x, y, test=test, level = 'col')

    # get adjusted p-values
    adj_p_values = multipleTestsCorrection(p_values)

    # get target information            
    targets_df = adat.var[var_names].copy()

    # combine results into a dataframe
    result = pd.concat([
        pd.Series(scores, index=adat.var.index, name='score'),
        pd.Series(p_values, index=adat.var.index, name=f'{test} pvalue'),
        pd.Series(adj_p_values, index=adat.var.index, name='BH adj_pvalue'),
    ], axis=1)
    
    # add targets information 
    result = pd.concat([targets_df, result], axis=1)

    return result


def compareByTargetGroup(adata, df_cond_ref, df_cond_test, keep_top_n, var_names='target', test='ttest', ctrl_label='negative_control', growth_rate=1, filter_type='mean', filter_threshold=40):
    """Calculate phenotype score and p-values comparing `cond_test` vs `cond_ref`.

    In this function, the phenotype calculation is done by comparing groups of 
    guide elements (e.g. sgRNAs) that target the same gene or groups of pseudogene (i.e.
    subsampled groups of non-targeting control elements) between `cond_test` vs `cond_ref`.

    Args:
        adata (AnnData): AnnData object
        df_cond_ref (pd.DataFrame): dataframe of condition reference
        df_cond_test (pd.DataFrame): dataframe of condition test
        keep_top_n (int): number of top guide elements to keep
        var_names (str): variable names to use as index in the result dataframe
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        ctrl_label (str): control label, default is 'negative_control'
        growth_rate (int): growth rate
        filter_type (str): filter type to apply to low counts ('mean', 'both', 'either')
        filter_threshold (int): filter threshold for low counts (default is 40)
        
    Returns:
        pd.DataFrame: result dataframe
    """

    adat = adata.copy()

    # apply NA to low counts
    df_cond_ref, df_cond_test = applyNAtoLowCounts(
        df_cond_ref=df_cond_ref, df_cond_test=df_cond_test,
        filter_type=filter_type, 
        filter_threshold=filter_threshold
    )

    # get control values
    x_ctrl = df_cond_ref[adat.var.targetType.eq(ctrl_label)].dropna().to_numpy()
    y_ctrl = df_cond_test[adat.var.targetType.eq(ctrl_label)].dropna().to_numpy()

    targets = []
    scores = []
    p_values = []
    target_sizes = []

    # group by target genes or pseudogenes to aggregate counts for score calculation
    for target_name, target_group in adat.var.groupby(var_names):

        # calculate phenotype scores and p-values for each target group
        target_score, target_p_value, target_size  = scoreTargetGroup(
            target_group=target_group, 
            df_cond_ref=df_cond_ref, 
            df_cond_test=df_cond_test,
            x_ctrl=x_ctrl, y_ctrl=y_ctrl, 
            test=test, growth_rate=growth_rate, 
            keep_top_n=keep_top_n
        )
        
        scores.append(target_score)
        p_values.append(target_p_value)
        targets.append(target_name)
        target_sizes.append(target_size)

    # average scores across replicates
    scores = [np.mean(s) for s in scores]

    # get adjusted p-values
    adj_p_values = multipleTestsCorrection(np.array(p_values))

    # get target information
    if type(var_names) == str:
        targets_df = pd.DataFrame(targets, columns=[var_names])
    elif type(var_names) == list:
        targets_df = pd.DataFrame(targets, columns=var_names)

    # combine results into a dataframe
    result = pd.concat([
        pd.Series(scores, name='score'),
        pd.Series(p_values, name=f'{test} pvalue'),
        pd.Series(adj_p_values, name='BH adj_pvalue'),
        pd.Series(target_sizes, name='number_of_guide_elements'),
    ], axis=1)

    # add targets information
    result = pd.concat([targets_df, result], axis=1)

    # set index to var_names 
    if type(var_names) == list and len(var_names) > 1:
        result.index = result[var_names].agg('-'.join, axis=1)
    else:
        result.index = result[var_names]

    return result


def getPhenotypeData(adata, score_tag, cond_ref, cond_test, growth_rate_reps=None, ctrl_label='negative_control'):
    """Calculate phenotype score for each pair of replicates

    Args:
        adata (AnnData): AnnData object
        score_tag (str): score tag. e.g. 'delta', 'gamma', 'tau', 'rho'.
        cond_ref (str): condition reference
        cond_test (str): condition test
        growth_rate_reps (dict): growth rate for each replicate. Key is replicate number, value is growth rate.
        ctrl_label (str): control label, default is 'negative_control'
    """
    score_name = f'{score_tag}:{cond_test}_vs_{cond_ref}'

    adat = adata.copy()

    adat_ctrl = adat[:, adat.var.targetType.eq(ctrl_label)].copy()
    
    results = {}
    
    if growth_rate_reps is None:
        growth_rate_reps = dict([(replicate, 1) for replicate in adat.obs.replicate.unique()])
    
    for replicate in adat.obs.replicate.unique():
        x=adat[adat.obs.query(f'condition == "{cond_ref}" & replicate == {str(replicate)}').index].X.T
        y=adat[adat.obs.query(f'condition == "{cond_test}" & replicate == {str(replicate)}').index].X.T
    
        x_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{cond_ref}" & replicate == {str(replicate)}').index].X.T
        y_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{cond_test}" & replicate == {str(replicate)}').index].X.T
    
        res = calculateDelta(
            x=x,y=y,x_ctrl=x_ctrl,y_ctrl=y_ctrl,
            growth_rate=growth_rate_reps[replicate],
        )
    
        results.update({f'{score_name}::replicate_{replicate}': res.reshape(-1)})
    
    out = ad.AnnData(
        pd.DataFrame(results, index=adat.var.index).T
    )

    return out


def calculateDelta(x, y, x_ctrl, y_ctrl, growth_rate):
    """Calculate phenotype score normalized by negative control and growth rate.
    
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
    
    Returns:
        np.array: array of scores
    """
    # calculate control median and std 
    ctrl_median = np.median(
        calculateLog2e(x=x_ctrl, y=y_ctrl), 
        axis=0 # for each individual sample (i.e. replicate)
    )

    # calculate log2e (i.e. log2 fold change enrichment y / x)
    log2e = calculateLog2e(x=x, y=y)

    # calculate delta score normalized by control median and growth rate
    delta = (log2e - ctrl_median) / growth_rate
    
    return delta


### Utility functions

def calculateLog2e(x, y):
    return np.log2(y) - np.log2(x)


def averageBestN(scores, numToAverage):
    # Sort and find top n guide per target, see #18
    return np.mean(sorted(scores, key=abs, reverse=True)[:numToAverage])


def getBestTargetByTSS(score_df,target_col,pvalue_col):
    """
    collapse the gene-transcript indices into a single score for a gene by best p-value
    """
    return score_df.dropna().groupby(target_col).apply(
        lambda x: x.loc[x[pvalue_col].idxmin()]
    )


def scoreTargetGroup(target_group, df_cond_ref, df_cond_test, x_ctrl, y_ctrl, test='ttest', growth_rate=1, keep_top_n=None):
    # select target group and convert to numpy arrays
    x = df_cond_ref.loc[target_group.index,:].dropna().to_numpy()
    y = df_cond_test.loc[target_group.index,:].dropna().to_numpy()

    # calculate phenotype scores
    target_scores = calculateDelta(
        x = x, y = y, 
        x_ctrl = x_ctrl, y_ctrl = y_ctrl, 
        growth_rate = growth_rate,
    )

    # get target size
    target_size = target_scores.shape[0] # number of guide elements in the target group
    
    if target_size == 0:
        target_score = np.full(target_scores.shape[1], np.nan)
    
    elif (keep_top_n is None or keep_top_n is False) or target_size <= keep_top_n:
        # average scores across guides
        target_score = np.mean(target_scores, axis=0)

    elif keep_top_n > 0 or target_size > keep_top_n:
        # get top n scores per target
        target_score = np.apply_along_axis(averageBestN, axis=0, arr=target_scores, numToAverage=keep_top_n)
        target_size = keep_top_n # update target size to keep_top_n
    
    else:
        raise ValueError(f'Invalid value for keep_top_n: {keep_top_n}')

    # compute p-value
    target_p_value = matrixStat(x, y, test=test, level='all')
    
    return target_score, target_p_value, target_size


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
    #TODO: check for `target` and `targetType` columns in adata.var
    #TODO: add input arg to replace "target" and name it `target_col` or something similar
    
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


def applyNAtoLowCounts(df_cond_ref, df_cond_test, filter_type, filter_threshold):
    # more flexible read filtering by adding np.nan scores and/or pvalues to low count rows
    # keep row if either both/all columns are above threshold, or if either/any column is
    # in other words, mask if any column is below threshold or only if all columns are below
    # source https://github.com/mhorlbeck/ScreenProcessing/blob/master/process_experiments.py#L464C1-L478C1

    df = pd.concat({'ref':df_cond_ref, 'test':df_cond_test},axis=1)
    
    if filter_type == 'mean':
        filter = df.apply(
            lambda row: np.mean(row) < filter_threshold, axis=1)
    elif filter_type == 'both' or filter_type == 'all':
        filter = df.apply(
            lambda row: min(row) < filter_threshold, axis=1)
    elif filter_type == 'either' or filter_type == 'any':
        filter = df.apply(
            lambda row: max(row) < filter_threshold, axis=1)
    else:
        raise ValueError('filter type not recognized or not implemented')

    df.loc[filter, :] = np.nan

    return df['ref'], df['test']
