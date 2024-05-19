"""
phenoscore module
"""

import numpy as np
import anndata as ad
import pandas as pd
from .phenostat import matrixStat, getFDR


def calculateDelta(x, y, transformation, level):
    """Calculate log ratio of y / x.
    `level` == 'all' (i.e. averaged across all values, sgRNA library elements and replicates)
    `level` == 'col' (i.e. averaged across columns, replicates)

    Args:
        x (np.array): array of values
        y (np.array): array of values
        transformation (str): transformation to use for calculating score
        level (str): level to use for calculating score
    
    Returns:
        np.array: array of log ratio values
    """
    # check if transformation is implemented
    if transformation not in ['log2', 'log2(x+1)', 'log10', 'log1p']:
        raise ValueError(f'transformation "{transformation}" not recognized')
    
    if level == 'all':
        # average across all values
        if transformation == 'log2':
            return np.log2(y) - np.log2(x)
        elif transformation == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1))
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x))
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x))
    elif level == 'row':
        # average across rows
        if transformation == 'log2':
            return np.log2(y) - np.log2(x)
        elif transformation == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=0)
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=0)
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=0)
    elif level == 'col':
        # average across columns
        if transformation == 'log2':
            return np.mean(np.log2(y) - np.log2(x), axis=1)
        elif transformation == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=1)
        elif transformation == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=1)
        elif transformation == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=1)


def calculatePhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, transformation, level):
    """Calculate phenotype score normalized by negative control and growth rate.
    
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
        transformation (str): transformation to use for calculating score
        level (str): level to use for calculating score
    
    Returns:
        np.array: array of scores
    """
    # calculate control median and std
    ctrl_median = np.median(calculateDelta(x=x_ctrl, y=y_ctrl, transformation=transformation, level=level))

    # calculate delta
    delta = calculateDelta(x=x, y=y, transformation=transformation, level=level)

    # calculate score
    return (delta - ctrl_median) / growth_rate


def matrixTest(x, y, x_ctrl, y_ctrl, transformation, level, test = 'ttest', growth_rate = 1):
    """Calculate phenotype score and p-values comparing `y` vs `x` matrices.

    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        transformation (str): transformation to use for calculating score
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
        growth_rate = growth_rate, transformation = transformation,
        level = level
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, level = level)

    return scores, p_values


def generatePseudoGeneAnnData(adata, num_pseudogenes='auto', pseudogene_size='auto', ctrl_label='negCtrl'):
    """Generate pseudogenes from negative control elements in the library.

    Args:
        adata (AnnData): AnnData object
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        ctrl_label (str): control label, default is 'negCtrl'
    
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


def calculateZScorePhenotypeScore(score_df,ctrl_label='negCtrl'):
    """Calculate z-score normalized phenotype score.
    
    Args:
        score_df (pd.DataFrame): dataframe of scores that includes `score` and `targetType` columns
        ctrl_label (str): control label, default is 'negCtrl'
    
    Returns:
        pd.Series: z-score normalized phenotype score
    """
    # calculate control median and std
    ctrl_std = score_df[score_df.targetType.eq(ctrl_label)].score.std()
    # z-score normalization
    out = score_df.score / ctrl_std

    return out


def runPhenoScore(adata, cond1, cond2, transformation, score_level, test,
                  growth_rate=1, n_reps=2, keep_top_n = None,num_pseudogenes='auto', pseudogene_size='auto',
                  count_layer=None, get_z_score=False, ctrl_label='negCtrl'):
    """Calculate phenotype score and p-values when comparing `cond2` vs `cond1`.

    Args:
        adata (AnnData): AnnData object
        cond1 (str): condition 1
        cond2 (str): condition 2
        transformation (str): transformation to use for calculating score
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        score_level (str): score level
        growth_rate (int): growth rate
        n_reps (int): number of replicates
        keep_top_n (int): number of top guides to keep per target
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        count_layer (str): count layer to use for calculating score, default is None (use default count layer in adata.X)
        get_z_score (bool): boolean to calculate z-score normalized phenotype score and add as a new column (default is False)
        ctrl_label (str): control label
    
    Returns:
        str: result name
        pd.DataFrame: result dataframe
    """
    # format result name
    result_name = f'{cond2}_vs_{cond1}'
    print(f'\t{cond2} vs {cond1}')

    # check if count_layer exists
    if count_layer is None:
        pass
    elif count_layer not in adata.layers.keys():
        raise ValueError(f"Layer '{count_layer}' not found in adata.layers.keys().")
    elif count_layer in adata.layers.keys():
        adata.X = adata.layers[count_layer].copy()
    
    # evaluate library table to get targets and riase error if not present
    required_columns = ['target', 'sequence']
    missing_columns = list(set(required_columns) - set(adata.var.columns))
    if len(missing_columns) > 0:
        raise ValueError(f"Missing required columns in library table: {missing_columns}")

    # calc phenotype score and p-value
    if score_level in ['compare_reps']:
        # prep counts for phenoScore calculation
        df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index[:n_reps],].to_df(count_layer).T
        df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index[:n_reps],].to_df(count_layer).T

        # convert to numpy arrays
        x = df_cond1.to_numpy()
        y = df_cond2.to_numpy()

        # get control values
        x_ctrl = df_cond1[adata.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond2[adata.var.targetType.eq(ctrl_label)].to_numpy()

        # calculate growth score and p_value
        scores, p_values = matrixTest(
            x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
            transformation=transformation, level='col', test=test, 
            growth_rate=growth_rate
        )
        # get adjusted p-values
        adj_p_values = getFDR(p_values)
                
        # get targets
        targets = adata.var['target'].to_list()

        # combine results into a dataframe
        result = pd.concat([
            pd.Series(targets, index=adata.var.index, name='target'),
            pd.Series(scores, index=adata.var.index, name='score'),
        ], axis=1)
        
        if get_z_score:
            # z-score normalization
            result['z_score'] = calculateZScorePhenotypeScore(result, ctrl_label=ctrl_label)
        # add p-values
        result[f'{test} pvalue'] = p_values
        result['BH adj_pvalue'] = adj_p_values
    
    elif score_level in ['compare_guides']:
        if n_reps == 2:
            pass
        else:
            raise ValueError('Currently, only n_reps=2 is supported for score_level="compare_guides"')
        # keep original adata for later use
        adata0 = adata.copy()

        # prep counts for phenoScore calculation
        df_cond1 = adata0[adata0.obs.query(f'condition=="{cond1}"').index].to_df().T
        df_cond2 = adata0[adata0.obs.query(f'condition=="{cond2}"').index].to_df().T
        # get control values
        x_ctrl = df_cond1[adata0.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond2[adata0.var.targetType.eq(ctrl_label)].to_numpy()
        del df_cond1, df_cond2
        
        adata_pseudo = generatePseudoGeneAnnData(adata0, num_pseudogenes=num_pseudogenes, pseudogene_size=pseudogene_size, ctrl_label=ctrl_label)
        adata = ad.concat([adata0[:,~adata0.var.targetType.eq(ctrl_label)], adata_pseudo], axis=1)
        adata.obs = adata0.obs.copy()

        # prep counts for phenoScore calculation
        df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index].to_df().T
        df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index].to_df().T
        
        targets = []
        scores = []
        p_values = []
        
        # group by target genes or pseudogenes to aggregate counts for score calculation
        for target_name, target_group in adata.var.groupby('target'):
            # convert to numpy arrays
            x = df_cond1.loc[target_group.index,:]
            y = df_cond2.loc[target_group.index,:]
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
        
        # combine results into a dataframe
        result = pd.concat([
            pd.Series([np.mean(s) for s in scores], index=targets, name='score'),
            pd.Series([np.mean(p) for p in p_values], index=targets, name=f'{test} pvalue'),
        ], axis=1)

        # # combine results into a dataframe
        # result = pd.concat({
        #     'replicate_1':pd.concat([
        #         pd.Series([s1 for s1,_ in scores], index=targets, name='score'),
        #         pd.Series([p1 for p1,_ in p_values], index=targets, name=f'{test} pvalue'),
                
        #     ],axis=1),
        #     'replicate_2':pd.concat([
        #         pd.Series([s2 for _,s2 in scores], index=targets, name='score'),
        #         pd.Series([p2 for _,p2 in p_values], index=targets, name=f'{test} pvalue'),
        #     ],axis=1),
        #     'replicate_ave':pd.concat([
        #         pd.Series([np.mean([s1,s2]) for s1,s2 in scores], index=targets, name='score'),
        #         pd.Series([np.mean([p1,p2]) for p1,p2 in p_values], index=targets, name=f'{test} pvalue'),
        #     ],axis=1)
        # },axis=1)
    
    else:
        raise ValueError(f'score_level "{score_level}" not recognized. Currently, "compare_reps" and "compare_guides" are supported.')
    
    
    return result_name, result


def runPhenoScoreForReplicate(adata, x_label, y_label, score, growth_factor_table, transformation, get_z_score=False, ctrl_label='negCtrl'):
    """Calculate phenotype score for each pair of replicates.

    Args:
        adata (AnnData): AnnData object
        x_label: name of the first condition in column `condition` of `screen.adata.obs`
        y_label: name of the second condition in column `condition` of `screen.adata.obs`
        score: score to use for calculating phenotype score, i.e. 'gamma', 'tau', or 'rho'
        growth_factor_table: dataframe of growth factors, i.e. output from `getGrowthFactors` function
        transformation (str): transformation to use for calculating score
        get_z_score: boolean to calculate z-score normalized phenotype score instead of regular score (default is False)
        ctrl_label: string to identify labels of negative control elements in sgRNA library (default is 'negCtrl')

    Returns:
        pd.DataFrame: dataframe of phenotype scores
    """
    adat = adata.copy()

    adat_ctrl = adat[:, adat.var.targetType.eq(ctrl_label)].copy()

    results = {}

    for replicate in adat.obs.replicate.unique():
        res = calculatePhenotypeScore(
            x=adat[adat.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y=adat[adat.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,

            x_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,

            growth_rate=growth_factor_table.query(f'score=="{score}" & replicate=={replicate}')['growth_factor'].values[0],
            transformation=transformation,
            level='row'  # there is only one column so `row` option here is equivalent to the value before averaging.
        )

        if get_z_score:
            res = calculateZScorePhenotypeScore(
                pd.DataFrame({'score':res,'targetType':adat.var.targetType},index=adat.var.index),
                ctrl_label=ctrl_label
            )

        results.update({f'replicate_{replicate}': res})

    out = pd.DataFrame(
        results,
        index=adat.var.index
    )

    return out
