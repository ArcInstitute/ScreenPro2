"""
phenoScore module
"""

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from .phenoStats import matrixStat, getFDR
from .utils import find_low_counts


def calculateDelta(x, y, math, level):
    """
    Calculate log ratio of y / x.
    `level` == 'all' (i.e. averaged across all values, oligo and replicates)
    `level` == 'col' (i.e. averaged across columns, replicates)
    Args:
        x (np.array): array of values
        y (np.array): array of values
        math (str): math to use for calculating score
        level (str): level to use for calculating score
    Returns:
        np.array: array of log ratio values
    """
    # check if math is implemented
    if math not in ['log2(x+1)', 'log10', 'log1p']:
        raise ValueError(f'math "{math}" not recognized')
    
    if level == 'all':
        # average across all values
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1))
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x))
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x))
    elif level == 'row':
        # average across rows
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=0)
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=0)
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=0)
    elif level == 'col':
        # average across columns
        if math == 'log2(x+1)':
            return np.mean(np.log2(y+1) - np.log2(x+1), axis=1)
        elif math == 'log10':
            return np.mean(np.log10(y) - np.log10(x), axis=1)
        elif math == 'log1p':
            return np.mean(np.log1p(y) - np.log1p(x), axis=1)


def calculatePhenotypeScore(x, y, x_ctrl, y_ctrl, growth_rate, math, level):
    """
    Calculate phenotype score normalized by negative control and growth rate.
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
        math (str): math to use for calculating score
        level (str): level to use for calculating score
    Returns:
        np.array: array of scores
    """
    # calculate control median and std
    ctrl_median = np.median(calculateDelta(x=x_ctrl, y=y_ctrl, math=math, level=level))

    # calculate delta
    delta = calculateDelta(x=x, y=y, math=math, level=level)

    # calculate score
    return (delta - ctrl_median) / growth_rate


def matrixTest(x, y, x_ctrl, y_ctrl, math, level, test = 'ttest', growth_rate = 1):
    """
    Calculate phenotype score and p-values comparing `y` vs `x` matrices.
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        math (str): math to use for calculating score
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
        growth_rate = growth_rate, math = math,
        level = level
    )

    # compute p-value
    p_values = matrixStat(x, y, test=test, level = level)

    return scores, p_values


def runPhenoScore(adata, cond1, cond2, math, score_level, test,
                  growth_rate=1, n_reps=2, keep_top_n = None,num_pseudogenes=None,
                  get_z_score=False,ctrl_label='negCtrl'):
    """
    Calculate phenotype score and p-values when comparing `cond2` vs `cond1`.
    Args:
        adata (AnnData): AnnData object
        cond1 (str): condition 1
        cond2 (str): condition 2
        math (str): math to use for calculating score
        test (str): test to use for calculating p-value ('MW': Mann-Whitney U rank; 'ttest' : t-test)
        score_level (str): score level
        growth_rate (int): growth rate
        n_reps (int): number of replicates
        get_z_score (bool): boolean to calculate z-score normalized phenotype score and add as a new column (default is False)
        ctrl_label (str): control label
    Returns:
        str: result name
        pd.DataFrame: result dataframe
    """
    # format result name
    result_name = f'{cond2}_vs_{cond1}'
    print(f'\t{cond2} vs {cond1}')

    count_layer = 'seq_depth_norm'
    # check if count_layer exists
    if 'seq_depth_norm' not in adata.layers.keys():
        seqDepthNormalization(adata)
    
    # evaluate library table to get targets and riase error if not present
    required_columns = ['target', 'sequence']
    missing_columns = list(set(required_columns) - set(adata.var.columns))
    if len(missing_columns) > 0:
        raise ValueError(f"Missing required columns in library table: {missing_columns}")

    # calc phenotype score and p-value
    if score_level == 'compare_reps':
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
            math=math, level='col', test=test, 
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
            raise ValueError(f'n_reps "{n_reps}" not recognized')

        find_low_counts(adata)
        adata = adata[:,~adata.var.low_count].copy()

        adata.var.drop(columns='low_count',inplace=True)

        # generate pseudo gene labels
        generatePseudoGeneLabels(adata, num_pseudogenes=num_pseudogenes, ctrl_label=ctrl_label)
        # drop categories from target column!
        adata.var['target'] = adata.var['target'].to_list()
        # replace values in target columns for negative control with pseudogene label
        adata.var.loc[
            adata.var.targetType.eq(ctrl_label), 'target'
        ] = adata.var.loc[adata.var.targetType.eq(ctrl_label), 'pseudoLabel']
        # drop pseudoLabel column 
        adata.var.drop(columns='pseudoLabel', inplace=True)

        # prep counts for phenoScore calculation
        df_cond1 = adata[adata.obs.query(f'condition=="{cond1}"').index].to_df().T
        df_cond2 = adata[adata.obs.query(f'condition=="{cond2}"').index].to_df().T        # get control values

        # get control values
        x_ctrl = df_cond1[adata.var.targetType.eq(ctrl_label)].to_numpy()
        y_ctrl = df_cond2[adata.var.targetType.eq(ctrl_label)].to_numpy()

        targets = []
        scores = []
        p_values = []

        # group by target genes or pseudogenes to aggregate counts for score calculation
        for target_name, target_group in adata.var.groupby('target'):
            # convert to numpy arrays
            x = df_cond1.loc[target_group.index,:].to_numpy()
            y = df_cond2.loc[target_group.index,:].to_numpy()
            # Sort and find top n guide per target, see #18
            if keep_top_n:
                x.sort()
                y.sort()
                x = x[:keep_top_n]
                y = y[:keep_top_n]
            
            # calculate growth score and p_value
            target_scores, target_p_values = matrixTest(
                x=x, y=y, x_ctrl=x_ctrl, y_ctrl=y_ctrl,
                math=math, level='row', test=test,
                growth_rate=growth_rate
            )

            scores.append(target_scores)
            p_values.append(target_p_values)
            targets.append(target_name)
        
        # combine results into a dataframe
        result = pd.concat({
            'replicate_1':pd.concat([
                pd.Series([s1 for s1,_ in scores], index=targets, name='score'),
                pd.Series([p1 for p1,_ in p_values], index=targets, name=f'{test} pvalue'),
                
            ],axis=1),
            'replicate_2':pd.concat([
                pd.Series([s2 for _,s2 in scores], index=targets, name='score'),
                pd.Series([p2 for _,p2 in p_values], index=targets, name=f'{test} pvalue'),
            ],axis=1),
            'replicate_ave':pd.concat([
                pd.Series([np.mean([s1,s2]) for s1,s2 in scores], index=targets, name='score'),
                pd.Series([np.mean([p1,p2]) for p1,p2 in p_values], index=targets, name=f'{test} pvalue'),
            ],axis=1)
        },axis=1)
    
    else:
        raise ValueError(f'score_level "{score_level}" not recognized')
    
    
    return result_name, result


def seqDepthNormalization(adata):
    """
    Normalize counts by sequencing depth.
    Args:
        adata (AnnData): AnnData object
    """
    # normalize counts by sequencing depth
    norm_counts, size_factors = preprocessing.deseq2_norm(adata.X)
    # update adata object
    adata.obs['size_factors'] = size_factors
    adata.layers['seq_depth_norm'] = norm_counts


def addPseudoCount():
    pass
    # # pseudocount
    # if pseudocountBehavior == 'default' or pseudocountBehavior == 'zeros only':
    #     def defaultBehavior(row):
    #         return row if min(
    #             row) != 0 else row + pseudocountValue
    #
    #     combinedCountsPseudo = combinedCounts.apply(defaultBehavior, axis=1)
    # elif pseudocountBehavior == 'all values':
    #     combinedCountsPseudo = combinedCounts.apply(
    #         lambda row: row + pseudocountValue, axis=1)
    # elif pseudocountBehavior == 'filter out':
    #     combinedCountsPseudo = combinedCounts.copy()
    #     zeroRows = combinedCounts.apply(lambda row: min(row) <= 0, axis=1)
    #     combinedCountsPseudo.loc[zeroRows, :] = np.nan
    # else:
    #     raise ValueError(
    #         'Pseudocount behavior not recognized or not implemented')


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


def runPhenoScoreForReplicate(screen, x_label, y_label, score, growth_factor_table, get_z_score=False, ctrl_label='negCtrl'):
    """
    Calculate phenotype score for each pair of replicates.
    Args:
        screen: ScreenPro object
        x_label: name of the first condition in column `condition` of `screen.adata.obs`
        y_label: name of the second condition in column `condition` of `screen.adata.obs`
        score: score to use for calculating phenotype score, i.e. 'gamma', 'tau', or 'rho'
        growth_factor_table: dataframe of growth factors, i.e. output from `getGrowthFactors` function
        get_z_score: boolean to calculate z-score normalized phenotype score instead of regular score (default is False)
        ctrl_label: string to identify labels of negative control oligos

    Returns:
        pd.DataFrame: dataframe of phenotype scores
    """
    adat = screen.adata.copy()

    adat_ctrl = adat[:, adat.var.targetType.eq(ctrl_label)].copy()

    results = {}

    for replicate in adat.obs.replicate.unique():
        res = calculatePhenotypeScore(
            x=adat[adat.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y=adat[adat.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,

            x_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{x_label}" & replicate == {str(replicate)}').index].X,
            y_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{y_label}" & replicate == {str(replicate)}').index].X,

            growth_rate=growth_factor_table.query(f'score=="{score}" & replicate=={replicate}')['growth_factor'].values[0],
            math=screen.math,
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


def generatePseudoGeneLabels(adata, num_pseudogenes=None, ctrl_label='negCtrl'):
    """
    Generate new labels per `num_pseudogenes` randomly selected non targeting oligo in `adata.var`.
    Args:
        adata (AnnData): AnnData object
        num_pseudogenes (int): number of pseudogenes
        ctrl_label (str): control label
    """
    # check if num_pseudogenes is defined. 
    ## If not, set to half of the number of non-targeting oligos
    if num_pseudogenes is None:
        num_pseudogenes = len(adata.var[adata.var.targetType.eq(ctrl_label)]) // 2
    # Get non-targeting oligos
    ctrl_oligos = adata.var[adata.var.targetType.eq(ctrl_label)].index
    adata.var['pseudoLabel'] = ''
    # Check if there are more than 1 non-targeting oligos to label as pseudogenes
    if len(ctrl_oligos) / 2 <= num_pseudogenes:
        # raise error if `num_pseudogenes` is greater than (total number of non-targeting oligos) / 2
        raise TypeError(
            "Define `num_pseudogenes` to be less than (total number of non-targeting oligos) / 2"
        )
    else:
        # randomly select `num` non-targeting oligos
        while len(ctrl_oligos) > num_pseudogenes:
            # randomly select `num` non-targeting oligos
            pseudo_oligos = np.random.choice(ctrl_oligos, num_pseudogenes, replace=False)
            # generate new labels
            pseudo_labels = [f'pseudo_{i}' for i in range(num_pseudogenes)]
            # update adata.var
            adata.var.loc[pseudo_oligos, 'pseudoLabel'] = pseudo_labels
            # remove selected oligos from ctrl_oligos
            ctrl_oligos = ctrl_oligos.drop(pseudo_oligos)

        # label remaining non-targeting oligos as pseudogenes
        adata.var.loc[adata.var.targetType.eq('gene'), 'pseudoLabel'] = 'gene'
        adata.var.loc[adata.var.pseudoLabel.eq(''), 'pseudoLabel'] = np.nan

