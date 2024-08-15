"""
annotate module
"""

import numpy as np
import pandas as pd

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


def getCombinedScore(df_in, score_col='score', pvalue_col='pvalue', target_col='target', ctrl_label='negative_control'):
    """
    Calculate the combined score column based on the given phenotypic scores and p-values.
    Combined score is calculated as:

        $combined\_score = \frac{score}{pseudo\_sd} \times -\log_{10}(pvalue)$

    Parameters:
        df_in (pandas.DataFrame): The input DataFrame.
        score_col (str): The column name for the individual scores. Default is 'score'.
        pvalue_col (str): The column name for the p-values. Default is 'pvalue'.
        target_col (str): The column name for the target variable. Default is 'target'.
        combined_score_col (str): The column name for the combined scores. Default is 'combined_score'.
        ctrl_label (str): The label for the control group. Default is 'control'.
    
    Returns:
        pandas.Series: The calculated combined score column.
    """
    # make a copy of input dataframe
    df = df_in.copy()

    for col in [score_col, pvalue_col, target_col]:
        if col not in df.columns:
            raise ValueError(f'Column "{col}" not found in the input DataFrame.')
    
    # calculate pseudo_sd
    pseudo_sd = df[df[target_col].eq(ctrl_label)][score_col].tolist()
    pseudo_sd = np.std(pseudo_sd)

    # calculate combined score
    return df[score_col]/pseudo_sd * -np.log10(df[pvalue_col])


def annotateScoreTable(df_in, up_hit, down_hit, threshold, score_col='score', pvalue_col='pvalue', target_col='target', ctrl_label='negative_control'):
    """
    Annotate the given score tabel 
    

    Parameters:
        df_in (pd.DataFrame): score dataframe
        up_hit (str): up hit label
        down_hit (str): down hit label
        threshold (int): threshold value
        score_col (str): score column name. Default is 'score'.
        target_col (str): column name for the target variable. Default is 'target'.
        pvalue_col (str): pvalue column name. Default is 'pvalue'.
        ctrl_label (str): control label value. Default is 'negative_control'.
    
    Returns:
        pd.DataFrame: annotated score dataframe
    """
    # make a copy of input dataframe
    df = df_in.copy()

    for col in [score_col, pvalue_col, target_col]:
        if col not in df.columns:
            raise ValueError(f'Column "{col}" not found in the input DataFrame.')

    df[score_col] = df[score_col].astype(float)
    df[pvalue_col] = df[pvalue_col].astype(float)

    # add combined score column
    df['combined_score'] = getCombinedScore(df, score_col, pvalue_col, ctrl_label)

    # add label column
    df['label'] = '.'

    # annotate hits: up
    df.loc[
        (df[score_col] > 0) & (~df[target_col].eq(ctrl_label)) &
        (df['combined_score'] >= threshold), 'label'
    ] = up_hit

    # annotate hits: down
    df.loc[
        (df[score_col] < 0) & (~df[target_col].eq(ctrl_label)) &
        (df['combined_score'] <= -threshold), 'label'
    ] = down_hit

    # annotate control
    df.loc[df[target_col].eq(ctrl_label), 'label'] = ctrl_label

    # annotate non-hit
    df.loc[df['label'] == '.', 'label'] = 'target_non_hit'

    # reorder factors
    df['label'] = pd.Categorical(
        df['label'],
        categories=[down_hit, up_hit, ctrl_label, 'target_non_hit']
    )

    return df
