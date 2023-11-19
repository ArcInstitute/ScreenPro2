import pandas as pd
import numpy as np


def ann_score_df(df_in, up_hit='resistance_hit', down_hit='sensitivity_hit', ctrl_label='non-targeting', threshold=10):
    """
    Annotate score dataframe with hit labels using given `threshold`
    (i.e. `score/pseudo_sd * -np.log10(pvalue) >= threshold`).
    Args:
        df_in (pd.DataFrame): score dataframe
        up_hit (str): up hit label
        down_hit (str): down hit label
        ctrl_label (str): control label
        threshold (int): threshold
    Returns:
        pd.DataFrame: annotated score dataframe
    """
    # make a copy of input dataframe
    df = df_in.copy()

    # rename/reformat columns
    df.columns = ['target', 'score', 'pvalue']
    df['score'] = df['score'].astype(float)
    df['pvalue'] = df['pvalue'].astype(float)

    # calculate pseudo_sd
    pseudo_sd = df[df['target'].str.contains(ctrl_label)]['score'].tolist()
    pseudo_sd = np.std(pseudo_sd)

    df['label'] = '.'

    # annotate hits: up
    df.loc[
        (df['score'] > 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) >= threshold), 'label'
    ] = up_hit

    # annotate hits: down
    df.loc[
        (df['score'] < 0) & (~df['target'].str.contains(ctrl_label)) &
        (df['score']/pseudo_sd * -np.log10(df['pvalue']) <= -threshold), 'label'
    ] = down_hit

    # annotate control
    df.loc[df['target'].str.contains(ctrl_label), 'label'] = ctrl_label

    # annotate non-hit
    df.loc[df['label'] == '.', 'label'] = 'target_non_hit'

    # reorder factors
    df['label'] = pd.Categorical(
        df['label'],
        categories=[down_hit, up_hit, ctrl_label, 'target_non_hit']
    )

    return df
