import pandas as pd
import numpy as np


def find_low_counts(adata, filter_type='either', minimum_reads=50):
    """
    Label variables with low counts in either or all samples.

    Args:
        adata: AnnData object
        filter_type: either or all
        minimum_reads:

    Returns: it's adding a True/False column to `adata.var['low_count']`
        None
    """
    count_bin = adata.X >= minimum_reads
    if filter_type == 'either':
        out = adata[:, count_bin.any(axis=0)].copy()
    elif filter_type == 'all':
        out = adata[:, count_bin.all(axis=0)].copy()

    # print the number of removed variables
    n_removed = adata.shape[1] - out.shape[1]
    print(
        f"{n_removed} variables with less than {minimum_reads} reads in {filter_type} replicates / experiment"
    )

    adata.var['low_count'] = ~adata.var.index.isin(out.var.index.to_list())


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
