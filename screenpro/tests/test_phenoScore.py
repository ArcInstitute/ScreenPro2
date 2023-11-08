import numpy as np
import pandas as pd
import scanpy as sc
from screenpro.phenoScore import runPhenoScore

def test_runPhenoScore():
    # create test data
    adata = sc.datasets.pbmc68k_reduced()
    assert isinstance(adata, sc.AnnData)
    adata.obs['condition'] = ['A'] * int(adata.obs.shape[0] / 2) + ['B'] * int(adata.obs.shape[0] / 2)
    adata.var['targetType'] = ['gene'] * adata.var.shape[0]
    adata.layers['seq_depth_norm'] = np.random.randint(0, 100, size=(adata.n_obs, adata.n_vars))

    # run function
    # result_name, result = runPhenoScore(
    #     adata=adata,
    #     cond1='A',
    #     cond2='B',
    #     math='log2(x+1)',
    #     test='MW',
    #     score_level='compare_reps',
    #     growth_rate=1,
    #     n_reps=2,
    #     ctrl_label='negCtrl'
    # )

    # # check result name
    # assert result_name == 'B_vs_A'

    # # check result dataframe
    # assert isinstance(result, pd.DataFrame)
    # assert result.shape == (adata.n_vars, 4)
    # assert set(result.columns) == {'target', 'score', 't-test pvalue', 'BH adj_pvalue'}