import numpy as np
import anndata as ad
import pandas as pd
from screenpro.phenoscore import runPhenoScore


def test_runPhenoScore():
    # create test data
    cond_A = np.random.randint(0, 30, size=(3, 10))
    cond_B = np.random.randint(20, 100, size=(3, 10))

    adat = ad.AnnData(
        X=np.concatenate([cond_A, cond_B], axis=0),
        obs=pd.DataFrame(
            {
                'condition': ['A'] * 3 + ['B'] * 3
            },
            index=pd.Index(['sample_' + str(i) for i in range(6)], name='sample')
        ),
        var=pd.DataFrame(
            {
                'target': ['targetID_' + str(i) for i in range(10)],
                'sequence': [
                    'ATGCGTACATGTATGCGTG',
                    'ATGCGTATGCATATGCGTC',
                    'ATGCGTATGCGTCATCGTG',
                    'ATGCGTATGCGTATGCATC',
                    'CATCGTATGCGTATGCGTG',
                    'ATGCGTACATGTATGCGTG',
                    'ATGCGTATGCATATGCGTC',
                    'ATGCGTATGCGTCATCGTG',
                    'ATGCGTATGCGTATGCATC',
                    'CATCGTATGCGTATGCGTG'
                ],
                'targetType': ['gene'] * 8 + ['negCtrl'] * 2
            },
            index=pd.Index(['targetID_' + str(i) for i in range(10)], name='target')
        )
    )

    print(adat.to_df())

    assert isinstance(adat, ad.AnnData)

    # run function
    result_name, result = runPhenoScore(
        adata=adat,
        cond1='A',
        cond2='B',
        transformation='log2(x+1)',
        test='ttest',
        score_level='compare_reps',
        growth_rate=1,
        n_reps=2,
        ctrl_label='negCtrl'
    )

    # check result name
    assert result_name == 'B_vs_A'

    # check result dataframe
    assert isinstance(result, pd.DataFrame)