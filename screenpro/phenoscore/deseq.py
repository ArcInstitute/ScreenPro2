"""
deseq module: adapt pyDESeq2 for use in ScreenPro2 package
"""

import numpy as np
import pandas as pd
import anndata as ad
import os, contextlib

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def runDESeq(adata, design, n_cpus=8,quiet=False):

    inference = DefaultInference(n_cpus=n_cpus)

    print(f'\tcreating `dds` object...')
    
    dds = DeseqDataSet(
        counts=adata.to_df().astype(int),
        metadata=adata.obs,
        design_factors=design,  # compare samples based on the "condition"
        refit_cooks=True,
        inference=inference,
        quiet=quiet
    )

    dds.deseq2()

    return dds


def extractDESeqResults(dds, design, ref_level, tested_level, n_cpus=8, quiet=False):

    inference = DefaultInference(n_cpus=n_cpus)

    result_name = f'{tested_level}_vs_{ref_level}'

    print(f'\t{tested_level}_vs_{ref_level}')

    stat_res = DeseqStats(
        dds, 
        contrast=[design, tested_level, ref_level], 
        inference=inference,
        quiet=quiet
    )
    
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            stat_res.summary()

    results = stat_res.results_df

    return result_name, results
