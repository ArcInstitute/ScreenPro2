"""
deseq module: adapt pyDESeq2 for use in ScreenPro2 package
"""

import sys 
import numpy as np
import pandas as pd
import anndata as ad

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def runDESeq(adata, design, tested_level, ref_level, n_cpus=8,quiet=False):

    inference = DefaultInference(n_cpus=n_cpus)
    
    dds = DeseqDataSet(
        counts=adata.to_df().astype(int),
        metadata=adata.obs,
        design_factors=design,  # compare samples based on the "condition"
        refit_cooks=True,
        inference=inference,
        quiet=quiet
    )

    dds.deseq2()

    stat_res = DeseqStats(
        dds, 
        contrast=[design, tested_level, ref_level], 
        inference=inference,
        quiet=quiet
    )
    
    sys.stdout = open(stat_res.summary(), 'w')

    df = stat_res.results_df

    return df
