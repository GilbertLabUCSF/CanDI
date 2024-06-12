import numpy as np
import pandas as pd
import anndata as ad

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from adpbulk import ADPBulk


def pseudobulk_by_group(adt, groups, method="mean"):
    # initialize the object
    adpb = ADPBulk(adt, groupby=groups, method=method)

    # perform the pseudobulking
    pseudobulk_matrix = adpb.fit_transform()

    # retrieve the sample metadata (useful for easy incorporation with edgeR)
    sample_meta = adpb.get_meta()

    out = ad.AnnData(
        X=pseudobulk_matrix,
        obs=sample_meta.set_index('SampleName')
    )

    return out


def run_deseq(adata, design, tested_level, ref_level, n_cpus=8):

    inference = DefaultInference(n_cpus=n_cpus)
    
    dds = DeseqDataSet(
        counts=adata.to_df().astype(int),
        metadata=adata.obs,
        design_factors=design,  # compare samples based on the "condition"
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    stat_res = DeseqStats(
        dds, 
        contrast=[design, tested_level, ref_level], 
        inference=inference
    )
    stat_res.summary()

    df = stat_res.results_df

    return df