import numpy as np
import pandas as pd
import polars as pl
from CanDI import candi
from pathlib import Path


def coessentiality(pvalue_threshold = 10**-3, data_dir='auto'):
    if data_dir == 'auto':
        data_dir=str(Path(candi.__path__[0]).parent.absolute()) + '/setup/data/coessentiality'
    else:
        # check if the path exists and it contains the necessary files
        if not Path(data_dir).exists():
            raise ValueError(f"Path {data_dir} does not exist")
        if not Path(data_dir+'/genes.txt').exists():
            raise ValueError(f"Path {data_dir}/genes.txt does not exist")
        if not Path(data_dir+'/GLS_sign.npy').exists():
            raise ValueError(f"Path {data_dir}/GLS_sign.npy does not exist")
        if not Path(data_dir+'/GLS_p.npy').exists():
            raise ValueError(f"Path {data_dir}/GLS_p.npy does not exist")
    
    gene_names = pd.read_csv(f'{data_dir}/genes.txt',header=None,names=['gene_name'])['gene_name']
    
    GLS_sign = np.load(f'{data_dir}/GLS_sign.npy')
    GLS_p = np.load(f'{data_dir}/GLS_p.npy')
    
    coessentiality_mat = pd.DataFrame((-1*np.log10(GLS_p)) * GLS_sign, columns = gene_names, index = gene_names).reset_index()
    coessentiality_mat = pl.from_dataframe(coessentiality_mat)
    
    coessentiality_df = coessentiality_mat.melt('gene_name')
    coessentiality_df.columns = ['gene_1','gene_2','coessentiality']
    coessentiality_df = coessentiality_df.filter(~(pl.col('gene_1') == pl.col('gene_2')))
    coessentiality_df = coessentiality_df.filter(pl.col('coessentiality') > -np.log10(pvalue_threshold))
    
    out = coessentiality_df.to_pandas()

    return out
