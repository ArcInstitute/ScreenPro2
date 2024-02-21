'''Scripts to work with NGS data
'''

from time import time
import pandas as pd
import polars as pl
import biobear as bb
from biobear.compression import Compression


def fastq_to_dataframe(fastq_file_path: str,n_bp_from_5p=0, n_bp_from_3p=0) -> pl.DataFrame:
    """
    Reads a FASTQ file and returns a four columns Polars DataFrame
    - 'name': the name of the sequence
    - 'description': the description line 
    - 'sequence': the nucleotide sequence
    - 'quality_scores': the sequence quality scores

    Parameters
    ----------
    fastq_file_path : str
        The path to the FASTQ file
    """
    t0 = time()
    print('load FASTQ file as a Polars DataFrame')

    if '.gz' in fastq_file_path:
        df = bb.FastqReader(fastq_file_path,compression=Compression.GZIP).to_polars()
    else:
        df = bb.FastqReader(fastq_file_path).to_polars()

    print("done in %0.3fs" % (time() - t0))

    return df


def trim_ngs_reads(df: pl.DataFrame, n_bp_from_5p: int, n_bp_from_3p: int) -> pl.DataFrame:
    """
    Trims the sequences in a Polars DataFrame

    Parameters
    ----------
    df : pl.DataFrame
        The Polars DataFrame containing the sequences
    n_bp_from_5p : int
        The number of base pairs to trim from the 5' end
    n_bp_from_3p : int
        The number of base pairs to trim from the 3' end
    """
    t0 = time()
    print('Trim sequences')

    df['sequence'] = df['sequence'].str.slice(n_bp_from_5p, -n_bp_from_3p)

    print("done in %0.3fs" % (time() - t0))

    return df


def fastq_to_count_unique_seq(fastq_file_path: str) -> pl.DataFrame:
    df = fastq_to_dataframe(fastq_file_path)

    t0 = time()
    print('Count unique sequences')

    df_count = df.groupby('sequence').count()

    print("done in %0.3fs" % (time() - t0))

    return df_count


def map_sample_counts_to_library(library, sample):
    counts_df = library.copy()

    overlap = list(set(library.index.tolist()) & set(sample['sequence'].to_list()))

    counts_df['counts'] = 0
    counts_df.loc[overlap, 'counts'] = sample.set_index('sequence').loc[overlap, 'count']

    return counts_df.reset_index(drop=True).set_index('oligoname')