'''Scripts to work with NGS data
'''

from time import time
import pandas as pd
import polars as pl
import biobear as bb
from biobear.compression import Compression


def fastq_to_dataframe(fastq_file_path: str, num_threads: int, seq_only=True) -> pl.DataFrame:
    """
    Reads a FASTQ file and returns a Polars DataFrame with the following columns:
    - 'id': the sequence ID (e.g. "@SEQ_ID")
    - 'seq': the nucleotide sequence
    - 'qual': the sequence quality scores
    """
    t0 = time()
    print('load FASTQ file as a Polars DataFrame')

    if '.gz' in fastq_file_path:
        df = bb.FastqReader(fastq_file_path,compression=Compression.GZIP).to_polars()
    else:
        df = bb.FastqReader(fastq_file_path).to_polars()
    print(df.head())

    print("done in %0.3fs" % (time() - t0))

    # return df


def fastq_to_count_unique_seq(fastq_file_path: str, num_threads: int) -> pl.DataFrame:
    df = fastq_to_dataframe(fastq_file_path, num_threads)

    t0 = time()
    print('Count unique sequences')

    df_count = df.groupby('seq').count()

    print("done in %0.3fs" % (time() - t0))

    return df_count


def map_sample_counts_to_library(library, sample):
    counts_df = library.copy()

    ol = list(set(library.index.tolist()) & set(sample['seq'].to_list()))

    counts_df['counts'] = 0
    counts_df.loc[ol, 'counts'] = sample.to_pandas().set_index('seq').loc[ol, 'count']

    return counts_df.reset_index(drop=True).set_index('oligoname')