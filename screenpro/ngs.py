'''Scripts to work with NGS data
'''

import click
from time import time
import pandas as pd
import polars as pl
import biobear as bb


def load_fastq(fastq_file_path: str, verbose: bool=False):
    t0 = time()
    if verbose: print('load FASTQ file ...')

    # Read the FASTQ file using read_records function
    if '.gz' in fastq_file_path:
        out = bb.FastqReader(fastq_file_path,compression=bb.Compression.GZIP)
    else:
        out = bb.FastqReader(fastq_file_path)

    if verbose: print("done in %0.3fs" % (time() - t0))

    return out


def fastq_to_count_unique_seq(fastq_file_path:str, trim5p_start:int=None, trim5p_length:int=None, verbose: bool=False) -> pl.DataFrame:
    
    if verbose: ('count unique sequences ...')
    t0 = time()
    
    session = bb.connect()
    if trim5p_start and trim5p_length:
        df_count = session.sql(f"""
        SELECT substr(f.sequence, {trim5p_start}, {trim5p_length}) AS sequence, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY substr(f.sequence, {trim5p_start}, {trim5p_length})
        """
        ).to_polars()
    else:
        df_count = session.sql(f"""
        SELECT f.sequence AS sequence, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY f.sequence
        """
        ).to_polars()
    
    if verbose: print("done in %0.3fs" % (time() - t0))

    return df_count


def map_sample_counts_to_library(library, sample):
    counts_df = library.copy()

    overlap = list(set(library.index.tolist()) & set(sample['sequence'].to_list()))
    non_overlap = list(set(sample['sequence'].to_list()) - set(library.index.tolist()))
    
    # sum counts of overlapping sequences
    n_mapped_counts = sample.set_index('sequence').loc[overlap, 'count'].sum()
    n_unmapped_counts = sample.set_index('sequence').loc[non_overlap, 'count'].sum()

    # print number of overlapping and non-overlapping sequences
    print(f"% mapped sequences: {n_mapped_counts/sample['count'].sum():.2f}")
    print(f"% non-mapped sequences: {n_unmapped_counts/sample['count'].sum():.2f}")

    counts_df['counts'] = 0
    counts_df.loc[overlap, 'counts'] = sample.set_index('sequence').loc[overlap, 'count']

    return counts_df


@click.command()
@click.argument('fastq_file_path', type=click.Path(exists=True))
@click.option('--engine', default='biopython', help='Engine to use for reading FASTQ file')
@click.option('--slice-seq', nargs=2, type=int, help='Slice sequence range')
def main(fastq_file_path, engine, slice_seq):
    """
    Command line interface for working with NGS data.
    """
    fastq_to_count_unique_seq(fastq_file_path, engine=engine, slice_seq=slice_seq)


if __name__ == '__main__':
    main()
