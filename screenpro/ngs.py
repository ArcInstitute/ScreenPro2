## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

## This part of the software is conceptualized and developed by Abolfazl Arab (@abearab)
## with support from the Nick Youngblut (@nick-youngblut).

'''Scripts to work with NGS data

This module provides functions to process FASTQ files from screens with single or dual guide
libraries. In general, the algorithm is fairly simple:
1. Read the FASTQ file and extract the proper sequences
2. Count the exact number of occurrences for each unique sequence
3. Map the counted sequences to the reference sequence library
4. Return the counted mapped or unmapped events as a dataframe(s)

For single-guide screens, the sequences are counted as single protospacer
from a single-end read file (R1). Then, these sequences are mapped to the reference
library of protospacer sequences.

For dual-guide screens, the sequences are counted as pairs of protospacer A and B
from paired-end read files (R1 and R2). Then, sequences are mapped to the reference
library of protospacer A and B pairs.

Theoretically, the algorithm is able to detect any observed sequence since it is counting first
and then mapping. Therefore, the recombination events can be detected. In dual-guide design
protospacer A and B are not the same pairs as in the reference library. These events include:
- Protospacer A and B pairs are present in the reference library but paired differently
- Only one of the protospacer A and B is present in the reference library
- None of the protospacer A and B is present in the reference library
'''

from time import time
import pandas as pd
import polars as pl
import biobear as bb


def fastq_to_count_single_guide(
        fastq_file_path:str, 
        trim5p_start:int=None, trim5p_length:int=None, 
        verbose: bool=False) -> pl.DataFrame:
    
    if verbose: ('count unique sequences ...')
    t0 = time()
    
    session = bb.connect()
    
    if trim5p_start and trim5p_length:
        sql_cmd = f"""
        SELECT substr(f.sequence, {trim5p_start}, {trim5p_length}) AS sequence, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY sequence
        """
    else:
        sql_cmd = f"""
        SELECT f.sequence AS sequence, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY sequence
        """
    
    df_count = session.sql(sql_cmd).to_polars()

    if verbose: print("done in %0.3fs" % (time() - t0))

    return df_count


def fastq_to_count_dual_guide(
        R1_fastq_file_path:str, R2_fastq_file_path:str,
        trim5p_pos1_start:int=None, trim5p_pos1_length:int=None,
        trim5p_pos2_start:int=None, trim5p_pos2_length:int=None,
        verbose: bool=False) -> pl.DataFrame:
    
    if verbose: ('count unique sequences ...')
    t0 = time()

    session = bb.connect()

    if trim5p_pos1_start and trim5p_pos1_length and trim5p_pos2_start and trim5p_pos2_length:
        sql_cmd = f"""
            WITH pos1 AS (
                SELECT REPLACE(name, '_R1', '') trimmed_name, *
                FROM fastq_scan('{R1_fastq_file_path}')
            ), pos2 AS (
                SELECT REPLACE(name, '_R2', '') trimmed_name, *
                FROM fastq_scan('{R2_fastq_file_path}')
            )
            SELECT substr(pos1.sequence, {trim5p_pos1_start}, {trim5p_pos1_length}) protospacer_A, substr(reverse_complement(pos2.sequence), {trim5p_pos2_start}, {trim5p_pos2_length}) protospacer_B, COUNT(*) count
            FROM pos1
            JOIN pos2
                ON pos1.name = pos2.name
            GROUP BY protospacer_A, protospacer_B
        """
    elif trim5p_pos1_start==None and trim5p_pos1_length==None and trim5p_pos2_start==None and trim5p_pos2_length==None:
        sql_cmd = f"""
            WITH pos1 AS (
                SELECT REPLACE(name, '_R1', '') trimmed_name, *
                FROM fastq_scan('{R1_fastq_file_path}')
            ), pos2 AS (
                SELECT REPLACE(name, '_R2', '') trimmed_name, *
                FROM fastq_scan('{R2_fastq_file_path}')
            )
            SELECT pos1.sequence protospacer_A, reverse_complement(pos2.sequence) protospacer_B, COUNT(*) count
            FROM pos1
            JOIN pos2
                ON pos1.name = pos2.name
            GROUP BY protospacer_A, protospacer_B
        """
    else:
        raise ValueError("trim5p_pos1_start, trim5p_pos1_length, \
                         trim5p_pos2_start, and trim5p_pos2_length \
                         must be provided concurrently!")

    df_count = session.sql(sql_cmd).to_polars()

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