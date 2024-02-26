'''Scripts to work with NGS data
'''

import gzip
from time import time
import pandas as pd
import polars as pl
import biobear as bb
from Bio import SeqIO	
from biobear.compression import Compression


def read_records(file_path):	
    if file_path.endswith('.gz'):	
        handle = gzip.open(file_path, "rt")	
    else:	
        handle = open(file_path, "rt")	

    records = []	
    for record in SeqIO.parse(handle, "fastq"): 
        records.append((str(record.id), str(record.description), str(record.seq), str(record.letter_annotations["phred_quality"])))
    handle.close()	

    return records	


def fastq_to_dataframe(fastq_file_path: str, engine='biopython') -> pl.DataFrame: 
    t0 = time()
    print('load FASTQ file as a Polars DataFrame')

    # Read the FASTQ file using read_records function
    records = read_records(fastq_file_path)
    
    # Create a Polars DataFrame from the list of tuples
    if engine == 'biopython':
        df = pd.DataFrame(records, columns=['name', 'description', 'sequence', 'quality_scores'], dtype='str')
        df = pl.from_pandas(df)	
    elif engine == 'biobear':
        if '.gz' in fastq_file_path:
            df = bb.FastqReader(fastq_file_path,compression=Compression.GZIP).to_polars()
        else:
            df = bb.FastqReader(fastq_file_path).to_polars()

    print("done in %0.3fs" % (time() - t0))


def fastq_to_count_unique_seq(fastq_file_path: str, engine: str='biopython', slice_seq: list=None) -> pl.DataFrame:
    df = fastq_to_dataframe(fastq_file_path, engine=engine)

    t0 = time()
    print('Count unique sequences')

    # keep full sequence or slice it
    if slice_seq:
        # make a copy of the original sequence column into a new column called 'fullsequence'
        df = df.rename({"sequence":"fullsequence"})
        
        df = df.with_columns(
            sequence = df.get_column('fullsequence').str.slice(slice_seq[0], slice_seq[1])
        )
        
        # drop fullsequence column
        df = df.drop('fullsequence')
        
    df = df.drop(['name','description','quality_scores'])
    
    df_count = df.group_by('sequence').len().rename({"len":"count"})

    print("done in %0.3fs" % (time() - t0))

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

