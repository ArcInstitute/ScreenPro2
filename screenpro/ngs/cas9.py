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


def map_to_library_dual_guide(df_count, library, get_recombinant=False, return_type='all', verbose=False):
    # get counts for given input
    res = df_count.copy()
    res = res.sort('count', descending=True)
    res = res.with_columns(
        pl.concat_str(
            [
                pl.col("protospacer_a"),
                pl.col("protospacer_b")
            ],
            separator=";",
        ).alias("sequence"),
    )

    res_map = pl.DataFrame(library).join(
            res, on="sequence", how="left"
        )
    if verbose:
        print("% mapped reads",
            100 * \
            res_map.to_pandas()['count'].fillna(0).sum() / \
            int(res.select(pl.sum("count")).to_pandas()['count'])
        )
    
    if get_recombinant:
        res_unmap = res.join(
            pl.DataFrame(library), on="sequence", how="anti"
        )

        if verbose:
            print("% unmapped reads",
                100 * \
                res_unmap.to_pandas()['count'].fillna(0).sum() / \
                int(res.select(pl.sum("count")).to_pandas()['count'])
            )
        
        res_unmap_remapped_a = res_unmap.join(
            pl.DataFrame(library[['sgID_A','protospacer_a']]), on=["protospacer_a"], how="left"
        )

        res_recomb_events = res_unmap_remapped_a.join(
            pl.DataFrame(library[['sgID_B','protospacer_b']]), 
            on=["protospacer_b"], how="left"
        )
        if verbose:            
            print("% fully remapped recombination events",
                100 * \
                res_recomb_events.drop_nulls().to_pandas()['count'].fillna(0).sum() / \
                int(res.select(pl.sum("count")).to_pandas()['count'])
            )
    
    if get_recombinant and return_type == 'all':
        sample_count = {'full': res,'mapped': res_map,'recomb': res_recomb_events}
        return sample_count
    elif get_recombinant and return_type == 'recomb':
        return res_recomb_events
    elif not get_recombinant and return_type == 'all':
        sample_count = {'full': res,'mapped': res_map}
        return sample_count
    elif not get_recombinant and return_type == 'mapped':
        return res_map
    else:
        raise ValueError("return_type must be either 'all' or 'mapped' or 'recomb'")
