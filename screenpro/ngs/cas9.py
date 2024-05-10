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
        SELECT substr(f.sequence, {trim5p_start}, {trim5p_length}) AS protospacer, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY protospacer
        """
    else:
        sql_cmd = f"""
        SELECT f.sequence AS protospacer, COUNT(*) as count
        FROM fastq_scan('{fastq_file_path}') f
        GROUP BY protospacer
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
            SELECT substr(pos1.sequence, {trim5p_pos1_start}, {trim5p_pos1_length}) protospacer_A, reverse_complement(substr(pos2.sequence, {trim5p_pos2_start}, {trim5p_pos2_length})) protospacer_B, COUNT(*) count
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


def map_to_library_single_guide(df_count, library, return_type='all', verbose=False):
    # get counts for given input
    res = df_count.clone() #cheap deepcopy/clone
    res = res.sort('count', descending=True)

    res = res.with_columns(
        pl.col("protospacer").alias("sequence"),
    )

    res_map = pl.DataFrame(library).join(
        res, on="sequence", how="left"
    )

    if return_type == 'unmapped' or return_type == 'all':
        res_unmap = res.join(
            pl.DataFrame(library), on="sequence", how="anti"
        )

    if verbose:
        print("% mapped reads",
            100 * \
            res_map.to_pandas()['count'].fillna(0).sum() / \
            int(res.select(pl.sum("count")).to_pandas()['count'])
        )
    
    if return_type == 'unmapped':
        return res_unmap
    elif return_type == 'mapped':
        return res_map
    elif return_type == 'all':
        return {'full': res, 'mapped': res_map, 'unmapped': res_unmap}
    else:
        raise ValueError("return_type must be either 'unmapped', 'mapped', or 'all'")


def map_to_library_dual_guide(df_count, library, get_recombinant=False, return_type='all', verbose=False):
    # get counts for given input
    res = df_count.clone() #cheap deepcopy/clone
    res = res.rename(
        {'protospacer_a':'protospacer_A','protospacer_b':'protospacer_B'}
    )
    res = res.sort('count', descending=True)
    res = res.with_columns(
        pl.concat_str(
            [
                pl.col("protospacer_A"),
                pl.col("protospacer_B")
            ],
            separator=";",
        ).alias("sequence"),
    )

    res_map = pl.DataFrame(library).join(
            res, on="sequence", how="left"
        )

    if get_recombinant or return_type == 'unmapped' or return_type == 'all':
        res_unmap = res.join(
            pl.DataFrame(library), on="sequence", how="anti"
        )

    if verbose:
        print("% mapped reads",
            100 * \
            res_map.to_pandas()['count'].fillna(0).sum() / \
            int(res.select(pl.sum("count")).to_pandas()['count'])
        )
    
    if get_recombinant:

        if verbose:
            print("% unmapped reads",
                100 * \
                res_unmap.to_pandas()['count'].fillna(0).sum() / \
                int(res.select(pl.sum("count")).to_pandas()['count'])
            )
        
        res_unmap_remapped_a = res_unmap.join(
            pl.DataFrame(library[['sgID_A','protospacer_A']]), on=["protospacer_A"], how="left"
        )

        res_recomb_events = res_unmap_remapped_a.join(
            pl.DataFrame(library[['sgID_B','protospacer_B']]), 
            on=["protospacer_B"], how="left"
        )
        if verbose:            
            print("% fully remapped recombination events",
                100 * \
                res_recomb_events.drop_nulls().to_pandas()['count'].fillna(0).sum() / \
                int(res.select(pl.sum("count")).to_pandas()['count'])
            )
    
    if return_type == 'unmapped':
        return res_unmap
    elif return_type == 'mapped':
        return res_map
    elif return_type == 'recombinant':
        if get_recombinant:
            return res_recomb_events
        else:
            raise ValueError("get_recombinant must be set to True to calculate recombinant events")
    elif return_type == 'all':
        if get_recombinant:
            return {'full': res,'mapped': res_map,'recombinant': res_recomb_events, 'unmapped': res_unmap}
        else:
            return {'full': res,'mapped': res_map, 'unmapped': res_unmap}
    else:
        raise ValueError("return_type must be either 'unmapped', 'mapped', 'recombinant', or 'all'")
