import biobear as bb
import polars as pl
from time import time


def fastq_to_count_merged_reads(
        fastq_file_path:str, 
        verbose: bool=False) -> pl.DataFrame:
    if verbose: ('count unique sequences ...')
    t0 = time()
    
    session = bb.connect()
    
    sql_cmd = f"""
    SELECT f.sequence AS sequence, COUNT(*) as count
    FROM fastq_scan('{fastq_file_path}') f
    GROUP BY sequence
    """
    
    df_count = session.sql(sql_cmd).to_polars()

    if verbose: print("done in %0.3fs" % (time() - t0))

    return df_count
        

def get_spacers_cas12(df_count,DRref = DRref):

    df_count = df_count.with_columns(
        DR1_loc = df_count['sequence'].str.find(DRref['DR-1']),
    ).with_columns(
        pl.col('sequence').str.slice(
            pl.col("DR1_loc")-23, 23
        ).alias("SP1_sequence"),

        pl.col('sequence').str.slice(
            pl.col("DR1_loc")+19, 23
        ).alias("SP2_sequence")
    )
    
    out = df_count.select([
        'SP1_sequence','SP2_sequence','count'
    ]).group_by(['SP1_sequence','SP2_sequence']).sum()

    return df_count, out


def map_to_cas12_pairs_library(df_count,library,verbose=False):
    
    t0 = time()
    
    df_count, df_count_split = get_spacers_cas12(df_count)
    
    if verbose:
        perc_DR1 = df_count.with_columns(
            pl.col('DR1_loc').fill_null(0).gt(0)
        ).filter(
            pl.col('DR1_loc')
        ).get_column('count').sum() / df_count['count'].sum() * 100
        
        print(f"% counts with DR1: {perc_DR1}")

    df_res = pl.DataFrame(library[['SP1_sequence','SP2_sequence']].reset_index()).join(
        df_count_split, 
        on=["SP1_sequence","SP2_sequence"], how="left"
    )
    
    if verbose:
        perc_mapped = df_res['count'].drop_nulls().sum() / df_count['count'].sum() * 100
        
        print(f"% counts mapped to library: {perc_mapped}")
        

    df_res_unmap = df_count_split.join(
        pl.DataFrame(
            library[['SP1_sequence','SP2_sequence']].reset_index()
        ), on=["SP1_sequence","SP2_sequence"], how="anti"
    )

    df_res_sp1_mapped = pl.DataFrame(
        library[['SP1_name','SP1_id','SP1_sequence','SP2_name','SP2_id','SP2_sequence']].reset_index()
    ).join(
        df_res_unmap,
        on=["SP1_sequence"], how="anti"
        )

    df_res_unmap_remapped_sp1 = df_res_unmap.join(
        pl.DataFrame(library[['SP1_name','SP1_id','SP1_sequence']]), on=["SP1_sequence"], how="left"
    )

    df_res_unmap_remapped_sp1_and_sp2 = df_res_unmap_remapped_sp1.join(
        pl.DataFrame(library[['SP2_name','SP2_id','SP2_sequence']]),
        on=["SP2_sequence"], how="left"
    ).drop_nulls().unique().with_columns(
        oligoname=pl.concat_str(
            [
                pl.col("SP1_name"),
                pl.col("SP1_id"),
                pl.col("SP2_name"),
                pl.col("SP2_id"),
            ],
            separator="_"
        )
    ).select(['oligoname','SP1_sequence','SP2_sequence','count'])
        
    if verbose:
        perc_remapped = df_res_unmap_remapped_sp1_and_sp2['count'].drop_nulls().sum() / df_count['count'].sum() * 100
        
        print(f"% counts remapped to library: {perc_remapped} [fully remapped recombination events]")
        
    if verbose: print("done in %0.3fs" % (time() - t0))
    
    return df_res, df_res_unmap_remapped_sp1_and_sp2
