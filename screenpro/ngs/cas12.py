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
        

def get_spacers_cas12(df_count,DRref):

    # get 1st spacer sequence
    df_count = df_count.with_columns(
            DR1_loc = df_count['sequence'].str.find(DRref['DR-1']),
        ).with_columns(
            pl.col('sequence').str.slice(
                pl.col("DR1_loc")-23, 23
            ).alias("SP1_sequence")
        )

    # get 2ed+ spacer sequence
    
    for DR_key, DR_seq in DRref.items():
        DR_n = int(DR_key[-1]) # get the number of the DR
        spacer_i = DR_n+1 # get the spacer number, 3' of the DR sequence
        df_count = df_count.with_columns(
            df_count['sequence'].str.find(DR_seq).alias(f"DR{DR_n}_loc")
        ).with_columns(
            pl.col('sequence').str.slice(
                pl.col(f"DR{DR_n}_loc")+19, 23
            ).alias(f"SP{spacer_i}_sequence")
        )

    out = df_count.select([
        f'SP{i}_sequence' for i in range(1,len(DRref)+2)
    ]+ ['count']).group_by([
        f'SP{i}_sequence' for i in range(1,len(DRref)+2)
    ]).sum()

    return df_count, out


def map_to_cas12_pairs_library(df_count,library,DR1_seq, get_recombinant=False, verbose=False):
    
    t0 = time()
    
    df_count, df_count_split = get_spacers_cas12(df_count, DRref = {'DR-1': DR1_seq})
    
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
        

    if get_recombinant:
        df_res_unmap = df_count_split.join(
            pl.DataFrame(
                library[['SP1_sequence','SP2_sequence']].reset_index()
            ), on=["SP1_sequence","SP2_sequence"], how="anti"
        )

        df_res_unmap_remapped_sp1 = df_res_unmap.join(
            pl.DataFrame(library[['SP1_name','SP1_id','SP1_sequence']]), on=["SP1_sequence"], how="left"
        )

        df_res_unmap_remapped_sp1_sp2 = df_res_unmap_remapped_sp1.join(
            pl.DataFrame(library[['SP2_name','SP2_id','SP2_sequence']]),
            on=["SP2_sequence"], how="left"
        ).drop_nulls().unique().with_columns(
            recombinant_name=pl.concat_str(
                [
                    pl.col("SP1_name"),
                    pl.col("SP1_id"),
                    pl.col("SP2_name"),
                    pl.col("SP2_id"),
                ],
                separator="_"
            )
        ).select(['recombinant_name','SP1_sequence','SP2_sequence','count'])
            
        if verbose:
            perc_remapped = df_res_unmap_remapped_sp1_sp2['count'].drop_nulls().sum() / df_count['count'].sum() * 100
            
            print(f"% counts remapped to library: {perc_remapped} [fully remapped recombination events]")
            
        if verbose: print("done in %0.3fs" % (time() - t0))
        
        return df_res, df_res_unmap_remapped_sp1_sp2
    
    else:
        if verbose: print("done in %0.3fs" % (time() - t0))

        return df_res


def map_to_cas12_triplets_library(df_count,library,DR1_seq, DR2_seq, get_recombinant=False, verbose=False):
    
    t0 = time()
    
    df_count, df_count_split = get_spacers_cas12(df_count, DRref = {'DR-1': DR1_seq, 'DR-2': DR2_seq})
    
    if verbose:
        perc_DR1 = df_count.with_columns(
            pl.col('DR1_loc').fill_null(0).gt(0)
        ).filter(
            pl.col('DR1_loc')
        ).get_column('count').sum() / df_count['count'].sum() * 100
        
        print(f"% counts with DR1: {perc_DR1}")

        perc_DR2 = df_count.with_columns(
            pl.col('DR2_loc').fill_null(0).gt(0)
        ).filter(
            pl.col('DR2_loc')
        ).get_column('count').sum() / df_count['count'].sum() * 100

        print(f"% counts with DR2: {perc_DR2}")
    
    df_res = pl.DataFrame(library[['SP1_sequence','SP2_sequence','SP3_sequence']].reset_index()).join(
        df_count_split, 
        on=["SP1_sequence","SP2_sequence","SP3_sequence"], how="left"
    )

    if verbose:
        perc_mapped = df_res['count'].drop_nulls().sum() / df_count['count'].sum() * 100
        
        print(f"% counts mapped to library: {perc_mapped}")

    if get_recombinant:
        df_res_unmap = df_count_split.join(
            pl.DataFrame(
                library[['SP1_sequence','SP2_sequence','SP3_sequence']].reset_index()
            ), on=["SP1_sequence","SP2_sequence","SP3_sequence"], how="anti"
        )

        df_res_unmap_remapped_sp1 = df_res_unmap.join(
            pl.DataFrame(library[['SP1_name','SP1_id','SP1_sequence']]), on=["SP1_sequence"], how="left"
        )

        df_res_unmap_remapped_sp1_sp2 = df_res_unmap_remapped_sp1.join(
            pl.DataFrame(library[['SP2_name','SP2_id','SP2_sequence']]),
            on=["SP2_sequence"], how="left"
        )

        df_res_unmap_remapped_sp1_sp2_sp3 = df_res_unmap_remapped_sp1_sp2.join(
            pl.DataFrame(library[['SP3_name','SP3_id','SP3_sequence']]),
            on=["SP3_sequence"], how="left"
        ).drop_nulls().unique().with_columns(
            recombinant_name=pl.concat_str(
                [
                    pl.col("SP1_name"),
                    pl.col("SP1_id"),
                    pl.col("SP2_name"),
                    pl.col("SP2_id"),
                    pl.col("SP3_name"),
                    pl.col("SP3_id"),
                ],
                separator="_"
            )
        ).select(['recombinant_name','SP1_sequence','SP2_sequence','SP3_sequence','count'])

        if verbose:
            perc_remapped = df_res_unmap_remapped_sp1_sp2_sp3['count'].drop_nulls().sum() / df_count['count'].sum() * 100
            
            print(f"% counts remapped to library: {perc_remapped} [fully remapped recombination events]")

        if verbose: print("done in %0.3fs" % (time() - t0))

        return df_res, df_res_unmap_remapped_sp1_sp2_sp3
    
    else:
        
        if verbose: print("done in %0.3fs" % (time() - t0))

        return df_res
