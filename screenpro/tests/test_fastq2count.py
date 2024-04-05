import polars as pl
from screenpro import ngs


def test_cas9_single_guide():
    df_count = ngs.cas9.fastq_to_count_single_guide(
        fastq_file_path='demo/step1_process_fastq_files/example_crispri_v2_sample.fastq.gz',
        trim5p_start=2,
        trim5p_length=19,
        verbose=True
    )

    assert isinstance(df_count, pl.DataFrame)


def test_cas9_dual_guide():
    df_count = ngs.cas9.fastq_to_count_dual_guide(
        R1_fastq_file_path='demo/step1_process_fastq_files/example_crispri_v3_sample_R1.fastq.gz',
        R2_fastq_file_path='demo/step1_process_fastq_files/example_crispri_v3_sample_R2.fastq.gz',
        trim5p_pos1_start=2,
        trim5p_pos1_length=19,
        trim5p_pos2_start=2,
        trim5p_pos2_length=19,
        verbose=True
    )

    assert isinstance(df_count, pl.DataFrame)
