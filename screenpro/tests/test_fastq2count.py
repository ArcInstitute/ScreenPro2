from .screenpro import ngs


def test_single_guide():
    ngs.fastq_to_count_single_guide(
        fastq_file_path='demo/step1_process_fastq_files/example_crispri_v2_sample.fastq',
        trim5p_start=2,
        trim5p_length=19,
        verbose=True
    )


def test_dual_guide():
    ngs.fastq_to_count_dual_guide(
        R1_fastq_file_path='demo/step1_process_fastq_files/example_crispri_v3_sample_R1.fastq',
        R2_fastq_file_path='demo/step1_process_fastq_files/example_crispri_v3_sample_R2.fastq',
        trim5p_pos1_start=2,
        trim5p_pos1_length=19,
        trim5p_pos2_start=2,
        trim5p_pos2_length=19,
        verbose=True
    )