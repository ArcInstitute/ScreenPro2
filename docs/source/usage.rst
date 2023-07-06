Detailed documentation
======================

Step 1
______
counting sgRNAs in raw sequencing files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: bash

  python screenpro/processing/counts/fq2count_v2.py -h

::

   usage: fq2count_v2.py [-h] [-p PROCESSORS] [--trim_start TRIM_START] [--trim_end TRIM_END] [--test] Library_Fasta Out_File_Path Seq_File_Names [Seq_File_Names ...]

   Process raw sequencing data from screens to counts files in parallel

   positional arguments:
     Library_Fasta         Fasta file of expected library reads.
     Out_File_Path         Directory where output files should be written.
     Seq_File_Names        Name(s) of sequencing file(s). Unix wildcards can be used to select multiple files at once. The script will search for all *.fastq.gz,
                           *.fastq, and *.fa(/fasta/fna) files with the given wildcard name.

   optional arguments:
     -h, --help            show this help message and exit
     -p PROCESSORS, --processors PROCESSORS
     --trim_start TRIM_START
     --trim_end TRIM_END
     --test                Run the entire script on only the first 10000 reads of each file. Be sure to delete or move all test files before re-running script as they
                           will not be overwritten.

--------------

.. code:: bash

   python screenpro/processing/counts/fq2count_v3.py -h

::

    usage: fq2count_v3.py [-h] [–test] Guide_Table UMI_Table
    Out_File_Path Seq_File_Names [Seq_File_Names …]

    Process raw sequencing data from screens to counts files in parallel.

    positional arguments: Guide_Table Table of sgRNA pairs in the library.
    UMI_Table Table of sgRNA pairs in the library. Out_File_Path Directory
    where output files should be written. Seq_File_Names Name(s) of
    sequencing file(s). Unix wildcards can be used to select multiple files
    at once. The script will search for all *.fastq.gz,*.fastq, and
    \*.fa(/fasta/fna) files with the given wildcard name.

    optional arguments: -h, –help show this help message and exit –test Run
    the entire script on only the first 100000 reads of each file. Be sure
    to delete or move all test files before re-running script as they will
    not be overwritten.

--------------

Step 2
______
calculating sgRNA-level and gene-level phenotypes and p-values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: bash

   python -m screenpro.processing.scores -h

::

   usage: scores.py [-h] [--plot_extension PLOT_EXTENSION] Config_File Library_File_Directory

   Calculate sgRNA- and gene-level phenotypes based on sequencing read counts, as specified by the experiment config file.

   positional arguments:
     Config_File           Experiment config file specifying screen analysis settings (see accomapnying BLANK and DEMO files).
     Library_File_Directory
                           Directory containing reference library tables and the library_config.txt file.

   optional arguments:
     -h, --help            show this help message and exit
     --plot_extension PLOT_EXTENSION
                           Image extension for plot files, or "off". Default is png.

--------------

Step 3
______

making custom graphs
~~~~~~~~~~~~~~~~~~~~

COMING SOON.
