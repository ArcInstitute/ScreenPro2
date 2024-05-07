=======
History
=======

0.3.0 (coming soon)
~~~~~~~~~~~~~~~~
* add command line interface

0.2.11 (May 2024)
~~~~~~~~~~~~~~~~
* introduce `counter` module
* improve core functionalities for CLI

0.2.7 - 0.2.10 (Mar 2024 - Apr 2024)
~~~~~~~~~~~~~~~~
* major reformatting and introduce alpha version of CLI
* introduce `ngs` module – process fastq files from single, dual, or multiplexed Cas9 and Cas12 screens

  * Support both single-guide and dual-guide Cas9 library design #37
  (i.e. V2 and V3 genome-scale dCas9 screen platforms)
  * Support higher-order Cas12 screens #39

* phenoScore and phenoStats modules

  * add missing features to support single-guide-design screens (i.e. V2 CRISPRi/a screens)

0.2.5 – 0.2.6 (Dec 2023 - Feb 2024)
~~~~~~~~~~~~~~~~
* public release on Arc's website https://arcinstitute.org/tools/screenpro2
* improve the conda environment and docker config files
* improve documentation
* split phenoScore and phenoStats modules
* fix bugs in phenotype score calculation and growth rate normalization

0.2.1 - 0.2.4 (July 2023 - Nov 2023)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Included PyPi package installation.
* Introduced `ScreenPro` class.
* Introduced a method in `ScreenPro` class for common drug screen analysis.

0.2.0 (May 2022)
~~~~~~~~~~~~~~~~
* Stable python 3 version.
* Reorganized `ScreenProcessing`_ pipeline into python modules
* Include ReadTheDocs documentations – `screenpro2.rtfd.io`_
* Enable exploratory data analysis

.. _ScreenProcessing: https://github.com/mhorlbeck/ScreenProcessing
.. _screenpro2.rtfd.io: https://screenpro2.rtfd.io
