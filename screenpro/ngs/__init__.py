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

from . import cas12
from . import cas9
from .counter import Counter