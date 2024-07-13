## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

'''ScreenPro2: A Python package for pooled CRISPR screens analysis

This package contains several modules, including:

**Main modules:**
- ngs: tools for generating counts from NGS data
- phenoscore: tools for calculating phenoscores
- assays: wrappers for analyzing CRISPR screens data from standard assays

**Additional modules:**
- load: tools for loading and saving data
- visualize: tools for visualizing data
- datasets: API for accessing pre-processed datasets
'''

from . import phenoscore as ps

from . import ngs
from . import assays
from . import load
from . import visualize as viz

from .ngs import GuideCounter
from .assays import PooledScreens, GImaps

__version__ = "0.4.1"
__author__ = "Abe Arab"
__email__ = 'abea@arcinstitute.org' # "abarbiology@gmail.com"
