## Copyright (c) 2022-2025 ScreenPro2 Development Team.
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

from . import ngs
from . import load
from . import preprocessing as pp
from . import phenoscore as ps
from . import assays
from . import plotting as pl
from . import dashboard

from .ngs import GuideCounter
from .assays import PooledScreens, GImaps
from .dashboard import DrugScreenDashboard


def _get_version():

    import os

    pyproject_path = os.path.join(os.path.dirname(__file__), "..", "pyproject.toml")

    with open(pyproject_path, "r") as pyproject_file:
        for line in pyproject_file.readlines():
            if "version" in line:
                return line.split("=")[1].strip().strip('"')


try:
    __version__ = _get_version()
except Exception:
    __version__ = "Unknown"
