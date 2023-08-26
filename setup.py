from setuptools import setup, find_packages
from screenpro import __version__
from pathlib import Path

this_directory = Path(__file__)
long_description = (this_directory / "README.rst").read_text()

setup(
    python_requires='>=3.8',
    name='ScreenPro2',
    description="Analyze pooled CRISPR screens",
    packages=find_packages(include=['screenpro', 'screenpro.*']),
    # install_requires=[
    #     "anndata",
    #     "numpy",
    #     "pandas",
    #     "scipy",
    #     "biopython",
    #     "matplotlib",
    #     "seaborn",
    #     "pandas",
    #     "numpy",
    #     "anndata>=0.8.0",
    #     "pydeseq2"
    # ],
    version=__version__,
    # other arguments omitted
    long_description=long_description,
    long_description_content_type='text/x-rst',
)
