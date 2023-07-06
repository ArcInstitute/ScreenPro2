from setuptools import setup, find_packages
from screenpro import __version__

setup(
    python_requires='>=3.6',
    name='ScreenPro2',
    description="Analyze pooled CRISPR screens",
    packages=find_packages(include=['screenpro', 'screenpro.*']),
    install_requires=[
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "biopython",
        "matplotlib",
        "seaborn",
        "pandas",
        "numpy",
        "anndata>=0.8.0",
        "pydeseq2"
    ],
    version=__version__
)
