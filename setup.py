from setuptools import setup, find_packages
from screenpro import __version__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    python_requires='>=3.9',
    name='ScreenPro2',
    description="Analyze pooled CRISPR screens",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Abe Arab',
    author_email='abea@arcinstitute.org',
    url='https://github.com/abearab/ScreenPro2',
    packages=find_packages(include=['screenpro', 'screenpro.*']),
    version=__version__,
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ]
)
