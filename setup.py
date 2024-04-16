from setuptools import setup, find_packages
from screenpro.__init__ import __version__, __author__, __email__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    python_requires='>=3.9',
    name='ScreenPro2',
    description="Flexible analysis of high-content CRISPR screening",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="MIT License",

    version=__version__,
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,

    url='https://github.com/ArcInstitute/ScreenPro2',
    packages=find_packages(include=['screenpro', 'screenpro.*']),
    entry_points={
        "console_scripts": ["screenpro=screenpro.main:main"],
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ]
)
