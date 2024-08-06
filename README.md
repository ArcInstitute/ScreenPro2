[![website](https://img.shields.io/badge/website-live-brightgreen)](https://arcinstitute.org/tools/screenpro2)
[![PyPI version](https://badge.fury.io/py/ScreenPro2.svg)](https://badge.fury.io/py/ScreenPro2)
[![Documentation Status](https://readthedocs.org/projects/screenpro2/badge/?version=latest)](https://screenpro2.readthedocs.io/en/latest/?version=latest)
[![Downloads](https://static.pepy.tech/badge/screenpro2)](https://pepy.tech/project/screenpro2)
[![Downloads](https://static.pepy.tech/badge/screenpro2/month)](https://pepy.tech/project/screenpro2)
[![CodeQL](https://github.com/ArcInstitute/ScreenPro2/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/ArcInstitute/ScreenPro2/actions/workflows/github-code-scanning/codeql)
# ScreenPro2

## Introduction

### TL;DR

[**ReadTheDocs**](https://screenpro2.readthedocs.io) |
[**PyPI**](https://pypi.org/project/ScreenPro2)

ScreenPro2 enables perform flexible analysis on high-content CRISPR screening datasets. It has functionalities to process data from diverse CRISPR screen platforms and is designed to be modular to enable easy extension to custom CRISPR screen platforms or other commonly used platforms in addition to the ones currently implemented.

___
<details>
  <summary>Background</summary>

  Functional genomics field is evolving rapidly and many more CRISPR screen platforms are now developed. Therefore, 
  it's important to have a standardized workflow to analyze the data from these screens. ScreenPro2 is provided to 
  enable researchers to easily process and analyze data from CRISPR screens. Currently, you need to have a basic background in programming (especially Python) to use ScreenPro2.

  ScreenPro2 is conceptually similar to the [**ScreenProcessing**](https://github.com/mhorlbeck/ScreenProcessing) pipeline but **ScreenPro2** is designed to be more modular, flexible, and extensible. Common CRISPR screen methods that we have implemented here are illustrated in a recent review paper:

  > From: [A new era in functional genomics screens](https://www.nature.com/articles/s41576-021-00409-w)

  > Fig. 1: Common types of CRISPR screening modalities indicating advances in CRISPR methods.

  > <img width="1000" alt="image" src="https://github.com/GilbertLabUCSF/ScreenPro2/assets/53412130/a39400ad-b24f-4859-b6e7-b4d5f269119c">

</details>

___

<details>
  <summary>Benchmarking</summary>
  <div style="margin-top: 20px;"></div>
  
  Benchmarking ScreenPro2 with other CRISPR screen analysis tools

  ### More thoughtful NGS read trimming recovers more sgRNA counts

  ### ScreenPro2 statistical analysis is more accurate than ScreenProcessing

  ### ScreenPro2 is more flexible than ScreenProcessing

  Not only does ScreenPro2 have more features than ScreenProcessing, but it is also more flexible. ScreenPro2 can process data from diverse CRISPR screen platforms and is designed to be modular to enable easy extension to custom CRISPR screen platforms or other commonly used platforms in addition to the ones currently implemented.

  ### ScreenPro2 is faster than ScreenProcessing

  Last but not least, ScreenPro2 runs faster than ScreenProcessing (thanks to [biobear](https://github.com/wheretrue/biobear)) for processing FASTQ files.

</details>

___

## Installation
ScreenPro2 is available on [PyPI](https://pypi.org/project/ScreenPro2/) and can be installed with pip:
```bash
pip install ScreenPro2
```
___
For the latest version (development version) install from GitHub:
```bash
pip install git+https://github.com/ArcInstitute/ScreenPro2.git
```

## Usage

### Command Line Interface (CLI)
ScreenPro2 has a built-in command line interface (CLI). You can access the CLI by running the following command in your terminal:

```bash
screenpro --help
```

### Python Package Usage
First, import the ScreenPro2 package:

```python
import screenpro as scp
```

## Analysis Workflow

Data analysis for CRISPR screens with NGS readouts can be broken down into three main steps:

- [Step 1: FASTQ processing](#step-1-fastq-processing)
- [Step 2: Phenotype calculation](#step-2-phenotype-calculation)
- [Step 3: Data visualization](#step-3-data-visualization)

### Step 1: FASTQ processing

ScreenPro2 has a built-in command line interface (CLI) to process FASTQ files and generate counts.

```bash
screenpro guidecounter --help
```

A draft code to process FASTQ files and generate counts for [CRISPRa/i-single-sgRNA-screens](#dcas9-crisprai-single-sgrna-screens) dataset:

```bash
screenpro guidecounter
  --cas-type dCas9
  --single-guide-design
  -l <path-to-CRISPR-library-table>
  -p <path-to-fastq-directory>
  -s <sample-id-1>,<sample2-id>       # comma-separated list of sample ids, i.e. `<sample_id>.fastq.gz` for single sgRNA screens
  -o <output-directory>
  --write-count-matrix
```

A draft code to process FASTQ files and generate counts for [CRISPRa/i-dual-sgRNA-screens](#dcas9-crisprai-dual-sgrna-screens) dataset:
  
```bash
screenpro guidecounter
  --cas-type dCas9
  --dual-guide-design
  -l <path-to-CRISPR-library-table>
  -p <path-to-fastq-directory>
  -s <sample-id-1>,<sample2-id>       # comma-separated list of sample ids, i.e. `<sample_id>_R[1,2].fastq.gz` for dual sgRNA screens
  -o <output-directory>
  --write-count-matrix
```

___

In addition to the CLI, ScreenPro2 has a built-in method to process FASTQ files and generate counts in Python.
This method is implemented in the `ngs` module and relvent submodules. 
A minor novelty here has enabled processing single, dual, or multiple sgRNA 
CRISPR screens. Also, this approach can retain recombination events which can
occur in dual or higher order sgRNA CRISPR screens.

Currently, `GuideCounter` class from the `ngs` module can process FASTQ files and generate counts for standard 
CRISPR screens with [single](#dcas9-crisprai-single-sgrna-screens) or [dual](#dcas9-crisprai-dual-sgrna-screens) 
guide design. 

Here is a draft code to process FASTQ files and generate counts for an experiment with [CRISPRa/i-dual-sgRNA-screens](#dcas9-crisprai-dual-sgrna-screens):

```python
# Initialize the GuideCounter object
counter = scp.GuideCounter(cas_type = 'cas9', library_type = 'single_guide_design')

# Load the reference library
counter.load_library("<path-to-CRISPR-library-table>", sep = '\t', verbose = True, index_col=None)

# Define the samples
samples = [] 
## `samples` is a list of sample ids in the experiment. 
## Each sample id should match the sample name in the FASTQ files, i.e. <sample_id>.fastq.gz

# Process the FASTQ files and generate counts
counter.get_counts_matrix(
    fastq_dir = '<path-to-fastq-directory>',
    samples = samples,
    verbose = True
)
```

Here is a draft code to process FASTQ files and generate counts for an experiment with [CRISPRa/i-dual-sgRNA-screens](#crispri-dual-sgrna-screens):


```python
# Initialize the Counter object
counter = scp.GuideCounter(cas_type = 'dCas9', library_type = 'dual_guide_design')

# Load the reference library
counter.load_library("<path-to-CRISPR-library-table>", sep = '\t', verbose = True, index_col=None)

# Define the samples
samples = []
## `samples` is a list of sample ids in the experiment.
## Each sample id should match the sample name in the FASTQ files, i.e. <sample_id>_R[1,2].fastq.gz

# Process the FASTQ files and generate counts
counter.get_counts_matrix(
    fastq_dir = '<path-to-fastq-directory>',
    samples = samples,
    verbose = True
)
```

After this, you have `.counts_mat` calculated in the `GuideCounter` object.

___

To proceed, you need to create an `AnnData` object from the counts matrix and metadata. You can use the following code to create an `AnnData` object:

```python
adata = counter.build_counts_anndata()
```

### Step 2: Phenotype calculation

Once you have the counts, you can use ScreenPro2 `phenoscore` and `phenostats` modules to calculate the phenotype scores and statistics between screen arms.

#### Load Data
First, load your data into an `AnnData` object (see [anndata](https://anndata.readthedocs.io/en/latest/index.html) for more information).

The `AnnData` object must have the following contents:
- `adata.X` – counts matrix (samples x targets) where each value represents the sequencing count from NGS data.
- `adata.obs` – a pandas dataframe of samples metadata including "condition" and "replicate" columns.
  - "condition": the condition for each sample in the experiment.
  - "replicate": the replicate number for each sample in the experiment.
- `adata.var` – a pandas dataframe of targets in sgRNA library including "target" and "targetType" columns.
  - "target": the target for each entry in reference sgRNA library. For single sgRNA libraries, this column can be 
    used to store gene names. For dual or multiple targeting sgRNA libraries, this column can be used to store gene pairs
    or any other relevant information about the target.
  - "targetType": the type of target for each entry in reference sgRNA library. Note that this column is used to 
    distinguish between different types of sgRNAs in the library and negative control sgRNAs can be defined as `"targetType" == "negative_control"`.
    This is important for the phenotype calculation step.


ScreenPro2 has a built-in class for different types of CRISPR screen assays. Currently, there is a class called `PooledScreens` 
that can be used to process data from pooled CRISPR screens. To create a `PooledScreens` object from an `AnnData` object, 
you can use the following example code:

```python
import pandas as pd
import anndata as ad
from screenpro.assays import PooledScreens

adata = ad.AnnData(
    X   = counts_df, # pandas dataframe of counts (samples x targets)
    obs = meta_df,   # pandas dataframe of samples metadata including "condition" and "replicate" columns
    var = target_df  # pandas dataframe of targets metadata including "target" and "targetType" columns
)

screen = PooledScreens(adata)
```

<img width="600" alt="image" src="https://github.com/ArcInstitute/ScreenPro2/assets/53412130/bb38d119-8f24-44fa-98ab-7ef4457ef8d2">

#### Perform Screen Processing Analysis
Once the screen object is created, you can use several available workflows to calculate the phenotype scores and statisitics by comparing each entry in reference sgRNA library between screen arms. Then, these scores and statistics are used to nominate hits.

##### Drug Screen Workflow: calculate `gamma`, `rho`, and `tau` scores
`.calculateDrugScreen` method can be used to calculate the enrichment of each gene between screen arms for a drug 
screen experiment. This method calculates `gamma`, `rho`, and `tau` scores for each gene and adds them to the 
`.phenotypes` attribute of the `PooledScreens` object.

Here is an example for running the workflow on a [CRISPRi-dual-sgRNA-screens](#dcas9-crisprai-dual-sgrna-screens) dataset:

```python
# Run the ScreenPro2 workflow for CRISPRi-dual-sgRNA-screens
screen.calculateDrugScreen(
  t0='T0',
  untreated='DMSO',  # replace with the label for untreated condition
  treated='Drug',    # replace with the label for treated condition
  score_level='compare_reps'
)
```
___
For example, in a Decitabine CRISPRi drug screen (see Figure 1B-C in [this bioRxiv paper](https://www.biorxiv.org/content/10.1101/2022.12.14.518457v2.full)), each phenotype score represents a comparison between different arms of the screen and `rho` scores shows the main drug phenotype as illustrated here:
<img width="800" alt="image" src="https://github.com/abearab/ScreenPro2/assets/53412130/b84b3e1f-e049-4da6-b63d-d4c72bc97cda">

##### Flow cytometry based screen workflow: calculate phenotype score to compare high and low bins
`.calculateFlowBasedScreen` method can be used to calculate the enrichment of each target between high bin vs. low bin 
of a flow cytometry-based screen experiment. This method calculates `PhenoScore` for each target and adds them to the 
`.phenotypes` attribute of the `ScreenPro` object.

```python
# Run the ScreenPro2 workflow for CRISPRi-dual-sgRNA-screens
screen.calculateFlowBasedScreen(
  low_bin='low_bin', high_bin='high_bin',
  score_level='compare_reps'
)
```

### Step 3: Data visualization

Once the phenotypes are calculated, you can extract and explore the results using the `.phenotypes` attribute of the `ScreenPro` object. Currently, there are very limited functionalities built-in to visualize the results, but we are working on adding more features to make it easier for users. However, you can easily extract the results and use other libraries like `seaborn` and `matplotlib` in Python or `ggplot2` in R to visualize the results.

___

## Supported CRISPR Screen Platforms
One of the main goals of ScreenPro2 is to make it easy to process data from commonly used CRISPR screen platforms.
Also, it is designed to be modular to enable easy extension to custom CRISPR screen platforms or other commonly used
platforms in addition to the ones currently implemented.

___
Currently, ScreenPro2 has easy-to-use workflows for the following CRISPR screen platforms:
### dCas9 CRISPRa/i single-sgRNA screens
[Horlbeck et al., _eLife_ (2016)](http://dx.doi.org/10.7554/eLife.19760)

Horlbeck et al. developed a CRISPR interference (CRISPRi) and CRISPR activation (CRISPRa) screening platform that uses a single sgRNA within a single plasmid and then there are up to 10 sgRNAs per gene. The multiple sgRNAs per gene can be used to perfrom statistical comparisons in guide-level or gene-level between screen arms. [ScreenProcessing](https://github.com/mhorlbeck/ScreenProcessing) has been developed to process data from this type of screen. We reimplemented the same workflow in ScreenPro2 and it has all the necessary tools to process data from this type of screen.

<!-- TODO: Add link to example / tutorial -->

### dCas9 CRISPRa/i dual-sgRNA screens
[Replogle et al., _eLife_ (2022)](https://elifesciences.org/articles/81856)

Replogle et al. developed a CRISPR interference (CRISPRi) and CRISPR activation (CRISPRa) screening platform that uses two sgRNAs per gene within a single plasmid, and it has been used to perform genome-scale CRISPRi screens. ScreenPro2 has all the necessary tools to process data from this type of screen.

<!-- TODO: Add link to example / tutorial -->

<!-- ### multiCas12a CRISPRi screens -->

## License
ScreenPro2 is licensed under the terms of the MIT license (see [LICENSE](LICENSE) for more information) and developed 
by Abolfazl (Abe) Arab ([@abearab](https://github.com/abearab)) as a Research Associate in the Gilbert lab at UCSF and Arc Institute.  

## Citation
If you use ScreenPro2 in your research, please cite the following paper.

  Coming soon...
