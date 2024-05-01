[![website](https://img.shields.io/badge/website-live-brightgreen)](https://arcinstitute.org/tools/screenpro2)
![image](https://img.shields.io/pypi/v/screenpro2.svg) [![Documentation Status](https://readthedocs.org/projects/screenpro2/badge/?version=latest)](https://screenpro2.readthedocs.io/en/latest/?version=latest) ![Downloads](https://static.pepy.tech/badge/screenpro2)
[![Downloads](https://static.pepy.tech/badge/screenpro2/month)](https://pepy.tech/project/screenpro2)
[![CodeQL](https://github.com/ArcInstitute/ScreenPro2/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/ArcInstitute/ScreenPro2/actions/workflows/github-code-scanning/codeql)
# ScreenPro2


The complete docs are available at [screenpro2.rtfd.io](https://screenpro2.readthedocs.io).

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  * [Step 1: FASTQ to counts](#step-1-fastq-to-counts)
  * [Step 2: Phenotype calculation](#step-2-phenotype-calculation)
    + [Load Data](#load-data)
    + [Perform Screen Processing Analysis](#perform-screen-processing-analysis)
      - [Drug Screen Workflow: calculate `gamma`, `rho`, and `tau` scores](#drug-screen-workflow-calculate-gamma-rho-and-tau-scores)
      - [Flow cytometry based screen workflow: calculate phenotype score to compare high and low bins](#flow-cytometry-based-screen-workflow-calculate-phenotype-score-to-compare-high-and-low-bins)
  * [Step 3: Explore results and QC reports](#step-3-explore-results-and-qc-reports)

- [Supported CRISPR Screen Platforms](#supported-crispr-screen-platforms)
  * [dCas9 CRISPRa/i single-sgRNA screens](#dcas9-crispra/i-single-sgrna-screens)
  * [dCas9 CRISPRa/i dual-sgRNA screens](#dcas9-crispra/i-dual-sgrna-screens)
  <!-- * [multiCas12a CRISPRi screens](#multicas12a-crispri-screens) -->
- [License](#license)
- [Citation](#citation)

## Introduction
Functional genomics field is evolving rapidly and many more CRISPR screen platforms are now developed. Therefore, 
it's important to have a standardized workflow to analyze the data from these screens. ScreenPro2 is provided to 
enable researchers to easily process and analyze data from CRISPR screens. Currently, you need to have a basic background in programming (especially Python) to use ScreenPro2.

ScreenPro2 is conceptually similar to the [**ScreenProcessing**](https://github.com/mhorlbeck/ScreenProcessing) pipeline but **ScreenPro2** is designed to be more modular, flexible, and extensible. Common CRISPR screen methods that we have implemented here are illustrated in a recent review paper:

> From: [A new era in functional genomics screens](https://www.nature.com/articles/s41576-021-00409-w)

> Fig. 1: Common types of CRISPR screening modalities indicating advances in CRISPR methods.

> <img width="1000" alt="image" src="https://github.com/GilbertLabUCSF/ScreenPro2/assets/53412130/a39400ad-b24f-4859-b6e7-b4d5f269119c">

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
Data analysis for CRISPR screens with NGS readouts can be broken down into three main steps:

- [Step 1: FASTQ to counts](#step-1-fastq-to-counts)
- [Step 2: Phenotype calculation](#step-2-phenotype-calculation)
- [Step 3: Explore results and QC reports](#step-3-explore-results-and-qc-reports)

### Step 1: FASTQ to counts

Since version 0.2.7, ScreenPro2 has a built-in method to process FASTQ files and generate counts. This method is implemented in the `ngs` module 
and relvent submodules. A minor novelty here has enabled processing single, dual, or multiple sgRNA CRISPR screens. Also, this approach can retain 
recombination events which can occur in dual or higher order sgRNA CRISPR screens.

There is no example code for this step yet, but a command line interface (CLI) will be available soon. 

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
    distinguish between different types of sgRNAs in the library and negative control sgRNAs can be defined as `"targetType" == "negCtrl"`.
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
<img width="600" alt="image" src="https://github.com/abearab/ScreenPro2/assets/53412130/d1c8c3ad-3668-4390-8b1d-bf72b591a927">

#### Perform Screen Processing Analysis
Once the screen object is created, you can use several available workflows to calculate the phenotype scores and statisitics by comparing each entry in reference sgRNA library between screen arms. Then, these scores and statistics are used to nominate hits.

##### Drug Screen Workflow: calculate `gamma`, `rho`, and `tau` scores
`.calculateDrugScreen` method can be used to calculate the enrichment of each gene between screen arms for a drug 
screen experiment. This method calculates `gamma`, `rho`, and `tau` scores for each gene and adds them to the 
`.phenotypes` attribute of the `PooledScreens` object.

Here is an example for running the workflow on a [CRISPRi-dual-sgRNA-screens](#crispri-dual-sgrna-screens) dataset:

```python
# Run the ScreenPro2 workflow for CRISPRi-dual-sgRNA-screens
screen.calculateDrugScreen(
  t0='T0',
  untreated='DMSO',  # replace with the label for untreated condition
  treated='Drug',    # replace with the label for treated condition
  db_untreated=1,    # replace with doubling rate of untreated condition
  db_treated=1,      # replace with doubling rate of treated condition
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

### Step 3: Explore results and QC reports

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

Horlbeck et al. developed a CRISPR interference (CRISPRi) and CRISPR activation (CRISPRa) screening platform that uses a single sgRNA within a single plasmid and then there are up to 10 sgRNAs per gene. The multiple sgRNAs per gene can be used to perfrom statistical comparisons in guide-level or gene-level between screen arms. [ScreenProcessing](https://github.com/mhorlbeck/ScreenProcessing) has been developed to process data from this type of screen. We reimplemented the same workflow in ScreenPro2 and it has all the necessary tools to process data from this type of screen. An automated workflow / pipeline will be available soon.

### dCas9 CRISPRa/i dual-sgRNA screens
[Replogle et al., _eLife_ (2022)](https://elifesciences.org/articles/81856)

Replogle et al. developed a CRISPR interference (CRISPRi) and CRISPR activation (CRISPRa) screening platform that uses two sgRNAs per gene within a single plasmid, and it has been used to perform genome-scale CRISPRi screens. ScreenPro2 has all the necessary tools to process data from this type of screen. An automated workflow / pipeline will be available soon.

<!-- ### multiCas12a CRISPRi screens -->

## License
ScreenPro2 is licensed under the terms of the MIT license (see [LICENSE](LICENSE) for more information) and developed 
by Abolfazl (Abe) Arab ([@abearab](https://github.com/abearab)) as a Research Associate in the Gilbert lab at UCSF and Arc Institute.  

## Citation
If you use ScreenPro2 in your research, please cite the following paper.

  Coming soon...
