![image](https://img.shields.io/pypi/v/screenpro2.svg) [![Documentation Status](https://readthedocs.org/projects/screenpro2/badge/?version=latest)](https://screenpro2.readthedocs.io/en/latest/?version=latest) ![Downloads](https://static.pepy.tech/badge/screenpro2)


ScreenPro2
==========

<!-- The docs are available at [https://screenpro2.readthedocs.io](https://screenpro2.readthedocs.io) -->

This package is conceptually similar to the [**ScreenProcessing**](https://github.com/mhorlbeck/ScreenProcessing) 
pipeline but **ScreenPro2** is designed to be more modular, flexible, and extensible as the field of Functional 
Genomics evolves and newer CRISPR screen platforms are developed.

For more information about the statistical methods used in ScreenPro2, please refer to the detailed documentation
about [PhenoScore](https://screenpro2.readthedocs.io/en/latest/PhenoScore.html) module.

Note that ScreenPro2 starts with a counts matrix of oligo counts (samples x oligos) so you will need to process your 
raw sequencing data into a counts matrix before using ScreenPro2.

## Installation
ScreenPro2 is available on PyPI and can be installed with pip:
```bash
pip install ScreenPro2
```
___
For the latest version (development version) install from GitHub:
```bash
pip install git+https://github.com/abearab/ScreenPro2.git
```

## Usage

### Load Data
First, load your data into an `AnnData` object (see [anndata](https://anndata.readthedocs.io/en/latest/index.html) for 
more information).

The `AnnData` object should have the following structure:
- `adata.X` should be a pandas dataframe of counts (samples x oligos)
- `adata.obs` should be a pandas dataframe of sample metadata including "condition" and "replicate" columns
- `adata.var` should be a pandas dataframe of oligo metadata including "target" and "targetType" columns
  - "target" column should be the gene name or other identifier for the reference oligo
  - "targetType" column should be the type of reference oligo. Currently, negative control oligos should have
    `"targetType" == "negCtrl"`

Then you need create a `ScreenPro` object. Here is an example code making a `ScreenPro` object from an `AnnData` object:

```python
import pandas as pd
import anndata as ad
import screenpro as scp

adata = ad.AnnData(
    X   = counts_df, # pandas dataframe of counts (samples x oligos)
    obs = meta_df,   # pandas dataframe of sample metadata including "condition" and "replicate" columns
    var = target_df  # pandas dataframe of oligo metadata including "target" and "targetType" columns
)

screen = scp.ScreenPro(adata)
```

### Perform Screen Processing Analysis
Once the `ScreenPro` object is created, you can use several available workflows to calculate the enrichment of each oligo 
between screen arms. 

#### Drug Screen Workflow: calculate `gamma`, `rho`, and `tau` scores
`.calculateDrugScreen` method can be used to calculate the enrichment of each gene between screen arms for a drug 
screen experiment. This method calculates `gamma`, `rho`, and `tau` scores for each gene and adds them to the 
`.pheno

Here is an example for running the workflow on a [CRISPRi-dual-sgRNA-screens](#crispri-dual-sgrna-screens) dataset:

```python
# Run the ScreenPro2 workflow for CRISPRi-dual-sgRNA-screens
screen.calculateDrugScreen(
  t0='T0', untreated='DMSO', treated='Drug', 
  growth_rate=1, # can be replaced by the growth values (population doublings/doubling differences)
  score_level='compare_reps'
)
```

## Supported CRISPR Screen Platforms
One of the main goals of ScreenPro2 is to make it easy to process data from commonly used CRISPR screen platforms.
Also, it is designed to be modular to enable easy extension to custom CRISPR screen platforms or other commonly used
platforms in addition to the ones currently implemented.

___
Currently, ScreenPro2 has easy-to-use workflows for the following CRISPR screen platforms:
#### CRISPRi-dual-sgRNA-screens
[Replogle et al., _eLife_ (2022)](https://elifesciences.org/articles/81856)

Replogle et al. developed a CRISPRi screening platform that uses two sgRNAs per gene within a single plasmid and it has
been used to perform genome-scale CRISPRi screens. If you follow the codes in the provided [GitHub](https://github.com/josephreplogle/CRISPRi-dual-sgRNA-screens) repository, you 
will end up with oligo counts and once you make `ScreenPro` object, you can use the ScreenPro2 workflow for this
platform to calculate the enrichment of each gene between screen arms.

## License
ScreenPro2 is licensed under the terms of the MIT license (see [LICENSE](LICENSE) for more information) and developed 
by Abolfazl (Abe) Arab ([@abearab](https://github.com/abearab)), a Research Associate in the Gilbert lab at UCSF and Arc Institute.  
