# NGS screen processing module

`GuideCounter` class is a wrapper to run the functions for a
CRISPR screen experiment.

This module contains a set of python functions to process and analyze
NGS files from CRISPR screens. Based on the type of CRISPR-Cas system
used for the screen, the functions are divided into two classes:
`Cas9` and `Cas12`.

------------------------------------------------------------------------
```{eval-rst}  
.. automodule:: screenpro.ngs
   :members:
   :show-inheritance:
```

### Cas9 CRISPR-Cas system (single or dual sgRNA libraries)
```{eval-rst}  
.. automodule:: screenpro.cas9
   :members:
   :show-inheritance:
```

### Cas12 CRISPR-Cas system (multiplexed crRNA libraries)
```{eval-rst}  
.. automodule:: screenpro.cas12
   :members:
   :show-inheritance:
```
