PhenoScore module
=======================

Log ratio of :math:`y` vs :math:`x`:

.. math:: \Delta=\log(\frac{\begin{bmatrix}{N_{y}}\end{bmatrix}_{(a,b)} + 1}{\begin{bmatrix}{N_{x}}\end{bmatrix}_{(a,b)} + 1})

- :math:`y \rightarrow` condition :math:`x` (e.g. treated samples)
- :math:`x \rightarrow` condition :math:`y` (e.g. :math:`t_{0}` samples)
- :math:`a \rightarrow` number of oligo constructs with sgRNAs targeting :math:`T`
- :math:`b \rightarrow` number of biological replicates, :math:`R` (e.g. 2 or 3)
- :math:`N_{x}` | :math:`N_{y} \rightarrow` read fqcounter normalized for sequencing depth in condition :math:`x` or :math:`y`


Here is a formula for V3 library with single oligo construct per gene (i.e. 2 sgRNA in one oligo targeting same gene).

Phenotype score for each :math:`T` comparing :math:`y` vs :math:`x`:

.. math::
    \text{PhenoScore}(T,x,y) =
    \left(
    \frac{
        \overline{\Delta_{(x,y)}} - \text{median}({\overline{\Delta_{(x_{ctrl},y_{ctrl})}}})
    }{
        \sigma(\overline{\Delta_{(x_{ctrl},y_{ctrl})}})
    }
    \right)
    \times \frac{ 1 }{d_{growth}}

- :math:`\overline{\Delta(x,y)} \rightarrow` log ratio averaged across replicates
- :math:`T \rightarrow` oligo constructs with sgRNAs targeting :math:`T`
- :math:`\sigma(\text{...}) \rightarrow` standard deviation
- :math:`d_{growth} \rightarrow` growth rate...

-----------------------

.. automodule:: screenpro.phenoScore
   :members:
   :undoc-members:
   :show-inheritance:

