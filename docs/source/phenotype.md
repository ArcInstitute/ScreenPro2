# Phenotype calculation modules

## Phenotype score calculation

Log ratio of {math}`y` vs {math}`x`:

```{math}
\Delta =
\log(\frac
    {\begin{bmatrix}{N_{y}}\end{bmatrix}_{(a,b)}}
    {\begin{bmatrix}{N_{x}}\end{bmatrix}_{(a,b)}}
)
```

-   {math}`y \rightarrow` condition {math}`x` (e.g. treated samples)
-   {math}`x \rightarrow` condition {math}`y` (e.g. {math}`t_{0}` samples, or untreated samples)
-   {math}`a \rightarrow` number of library elements with sgRNAs targeting {math}`T`
-   {math}`b \rightarrow` number of biological replicates, {math}`R` (e.g. 2 or 3)
-   {math}`N_{x}` \| {math}`N_{y} \rightarrow` read counts normalized for sequencing
    depth in condition {math}`x` or {math}`y`

___

Here is a formula for V3 library with single library element per gene
(i.e. dual sgRNAs in one construct targeting same gene).

Phenotype score for each {math}`T` comparing {math}`y` vs {math}`x`:

```{math}
\text{PhenoScore}(T,x,y) =
\left(
\frac{
\overline{\Delta_{(x,y)}}
}{
\text{median}( {\overline{\Delta_{(x_{ctrl},y_{ctrl})}}} )
}
\right)
\times \frac{ 1 }{d_{growth}}
```

-   {math}`\overline{\Delta(x,y)} \rightarrow` log ratio averaged across
    replicates
-   {math}`T \rightarrow` library elements with sgRNAs targeting {math}`T`
-   {math}`d_{growth} \rightarrow` growth factor to normalize the phenotype
    score.

## Phenotype statistics calculation

Statistical test comparing {math}`y` vs {math}`x` per each target, {math}`T`:

```{math}
\text{p-value}(T,x,y) = \text{t-test} \left(
\begin{bmatrix}{N_{x}}\end{bmatrix}_{(a,b)},
\begin{bmatrix}{N_{y}}\end{bmatrix}_{(a,b)}
\right)
```

(see this wikipedia page: [Dependent t-test for paired
samples](https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples))

(see the link to the implemented tool: [ttest_rel, a scipy
module](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html))

> This is a test for the null hypothesis that two related or repeated
> samples have identical average (expected) values).

## Combined score calculation

```{math}
\text{combined score} = \left( \dfrac{T_{\text{phenotype score}}}{\sigma{\text{(negative controls)}}} \right) \times -\log_{10}(\text{pvalue})
```

___

```{eval-rst}
.. automodule:: screenpro.phenoscore
   :members:
   :show-inheritance:

```

## Other related modules and functions

```{eval-rst}

.. automodule:: screenpro.phenoscore.phenostat
   :members:
   :show-inheritance:

.. automodule:: screenpro.phenoscore.delta
    :members:
    :show-inheritance:

.. automodule:: screenpro.phenoscore.deseq
    :members:
    :show-inheritance:

.. automodule:: screenpro.phenoscore.annotate
    :members:
    :show-inheritance:

```
