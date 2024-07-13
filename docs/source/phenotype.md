# Phenotype calculation modules

Log ratio of $y$ vs $x$:

$$\Delta =
\log(\frac
    {\begin{bmatrix}{N_{y}}\end{bmatrix}_{(a,b)} + 1}
    {\begin{bmatrix}{N_{x}}\end{bmatrix}_{(a,b)} + 1}
)$$

-   $y \rightarrow$ condition $x$ (e.g. treated samples)
-   $x \rightarrow$ condition $y$ (e.g. $t_{0}$ samples)
-   $a \rightarrow$ number of library elements with sgRNAs targeting $T$
-   $b \rightarrow$ number of biological replicates, $R$ (e.g. 2 or 3)
-   $N_{x}$ \| $N_{y} \rightarrow$ read counts normalized for sequencing
    depth in condition $x$ or $y$

Here is a formula for V3 library with single library element per gene
(i.e. dual sgRNAs in one construct targeting same gene).

Phenotype score for each $T$ comparing $y$ vs $x$:

$$\text{PhenoScore}(T,x,y) =
\left(
\frac{
\overline{\Delta_{(x,y)}}
}{
\text{median}( {\overline{\Delta_{(x_{ctrl},y_{ctrl})}}} )
}
\right)
\times \frac{ 1 }{d_{growth}}$$

-   $\overline{\Delta(x,y)} \rightarrow$ log ratio averaged across
    replicates
-   $T \rightarrow$ library elements with sgRNAs targeting $T$
-   $d_{growth} \rightarrow$ growth factor to normalize the phenotype
    score.

Statistical test comparing $y$ vs $x$ per each target, $T$:

$$\text{p-value}(T,x,y) = \text{t-test} \left(
\begin{bmatrix}{N_{x}}\end{bmatrix}_{(a,b)},
\begin{bmatrix}{N_{y}}\end{bmatrix}_{(a,b)}
\right)$$

(see this wikipedia page: [Dependent t-test for paired
samples](https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples))

(see the link to the implemented tool: [ttest_rel, a scipy
module](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html))

> This is a test for the null hypothesis that two related or repeated
> samples have identical average (expected) values.

___

```{eval-rst}  
.. automodule:: screenpro.phenoscore
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: screenpro.phenostats
   :members:
   :undoc-members:
   :show-inheritance:
```