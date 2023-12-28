# Analysis of the differential gene expression using `DESeq2`

![](images/lamp.png)

----

DESeq2 is a great tool for Differential Gene Expression (DGE) analysis.
It takes read counts and combines them into a table (with genes in the rows and samples in the columns).
Importantly, it applies size factor normalization by:

- Computing for each gene the geometric mean of read counts across all samples
- Dividing every gene count by the geometric mean accross samples
- Using the median of these ratios as a sample’s size factor for normalization

Multiple factors with several levels can then be incorporated in the analysis.
After normalization we can compare the response of the expression of any gene to
the presence of different levels of a factor in a statistically reliable way.

In our example, we have samples with two varying factors that can contribute to
differences in gene expression:

- Treatment (either treated or untreated)
- Sequencing type (paired-end or single-end)

Here, treatment is the primary factor that we are interested in.

The sequencing type is further information we know about the data that might affect
the analysis. Multi-factor analysis allows us to assess the effect of the treatment,
while taking the sequencing type into account too.

```
We recommend that you add as many factors as you think may affect gene expression in
your experiment. It can be the sequencing type like here, but it can also be the
manipulation (if different persons are involved in the library preparation),
other batch effects, etc…
```
