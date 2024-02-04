# ![](images/lamp.png){align="absbottom} Statistical analysis of Differential Gene Expression

## The basic ideas of Differential expression analysis

### Gene Counts are Observations of Variables
We consider that each gene is a ==variable==. Thus DE analysis is dealing with
==multiple== variables (10 to several tens of thousands).

Each read count is an ==observation== of these variables. Thus, for instance, if your
experiment is based on biological triplicates in two conditions, you have three
observation of your gene variables under 2 different conditions (6 observations in total).

### Testing
When we say that we are testing for differential expression, we mean that we are
performing multiple statistical tests, ==one for each== gene variable. These tests are
well established mathematical treatments such as the Student Test (t-test), the
Mann-Whitney U test, the Wilcoxon test or the exact Fisher test, for the most used tests.
However, note that not all these tests are suitable for discrete count variables.

### Conditions for testing
As you probably learned during your university studies each of these tests have underlying
assumptions. The parametric tests require the a-priori knowledge of parameters (mean,
variance, etc.), or that the distribution of the tested variable follows a specific law
(normality, continuity, etc.). For instance the parametric student test requires that the
means of observation is normally distributed, which is the case if the number of
observation is > 30 (Central Limit Theorem) or if the variable is normally distributed
(which is not the case for a read count variable !)

### Read counts from NGS sequencing are not normally distributed

In contrast to the intensity of a probe in a microarray, a read count variable is does not
follow and gaussian (normal) distribution ! Statistician showed that read counts
variables, which are discrete variable, follow a generalized Poisson distribution, which
can be approximated by an inverse binomial distribution (also referred to as negative
binomial, NB) when the number of observations is low.

From this, it comes that the main tools for NGS DE analysis model read count variables with
a Poisson (Limma) or a NB law, and use specific tests for differential expression.

Note that these tests are ==parametric test==, which implies that the Mean and the Variance
must be approximated before the tests. In the case of the Negative Binomial distribution,
Mean == Variance.

### Shared information between read count variables

Although differential Expression analysis is based on the assumption that gene expression
variables are **independent**, it happens that these variables share information which can be
used for better modeling of test parameters *for each test*.

Thus, the main benefit from using Limma, DESeq or EdgeR packages is this modeling operation
which improves the accuracy of the statistical tests for differential expression.

### Normalization

Last, but certainly not the least, to test a read count variable for differential expression,
a Normalization operation must be performed, since different sequencing depth lead to different
estimation of gene expression !

This Normalisation is performed differently by the Limma, DESeq or EdgeR packages, which is
responsible a part of the differences between the package.

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
