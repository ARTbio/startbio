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

Thus, the main benefit from using Limma, DESeq or edgeR packages is this modeling operation
which improves the accuracy of the statistical tests for differential expression.

### Multiple Testing

Each test for DE gives rise to a p-value, which is the probability of wrongly rejecting
hypothesis H0 which is, remember, that there is no difference in gene expression in view
of the observations of the number of reads.

For instance, when you read p = 0.05, this implies that the gene is differentially
expressed, with the probability that this conclusion is false being < 0.05. This is the
type I error.

However, transcriptome analysis implies thousands of tests. It happens that these tests
also follow a statistical law ! Even if a given test returns a p-value < 0.05, there is,
in addition a probability that this test was wrong !

Thus, when you perform thousands of test, you know that a proportion of these tests will
return wrong p-values.

The adjustment of the p-values seeds from this situation: you **must** correct your
p-values for multiple testing, because in the context of dozen of thousands tests,
p-values are poor indicators and does not allow to control the False Discovery Rate (type
I error)

Several methods exist for this correction. The Bonferroni correction is popular and
relatively conservative, whereas the Benjamini and Hochberg correction, which controls a
priori the False Discovery Rate is considered as less stringent. We can also cite other
methods that are not widely used in Biology such as the  Bonferroni Step-down (Holm)
correction or the Westfall and Young Permutation.

### Normalization

Last, but certainly not the least, to test a read count variable for differential
expression, a Normalization operation must be performed, since different sequencing depth
lead to different estimation of gene expression !

This Normalisation is performed differently by the Limma, DESeq or edgeR packages, which
is responsible a significant part of the differences between the packages.

### R packages used in this companionship

We are going to use most popular R packages DESeq2 and edgeR and it will be interesting to
compare the results returned by both packages. We will also try to give a shot to Limma
which was a very popular packages for analysing microarray results. Interestingly, Limma
has evolved and incorpores several methods to adapt to the more recent NGS RNAseq results.

### References

- **DESeq2**: Anders and Huber, Genome Biology 2010, 11:R106
  [DOI](https://doi.org/10.1186/gb-2010-11-10-r106)
- **edgeR**: Robinson, McCarthy and Smyth,
  Bioinformatics 2010, 26 p 139 [DOI](https://doi.org/10.1093/bioinformatics/btp616)
- **Limma**: Ritchie, Phipson, Wu, Hu, Law, Shi, et al., Nucleic Acids Res. 2015;43: e47.
  [DOI](https://doi.org/10.1093/nar/gkv007)


---
