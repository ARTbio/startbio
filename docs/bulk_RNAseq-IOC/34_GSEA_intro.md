# Gene Set Enrichment Analysis

## Definition and Rationale Behind Gene Set Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) is a powerful computational method used in
bioinformatics to interpret gene expression data in the context of biological pathways,
processes, or sets of functionally related genes.

Unlike traditional methods that focus on individual genes, GSEA evaluates the **coordinated
expression changes** of predefined gene sets, providing a more holistic view of molecular
mechanisms underlying experimental conditions or phenotypes.

### Definition of GSEA
  - [x] GSEA assesses whether predefined sets of genes show statistically significant,
  concordant differences between two biological states (e.g., treatment vs. control,
  diseased vs. healthy).
  - [x] Rather than focusing on individual genes, GSEA operates on gene sets, which can
  represent pathways, molecular functions, cellular processes, or ==other biologically
  relevant groups of genes==. This last case actually represents the most common use of
  GSEA. Many groups of genes, sometimes also called "gene signatures" or "molecular
  signatures" are available in databases or published articles.
  - [x] It ranks all genes based on their expression changes between experimental
  conditions and then tests whether genes within a gene set tend to appear towards the top
  (or bottom) of the ranked list more than expected by chance.

### Rationale Behind GSEA
  - [x] **Biological Context**: GSEA acknowledges that genes rarely act in isolation but
  rather function in coordinated networks and pathways. Analyzing gene sets helps
  contextualize gene expression changes within the framework of biological processes.
  - [x] **Statistical Power**: By aggregating signals from groups of genes, GSEA enhances
  statistical power to detect subtle but coordinated changes in gene expression that might
  be missed when analyzing individual genes.
  - [x] **Reduction of Multiple Testing Burden**: GSEA reduces the multiple testing burden
  associated with examining thousands of individual genes by focusing on predefined gene
  sets. This reduces the risk of false positives and improves the reliability of results.
  - [x] **Interpretability**: GSEA provides interpretable results by associating gene
  expression changes with known biological pathways or processes, enabling researchers to
  generate hypotheses and gain insights into the underlying biology.
  - [x] **Robustness Across Platforms**: GSEA is platform-independent and can be applied
  to various types of gene expression data, including microarray and RNA sequencing
  (RNA-seq) data, making it widely applicable across different experimental settings and
  datasets.

### Key Features of GSEA
  - [x] **Enrichment Score**: Measures the degree to which a gene set is overrepresented
  at the top or bottom of the ranked gene list.
  - [x] **Normalized Enrichment Score (NES)**: Corrects for gene set size and data set
  size, facilitating comparison of results across different datasets.
  - [x] **False Discovery Rate (FDR)**: Estimates the proportion of false positive results
  among significant findings, controlling for multiple testing.

In summary, GSEA offers a systematic and biologically meaningful approach to analyze gene expression data, enabling researchers to uncover key molecular pathways and processes associated with different experimental conditions or phenotypes. Its ability to integrate complex genomic data with prior biological knowledge makes it a valuable tool in deciphering the mechanisms underlying biological phenomena and disease states.

## How to perform GSEA ?

### A video presentation by Katherine West (University of Glasgow) 
We strongly advise you to look at the excellent
[presentation by Katherine West](https://youtu.be/KY6SS4vRchY?si=cxbHjHdXdjc7uE-4){:target="_blank"}.
The aspects that have been presented above are all taken up and illustrated with graphics
in a very educational way.

### Practical focus: computation of Enrichment Score (ES)

The Enrichment Score (ES) is central in Gene Set Enrichment Analysis (GSEA) since it
quantifies the degree to which a gene set is overrepresented at the top or bottom of a
ranked list of genes based on their expression changes between two biological conditions.

The computation of the Enrichment Score involves several steps:

1. **Ranking Genes**: The first step is to rank all genes in the dataset based on a metric
that reflects their differential expression between the two biological conditions. This
metric is most often fold change, but could be t-statistic, or any other relevant
statistical measure.

2. **Cumulative Sum Calculation**: The Enrichment Score is calculated by walking down the
      ranked list of genes, accumulating a running sum statistic. At each step, the running sum
      is increased when a gene belongs to the gene set being evaluated and decreased otherwise.
      The running sum captures the degree of enrichment of the gene set at that point in the
      ranked list.
      
      The way the running sum is increased when a gene belongs to the gene set being evaluated
      and decreased otherwise varies depending on the GSEA implementation (there are several).
      What you need to remember is that increment and decrement are never calculated
      symmetrically.
      
      A simple example of running sum calculation is to add the GSEA metric (eg fold change)
      when a gene belongs to the gene set being evaluated and to remove a fixed value that
      depends on the total number of genes in the ranked gene list (eg 1 / N). This fixed value
      is typically referred to as the "penalty" or "decay" factor.
      
      The rationale behind using a penalty or decay factor is to adjust the running sum to
      account for the fact that genes not belonging to the gene set can still contribute to the
      overall distribution of scores. This adjustment helps to prevent the running sum from
      being overly biased by the presence or absence of genes in the gene set.
      
3. **Peak Enrichment Score**: The Enrichment Score reaches its maximum (peak) value when
the cumulative sum reaches its maximum deviation from zero. This peak reflects the
enrichment of the gene set at a particular position in the ranked list.

4. **Normalization of Enrichment Score**: To make Enrichment Scores comparable across
different gene sets and datasets, the Enrichment Score is normalized. This normalization
accounts for differences in gene set size and dataset size. One common normalization
method is to divide the Enrichment Score by the mean enrichment score from permuted
datasets.

5. **Estimation of Significance**: The significance of the Enrichment Score is assessed
through permutation testing. This involves repeatedly permuting the gene labels to
generate a null distribution of Enrichment Scores (ie computing many NES from gene sets
randomly sampled from the total gene list). The observed Enrichment Score is then compared
to the null distribution to determine its statistical significance, typically reported as
a nominal p-value or false discovery rate (FDR).


Overall, the Enrichment Score provides a quantitative measure of the degree to which a
predefined gene set is enriched towards the top or bottom of a ranked list of genes,
indicating the collective expression behavior of genes within that set under different
experimental conditions. It enables the identification of biologically relevant gene sets
associated with specific phenotypes or experimental treatments in gene expression studies.

## The main resource for GSEA

* GSEA software: [https://www.gsea-msigdb.org/gsea/msigdb](https://www.gsea-msigdb.org/gsea/msigdb)
    * Provides a user-friendly platform for performing GSEA analysis.
    * Provides access to a large database of curated gene sets in various format, including
    the GMT format (.gmt files) which is the format that we are going to use in this IOC.