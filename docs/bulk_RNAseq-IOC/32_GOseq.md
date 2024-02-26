# Analysis of Gene Ontology enrichments using the Galaxy Tool `goseq`

As mentioned previously NGS read count data are inherently biased by the length of
transcripts.

How does this bias translate in DE tables ?

==Not== in biased Fold Changes ! Fold changes are count ratios for a given gene, so bias
is compensated for.

In contrast, genes with a high number of counts (whether due to a high level of expression
or a large transcript) tend to be ==over-**detected**== since the p-value returned by the
statistical test depends on the count base mean. This becomes a major issue in GO-based
Gene Set Analyses since gene sets are selected on the basis of their p-values/p-adjusted-values.

The `goseq` tool provides methods for performing GO analysis of RNA-seq data, taking
length bias into account.

## `goseq` analysis of the DESeq2 DE tables in the use-case PRJNA630433

### Gather needed data in a new history
- [x] Although we are going to focus on the DESeq2 tables, the approach would be the same
with other DE callers: feel free to test it latter on.
- [x] to perform the analysis, we will need a gene length list. This list was generated
[previously by the featureCounts tool](../21_FeatureCounts) and is therefore available in
the history `PRJNA630433 FeatureCounts Counting on HISAT2 bam alignments`.
  
  --> Copy one of the collection among the `Dc`, `Mo` or `Oc FeatureCounts Feature Length`
  collections in this history in a new history which you will name `PRJNA630433 GOseq
  analysis`. Note that Feature lenght datasets are all identical, we will just need to
  extract one for the goseq analysis.
  
- [x] Finally, copy the collection of DEseq2 tables from the history `PRJNA630433 DESeq2
  analysis` in the new history `PRJNA630433 FeatureCounts Counting on HISAT2 bam
  alignments`. This collection should be named `DESeq2 Results Tables`. 
  
### Prepare the Gene Set(s)

Since we are going to perform a GO-based Gene Set Analysis we first need to define the
Gene Sets ! We will keept the advantage of collections and treat in parallel the three
comparisons (Mo vs Dc, Oc vs Dc and Oc vs Mo).

#### Clean up tables from `NA`s
We first clean up the Tables by removing all lines that contain `NA` values using the Tool
`Select lines that match an expression`

This is an important step because it allows to reduce the "gene space". As mentioned in
introduction, choosing a relevant background list representing the genes ==not
differentially expressed== is crucial for accurate enrichment analysis. By removing the
lines containing NAs, we shrink our datasets and retains only genes that are accessible to
tests, ie, that are at least significantly ==expressed== in the considered tissue/experience.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Select lines that match an expression` settings"
    - Select lines from
        
        --> Click the collection icon and select `DESeq2 Results Tables`
    - that
        
        --> `Not matching` (be careful, ==not== matching)
    - the pattern
        
        --> `\tNA$` (we seach for lines that contain any tabulation mark, followed by `NA` and
        this `NA` is a word, ie, it is followed by a space or another tabulation mark)
    - Keep header line
        
        --> `Yes`
    - `Run Tool`

#### Add a boolean column to tag the gene Set
Now we are going to select our gene sets from DE table collection using a treatment which
evaluates a each line of the table and returns a boolean value in a new column. Genes with
a `True` value will been considered as part of the gene set, whereas genes with a `False`
value will be considered as part of the background list for the GSA.

Let's do this, using a stringent evaluation, ie we will tag genes for which both `p-adjust
< 0.001` and `|log2FC| > 2`. In your own analyses you may have to use less stringent
filtering. It depends mostly of your datasets and whether you get a sufficiently large enough
gene set after filtering. Considere that the gene set size should be between 1 and 10 % of
your background genes.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Compute on rows` settings"
    - Input file
        
        --> Click the collection icon and select `Select on collection 9`
    - Input has a header line with column names?
        
        --> `Yes`
    - Add expression
        
        --> `c7 < 0.001 and abs(c3) > 2`
    - Mode of the operation
    
        --> `Append`
    - The new column name
    
        --> `boolean`
    - `Run Tool`

#### Cut columns 1 and 8 to get the final gene sets and backgrounds

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Advanced Cut columns from a table` settings"
    - File to cut
        
        --> Click the collection icon and select `Compute on collection 10`
    - Operation
        
        --> `Keep`
    - Delimited by
        
        --> `Tab`
    - Cut by
    
        --> `fields`
    - List of Fields
    
        --> `boolean`
    - `Run Tool`

:warning: As this is the last step of the construction of the gene set lists, you should
rename for clarity the returned collection as `Gene Lists`.

### Extract a single gene length table dataset
A single gene length table can be extract from the collection which we initially copied from
the history `PRJNA630433 FeatureCounts Counting on HISAT2 bam alignments`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Extract dataset` settings"
    - Input List
        
        --> Click the collection icon and select `Mo featureCounts: Feature lengths`
    - How should a dataset be selected?
        
        --> `The first dataset`
        
        :warning: it would work equally well with `Select by element identifier` or
        `Select by index`
    - `Run Tool`

:warning: For clarity, also rename this single dataset as `Gene lengths`

### Perform goseq analysis

Our two inputs are now ready for goseq:
- The list of tagged genes. `True` for the gene set, `False` for the background genes
- The length of genes (aggregated size of exons) to correct for the length bias of detection.

goseq proposes 3 types of GSA methods

- [x] The **Wallenius method** approximates the true distribution of numbers of members
  of a category amongst DE genes by the Wallenius non-central hypergeometric distribution.
  This distribution assumes that within a category all genes have the same probability of
  being chosen. Therefore, this approximation works best when the range in probabilities
  obtained by the probability weighting function (the regression over gene lenght) is small.
  This is the method specifically developed in the goseq package as well as the one we are
  going to choose in this IOC.

- [x] The **Sampling method** uses random sampling to approximate the true distribution
  and uses it to calculate the p-values for over (and under) representation of categories.
  Although this is the most accurate method given a high enough value of sampling number,
  its use quickly becomes computationally prohibitive. It may sometimes be desirable to
  use random sampling to generate the null distribution for category membership. For
  example, to check consistency against results from the Wallenius approximation. This is
  easily accomplished by using the method option to additionally specify sampling and the
  number of samples to generate.

- [x] The **Hypergeometric Method** assumes ==there is no bias== in power to detect
  differential expression at all and calculates the p-values using a standard hypergeometric
  distribution (no length bias correction is performed). Useful if you wish to test the
  effect of length bias on your results. :warning: Hypergeometric method should NEVER be
  used for producing results for biological interpretation of RNA-seq data. Indeed, if
  length bias is truly not present in your data, goseq will produce a nearly flat PWF plot,
  no length bias correction will be applied to your data, and all methods will produce the
  same results.


!!! info "![](images/tool_small.png){width="25" align="absbottom"} `goseq tests` settings"
    - Differentially expressed genes file
        
        --> Click the collection icon and select `Gene Lists` (the renamed collection)
    - Gene lengths file
        
        --> `Gene lengths` (the renamed single dataset)
    - Gene categories
        
        Because we are working here with Mus musculus, we can rely on the tool to fetch
        directly the GO categories from a remote database. :warning: For *Nesseria gon* As
        well as *Apis mel* you will have to construct and use your own Gene categorie file.
        And... Yes, the section title `Gene categories` is not sufficiently explicit and
        should rather be `GO categories`
        
        --> `Get categories`
    - Select a genome to use
        
        --> `Mouse (mm10)`
    - Select Gene ID format
        
        --> `Ensembl Gene ID`
    - Select one or more categories
        
        --> Check only `GO: Biological Process` (for simplicity in this first run)
    - Method Options
        
        --> Use Wallenius method `Yes`
        
        --> Use Hypergeometric method `No`
        
        --> Sampling number `5000` :warning: We also run this method to compare it with
        the Wallenius method.
    - Advanced Options 
        
        --> Select a method for multiple hypothesis testing correction `Benjamini-Hochberg [FDR]`
        
        --> Count genes without any category? `Yes`
    - Output Options
        
        --> Output Top GO terms plot? `Yes`
        
        --> Produce diagnostic plots? `Yes`
        
        --> Extract the DE genes for the categories (GO/KEGG terms)? `Yes`
    - `Run Tool`

`goseq` generates a big table with the following columns for each GO term:

|Column                       | Description                                                                                           |
|-----------------------------|-------------------------------------------------------------------------------------------------------|
|category                     | GO category                                                                                           |
|over_rep_pval                |p-value for over representation of the term in the differentially expressed genes                      |
|under_rep_pval               |p-value for under representation of the term in the differentially expressed genes                     |
|numDEInCat                   |number of differentially expressed genes in this category                                              |
|numInCat                     |number of genes in this category                                                                       |
|term                         |detail of the term                                                                                     |
|ontology                     |MF (Molecular Function - molecular activities of gene products), CC (Cellular Component - where gene products are active), BP (Biological Process - pathways and larger processes made up of the activities of multiple gene products)|
|p.adjust.over_represented    |p-value for over representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure |
|p.adjust.under_represented   |p-value for over representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure |

To identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to use the adjusted p-value.

??? question "How many GO terms are **over**-represented at adjusted P value < 0.05?"
    
??? question "How many GO terms are **under**-represented at adjusted P value < 0.05?"
    

