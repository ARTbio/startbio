# Single Cell RNAseq via Seurat R package

A classical single-cell RNA seq analysis consists in identifying populations
of cells and the associated marker genes. It can also look for the effect of
a treatment or a condition using differential analysis methods.

The dataset used for the analysis is composed of peripheral blood mononuclear
cells (PBMC). It is available on the 10X and Seurat website which can be found
on the [Seurat tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

## Libraries

Here is R packages needed for the analysis.

``` r
## Import des packages
library(Seurat)
library(ggraph)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(biomaRt)
library(plyr)
library(dplyr)
library(magrittr)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(knitr)
library(rmarkdown)
library(msigdbr)
library(vroom)
```

## Functions

In a R script, good programming practices encourage that at the beginning of
each script we put the functions that we develop after the call of the
libraries. In the next *chunk* of code, we will find the functions that will
be useful during the analysis.

Here it is a function that will allow to annotate the plots in order to have
the name of the gene in addition to its identifier.

``` r
add_title_gene_name <- function(plot,
                                   gene_format,
                                   from = "ensembl_gene_id",
                                   to = "external_gene_name"){
  ## Add gene name as title and leave gene ID as subtitle of a plot
  ### Inputs
  ## - plot (data) : ggplot to modify
  ## - gene_format (data) : dataframe that contains at least the type of gene id present in the plot (from) and the gene label type to use instead (to)
  ## - from (chr) : label gene type present in the plot (must be the column name of the annotated dataframe)
  ## - to (chr) : label gene type to use
  ### Output
  ## - Plot with a new title + subtitle

  ##check if a modification is possible
  ## check if one of the column data is a gene of the "from" column of gene_format
  test_matching <- colnames(plot$data) %in% gene_format[, from]
  if(sum(test_matching)){
    gene_to_rename <- colnames(plot$data)[test_matching]
    plot <- plot +
      ggtitle(gene_format[gene_format[, from] == gene_to_rename, to],
              gene_to_rename)
  }else{
    stop(paste("No matching between plot metadata and", from, "column of `gene_format` dataframe.\n", "Please check your parameters."))
  }
  return(plot)
}
```
