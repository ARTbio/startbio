# Markers Annotation

Each cluster is associated with a list of marker genes that we now need to
annotate in order to biologically identify the cell clusters. We will start
by annotating the gene identifiers with Biomart and then use functional
enrichment methods to find the functions shared by the marker genes.

## Via Biomart

We will use the dataframe we generated at the beginning of the analysis
to allow us to add the gene name and a description for each gene ID set.

``` r
## Merge markers results with biomart annotation
pbmc_markers_annotated <- merge(x = pbmc_markers,         #First df to merge
                                y = annotated_hg19,       #Second df to merge
                                by.x = "gene",            #Column name of first df used for matching lines
                                by.y = "ensembl_gene_id", #Column name of second df used for matching lines
                                all.x = TRUE)             #Keep all lines from first df even if there is no match with second df
```

We will now remove all markers where the adjusted p-value is greater than
a 5% threshold and with the absolute value of the mean log(Fold Change)
less than 0.25 in order to obtain the list of markers with a significant
expression differential.

``` r
pbmc_markers_signif <- subset(pbmc_markers_annotated,
                              p_val_adj < 0.05 &
                                abs(avg_log2FC) >= 0.25)       #Filter dataframe based on p_val_adj column

## Number of significative DEG per cluster
table(pbmc_markers_signif$cluster)
```

    ##
    ##   0   1   2   3   4   5   6   7   8
    ## 321 518 166 337 173 436 384 127 233

``` r
## Sorting results by cluster and by average log2(Fold Change)
pbmc_markers_signif <- pbmc_markers_signif %>%                 #Rearrange df with dplyr package
  group_by(cluster) %>%                                        #Group df based on cluster column
  arrange(desc(avg_log2FC), .by_group = TRUE)                  #Sort lines by descending the column avg_log2FC and by group

## Most DE gene marker for each cluster
kable(top_n(x= pbmc_markers_signif, n = 3, wt = avg_log2FC))
```
??? abstract "First annotated markers for each cluster"
    | gene            | p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | external_gene_name | description                                                                                                  | gene_biotype   | chromosome_name |
    |:----------------|------:|-----------:|------:|------:|----------:|:--------|:-------------------|:-------------------------------------------------------------------------------------------------------------|:---------------|:----------------|
    | ENSG00000111716 |     0 |  0.7920556 | 0.929 | 0.592 |         0 | 0       | LDHB               | lactate dehydrogenase B \[Source:HGNC Symbol;Acc:6541\]                                                      | protein_coding | 12              |
    | ENSG00000145425 |     0 |  0.7576946 | 0.997 | 0.976 |         0 | 0       | RPS3A              | ribosomal protein S3A \[Source:HGNC Symbol;Acc:10421\]                                                       | protein_coding | 4               |
    | ENSG00000118181 |     0 |  0.7149172 | 1.000 | 0.978 |         0 | 0       | RPS25              | ribosomal protein S25 \[Source:HGNC Symbol;Acc:10413\]                                                       | protein_coding | 11              |
    | ENSG00000163220 |     0 |  4.4535241 | 0.998 | 0.212 |         0 | 1       | S100A9             | S100 calcium binding protein A9 \[Source:HGNC Symbol;Acc:10499\]                                             | protein_coding | 1               |
    | ENSG00000090382 |     0 |  4.0594058 | 1.000 | 0.517 |         0 | 1       | LYZ                | lysozyme \[Source:HGNC Symbol;Acc:6740\]                                                                     | protein_coding | 12              |
    | ENSG00000143546 |     0 |  3.9631418 | 0.971 | 0.119 |         0 | 1       | S100A8             | S100 calcium binding protein A8 \[Source:HGNC Symbol;Acc:10498\]                                             | protein_coding | 1               |
    | ENSG00000227507 |     0 |  1.1123793 | 0.979 | 0.649 |         0 | 2       | LTB                | lymphotoxin beta (TNF superfamily, member 3) \[Source:HGNC Symbol;Acc:6711\]                                 | protein_coding | 6               |
    | ENSG00000008517 |     0 |  0.9990429 | 0.949 | 0.473 |         0 | 2       | IL32               | interleukin 32 \[Source:HGNC Symbol;Acc:16830\]                                                              | protein_coding | 16              |
    | ENSG00000168685 |     0 |  0.7342628 | 0.733 | 0.336 |         0 | 2       | IL7R               | interleukin 7 receptor \[Source:HGNC Symbol;Acc:6024\]                                                       | protein_coding | 5               |
    | ENSG00000019582 |     0 |  2.7304144 | 1.000 | 0.820 |         0 | 3       | CD74               | CD74 molecule, major histocompatibility complex, class II invariant chain \[Source:HGNC Symbol;Acc:1697\]    | protein_coding | 5               |
    | ENSG00000105369 |     0 |  2.5281075 | 0.939 | 0.041 |         0 | 3       | CD79A              | CD79a molecule, immunoglobulin-associated alpha \[Source:HGNC Symbol;Acc:1698\]                              | protein_coding | 19              |
    | ENSG00000204287 |     0 |  2.5010512 | 1.000 | 0.492 |         0 | 3       | HLA-DRA            | major histocompatibility complex, class II, DR alpha \[Source:HGNC Symbol;Acc:4947\]                         | protein_coding | 6               |
    | ENSG00000161570 |     0 |  2.6006662 | 0.968 | 0.223 |         0 | 4       | CCL5               | chemokine (C-C motif) ligand 5 \[Source:HGNC Symbol;Acc:10632\]                                              | protein_coding | 17              |
    | ENSG00000105374 |     0 |  1.8826628 | 0.906 | 0.215 |         0 | 4       | NKG7               | natural killer cell group 7 sequence \[Source:HGNC Symbol;Acc:7830\]                                         | protein_coding | 19              |
    | ENSG00000113088 |     0 |  1.6720773 | 0.582 | 0.050 |         0 | 4       | GZMK               | granzyme K (granzyme 3; tryptase II) \[Source:HGNC Symbol;Acc:4711\]                                         | protein_coding | 5               |
    | ENSG00000204482 |     0 |  2.5414949 | 1.000 | 0.315 |         0 | 5       | LST1               | leukocyte specific transcript 1 \[Source:HGNC Symbol;Acc:14189\]                                             | protein_coding | 6               |
    | ENSG00000203747 |     0 |  2.2148756 | 0.975 | 0.135 |         0 | 5       | FCGR3A             | Fc fragment of IgG, low affinity IIIa, receptor (CD16a) \[Source:HGNC Symbol;Acc:3619\]                      | protein_coding | 1               |
    | ENSG00000158869 |     0 |  2.2070349 | 0.994 | 0.317 |         0 | 5       | FCER1G             | Fc fragment of IgE, high affinity I, receptor for; gamma polypeptide \[Source:HGNC Symbol;Acc:3611\]         | protein_coding | 1               |
    | ENSG00000115523 |     0 |  4.2579309 | 0.968 | 0.131 |         0 | 6       | GNLY               | granulysin \[Source:HGNC Symbol;Acc:4414\]                                                                   | protein_coding | 2               |
    | ENSG00000105374 |     0 |  3.4664868 | 1.000 | 0.260 |         0 | 6       | NKG7               | natural killer cell group 7 sequence \[Source:HGNC Symbol;Acc:7830\]                                         | protein_coding | 19              |
    | ENSG00000100453 |     0 |  3.3326422 | 0.955 | 0.068 |         0 | 6       | GZMB               | granzyme B (granzyme 2, cytotoxic T-lymphocyte-associated serine esterase 1) \[Source:HGNC Symbol;Acc:4709\] | protein_coding | 14              |
    | ENSG00000223865 |     0 |  2.5518720 | 1.000 | 0.510 |         0 | 7       | HLA-DPB1           | major histocompatibility complex, class II, DP beta 1 \[Source:HGNC Symbol;Acc:4940\]                        | protein_coding | 6               |
    | ENSG00000204287 |     0 |  2.3823825 | 1.000 | 0.552 |         0 | 7       | HLA-DRA            | major histocompatibility complex, class II, DR alpha \[Source:HGNC Symbol;Acc:4947\]                         | protein_coding | 6               |
    | ENSG00000101439 |     0 |  2.3027895 | 1.000 | 0.389 |         0 | 7       | CST3               | cystatin C \[Source:HGNC Symbol;Acc:2475\]                                                                   | protein_coding | 20              |
    | ENSG00000163736 |     0 |  6.4755064 | 1.000 | 0.025 |         0 | 8       | PPBP               | pro-platelet basic protein (chemokine (C-X-C motif) ligand 7) \[Source:HGNC Symbol;Acc:9240\]                | protein_coding | 4               |
    | ENSG00000163737 |     0 |  5.4177748 | 1.000 | 0.011 |         0 | 8       | PF4                | platelet factor 4 \[Source:HGNC Symbol;Acc:8861\]                                                            | protein_coding | 4               |
    | ENSG00000120885 |     0 |  4.5424885 | 0.857 | 0.015 |         0 | 8       | CLU                | clusterin \[Source:HGNC Symbol;Acc:2095\]                                                                    | protein_coding | 8               |


``` r
## Generate feature plots and stock them into variable
plots <- FeaturePlot(pbmc_small,                                                                #SeuratObject
                     features = top_n(x= pbmc_markers_signif, n = 1, wt = avg_log2FC)$gene,     #Vector of genes to plot
                     cols = c("yellow", "red"),                                                 #Change color
                     label = TRUE,                                                              #Plot ident position
                     combine = FALSE,                                                           #Return list of plot objets instead of a combined plot (easier to process)
                     repel = TRUE)                                                              #Avoid label overlap

## Add gene name as title
plots <- lapply(plots,                                                                          #List of plots
                add_title_gene_name,                                                            #Function to apply to the list
                gene_format = annotated_hg19)                                                   #Fill in the dataframe parameter

## Plot list of plots
grid.arrange(grobs = plots)
```

<img src="../images/FilterMarkers-1.png" style="display: block; margin: auto;" />

``` r
## Generate violin plots and stock them into variable
vln_plots <- VlnPlot(pbmc_small,                                                                #SeuratObject
                     features = top_n(x= pbmc_markers_signif, n = 1, wt = avg_log2FC)$gene,     #Vector of genes to plot
                     combine = FALSE)                                                           #Return list of plot objets instead of a combined plot (easier to process)

## Add gene name as title
vln_plots <- lapply(vln_plots,                                                                  #List of plots
                    add_title_gene_name,                                                        #Function to apply to the list
                    gene_format = annotated_hg19)                                               #Fill in the dataframe parameter

## Remove unecessary legend
vln_plots <- lapply(vln_plots, function(plot){
  plot <- plot + theme(legend.position = "none")
  return(plot)
})

## Plot list of plots
grid.arrange(grobs = vln_plots)
```

<img src="../images/FilterMarkers-2.png" style="display: block; margin: auto;" />

We now have each Ensembl gene ID set associated with a gene name and
a description to help us identify the gene lists. This is easy if you
know the theory of the biology of your system but if you don't know
enough about the genes identified as markers, the enrichment methods
will help you.

## Via ClusterProfiler

To understand the relationship between genes specific to our clusters
we can use functional enrichment methods. There are two types of
functional enrichment methods:

- Over-Representation Analysis methods which are based on a ratio between
  the number of marker genes present in a functional gene set
  and the total number of genes present in this gene set.
- Gene Set Enrichment Analysis (GSEA) methods which calculate an enrichment
  rate from a ranking of genes.

An R package will allow us to perform these different analyses using several
databases. It is quite complete and I advise you to take the time to look at
the [documentation](https://yulab-smu.top/biomedical-knowledge-mining-book/index.html) because here we will only see a small overview.
