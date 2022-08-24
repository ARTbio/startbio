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
paged_table(top_n(x= pbmc_markers_signif, n = 3, wt = avg_log2FC))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["gene"],"name":[1],"type":["chr"],"align":["left"]},{"label":["p_val"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["avg_log2FC"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["pct.1"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["pct.2"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["p_val_adj"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["cluster"],"name":[7],"type":["fct"],"align":["left"]},{"label":["external_gene_name"],"name":[8],"type":["chr"],"align":["left"]},{"label":["description"],"name":[9],"type":["chr"],"align":["left"]},{"label":["gene_biotype"],"name":[10],"type":["chr"],"align":["left"]},{"label":["chromosome_name"],"name":[11],"type":["chr"],"align":["left"]}],"data":[{"1":"ENSG00000111716","2":"3.970716e-116","3":"0.7920556","4":"0.929","5":"0.592","6":"1.299933e-111","7":"0","8":"LDHB","9":"lactate dehydrogenase B [Source:HGNC Symbol;Acc:6541]","10":"protein_coding","11":"12"},{"1":"ENSG00000145425","2":"1.049161e-115","3":"0.7576946","4":"0.997","5":"0.976","6":"3.434742e-111","7":"0","8":"RPS3A","9":"ribosomal protein S3A [Source:HGNC Symbol;Acc:10421]","10":"protein_coding","11":"4"},{"1":"ENSG00000118181","2":"8.077812e-134","3":"0.7149172","4":"1.000","5":"0.978","6":"2.644514e-129","7":"0","8":"RPS25","9":"ribosomal protein S25 [Source:HGNC Symbol;Acc:10413]","10":"protein_coding","11":"11"},{"1":"ENSG00000163220","2":"0.000000e+00","3":"4.4535241","4":"0.998","5":"0.212","6":"0.000000e+00","7":"1","8":"S100A9","9":"S100 calcium binding protein A9 [Source:HGNC Symbol;Acc:10499]","10":"protein_coding","11":"1"},{"1":"ENSG00000090382","2":"5.690667e-272","3":"4.0594058","4":"1.000","5":"0.517","6":"1.863010e-267","7":"1","8":"LYZ","9":"lysozyme [Source:HGNC Symbol;Acc:6740]","10":"protein_coding","11":"12"},{"1":"ENSG00000143546","2":"0.000000e+00","3":"3.9631418","4":"0.971","5":"0.119","6":"0.000000e+00","7":"1","8":"S100A8","9":"S100 calcium binding protein A8 [Source:HGNC Symbol;Acc:10498]","10":"protein_coding","11":"1"},{"1":"ENSG00000227507","2":"4.202245e-87","3":"1.1123793","4":"0.979","5":"0.649","6":"1.375731e-82","7":"2","8":"LTB","9":"lymphotoxin beta (TNF superfamily, member 3) [Source:HGNC Symbol;Acc:6711]","10":"protein_coding","11":"6"},{"1":"ENSG00000008517","2":"4.821630e-93","3":"0.9990429","4":"0.949","5":"0.473","6":"1.578505e-88","7":"2","8":"IL32","9":"interleukin 32 [Source:HGNC Symbol;Acc:16830]","10":"protein_coding","11":"16"},{"1":"ENSG00000168685","2":"6.235001e-61","3":"0.7342628","4":"0.733","5":"0.336","6":"2.041215e-56","7":"2","8":"IL7R","9":"interleukin 7 receptor [Source:HGNC Symbol;Acc:6024]","10":"protein_coding","11":"5"},{"1":"ENSG00000019582","2":"1.164791e-185","3":"2.7304144","4":"1.000","5":"0.820","6":"3.813294e-181","7":"3","8":"CD74","9":"CD74 molecule, major histocompatibility complex, class II invariant chain [Source:HGNC Symbol;Acc:1697]","10":"protein_coding","11":"5"},{"1":"ENSG00000105369","2":"0.000000e+00","3":"2.5281075","4":"0.939","5":"0.041","6":"0.000000e+00","7":"3","8":"CD79A","9":"CD79a molecule, immunoglobulin-associated alpha [Source:HGNC Symbol;Acc:1698]","10":"protein_coding","11":"19"},{"1":"ENSG00000204287","2":"3.209089e-184","3":"2.5010512","4":"1.000","5":"0.492","6":"1.050592e-179","7":"3","8":"HLA-DRA","9":"major histocompatibility complex, class II, DR alpha [Source:HGNC Symbol;Acc:4947]","10":"protein_coding","11":"6"},{"1":"ENSG00000161570","2":"5.543749e-227","3":"2.6006662","4":"0.968","5":"0.223","6":"1.814912e-222","7":"4","8":"CCL5","9":"chemokine (C-C motif) ligand 5 [Source:HGNC Symbol;Acc:10632]","10":"protein_coding","11":"17"},{"1":"ENSG00000105374","2":"1.926967e-174","3":"1.8826628","4":"0.906","5":"0.215","6":"6.308504e-170","7":"4","8":"NKG7","9":"natural killer cell group 7 sequence [Source:HGNC Symbol;Acc:7830]","10":"protein_coding","11":"19"},{"1":"ENSG00000113088","2":"1.462695e-187","3":"1.6720773","4":"0.582","5":"0.050","6":"4.788572e-183","7":"4","8":"GZMK","9":"granzyme K (granzyme 3; tryptase II) [Source:HGNC Symbol;Acc:4711]","10":"protein_coding","11":"5"},{"1":"ENSG00000204482","2":"4.307394e-124","3":"2.5414949","4":"1.000","5":"0.315","6":"1.410155e-119","7":"5","8":"LST1","9":"leukocyte specific transcript 1 [Source:HGNC Symbol;Acc:14189]","10":"protein_coding","11":"6"},{"1":"ENSG00000203747","2":"8.080338e-183","3":"2.2148756","4":"0.975","5":"0.135","6":"2.645341e-178","7":"5","8":"FCGR3A","9":"Fc fragment of IgG, low affinity IIIa, receptor (CD16a) [Source:HGNC Symbol;Acc:3619]","10":"protein_coding","11":"1"},{"1":"ENSG00000158869","2":"6.999917e-116","3":"2.2070349","4":"0.994","5":"0.317","6":"2.291633e-111","7":"5","8":"FCER1G","9":"Fc fragment of IgE, high affinity I, receptor for; gamma polypeptide [Source:HGNC Symbol;Acc:3611]","10":"protein_coding","11":"1"},{"1":"ENSG00000115523","2":"2.486700e-196","3":"4.2579309","4":"0.968","5":"0.131","6":"8.140959e-192","7":"6","8":"GNLY","9":"granulysin [Source:HGNC Symbol;Acc:4414]","10":"protein_coding","11":"2"},{"1":"ENSG00000105374","2":"3.209039e-138","3":"3.4664868","4":"1.000","5":"0.260","6":"1.050575e-133","7":"6","8":"NKG7","9":"natural killer cell group 7 sequence [Source:HGNC Symbol;Acc:7830]","10":"protein_coding","11":"19"},{"1":"ENSG00000100453","2":"5.211408e-268","3":"3.3326422","4":"0.955","5":"0.068","6":"1.706111e-263","7":"6","8":"GZMB","9":"granzyme B (granzyme 2, cytotoxic T-lymphocyte-associated serine esterase 1) [Source:HGNC Symbol;Acc:4709]","10":"protein_coding","11":"14"},{"1":"ENSG00000223865","2":"1.415112e-21","3":"2.5518720","4":"1.000","5":"0.510","6":"4.632794e-17","7":"7","8":"HLA-DPB1","9":"major histocompatibility complex, class II, DP beta 1 [Source:HGNC Symbol;Acc:4940]","10":"protein_coding","11":"6"},{"1":"ENSG00000204287","2":"6.910937e-20","3":"2.3823825","4":"1.000","5":"0.552","6":"2.262502e-15","7":"7","8":"HLA-DRA","9":"major histocompatibility complex, class II, DR alpha [Source:HGNC Symbol;Acc:4947]","10":"protein_coding","11":"6"},{"1":"ENSG00000101439","2":"1.089086e-21","3":"2.3027895","4":"1.000","5":"0.389","6":"3.565450e-17","7":"7","8":"CST3","9":"cystatin C [Source:HGNC Symbol;Acc:2475]","10":"protein_coding","11":"20"},{"1":"ENSG00000163736","2":"1.575953e-59","3":"6.4755064","4":"1.000","5":"0.025","6":"5.159355e-55","7":"8","8":"PPBP","9":"pro-platelet basic protein (chemokine (C-X-C motif) ligand 7) [Source:HGNC Symbol;Acc:9240]","10":"protein_coding","11":"4"},{"1":"ENSG00000163737","2":"1.324207e-112","3":"5.4177748","4":"1.000","5":"0.011","6":"4.335188e-108","7":"8","8":"PF4","9":"platelet factor 4 [Source:HGNC Symbol;Acc:8861]","10":"protein_coding","11":"4"},{"1":"ENSG00000120885","2":"3.260924e-67","3":"4.5424885","4":"0.857","5":"0.015","6":"1.067561e-62","7":"8","8":"CLU","9":"clusterin [Source:HGNC Symbol;Acc:2095]","10":"protein_coding","11":"8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

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
