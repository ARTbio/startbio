# Over-representation analysis

## Gene Ontholgy

We will start by analysing our gene lists with the [Gene Ontology annotation](http://geneontology.org/)
which classifies genes into gene sets according to three main types of
information:

- Molecular Function (MF): protein activity of the gene product
- Biological Process (BP): set of protein activities leading to a common task
- Cellular Component (CC): location of the gene product

Each set of genes is called a gene set and is grouped according to Gene
Ontology terms (referred to here as GO terms). The GO terms are classified
in a tree structure with general gene sets that become more specific as we
go along. The relationships between GO terms can reflect different cases:

- `is a` : A GO term is a subtype of the GO term B
- `part of` : the term GO A is part of the term GO B, so if the term GO A is
  present then so is the term GO B
- `has part` : the term GO A necessarily contains the term GO B, but if there
  is the term GO B there is not necessarily the term GO A
- `regulates` : the term GO A necessarily impacts the term GO B, but the
  latter is not necessarily impacted by A

In the `ClusterProfiler` package we will use the `enrichGO` function which
calculates for each GO term the over-representation of the genes in the
analysed cluster among those in the term.

To use the Gene Ontology database, we use an R package containing all the
annotation of the desired organism. These packages are called `OrgDb` which
can be found in the form `org.[Genome Initials].eg.db` (for the human
annotation, the OrgDb is `org.Hs.eg.db`). They contain gene identifiers
from different resources (NCBI, Ensembl, RefSeq, etc...) with annotation
from different databases (GO, OMIM, PMID, Uniprot, etc...). The default
gene identifiers are in the format `entrez` (which is a sequence of
numbers). Since we only have the Ensembl identifiers or the gene name
available to us, we will use the `keyType` parameter of the `enrichGO`
function to use the gene names instead of the entrez identifiers and thus
use the different information contained in the `Org.Db`.

``` r
## What's inside an organism database ?
ls("package:org.Hs.eg.db")

## You must see the help section of category you want to know about and don't hesite to test the examples to understand the architecture
?org.Hs.egENSEMBL    #Link between ensembl ID and entrez ID
?org.Hs.egSYMBOL2EG  #Link between entrez ID and gene name
?org.Hs.egSYMBOL     #Link between gene name and entrez ID
?org.Hs.egGENENAME   #Beware ! It concern gene description and not gene name as we know

## Different available mapping variable name
columns(org.Hs.eg.db)
```

We will therefore apply the `enrichGO` function for each clusters through
a `lapply`. For each cluster :

- We filter the marker result dataframe to retrieve only the rows concerning
  the genes overexpressed by the cluster cells.
- We run the `enrichGO` function to obtain the GO terms (BP, CC and MF)
  enriched in the marker gene set.
- Add the cluster name in a new column to the resulting dataframe
- We visualize the results with the `dotplot` function of the `enrichR`
  package which allows to visualize the first 3 GO terms for each ontology type

The result of the lapply is a list where each element of the list is a
resultant dataframe for each cluster. We then assemble the results with
the `do.call` function which applies a function to the whole list. The
`rbind` will concatenate the rows of all the elements in the list so that
there is only one dataframe with the results for each cluster.

``` r
## GO enrichment for all clusters
enrich_go_list <- lapply(levels(pbmc_markers_signif$cluster), function(cluster_name){

  res_markers <- subset(pbmc_markers_signif, cluster == cluster_name & avg_log2FC > 0) #Filter markers to retrieve only positive markers for the specific cluster

  ego <- enrichGO(gene = res_markers$external_gene_name,                 #Vector of target genes
                  universe = unique(annotated_hg19$external_gene_name),  #Vector of reference genes (all genes from the differential analysis)
                  OrgDb = "org.Hs.eg.db",                                #Organisme database
                  keyType = "SYMBOL",                                    #Column name of the OrgDB that convert gene format in `gene`parameter to entrez ID
                  ont = "ALL")                                           #What category of GO you want to analyse (BP, CC, MF or ALL)

  ego@result$cluster <- cluster_name                                     #Add cluster name in a new column named "cluster"

  ## visualisation
  ### don't forget to add print when inside a function/loop/lapply
  print(dotplot(ego,                                                     #enrichResult object
                split = "ONTOLOGY",                                      #Do separated plot for each ontology type (only valable fo GO results)
                showCategory = 3,                                        #Only show first three categories
                title = paste("Cluster", cluster_name)) +                #Add title
    facet_grid(ONTOLOGY~., scales = "free_y") +                          #Create subplot according to type used with `split = "ONTOLOGY"`
    theme(axis.text.y = element_text(size = 5),
                legend.key.size = unit(0.2, 'cm')))                      #Reduce ontology labels names
  return(ego@result)
})
```

<img src="../images/clusterGO-1.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-2.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-3.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-4.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-5.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-6.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-7.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-8.png" style="display: block; margin: auto;" /><img src="../images/clusterGO-9.png" style="display: block; margin: auto;" />

``` r
enrich_go <- do.call("rbind", enrich_go_list)                            #Bind all results together
enrich_go <- enrich_go %>%
  group_by(cluster, ONTOLOGY)                                            #Rearrange df according to cluster and ontology type


## show first results
paged_table(top_n(x= enrich_go, n = 3, wt = (Count))[,-9])               #Only remove list of genes for the visualisation
```

The result is a dataframe where the GO terms have been considered as enriched:

- `ONTOLOGY` :  ontology type (BP, CC, or MF)
- `ID` : Unique identifier of the GO term
- `Description` : Description of the GO term
- `GeneRatio` : Fraction representing the number of marker genes present
  in the GO term `nbr_marker_gene_in_GOterm / nbr_marker_gene`.
- `BgRatio` : Fraction representing the number of reference genes present
  in the GO term `nbr_total_gene_in_GOterm / nbr_total_gene`.
- `pvalue` : p-value of the enrichment test
- `p.adjust` : adjusted p-value of the Benjamini Hochberg test
- `qvalue` : q-value after FDR (False Discovery Rate) check
- `geneID` : String containing the list of marker genes present in the GO
  term (separated by `/`)
- `count` : Number of marker genes present in the GO term
- `cluster` : Column added before the `do.call("rbind", list)` to identify
  in which cluster the GO term has been considered as enriched.

A GO term has been considered as enriched if :

- the p-value \< 0.05
- the adjusted p-value \< 0.05
- the q-value \< 0.2

## Kyoto Encyclopedia of Genes and Genomes (KEGG)

[KEGG](http://www.genome.jp/kegg/) is a database that focuses on the molecular
annotation of the different metabolic pathways. It allows the description of
the biochemical reactions that compose the pathway. KEGG also allows the
visualization of these pathways through hand-drawn maps representing the
different reactions and the relationship between the genes. There are
several main categories of KEGG pathways:

- Metabolism
- Genetic information processing
- Environmental information processing
- Cellular processes
- Organ systems
- Human diseases
- Drug development

We will use the `enrichKEGG` function from the `ClusterProfiler` package.
This function doesn't work exactly like `enrichGO` because it calls directly
on the database and doesn't use an `Orgdb` package but it does require our
genes to be annotated with an entrez id. We will have to find another way
to convert our genes into the correct format.

To do this, we will use the `Orgdb` or `org.Hs.egSYMBOL` object to list
each enter id as its "*gene symbol*" (gene name). When we convert this object
to a dataframe we get a table with two columns (`gene_id` and `symbol`).
We now have the possibility to switch from a gene name to an entrez id.

We will therefore apply the `enrichKEGG` function for each clusters
through a `lapply`. For each cluster :

- We filter the marker result dataframe to retrieve only the rows
  concerning the genes overexpressed by the cluster cells.
- Run the `enrichGO` function to get the enriched KEGG terms in the marker
  gene set, filter our enter/symbol mapping table for the `gene` and
  `universe` parameters with the gene names contained in the marker and
  biomart annotation tables.
- Add the cluster name in a new column to the resulting dataframe
- We visualize the results with the `dotplot` function of the `enrichR`
  package which allows to visualize the first 5 KEGG terms

``` r
## Retrieve a corresponding table between entrez id and gene name (called gene symbol in org.db)
corresp_entrez <- as.data.frame(org.Hs.egSYMBOL)  #Change format to df

## Apply enrichKEGG for each cluster
enrich_kegg_list <- lapply(levels(pbmc_markers_signif$cluster), function(cluster_name){

  res_markers <- subset(pbmc_markers_signif, cluster == cluster_name & avg_log2FC > 0) #Filter markers dataframe based on cluster

  ## Perform enrichKEGG analysis
  ekegg <- enrichKEGG(gene = subset(corresp_entrez, symbol %in% res_markers$external_gene_name)$gene_id,                #Genes to analyse
                      universe = subset(corresp_entrez, symbol %in% unique(annotated_hg19$external_gene_name))$gene_id, #Background genes, here we take all genes from our expression matrix
                      organism = "hsa")

  ekegg@result$cluster <- cluster_name            #Add cluster name as column

  ## Add plot
  print(dotplot(ekegg,
                label_format = 30,
                font.size = 10,
                showCategory = 5,
                title = paste("Cluster", cluster_name)) +
    theme(axis.text.y = element_text(size = 10),
                legend.key.size = unit(0.2, 'cm')))

  return(ekegg@result)                            #Return dataframe result
})
```

<img src="../images/clusterKEGG-1.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-2.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-3.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-4.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-5.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-6.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-7.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-8.png" style="display: block; margin: auto;" /><img src="../images/clusterKEGG-9.png" style="display: block; margin: auto;" />

``` r
## Concatenate all results in one dataframe
enrich_kegg <- do.call("rbind", enrich_kegg_list)

## Group result by cluster (easier to manipulate with dplyr)
enrich_kegg <- enrich_kegg %>%
  group_by(cluster)

## Visualise first 3 KEGG categories for each cluster (removing the vector of genes just for the visualisation)
paged_table(top_n(x= enrich_kegg, n = 3, wt = (Count))[,-8])
```

The result is a dataframe where KEGG categories have been considered
as enriched with the following columns :

- `ID` : Unique identifier of the KEGG category
- `Description` : Description of the KEGG category
- `GeneRatio` : Fraction representing the number of marker genes present in
  the KEGG category $\GeneRatio = {frac{nbrMarkerGeneInKEGGcat}{nbrMarkerGene}$
- `BgRatio` : Fraction representing the number of genes of the reference
  present in the KEGG category
  $\BgRatio = \frac{nbrTotalGeneInKEGGcat}{nbrTotalGene}$
- `pvalue` : p-value of the enrichment test
- `p.adjust` : adjusted p-value of the Benjamini Hochberg test
- `qvalue` : q-value after FDR (False Discovery Rate) check
- `geneID` : String containing the list of marker genes present in the KEGG
  category (separated by `/`)
- `Count` : Number of marker genes present in the KEGG category
- `cluster` : Column added before the `do.call("rbind", list)` in order to
  identify in which cluster the KEGG category has been considered as enriched.

A KEGG category has been considered enriched if :

- the p-value \< 0.05
- the adjusted p-value \< 0.05
- q-value \< 0.2


These over-representation methods are highly dependent on the database
containing fairly generalized groups of genes. However, we can see that
the results of `enrichGO`, `enrichKEGG` and Biomart intersect for some
clusters where :

- Cluster 6 would represent the *Natural Killer* cells as well as cluster 4:
  knowing that they are very close on the UMAP it reflects their proximity
  of the transcriptomes of the cells that compose these two clusters. The GO
  analysis would however lean more towards T cells for cluster 4.
- Cluster 8 would be composed of platelet cells

Knowing when taking only the first 3 or 5 results (so very restricted),
we are extremely stringent in identifying clusters.
