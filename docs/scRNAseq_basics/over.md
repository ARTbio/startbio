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
kable(top_n(x= enrich_go, n = 3, wt = (Count))[,-9])               #Only remove list of genes for the visualisation
```

??? abstract "First Enriched GO terms for each cluster"

    | ONTOLOGY | ID           | Description                                                        | GeneRatio | BgRatio   |    pvalue |  p.adjust |    qvalue | Count | cluster |
    |:---------|:-------------|:-------------------------------------------------------------------|:----------|:----------|----------:|----------:|----------:|------:|:--------|
    | BP       | <GO:0006413> | translational initiation                                           | 69/111    | 191/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 0       |
    | BP       | <GO:0006402> | mRNA catabolic process                                             | 69/111    | 359/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 0       |
    | BP       | <GO:0006401> | RNA catabolic process                                              | 69/111    | 395/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 0       |
    | CC       | <GO:0022626> | cytosolic ribosome                                                 | 66/113    | 104/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | CC       | <GO:0044391> | ribosomal subunit                                                  | 66/113    | 179/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | CC       | <GO:0005840> | ribosome                                                           | 66/113    | 216/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | MF       | <GO:0003735> | structural constituent of ribosome                                 | 66/112    | 158/16533 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | MF       | <GO:0019843> | rRNA binding                                                       | 15/112    | 58/16533  | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 0       |
    | MF       | <GO:0003729> | mRNA binding                                                       | 14/112    | 284/16533 | 0.0000000 | 0.0000005 | 0.0000005 |    14 | 0       |
    | BP       | <GO:0043312> | neutrophil degranulation                                           | 69/274    | 468/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 1       |
    | BP       | <GO:0036230> | granulocyte activation                                             | 70/274    | 486/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    70 | 1       |
    | BP       | <GO:0002283> | neutrophil activation involved in immune response                  | 69/274    | 471/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 1       |
    | BP       | <GO:0042119> | neutrophil activation                                              | 69/274    | 481/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 1       |
    | BP       | <GO:0002446> | neutrophil mediated immunity                                       | 69/274    | 482/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    69 | 1       |
    | CC       | <GO:0034774> | secretory granule lumen                                            | 37/285    | 314/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    37 | 1       |
    | CC       | <GO:0060205> | cytoplasmic vesicle lumen                                          | 37/285    | 318/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    37 | 1       |
    | CC       | <GO:0031983> | vesicle lumen                                                      | 37/285    | 320/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    37 | 1       |
    | MF       | <GO:0009055> | electron transfer activity                                         | 18/278    | 121/16533 | 0.0000000 | 0.0000000 | 0.0000000 |    18 | 1       |
    | MF       | <GO:0033218> | amide binding                                                      | 17/278    | 363/16533 | 0.0001466 | 0.0041824 | 0.0035415 |    17 | 1       |
    | MF       | <GO:0004857> | enzyme inhibitor activity                                          | 16/278    | 362/16533 | 0.0004353 | 0.0090741 | 0.0076836 |    16 | 1       |
    | BP       | <GO:0042110> | T cell activation                                                  | 17/67     | 458/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    17 | 2       |
    | BP       | <GO:0006402> | mRNA catabolic process                                             | 15/67     | 359/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 2       |
    | BP       | <GO:0006401> | RNA catabolic process                                              | 15/67     | 395/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 2       |
    | BP       | <GO:1903131> | mononuclear cell differentiation                                   | 15/67     | 402/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 2       |
    | CC       | <GO:0022626> | cytosolic ribosome                                                 | 12/70     | 104/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    12 | 2       |
    | CC       | <GO:0044391> | ribosomal subunit                                                  | 12/70     | 179/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    12 | 2       |
    | CC       | <GO:0005840> | ribosome                                                           | 12/70     | 216/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    12 | 2       |
    | MF       | <GO:0003735> | structural constituent of ribosome                                 | 12/69     | 158/16533 | 0.0000000 | 0.0000000 | 0.0000000 |    12 | 2       |
    | MF       | <GO:0030674> | protein-macromolecule adaptor activity                             | 7/69      | 254/16533 | 0.0000891 | 0.0037620 | 0.0029300 |     7 | 2       |
    | MF       | <GO:0060090> | molecular adaptor activity                                         | 7/69      | 307/16533 | 0.0002863 | 0.0079818 | 0.0062167 |     7 | 2       |
    | MF       | <GO:0050839> | cell adhesion molecule binding                                     | 8/69      | 500/16533 | 0.0010961 | 0.0223564 | 0.0174124 |     8 | 2       |
    | BP       | <GO:0006413> | translational initiation                                           | 27/120    | 191/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | BP       | <GO:0006612> | protein targeting to membrane                                      | 27/120    | 197/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | BP       | <GO:0090150> | establishment of protein localization to membrane                  | 27/120    | 328/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | BP       | <GO:0006401> | RNA catabolic process                                              | 28/120    | 395/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    28 | 3       |
    | BP       | <GO:0006605> | protein targeting                                                  | 28/120    | 420/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    28 | 3       |
    | BP       | <GO:0002250> | adaptive immune response                                           | 27/120    | 401/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | CC       | <GO:0022626> | cytosolic ribosome                                                 | 25/125    | 104/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    25 | 3       |
    | CC       | <GO:0005840> | ribosome                                                           | 27/125    | 216/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | CC       | <GO:0044391> | ribosomal subunit                                                  | 25/125    | 179/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    25 | 3       |
    | MF       | <GO:0003735> | structural constituent of ribosome                                 | 26/123    | 158/16533 | 0.0000000 | 0.0000000 | 0.0000000 |    26 | 3       |
    | MF       | <GO:0019843> | rRNA binding                                                       | 11/123    | 58/16533  | 0.0000000 | 0.0000000 | 0.0000000 |    11 | 3       |
    | MF       | <GO:0140375> | immune receptor activity                                           | 11/123    | 131/16533 | 0.0000000 | 0.0000001 | 0.0000001 |    11 | 3       |
    | BP       | <GO:0002250> | adaptive immune response                                           | 20/66     | 401/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    20 | 4       |
    | BP       | <GO:0042110> | T cell activation                                                  | 17/66     | 458/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    17 | 4       |
    | BP       | <GO:0002694> | regulation of leukocyte activation                                 | 15/66     | 493/16175 | 0.0000000 | 0.0000001 | 0.0000001 |    15 | 4       |
    | CC       | <GO:0009897> | external side of plasma membrane                                   | 12/67     | 312/16891 | 0.0000000 | 0.0000002 | 0.0000001 |    12 | 4       |
    | CC       | <GO:0030135> | coated vesicle                                                     | 10/67     | 283/16891 | 0.0000002 | 0.0000027 | 0.0000021 |    10 | 4       |
    | CC       | <GO:0010008> | endosome membrane                                                  | 11/67     | 483/16891 | 0.0000028 | 0.0000330 | 0.0000251 |    11 | 4       |
    | MF       | <GO:0042608> | T cell receptor binding                                            | 5/67      | 11/16533  | 0.0000000 | 0.0000001 | 0.0000001 |     5 | 4       |
    | MF       | <GO:0003823> | antigen binding                                                    | 7/67      | 52/16533  | 0.0000000 | 0.0000002 | 0.0000001 |     7 | 4       |
    | MF       | <GO:0042605> | peptide antigen binding                                            | 5/67      | 23/16533  | 0.0000000 | 0.0000021 | 0.0000017 |     5 | 4       |
    | MF       | <GO:0004252> | serine-type endopeptidase activity                                 | 5/67      | 156/16533 | 0.0004229 | 0.0113647 | 0.0089582 |     5 | 4       |
    | MF       | <GO:0008236> | serine-type peptidase activity                                     | 5/67      | 174/16533 | 0.0006950 | 0.0146216 | 0.0115254 |     5 | 4       |
    | MF       | <GO:0017171> | serine hydrolase activity                                          | 5/67      | 176/16533 | 0.0007318 | 0.0146216 | 0.0115254 |     5 | 4       |
    | BP       | <GO:0002446> | neutrophil mediated immunity                                       | 59/274    | 482/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    59 | 5       |
    | BP       | <GO:0043312> | neutrophil degranulation                                           | 58/274    | 468/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    58 | 5       |
    | BP       | <GO:0002283> | neutrophil activation involved in immune response                  | 58/274    | 471/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    58 | 5       |
    | BP       | <GO:0042119> | neutrophil activation                                              | 58/274    | 481/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    58 | 5       |
    | BP       | <GO:0036230> | granulocyte activation                                             | 58/274    | 486/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    58 | 5       |
    | CC       | <GO:0060205> | cytoplasmic vesicle lumen                                          | 33/285    | 318/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    33 | 5       |
    | CC       | <GO:0031983> | vesicle lumen                                                      | 33/285    | 320/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    33 | 5       |
    | CC       | <GO:0034774> | secretory granule lumen                                            | 32/285    | 314/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    32 | 5       |
    | MF       | <GO:0003779> | actin binding                                                      | 26/277    | 424/16533 | 0.0000000 | 0.0000066 | 0.0000056 |    26 | 5       |
    | MF       | <GO:0044389> | ubiquitin-like protein ligase binding                              | 16/277    | 300/16533 | 0.0000487 | 0.0041852 | 0.0035090 |    16 | 5       |
    | MF       | <GO:0050839> | cell adhesion molecule binding                                     | 20/277    | 500/16533 | 0.0003057 | 0.0162987 | 0.0136653 |    20 | 5       |
    | BP       | <GO:0002768> | immune response-regulating cell surface receptor signaling pathway | 30/201    | 385/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    30 | 6       |
    | BP       | <GO:0002764> | immune response-regulating signaling pathway                       | 30/201    | 388/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    30 | 6       |
    | BP       | <GO:0002831> | regulation of response to biotic stimulus                          | 29/201    | 395/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    29 | 6       |
    | BP       | <GO:0002253> | activation of immune response                                      | 29/201    | 437/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    29 | 6       |
    | BP       | <GO:0032103> | positive regulation of response to external stimulus               | 29/201    | 479/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    29 | 6       |
    | BP       | <GO:0036230> | granulocyte activation                                             | 29/201    | 486/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    29 | 6       |
    | CC       | <GO:0005925> | focal adhesion                                                     | 28/212    | 407/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    28 | 6       |
    | CC       | <GO:0030055> | cell-substrate junction                                            | 28/212    | 413/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    28 | 6       |
    | CC       | <GO:0015629> | actin cytoskeleton                                                 | 22/212    | 478/16891 | 0.0000002 | 0.0000062 | 0.0000048 |    22 | 6       |
    | MF       | <GO:0003779> | actin binding                                                      | 15/211    | 424/16533 | 0.0003694 | 0.0138507 | 0.0117285 |    15 | 6       |
    | MF       | <GO:0050839> | cell adhesion molecule binding                                     | 16/211    | 500/16533 | 0.0007065 | 0.0244203 | 0.0206787 |    16 | 6       |
    | MF       | <GO:0004175> | endopeptidase activity                                             | 14/211    | 411/16533 | 0.0008361 | 0.0250827 | 0.0212396 |    14 | 6       |
    | BP       | <GO:0002764> | immune response-regulating signaling pathway                       | 26/112    | 388/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    26 | 7       |
    | BP       | <GO:0002768> | immune response-regulating cell surface receptor signaling pathway | 24/112    | 385/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 7       |
    | BP       | <GO:0002253> | activation of immune response                                      | 25/112    | 437/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    25 | 7       |
    | BP       | <GO:0002283> | neutrophil activation involved in immune response                  | 24/112    | 471/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 7       |
    | BP       | <GO:0042119> | neutrophil activation                                              | 24/112    | 481/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 7       |
    | BP       | <GO:0002446> | neutrophil mediated immunity                                       | 24/112    | 482/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 7       |
    | BP       | <GO:0036230> | granulocyte activation                                             | 24/112    | 486/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 7       |
    | CC       | <GO:0005765> | lysosomal membrane                                                 | 20/115    | 351/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    20 | 7       |
    | CC       | <GO:0098852> | lytic vacuole membrane                                             | 20/115    | 351/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    20 | 7       |
    | CC       | <GO:0005774> | vacuolar membrane                                                  | 21/115    | 402/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    21 | 7       |
    | MF       | <GO:0042277> | peptide binding                                                    | 14/113    | 293/16533 | 0.0000000 | 0.0000013 | 0.0000010 |    14 | 7       |
    | MF       | <GO:0033218> | amide binding                                                      | 15/113    | 363/16533 | 0.0000000 | 0.0000015 | 0.0000012 |    15 | 7       |
    | MF       | <GO:0004857> | enzyme inhibitor activity                                          | 11/113    | 362/16533 | 0.0000377 | 0.0012609 | 0.0010186 |    11 | 7       |
    | BP       | <GO:0042060> | wound healing                                                      | 31/178    | 497/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    31 | 8       |
    | BP       | <GO:0007596> | blood coagulation                                                  | 23/178    | 321/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    23 | 8       |
    | BP       | <GO:0007599> | hemostasis                                                         | 23/178    | 325/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    23 | 8       |
    | BP       | <GO:0050817> | coagulation                                                        | 23/178    | 325/16175 | 0.0000000 | 0.0000000 | 0.0000000 |    23 | 8       |
    | BP       | <GO:0050878> | regulation of body fluid levels                                    | 23/178    | 477/16175 | 0.0000000 | 0.0000009 | 0.0000008 |    23 | 8       |
    | CC       | <GO:0015629> | actin cytoskeleton                                                 | 25/181    | 478/16891 | 0.0000000 | 0.0000000 | 0.0000000 |    25 | 8       |
    | CC       | <GO:0030055> | cell-substrate junction                                            | 20/181    | 413/16891 | 0.0000000 | 0.0000012 | 0.0000010 |    20 | 8       |
    | CC       | <GO:0005925> | focal adhesion                                                     | 19/181    | 407/16891 | 0.0000001 | 0.0000024 | 0.0000020 |    19 | 8       |
    | MF       | <GO:0003779> | actin binding                                                      | 19/179    | 424/16533 | 0.0000002 | 0.0000735 | 0.0000652 |    19 | 8       |
    | MF       | <GO:0051015> | actin filament binding                                             | 10/179    | 201/16533 | 0.0000685 | 0.0090608 | 0.0080366 |    10 | 8       |
    | MF       | <GO:0003924> | GTPase activity                                                    | 11/179    | 291/16533 | 0.0003372 | 0.0278493 | 0.0247013 |    11 | 8       |


The result is a dataframe where the GO terms have been considered as enriched:

- `ONTOLOGY` :  ontology type (BP, CC, or MF)
- `ID` : Unique identifier of the GO term
- `Description` : Description of the GO term
- `GeneRatio` : Fraction representing the number of marker genes present
  in the GO term, $GeneRatio = \frac{nbrMarkerGeneInKEGGcat}{nbrMarkerGene}$.
- `BgRatio` : Fraction representing the number of reference genes present
  in the GO term, $BgRatio = \frac{nbrTotalGeneInKEGGcat}{nbrTotalGene}$.
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
kable(top_n(x= enrich_kegg, n = 3, wt = (Count))[,-8])
```

??? abstract "First Enriched KEGG categories for each cluster"
    | ID       | Description                                       | GeneRatio | BgRatio  |    pvalue |  p.adjust |    qvalue | Count | cluster |
    |:---------|:--------------------------------------------------|:----------|:---------|----------:|----------:|----------:|------:|:--------|
    | hsa03010 | Ribosome                                          | 66/94     | 134/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | hsa05171 | Coronavirus disease - COVID-19                    | 66/94     | 228/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    66 | 0       |
    | hsa04640 | Hematopoietic cell lineage                        | 6/94      | 94/7558  | 0.0010652 | 0.0242334 | 0.0241073 |     6 | 0       |
    | hsa05415 | Diabetic cardiomyopathy                           | 24/179    | 175/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    24 | 1       |
    | hsa05020 | Prion disease                                     | 23/179    | 246/7558 | 0.0000000 | 0.0000004 | 0.0000003 |    23 | 1       |
    | hsa05022 | Pathways of neurodegeneration - multiple diseases | 23/179    | 443/7558 | 0.0002945 | 0.0018710 | 0.0013802 |    23 | 1       |
    | hsa03010 | Ribosome                                          | 12/44     | 134/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    12 | 2       |
    | hsa05171 | Coronavirus disease - COVID-19                    | 14/44     | 228/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    14 | 2       |
    | hsa05170 | Human immunodeficiency virus 1 infection          | 9/44      | 204/7558 | 0.0000020 | 0.0000413 | 0.0000317 |     9 | 2       |
    | hsa03010 | Ribosome                                          | 27/81     | 134/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    27 | 3       |
    | hsa05171 | Coronavirus disease - COVID-19                    | 28/81     | 228/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    28 | 3       |
    | hsa04640 | Hematopoietic cell lineage                        | 17/81     | 94/7558  | 0.0000000 | 0.0000000 | 0.0000000 |    17 | 3       |
    | hsa04650 | Natural killer cell mediated cytotoxicity         | 10/47     | 124/7558 | 0.0000000 | 0.0000002 | 0.0000001 |    10 | 4       |
    | hsa04514 | Cell adhesion molecules                           | 10/47     | 147/7558 | 0.0000000 | 0.0000005 | 0.0000004 |    10 | 4       |
    | hsa05170 | Human immunodeficiency virus 1 infection          | 11/47     | 204/7558 | 0.0000000 | 0.0000008 | 0.0000006 |    11 | 4       |
    | hsa05166 | Human T-cell leukemia virus 1 infection           | 11/47     | 219/7558 | 0.0000001 | 0.0000011 | 0.0000008 |    11 | 4       |
    | hsa04145 | Phagosome                                         | 21/180    | 146/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    21 | 5       |
    | hsa05152 | Tuberculosis                                      | 21/180    | 176/7558 | 0.0000000 | 0.0000001 | 0.0000001 |    21 | 5       |
    | hsa05022 | Pathways of neurodegeneration - multiple diseases | 21/180    | 443/7558 | 0.0018748 | 0.0133316 | 0.0100316 |    21 | 5       |
    | hsa04650 | Natural killer cell mediated cytotoxicity         | 23/126    | 124/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    23 | 6       |
    | hsa05163 | Human cytomegalovirus infection                   | 17/126    | 222/7558 | 0.0000001 | 0.0000043 | 0.0000033 |    17 | 6       |
    | hsa05170 | Human immunodeficiency virus 1 infection          | 16/126    | 204/7558 | 0.0000002 | 0.0000064 | 0.0000049 |    16 | 6       |
    | hsa05152 | Tuberculosis                                      | 17/72     | 176/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    17 | 7       |
    | hsa04145 | Phagosome                                         | 15/72     | 146/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 7       |
    | hsa05164 | Influenza A                                       | 15/72     | 166/7558 | 0.0000000 | 0.0000000 | 0.0000000 |    15 | 7       |
    | hsa05168 | Herpes simplex virus 1 infection                  | 15/72     | 484/7558 | 0.0000389 | 0.0002724 | 0.0002048 |    15 | 7       |
    | hsa04611 | Platelet activation                               | 9/97      | 122/7558 | 0.0000240 | 0.0043154 | 0.0039621 |     9 | 8       |
    | hsa04530 | Tight junction                                    | 9/97      | 165/7558 | 0.0002499 | 0.0224924 | 0.0206509 |     9 | 8       |
    | hsa04510 | Focal adhesion                                    | 9/97      | 198/7558 | 0.0009458 | 0.0340474 | 0.0312599 |     9 | 8       |
    | hsa05022 | Pathways of neurodegeneration - multiple diseases | 9/97      | 443/7558 | 0.1141696 | 0.5182323 | 0.4758039 |     9 | 8       |


The result is a dataframe where KEGG categories have been considered
as enriched with the following columns :

- `ID` : Unique identifier of the KEGG category
- `Description` : Description of the KEGG category
- `GeneRatio` : Fraction representing the number of marker genes present in
  the KEGG category, $GeneRatio = \frac{nbrMarkerGeneInKEGGcat}{nbrMarkerGene}$
- `BgRatio` : Fraction representing the number of genes of the reference
  present in the KEGG category,
  $BgRatio = \frac{nbrTotalGeneInKEGGcat}{nbrTotalGene}$
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
