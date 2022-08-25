# GSEA : Gene Set Enrichment Analysis

Enrichment analyses [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) are
based on *ranking* the genes. In order to take into account both the
direction of the deregulation and its significance we order the genes
according to :

$$
x = sign(avg.log2FC) * -log10(pval)
$$

Knowing that the `sign` function in $R$ returns the sign of the average
log2(FC). That is, when the genes are under expressed the value returned
by the function is `-1`, `+1` when the FC is positive and `0` when the gene
is not deregulated.

Taking into account the significance of the deregulation sometimes leads
to complications. Indeed, it happens that the function `FindAllMarkers`
returns p-values so low that they become zero, which produces `Inf` values
that cannot be processed by GSEA. To overcome this problem and only for
this case, we replace $-log10(pval)$ by $-log10(1e-323)$ because `1e-323`
would be the smallest value that could be represented by a computer.

``` r
print(1e-323) ## equal 9.881313e-324
```

\[1\] 9.881313e-324

``` r
print(1e-324) ## equal 0
```

\[1\] 0

GSEA will calculate an enrichment score for each gene set analysed from
this vector containing genes ordered according to their significance.

For each signature (and for each cluster), GSEA will calculate the
enrichment score by running the vector of ordered genes, increasing
the score if it encounters a gene from the gene set being analysed and
decreasing it if it encounters a gene not in the gene set. The enrichment
score (ES) is the maximum value during the increment.

The results of the `GSEA` function will be a dataframe with the
following columns:

- `ID ` : Identifier of the gene group being analysed
- `setSize` : Size of the gene group
- `EnrichmentScore` (ES): Enrichment score representing the degree of
  presence of the gene set in the ordered list of genes
- `NES` (Normalized Enrichment Score) : Normalized Enrichment Score such that :
  $NES = \frac{actual ES}{mean(ESs Against All Permutations Of The Dataset)}$
- `p-value` : p-value of the enrichment test
- `p.adjust` : adjusted p-value of the Benjamini Hochberg test
- `qvalue` : q-value after FDR (False Discovery Rate) control
- `rank` : Position in the list of genes for which ES is reached
- `leading_edge` : Three statistics calculated during the analysis:
    - `Tags` : Percentage of genes in the gene set before or after the
      ES peak depending on whether it is positive or negative
    - `List` : Percentage of genes before or after the ES peak that are
      positive or negative
    - `Signal` : Strength of the enrichment signal calculated:
      $(Tag)(1 - List)(\displaystyle \frac{N}{N - Nh})$
- `cluster` : Name of the cluster for which the gene set has been
  identified as significant

The format of the gene sets used for the `GSEA` function must be a mapping
table (`TERM2GENE`) which associates the name of a signature with the genes
which compose it. However, this is not like the `gmt` format where one row
corresponds to one signature, here there is one row for each possible
(signature/gene) pair. We will use the gene sets from the [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) database where
we will focus on the molecular signatures of cell types (Category `C8`)
and [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/) another
database which provides us with lists of specific genes for cell types.

The results will be discussed once we have completed the analysis for
both databases.

## Molecular Signature Database (MSigDB)

The Broad Institute's MSigDB database contains several collections of
gene signatures:

- `H` or *HallMarks gene sets*: A set of gene sets that co-express in
  identified biological processes or states with respect to other
  collections
- `C1` or *Positional gene sets* : A set of gene sets based on their
  cytogenetic and chromosomal position
- `C2` or *Curated gene sets* : A set of gene sets found in databases
  and in the scientific literature
    - Biocarta
    - KEGG
    - PID
    - Reactome
    - WikiPathways
- `C3` or *Regulatory target gene sets* : Set of potential microRNA
  target genes (MIR) or transcription factors (TFT)
    - MIR: miRDB prediction
    - TFT: prediction based on the work of Kolmykov et al. 2021 and
      Xie et al. 2005
- `C4` or *Computational gene sets* : Gene set based on two rather
  cancer-oriented microarray papers (Subramanian, Tamayo et al. 2005
  and Segal et al. 2004) that generated over 800 gene sets.
- `C5` or *Ontology gene sets* : Gene sets based on ontology databases.
    - Gene Ontology (GO): MF, CC, BP
    - Human Phenotype Ontology (HPO)
- `C6` or *Oncogenic signature gene sets* : A set of genes based on
  microarray results mainly concerning pathways that are often
  deregulated in cancer
- `C7` or *Immunologic signature gene sets* : Gene sets based on databases
  of the immune system and its possible perturbations.
    - ImmuneSigDB (human + mouse)
    - VAX: cured by the Human Immunology Project Consortium (HIPC)
- `C8` or *Cell type signature gene sets* : A set of gene sets corresponding
  to cell type markers defined mainly in single cell analysis

The R package `msigdbr` allows to query the database directly. With the
`msigdbr` function we select the `C8` collection in order to obtain the
gene sets corresponding to human cell lines.

We will then apply the `GSEA` function to each cluster using a
`lapply`. For each cluster :

- We filter the dataframe of the unfiltered markers result to get only
  the lines concerning the genes deregulated by the cells of the cluster.
- Run the `GSEA` function to get the enriched signatures in the marker
  gene set
- Add the cluster name in a new column to the resulting dataframe
- We visualize the results with the function `gseaplot2` of the package
  `enrichR` which allows to visualize the first 3 signatures

``` r
## Retrieve MSigDB Database for human cell types signatures gene sets
C8_t2g <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, ensembl_gene)

  kable(head(C8_t2g, 10))
  ```

  | gs_name                          | ensembl_gene    |
  |:---------------------------------|:----------------|
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000123146 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000205336 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000106948 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000132965 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000244509 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000239713 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000163219 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000186517 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000180448 |
  | AIZARANI_LIVER_C1_NK_NKT_CELLS_1 | ENSG00000123329 |

``` r
## Apply GSEA for each cluster
GSEA_list <- lapply(levels(pbmc_markers_annotated$cluster), function(cluster_name){

  res_markers <- subset(pbmc_markers_annotated, cluster == cluster_name)                   #Filter markers dataframe by cluster

  ## Generate named vector of ranked gene mandatory for GSEA analysis that take into account DE and significativity
  geneList_byclus <- sign(res_markers$avg_log2FC) * -log10(ifelse(res_markers$p_val == 0,  #Deal with pval = 0
                                                                  1e-323,                  #Smallest interpretable number
                                                                  res_markers$p_val))
  names(geneList_byclus) <- res_markers$gene

  ## Order by avg log FC and significativity
  geneList_byclus <- sort(geneList_byclus, decreasing = TRUE)

  ## Perform GSEA analysis
  gseaC8 <- GSEA(geneList_byclus, TERM2GENE = C8_t2g)
  gseaC8@result$cluster <- cluster_name #add cluster name as column

  ## Add plot
  # print(ridgeplot(gseaC8,
  #                 showCategory = 3,
  #                 orderBy = "NES") +
  #         ggtitle(paste("Cluster", cluster_name)) +
  #         theme(axis.text.y = element_text(size = 10),
  #               legend.key.size = unit(0.2, 'cm')))

  print(gseaplot2(gseaC8,
                  geneSetID = rownames(gseaC8@result %>%
                                         arrange(desc(NES)))[1:ifelse(nrow(gseaC8) < 3,
                                                                      nrow(gseaC8),
                                                                      3)],
                  base_size = 8,
                  pvalue_table = TRUE,
                  subplots = 1:2,
                  title = paste("Cluster", cluster_name)))

  return(gseaC8@result) #Return dataframe result
})
```

<img src="../images/clusterGSEA-1.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-2.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-3.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-4.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-5.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-6.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-7.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-8.png" style="display: block; margin: auto;" /><img src="../images/clusterGSEA-9.png" style="display: block; margin: auto;" />

``` r
## Concatenate all results in one dataframe
GSEA_res <- do.call("rbind", GSEA_list)

## Group result by cluster (easier to manipulate with dplyr)
GSEA_res <- GSEA_res %>%
  group_by(cluster)

## Visualize first 3 signatures for each cluster (removing the vector of genes just for the visualization and the description that match ID column for this dataset MSigDB)
kable(top_n(x= GSEA_res, n = 3, wt = NES)[, -c(2,11)])
```

??? abstract "First Enriched MSigDB gene sets for each cluster"
    | ID                                                      | setSize | enrichmentScore |      NES |  pvalue |  p.adjust |   qvalues | rank | leading_edge                   | cluster |
    |:--------------------------------------------------------|--------:|----------------:|---------:|--------:|----------:|----------:|-----:|:-------------------------------|:--------|
    | HAY_BONE_MARROW_NAIVE_T\_CELL                           |     208 |       0.9209917 | 3.046758 | 0.0e+00 | 0.0000000 | 0.0000000 |  216 | tags=76%, list=17%, signal=75% | 0       |
    | RUBENSTEIN_SKELETAL_MUSCLE_T\_CELLS                     |     153 |       0.8553262 | 2.739346 | 0.0e+00 | 0.0000000 | 0.0000000 |  146 | tags=67%, list=11%, signal=67% | 0       |
    | TRAVAGLINI_LUNG_CD4_NAIVE_T\_CELL                       |     116 |       0.9155134 | 2.824719 | 0.0e+00 | 0.0000000 | 0.0000000 |  135 | tags=75%, list=11%, signal=74% | 0       |
    | AIZARANI_LIVER_C23_KUPFFER_CELLS_3                      |     142 |       0.8455834 | 2.003393 | 0.0e+00 | 0.0000000 | 0.0000000 |  257 | tags=68%, list=16%, signal=62% | 1       |
    | DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_DENDRITIC_CELLS |     101 |       0.8823161 | 2.065142 | 0.0e+00 | 0.0000000 | 0.0000000 |  192 | tags=72%, list=12%, signal=68% | 1       |
    | TRAVAGLINI_LUNG_CLASSICAL_MONOCYTE_CELL                 |     168 |       0.8365088 | 1.997979 | 0.0e+00 | 0.0000000 | 0.0000000 |  265 | tags=69%, list=16%, signal=64% | 1       |
    | HAY_BONE_MARROW_NAIVE_T\_CELL                           |     230 |       0.7672855 | 3.015408 | 0.0e+00 | 0.0000000 | 0.0000000 |  260 | tags=67%, list=25%, signal=64% | 2       |
    | RUBENSTEIN_SKELETAL_MUSCLE_T\_CELLS                     |     145 |       0.7263295 | 2.725064 | 0.0e+00 | 0.0000000 | 0.0000000 |  192 | tags=58%, list=19%, signal=55% | 2       |
    | TRAVAGLINI_LUNG_CD4_NAIVE_T\_CELL                       |     107 |       0.7627843 | 2.748389 | 0.0e+00 | 0.0000000 | 0.0000000 |  180 | tags=64%, list=17%, signal=59% | 2       |
    | AIZARANI_LIVER_C34_MHC_II_POS_B\_CELLS                  |      83 |       0.9034082 | 2.063760 | 0.0e+00 | 0.0000000 | 0.0000000 |  152 | tags=80%, list=13%, signal=75% | 3       |
    | FAN_EMBRYONIC_CTX_BRAIN_B\_CELL                         |      70 |       0.9140261 | 2.070592 | 0.0e+00 | 0.0000000 | 0.0000000 |  103 | tags=66%, list=9%, signal=64%  | 3       |
    | TRAVAGLINI_LUNG_B\_CELL                                 |     112 |       0.9092172 | 2.119503 | 0.0e+00 | 0.0000000 | 0.0000000 |  121 | tags=68%, list=10%, signal=67% | 3       |
    | AIZARANI_LIVER_C1_NK_NKT_CELLS_1                        |      70 |       0.8838737 | 2.284841 | 0.0e+00 | 0.0000000 | 0.0000000 |   90 | tags=64%, list=13%, signal=62% | 4       |
    | HAY_BONE_MARROW_NK_CELLS                                |     105 |       0.8450920 | 2.260853 | 0.0e+00 | 0.0000000 | 0.0000000 |  132 | tags=68%, list=18%, signal=65% | 4       |
    | TRAVAGLINI_LUNG_CD8_NAIVE_T\_CELL                       |      92 |       0.8842939 | 2.338596 | 0.0e+00 | 0.0000000 | 0.0000000 |   90 | tags=62%, list=13%, signal=62% | 4       |
    | AIZARANI_LIVER_C31_KUPFFER_CELLS_5                      |      53 |       0.8510682 | 1.948702 | 0.0e+00 | 0.0000000 | 0.0000000 |  191 | tags=68%, list=10%, signal=63% | 5       |
    | HAY_BONE_MARROW_MONOCYTE                                |     191 |       0.8332887 | 2.001325 | 0.0e+00 | 0.0000000 | 0.0000000 |  295 | tags=73%, list=16%, signal=68% | 5       |
    | TRAVAGLINI_LUNG_NONCLASSICAL_MONOCYTE_CELL              |     163 |       0.8581938 | 2.050371 | 0.0e+00 | 0.0000000 | 0.0000000 |  301 | tags=81%, list=16%, signal=75% | 5       |
    | DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_NK_CELLS        |      72 |       0.8853806 | 2.022495 | 0.0e+00 | 0.0000000 | 0.0000000 |  131 | tags=71%, list=13%, signal=66% | 6       |
    | HAY_BONE_MARROW_NK_CELLS                                |     236 |       0.8420024 | 1.990762 | 0.0e+00 | 0.0000000 | 0.0000000 |  196 | tags=56%, list=19%, signal=59% | 6       |
    | TRAVAGLINI_LUNG_NATURAL_KILLER_CELL                     |     108 |       0.8994932 | 2.104728 | 0.0e+00 | 0.0000000 | 0.0000000 |  124 | tags=69%, list=12%, signal=67% | 6       |
    | HAY_BONE_MARROW_DENDRITIC_CELL                          |      96 |       0.8483908 | 2.413478 | 0.0e+00 | 0.0000000 | 0.0000000 |   99 | tags=49%, list=6%, signal=49%  | 7       |
    | DESCARTES_FETAL_LUNG_MYELOID_CELLS                      |      62 |       0.8072714 | 2.240695 | 0.0e+00 | 0.0000000 | 0.0000000 |  145 | tags=55%, list=8%, signal=52%  | 7       |
    | TRAVAGLINI_LUNG_PLASMACYTOID_DENDRITIC_CELL             |      41 |       0.8240860 | 2.252436 | 8.0e-07 | 0.0000107 | 0.0000073 |  128 | tags=51%, list=7%, signal=49%  | 7       |
    | DESCARTES_FETAL_LIVER_MEGAKARYOCYTES                    |      55 |       0.8328324 | 1.392413 | 2.1e-06 | 0.0000795 | 0.0000704 |   83 | tags=56%, list=13%, signal=54% | 8       |
    | DESCARTES_FETAL_ADRENAL_MEGAKARYOCYTES                  |      43 |       0.8505604 | 1.416959 | 4.6e-06 | 0.0001235 | 0.0001092 |   74 | tags=63%, list=11%, signal=60% | 8       |
    | DESCARTES_FETAL_KIDNEY_MEGAKARYOCYTES                   |      32 |       0.8561275 | 1.417434 | 5.7e-05 | 0.0010710 | 0.0009475 |   69 | tags=66%, list=11%, signal=62% | 8       |


## CellMarkers

CellMarker is a hand-curated database of scientific literature and
other resources to describe over 400 cell types (human and mouse only).
Human cell markers will be retrieved directly from the Cell Marker site.

After manipulating the data to format a two column dataframe where the
first column is the term name and the second column is the gene name
associated with that term.

We will then apply the `GSEA` function for each of the clusters using
a `lapply`. For each cluster :

- We filter the unfiltered marker result dataframe to retrieve only the
  rows concerning the genes deregulated by the cluster cells.
- Run the `GSEA` function to get the enriched signatures in the marker
  gene set
- Add the cluster name in a new column to the resulting dataframe
- We visualize the results with the function `gseaplot2` of the package
  `enrichR` which allows to visualize the first 5 signatures

``` r
## Retrieve Cell Markers Database for human cell types signatures gene sets
cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt')

## Instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
    dplyr::select(cellName, geneSymbol) %>%                            #Select only the two columns
    dplyr::mutate(geneSymbol = strsplit(geneSymbol, ', ')) %>%         #Split gene names based on the comma
    tidyr::unnest()                                                    #Flatten gene vector in order to have a line for each gene in terme

## Remove [ and ] found in gene names due to the Cell Marker annotation
cells$geneSymbol <- gsub("\\[|\\]",
                         "",
                         cells$geneSymbol,
                         fixed = FALSE)

paged_table(head(cell_marker_data, 30))

## Apply GSEA for each cluster
GSEA_CM_list <- lapply(levels(pbmc_markers_annotated$cluster), function(cluster_name){

  res_markers <- subset(pbmc_markers_annotated, cluster == cluster_name)                     #Filter markers dataframe by cluster

  ## Generate named vector of ranked gene mandatory for GSEA analysis that take into account DE importance and significativity
  geneList_byclus <- sign(res_markers$avg_log2FC) * -log10(ifelse(res_markers$p_val == 0,    #Deal with pval = 0
                                                                  1e-323,                    #Smallest interpretable number
                                                                  res_markers$p_val))
  names(geneList_byclus) <- res_markers$external_gene_name

  ## Order by avg log FC and significativity
  geneList_byclus <- sort(geneList_byclus, decreasing = TRUE)

  ## Perform GSEA analysis
  gseaCM <- GSEA(geneList_byclus, TERM2GENE = cells)
  gseaCM@result$cluster <- cluster_name #add cluster name as column

  ## Add plot
  # print(ridgeplot(gseaCM,
  #                 showCategory = 5,
  #                 orderBy = "NES") +
  #         ggtitle(paste("Cluster", cluster_name)) +
  #         theme(axis.text.y = element_text(size = 10),
  #               legend.key.size = unit(0.2, 'cm')))

  print(gseaplot2(gseaCM,
                  geneSetID = rownames(gseaCM@result %>%
                                         arrange(desc(NES)))[1:ifelse(nrow(gseaCM) < 3,
                                                                      nrow(gseaCM),
                                                                      3)],
                  base_size = 8,
                  pvalue_table = TRUE,
                  subplots = 1:2,
                  title = paste("Cluster", cluster_name)))

  return(gseaCM@result) #Return dataframe result
})

## Concatenate all results in one dataframe
GSEA_CM_res <- do.call("rbind", GSEA_CM_list)

## Group result by cluster (easier to manipulate with dplyr)
GSEA_CM_res <- GSEA_CM_res %>%
  group_by(cluster)

## Visualise first 3 signatures for each cluster (removing the vector of genes just for the visualisation and the description that match ID column for this dataset MSigDB)
paged_table(top_n(x= GSEA_CM_res, n = 3, wt = NES)[, -c(2,11)])
```

## Analyse des resultats d’enrichissement

The different GSEA plots represent the first three signatures for which
we had the highest enrichment scores. The increment of the score can be
followed on the top panel. The position of the signature genes is shown
on the bottom panel. For each of the three signatures the p-value and the
adjusted p-value are also shown in a table.

Many of the results confirm what was already observed with the
over-representation analyses:

- The **cluster 0** would be composed of native T cells (either CD4+
  for MSigDB or CD8+ for CellMarker)
- **cluster 1** would represent monocytes as well as **cluster 5** which
  is close on the UMAP however the enrichment analyses cannot detect their
  difference, it would be necessary to investigate the entire
  ClusterProfiler results to differentiate them
- **cluster 2** would be composed of CD4+ T cells considered as native
  for MSigDB (again cluster 2 cells are very close on the UMAP to those
  of cluster 0)
- **cluster 3** is always related to B cells
- **cluster 4** and **cluster 6** (very close on the UMAP) are still
  considered to be composed of NK cells. However, both databases also
  detected that cluster 4 could be composed of CD8+ T cells
- **cluster 7** cells would be dentritic cells.
- **cluster 8** would be composed of megakaryocytes which are defined
  as cells at the origin of platelet formation which confirms the
  previous results. However, CellMarker does not allow us to rule out
  this cluster since only one signature is significant and the enrichment
  of this one is negative.

Enrichment analyses are highly dependent on the signature whose
enrichment or over-representation is being measured. If the signatures
are too general or based on something too far from our dataset (here we
were on blood cells) then it will be difficult to get results that make
sense the first time. You will have to rely on your own signatures or
marker genes to identify your population.
