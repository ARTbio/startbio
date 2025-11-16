## Script used to generated the tutorial scRNAseq for IOC in STARTbio

###############################################################################
##                               LOAD LIBRARIES                              ##
###############################################################################

## Import packages
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


###############################################################################
##                                ADD FUNCTIONS                              ##
###############################################################################

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

###############################################################################
##                                 IMPORT DATA                               ##
###############################################################################

## Import expression matrix
system("wget -P ./test-data/ https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
system("tar -zxvf ./test-data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C ./test-data/")
tenX_matrix <- Read10X(data.dir = "./test-data/filtered_gene_bc_matrices/hg19", #Path to the directory containing the three 10X files (matrix.mtx, genes.tsv (or features.tsv depending on the CellRanger version), barcodes.tsv)
                       gene.column = 1) #Choice of the genes.tsv column that will determine the names of the genes in the expression matrix

## Matrix dimension : first value is the number of row (genes) and second value is the number of column (sample/barcodes)
dim(tenX_matrix)

## Visualize first three rows and columns
tenX_matrix[1:3, 1:3]

## Import of the biomart database for hg19
ensembl_hg19 <- useEnsembl(biomart = "genes",                  #Import ensembl genes database
                           dataset = "hsapiens_gene_ensembl",  #Genome
                           GRCh = 37)                          #Genome version, only 38 or 37 are accepted for now

## If you don't know what to give to biomart parameter of useEnsembl just run :
#listEnsembl() #dataframe where the first column is the value to give to biomart parameter of useEnsembl

## If you don't know the available datasets :
#listDatasets(mart = useEnsembl(biomart = "genes")) #dataframe where the first column (dataset) contains the value to give to dataset parameter of useEnsembl, add GRCh parameter in the useEnsembl function to precise your preferred genome version

## Recovering attributes according to our gene list (genes present in our expression matrix)
annotated_hg19 <- getBM(attributes = c("ensembl_gene_id",
                                       "external_gene_name",
                                       "description",
                                       "gene_biotype",
                                       "chromosome_name"),    #Informations to retrieve
                           filters = "ensembl_gene_id",       #Which variable to choose for the filtering
                           values = rownames(tenX_matrix),    #Values that will filter the database
                           mart = ensembl_hg19)               #Database

## If you need to know all available attributes and their name
#listAttributes(mart = ensembl_hg19) #dataframe where the first column is the attribute' names needed for the attributes parameter of getBM

## Preview of the resulting dataframe
kable(head(annotated_hg19), "simple")

## Creation of the Seurat Object
pbmc_small <- CreateSeuratObject(tenX_matrix,                    #Expression matrix
                                 project = "PBMC analysis",      #Name of the project : something meaningful for your dataset
                                 assay = "RNA",                  #Name of the initial assay, (others can be created during the analysis), default is RNA
                                 names.field = 1,                #Associated with the names.delim, it allows to separate the name of the barcode when this one is composed of several information ex: BARCODE_CLUSTER_CELLTYPE and to choose after split on the names.delim which part we choose as the barcode name
                                 names.delim = "_",              #Character that allows to dissociate the different information contained in the barcode names
                                 meta.data = NULL,               #We can add the metadata on the transcriptomes with a dataframe where the barcodes are in line and the different information in column
                                 min.cells = 0,                  #Filtering genes that are not detected in less than min.cells
                                 min.features = 1)               #Filtering cells that do not detect at least min.features, here we filter all barcodes that detect no gene

pbmc_small #Small presentation of the Seurat object in R console

## Discovery of SeuratObject
str(pbmc_small)

dim(pbmc_small@assays$RNA@counts)
dim(pbmc_small@assays$RNA@data)

#Preview of the cell metadata
head(pbmc_small@meta.data)

###############################################################################
##                              QUALITY CONTROLS                             ##
###############################################################################

## Retrieve genes from the MT genome using biomart
genes_MT <- annotated_hg19$ensembl_gene_id[annotated_hg19$chromosome_name == "MT"]


pbmc_small <- PercentageFeatureSet(pbmc_small,
                                   features = genes_MT,        #Vectors of gene names present on the MT genome
                                   col.name = "percent_mito",  #Defines the name of the new column generated in the metadata of the Seurat object
                                   assay = "RNA")

## Violin plot of QC (Quality Controls)
VlnPlot(object = pbmc_small,
        features = c("nCount_RNA", "nFeature_RNA", "percent_mito"), #Parameters to plot (either gene expression or continuous variable in cell metadata)
        ncol = 3,                                                   #Number of columns if several figures are to be ploted
        pt.size = 0.01)                                             #Point size

## Graphical representation of QC
ggplot(pbmc_small@meta.data,
       aes(y = nCount_RNA,
           x = nFeature_RNA,
           color = percent_mito)) +
  geom_point() +
  geom_hline(yintercept = 650, linetype = 'dotted') +
  geom_vline(xintercept = 300, linetype = 'dotted') +
  scale_y_log10() +
  scale_color_gradient2(low = "green",
                        high = "red",
                        mid = "yellow",
                        midpoint = 20) +
  ggtitle("QC plot", "Number of detected genes in function of number of UMI")+
  labs(y = "Number of UMI per cell", x = "Number of detected genes by cell")

hist(pbmc_small$nCount_RNA,
     breaks = 100,
     xlab = "Number of UMI per cell",
     main = "")
abline(v = 650, col = "red")
abline(v = 10000, col = "red")

hist(pbmc_small$nFeature_RNA,
     breaks = 100,
     xlab = "Number of detected genes by cell",
     main = "")
abline(v = 300, col = "red")
abline(v = 2300, col = "red")

## Filtering SeuratObject
pbmc_small <- subset(pbmc_small,
                     percent_mito < 10 &
                       nCount_RNA > 650 &
                       nCount_RNA < 10000 &
                       nFeature_RNA > 300 &
                       nFeature_RNA < 2300)

## Plot
ggplot(pbmc_small@meta.data,
               aes(x = nCount_RNA,
                   y = nFeature_RNA,
                   color = percent_mito)) +
  geom_point() +
  scale_y_log10() +
  scale_color_gradient2(low = "green",
                        high = "red",
                        mid = "yellow",
                        midpoint = 20) +
  ggtitle("QC plot after filtering", "Number of detected genes in function of number of UMI")+
  labs(y = "Number of UMI per cell", x = "Number of detected genes by cell")

## Update object in R console
pbmc_small


###############################################################################
##                             CELL NORMALIZATION                            ##
###############################################################################

## Inter-cell normalization
pbmc_small <- NormalizeData(pbmc_small,                                   #SeuratObject
                            assay = "RNA",                                #Assay to use
                            normalization.method = "LogNormalize",        #Normalization method
                            scale.factor = median(pbmc_small$nCount_RNA), #Scale factor
                            verbose = TRUE)


###############################################################################
##                         HIGHLY VARIABLE GENES                             ##
###############################################################################

pbmc_small <- FindVariableFeatures(pbmc_small,                 #SeuratObject
                                   selection.method = "vst",   #Method
                                   nfeatures = 2000)           #Top HVG (Highly Variable Gene), default value
pbmc_small

## Plot
VariableFeaturePlot(pbmc_small)

###############################################################################
##                         REDUCTION OF DIMENSIONALITY                       ##
###############################################################################

## Scaling data
pbmc_small <- ScaleData(pbmc_small)

## PCA (Principal Component Analysis)
pbmc_small <- RunPCA(pbmc_small,                 #SeuratObject
                     reduction.name = "pca",     #Name of the reduction stored in the reduction slot
                     npcs = 50,                  #Total Number of PCs to compute and store (50 by default)
                     seed.use = 42,              #Set a random seed. By default, sets the seed to 42.
                     verbose = TRUE)

## Graphic representation of cells
PCAPlot(pbmc_small,                              #SeuratObject
        dims = c(1, 2))                          #Dimensions (PCs) to plot, default is the first two

## PC Selection

## JackStraw : Determine statistical significance of PCA scores
pbmc_small <- JackStraw(pbmc_small,          #SeuratObject
                        reduction = "pca",   #Reduction to analyse
                        dims = 50,           #Number of dimension to analyse
                        assay = "RNA")       #Assay to use

## Compute Jackstraw scores significance.
pbmc_small <- ScoreJackStraw(pbmc_small,     #SeuratObject
                             dims = 1:50)    #Number of dimension to analyse

## Representation JackStrow
JackStrawPlot(pbmc_small,                    #SeuratObject
              dims = 1:50)                   #Number of dimension to plot

## Elbow Plot
ElbowPlot(pbmc_small,           #SeuratObject
          ndims = 20,           #Number of dimension to analyse
          reduction = "pca")    #Reduction to analyse

## UMAP
pbmc_small <- RunUMAP(pbmc_small,               #SeuratObject
                      reduction = "pca",        #Reduction used to compute UMAP
                      reduction.key = "UMAP_",  #Dimension prefix
                      assay = "RNA",            #Assay to use
                      dims = 1:10)              #Number of PCs to keep (previously determined)

## Plot
UMAPPlot(pbmc_small)

###############################################################################
##                                 CLUSTERING                                ##
###############################################################################

## Create Graph
pbmc_small <- FindNeighbors(pbmc_small,          #SeuratObject
                            reduction = "pca",   #Reduction to used
                            k.param = 20,
                            dims = 1:10)         #Number of PCs to keep (previously determined)

## Create partition based on graph
pbmc_small <- FindClusters(pbmc_small,                                        #SeuratObject
                           resolution = seq(from = 0.2, to = 1.2, by = 0.2),  #Compute clustering with several resolutions (from 0.2 to 1.2 : values usually used)
                           verbose = FALSE)

#Preview of the cell metadata
head(pbmc_small@meta.data)

## Which resolution to choose ?
clustree(pbmc_small,               #SeuratObject
         prefix = "RNA_snn_res.")  #Prefix that retrieve all resolution to analyse in cell metadata slot

## Set the default resolution level
Idents(pbmc_small) <- "RNA_snn_res.0.4"

## Visualize clusters in UMAP coordinates
UMAPPlot(pbmc_small,          #SeuratObject
         label = TRUE,        #Plot label on the plot
         label.size = 4)      #Change label size


###############################################################################
##                           MARKERS IDENTIFICATION                          ##
###############################################################################

FeaturePlot(pbmc_small,                     #SeuratObject
            features = "ENSG00000156738",   #Value to plot, can be a vector of several variable
            reduction = "umap",             #Dimensional reduction to use
            label = TRUE,                   #Plot label on the plot
            label.size = 4) +               #Change label size
  ggtitle(annotated_hg19[annotated_hg19$ensembl_gene_id == "ENSG00000156738", "external_gene_name"],
          "ENSG00000156738")

VlnPlot(pbmc_small,                         #SeuratObject
        features = "ENSG00000156738") +     #Variable to plot
  ggtitle(annotated_hg19[annotated_hg19$ensembl_gene_id == "ENSG00000156738", "external_gene_name"],
          "ENSG00000156738")


###############################################################################
##                       DIFFERENTIAL EXPRESSION ANALYSIS                    ##
###############################################################################

pbmc_markers <- FindAllMarkers(pbmc_small,              #SeuratObject
                               only.pos = FALSE,        #Returns positive and negative gene markers
                               min.pct = 0.1,           #Take into account genes that are detected in at least 10% of the cells
                               logfc.threshold = 0,     #Return markers with a logFC superior to threshold
                               test.use = "wilcox",     #Method used
                               verbose = FALSE)

## Preview of the resulting dataframe
head(pbmc_markers)


###############################################################################
##                             MARKERS ANNOTATIONS                           ##
###############################################################################

## Via biomart
## Merge markers results with biomart annotation
pbmc_markers_annotated <- merge(x = pbmc_markers,         #First df to merge
                                y = annotated_hg19,       #Second df to merge
                                by.x = "gene",            #Column name of first df used for matching lines
                                by.y = "ensembl_gene_id", #Column name of second df used for matching lines
                                all.x = TRUE)             #Keep all lines from first df even if there is no match with second df

## Filter out non significant results
pbmc_markers_signif <- subset(pbmc_markers_annotated,
                              p_val_adj < 0.05 &
                                abs(avg_log2FC) >= 0.25)       #Filter dataframe based on p_val_adj column

## Number of significative DEG per cluster
table(pbmc_markers_signif$cluster)

## Sorting results by cluster and by average log2(Fold Change)
pbmc_markers_signif <- pbmc_markers_signif %>%                 #Rearrange df with dplyr package
  group_by(cluster) %>%                                        #Group df based on cluster column
  arrange(desc(avg_log2FC), .by_group = TRUE)                  #Sort lines by descending the column avg_log2FC and by group

## Most DE gene marker for each cluster
top_n(x= pbmc_markers_signif, n = 3, wt = avg_log2FC)

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

## Via ClusterProfilter
# ORA

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

enrich_go <- do.call("rbind", enrich_go_list)                            #Bind all results together
enrich_go <- enrich_go %>%
  group_by(cluster, ONTOLOGY)                                            #Rearrange df according to cluster and ontology type


## show first results
kable(top_n(x= enrich_go, n = 3, wt = (Count))[,-9])               #Only remove list of genes for the visualisation

## KEGG
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

## Concatenate all results in one dataframe
enrich_kegg <- do.call("rbind", enrich_kegg_list)

## Group result by cluster (easier to manipulate with dplyr)
enrich_kegg <- enrich_kegg %>%
  group_by(cluster)

## Visualise first 3 KEGG categories for each cluster (removing the vector of genes just for the visualisation)
top_n(x= enrich_kegg, n = 3, wt = (Count))[,-8]

# GSEA
## Retrieve MSigDB Database for human cell types signatures gene sets
C8_t2g <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, ensembl_gene)

  head(C8_t2g, 10)

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

## Concatenate all results in one dataframe
GSEA_res <- do.call("rbind", GSEA_list)

## Group result by cluster (easier to manipulate with dplyr)
GSEA_res <- GSEA_res %>%
  group_by(cluster)

## Visualize first 3 signatures for each cluster (removing the vector of genes just for the visualization and the description that match ID column for this dataset MSigDB)
top_n(x= GSEA_res, n = 3, wt = NES)[, -c(2,11)]

## CellMarkers
## Retrieve Cell Markers Database for human cell types signatures gene sets
cell_marker_data <- vroom::vroom('http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt')

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

head(cell_marker_data, 10)

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
top_n(x= GSEA_CM_res, n = 3, wt = NES)[, -c(2,11)]

###############################################################################
##                           CLUSTERS ANNOTATION                             ##
###############################################################################

data.frame(Marker = c("IL7R, CCR7", "CD14, LYZ", "IL7R, S100A4", "MS4A1", "CD8A", "FCGR3A, MS4A7", "GNLY, NKG7", "FCER1A, CST3", "PPBP"),
                 CellType = c("Naive CD4+ T", "CD14+ Mono", "Memory CD4+", "B cells", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"),
                 Cluster = 0:8)

## Retrieve specific markers based on their gene name
markers_pop <- subset(annotated_hg19,
                      subset = external_gene_name %in%
                        c("CCR7",                                #Naive CD4+ T
                          "CD14",                                #CD14+ Mono
                          "S100A4",                              #Memory CD4+
                          "MS4A1",                               #B
                          "CD8A",                                #CD8+ T
                          "FCGR3A",                              #FCGR3A+ Mono
                          "NKG7",                                #NK
                          "CST3",                                #DC
                          "PPBP"))                               #Platelet

## Generate violin plots and stock them into variable
vln_plots <- VlnPlot(pbmc_small,                                 #SeuratObject
                     features = markers_pop$ensembl_gene_id,     #Vector of genes to plot
                     combine = FALSE)                            #Return list of plot objets instead of a combined plot (easier to process)

## Add gene name as title
vln_plots <- lapply(vln_plots,                                   #List of plots
                    add_title_gene_name,                         #Function to apply to the list
                    gene_format = annotated_hg19)                #Fill in the dataframe parameter

## Remove unecessary legend
vln_plots <- lapply(vln_plots, function(plot){
  plot <- plot + theme(legend.position = "none")
  return(plot)
})

## Plot list of plots
grid.arrange(grobs = vln_plots)

## Vector of new cluster labels         #Correspond to cluster :
new_cluster_ids <- c("Naive CD4+ T",    #0
                     "CD14+ Mono",      #1
                     "Memory CD4+",     #2
                     "B",               #3
                     "CD8+ T",          #4
                     "FCGR3A+ Mono",    #5
                     "NK",              #6
                     "DC",              #7
                     "Platelet")        #8

## Create a named vector with the actual cell identifiers
names(new_cluster_ids) <- levels(pbmc_small)

## Renamed cell identities in the Seurat Object
pbmc_small <- RenameIdents(pbmc_small, new_cluster_ids)

## Plot
UMAPPlot(pbmc_small,
         label = TRUE,
         pt.size = 0.5)


###############################################################################
##                         COMPARING TWO POPULATIONS                         ##
###############################################################################

## Perform differential expression analysis
NK_CD8_diff_markers <- FindMarkers(pbmc_small,
                           ident.1 = "NK",
                           ident.2 = "CD8+ T")

## Merge markers results with biomart annotation
NK_CD8_diff_markers_annotated <- merge(x = NK_CD8_diff_markers,  #First df to merge
                                       y = annotated_hg19,       #Second df to merge
                                       by.x = 0,                 #Column name of first df used for matching lines, 0 for rownames
                                       by.y = "ensembl_gene_id", #Column name of second df used for matching lines
                                       all.x = TRUE)             #Keep all lines from first df even if there is no match with second df

## Filter dataset based on Fold change and p-value adjusted
NK_CD8_diff_markers_annotated_signif <- subset(NK_CD8_diff_markers_annotated,
                                               p_val_adj < 0.05 &
                                                 abs(avg_log2FC) >= 0.25)       #Filter dataframe based on p_val_adj column

## Sorting results by average log2(Fold Change)
NK_CD8_diff_markers_annotated_signif <- NK_CD8_diff_markers_annotated_signif %>%                 #Rearrange df with dplyr package
  arrange(desc(avg_log2FC))                  #Sort lines by descending the column avg_log2FC and by group

## Most DE gene marker for each cluster
NK_CD8_diff_markers_annotated_signif[(c(1:3, (nrow(NK_CD8_diff_markers_annotated_signif)-2):nrow(NK_CD8_diff_markers_annotated_signif))),]


###############################################################################
##                               VISUALIZATIONS                              ##
###############################################################################

## Visualize cells in UMAP coordinates where cells are colored by a certain clustering
UMAPPlot(pbmc_small,                                                           #SeuratObject
         group.by = "RNA_snn_res.0.4")                                         #Color cells based on different cell metadata

## Visualize cells in UMAP coordinates where cells are colored by two kind of variable separetely
UMAPPlot(pbmc_small,                                                           #SeuratObject
         group.by = c("RNA_snn_res.0.4", "RNA_snn_res.1.2"))                   #Color cells based on different cell metadata

## Visualize cells in UMAP coordinates where cells are splitted in different panels based on a variable
UMAPPlot(pbmc_small,                                                           #SeuratObject
         split.by = "RNA_snn_res.0.2")                                         #Separated cells based on a cell metadata variable

## Visualize cells in UMAP coordinates and adding cluster labels directly on the plot
UMAPPlot(pbmc_small,                                                           #SeuratObject
         label = TRUE,                                                         #Print cell identities directly on the plot
         repel = TRUE)                                                         #Avoid overlap of cell labels

## Dotplot to visualize target genes expression in the different cell identities
DotPlot(pbmc_small,                                                            #SeuratObject
        features = markers_pop$ensembl_gene_id,                                #Feature expression to plot
        cols = c("yellow", "red")) +                                           #Change expression color scale
  scale_x_discrete(labels = markers_pop$external_gene_name)                    #Change labels to print gene names instead of ensembl gene id

## Heatmap
DoHeatmap(pbmc_small,                                                          #SeuratObject
          features = markers_pop$ensembl_gene_id)                              #Feature expression to plot

## Visualize cells in UMAP coordinates where cells are colored by a continuous variable (here two expression genes)
FeaturePlot(pbmc_small,                                                        #SeuratObject
            features = c("ENSG00000126353", "ENSG00000168685"),                #Feature expression to plot
            cells = colnames(subset(pbmc_small, idents = "Naive CD4+ T")),     #Plot only Naive CD4+ T cells
            cols = c("white", "orange", "darkblue"),                           #Change color for the blend : first color : no expression, 2nd : expressed first gene, 3rd color : expressed gene 2
            blend = TRUE)                                                      #See the coexpression of the two genes


###############################################################################
##                                 SessionInfo                               ##
###############################################################################

sessionInfo()
