# Import Data

## Expression matrix

The first thing to import is the expression data. They are
generally in the form of 3 files:

- The name of the transcriptomes `barcodes.tsv` : file with one column with the
  transcriptome names (cells, barcodes)
- Gene names in `genes.tsv` format : three columns file
  (geneID of the used gtf annotation, gene name, "Gene Expression")
- the expression matrix `matrix.mtx` : three columns file
  (position of the transcriptome in the `barcodes.tsv` file, position of the of
  the gene in the file `genes.tsv`, associated expression level)

The `Read10X` function of Seurat allows to create from these files the
expression matrix in `dgCMatrix` format. It is a sparse matrix used for
efficient processing of large matrices.

!!! danger "Caution"
    A very important parameter of `Read10X` is `gene.column` which allows us
    to choose which column of the `genes.tsv` file to use. It is preferable
    to select the column composed only of gene ID during the analysis (column 1)
    and to use the gene names only during the annotation of the marker genes.
    Indeed column 2 of `genes.tsv` is in fact composed of gene name and
    *gene ID* (because there is not necessarily a gene name assigned to each
    unique identifier). Moreover, several *gene ID* can correspond to a single
    gene name making the annotation step very complex.
    If we choose column 2, duplicates in gene names are handled with the
    `make.unique` function which detects duplicates and then adds
    `.1, .2, ..., .n` after each occurrence (knowing that the first occurrence
    is not modified). With these modified gene names, it will be very difficult
    to use the different databases that will not be able to recognize textually
    these new gene names.

    We could think that it would be easy to find the genes whose names have
    been modified afterwards by searching all the genes that would have a dot
    in their names, unfortunately some gene names already contain dots making
    the search for patterns even more complex. To get rid of all these
    problems, I **strongly recommend** to use the first column and to use only
    gene names afterwards.


## Dataset test

For this usecase, the dataset used for the analysis is composed of peripheral 
blood mononuclear cells (PBMC). It is available on the 10X and Seurat website 
which can be found on the [Seurat tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). The first thing to do is to dowload the dataset directly from 10X website. 
We use the R function `system` to execute command lines that we normally run in a 
terminal. Then we use the function `Read10X` to read the bundle format that are the
direct output of CellRanger (with specific filenames). If you have another type of
bundle format, you might be interessed by the function `ReadMtx`.

``` r
## Download expression matrix from 10X website in my subfolder "test_data"
system("wget -P ./test-data/ https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
## Extract archive in my subfolder "test_data"
system("tar -zxvf ./test-data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C ./test-data/")

## Import expression matrix in R global environment
tenX_matrix <- Read10X(data.dir = "./test-data/filtered_gene_bc_matrices/hg19", #Path to the directory containing the three 10X files (matrix.mtx, genes.tsv (or features.tsv depending on the CellRanger version), barcodes.tsv)
                       gene.column = 1) #Choice of the genes.tsv column that will determine the names of the genes in the expression matrix

## Matrix dimension : first value is the number of row (genes) and second value is the number of column (sample/barcodes)
dim(tenX_matrix)
```

Description of the `dgCmatrix` object:

- `i`: row position of each non-zero value of your
    expression matrix (knowing that the index of the row starts at
    zero)
- `p` : vector of size `number of barcode + 1` and is described as
    representing the number of non-zero values in each column (`j`)
    such that: `tenX_matrix@p[j+2] - tenX_matrix@p[j+1]`
- `Dim` : vector containing the dimensions of your expression matrix
    expression matrix (number of genes, rows then number of barcodes,
    columns)
- `Dimnames` : list containing the names of the genes (rownames) and the
    barcode names (colnames)
- `x` : vector containing the values of non-null expressions
- `factors`.

The most important thing to remember about this complex object is the slots
corresponding to the number of dimensions and their names (to know the
annotation of your genes and your barcodes) and `x` which allows to visualize
the different expression values. Useful to know if what you have in your hands
are integer values, so probably raw expression levels, or if you have decimal
values which would reflect normalized expression levels.

A `dgCMatrix` is displayed in the console just like a dataframe where the
points correspond to zero:

``` r
## Visualize first three rows and columns
tenX_matrix[1:3, 1:3]
```

    ## 3 x 3 sparse Matrix of class "dgCMatrix"
    ##                 AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1
    ## ENSG00000243485                .                .                .
    ## ENSG00000237613                .                .                .
    ## ENSG00000186092                .                .                .

!!! question "What if I don't have a bundle format ?"
    If the format of the expression matrix is not a three-file format but a
    file containing a table (*n* genes x *n* barcodes), we can go directly to 
    the next step. You just need to import it in your R global environment as 
    a `matrix` or `data.frame` for example. Don't forget that the expression 
    table must have barcodes as columns and genes as rows.

## biomaRt Annotation

We will then import the biomart database which will be used to annotate the
genes of our analysis. This will allow us to go from a gene ID set to a
gene name. The data of Seurat being in hg19 we recover the database for
this annotation.  

The biomaRt package allows us to retrieve the databases that are available
via [Ensembl BioMart](https://www.ensembl.org/biomart/martview/aeb467155c9aeccaa44a70171cde4e15)
and to interact directly in the R environment.
The first function we will use is `useEnsembl` to connect to Ensembl. You can
see all the versions and genomes available through the function of the same
package `listEnsembl`. The result of `useEnsembl` is an S4 object of class
Mart which contains the database of the genome and the desired version.
A second step will allow us to obtain a dataframe with the annotation of the
genes. The `getBM` function will perform a query via the Mart object that we
have previously generated. It needs to be provided with a vector containing
different attributes (chromosome, gene name, homologs, orthologs, and many
others...). We can have the exhaustive list with the `listAttributes` function.

Here we will also filter to get the annotation only for our list of genes with
the parameters *filters* and *values*.  

``` r
## Import of the biomart database for hg19
ensembl_hg19 <- useEnsembl(biomart = "genes",                  #Import ensembl genes database
                           dataset = "hsapiens_gene_ensembl",  #Genome
                           GRCh = 37)                          #Genome version, only 37 or NULL for 38 are accepted for now

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
```

| ensembl_gene_id | external_gene_name | description                                                                           | gene_biotype   | chromosome_name |
|:----------------|:-------------------|:--------------------------------------------------------------------------------------|:---------------|:----------------|
| ENSG00000000457 | SCYL3              | SCY1-like 3 (S. cerevisiae) \[Source:HGNC Symbol;Acc:19285\]                          | protein_coding | 1               |
| ENSG00000000460 | C1orf112           | chromosome 1 open reading frame 112 \[Source:HGNC Symbol;Acc:25565\]                  | protein_coding | 1               |
| ENSG00000000938 | FGR                | feline Gardner-Rasheed sarcoma viral oncogene homolog \[Source:HGNC Symbol;Acc:3697\] | protein_coding | 1               |
| ENSG00000000971 | CFH                | complement factor H \[Source:HGNC Symbol;Acc:4883\]                                   | protein_coding | 1               |
| ENSG00000001460 | STPG1              | sperm-tail PG-rich repeat containing 1 \[Source:HGNC Symbol;Acc:28070\]               | protein_coding | 1               |
| ENSG00000001461 | NIPAL3             | NIPA-like domain containing 3 \[Source:HGNC Symbol;Acc:25233\]                        | protein_coding | 1               |


We obtained a dataframe with the information about the genes present
in our expression matrix.

!!! note
    The `kable` function is only useful when you are coding in RMarkdown notebooks. 
    It helps improve the table display. If you code in the R console or with a Rscript
    you just need to execute `head(annotated_hg19)`.

# Creation of the Seurat object

Then thanks to the `CreateSeuratObject` function we can import the expression
matrix in a Seurat object which will be the basis of our analysis. All the
different steps will be associated to this object which is an S4 object
(tree of other smaller objects, dataframe, vectors...), specific to Seurat
called a `SeuratObject`.  

There are three important parameters:

- `project`: the name of your project which will also be the primary identity
   of your cells. It is therefore very important to define a project name and
   not to leave the default value.
- `min.cells` : allows to filter out genes that are not detected in at least
  `min.cells`.
- `min.features` : allows to filter the cells which do not detect at least
  `min.features`. We will define `min.features = 1` to filter the barcodes
  (cells, nuclei) that do not contain any UMI. This also allows to reduce the weight
  of the Seurat object in our environment.

``` r
## Creation of the Seurat Object
pbmc_small <- CreateSeuratObject(tenX_matrix,                    #Expression matrix
                                 project = "PBMC analysis",      #Name of the project : something meaningful for your dataset
                                 assay = "RNA",                  #Name of the initial assay, (others can be created during the analysis), default is RNA
                                 names.field = 1,                #Associated with the names.delim, it allows to separate the name of the barcode when this one is composed of several information ex: BARCODE_CLUSTER_CELLTYPE and to choose after split on the names.delim which part we choose as the barcode name
                                 names.delim = "_",              #Character that allows to dissociate the different information contained in the barcode names
                                 meta.data = NULL,               #We can add the metadata on the transcriptomes with a dataframe where the barcodes are in line and the different information in column
                                 min.cells = 0,                  #Filtering genes that are not detected in less than min.cells
                                 min.features = 1)               #Filtering cells that do not detect at least min.features, here we filter all barcodes that detect no gene
```

#### Exploration of the SeuratObject

=== "SeuratObject Presentation"

    ``` r
    pbmc_small #Small presentation of the Seurat object in R console
    ```

        ## An object of class Seurat
        ## 32738 features across 2700 samples within 1 assay
        ## Active assay: RNA (32738 features, 0 variable features)

=== "Structure of the SeuratObject"

    ``` r
    ## Discovery of SeuratObject
    str(pbmc_small)
    ```

        ## Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ##   ..@ assays      :List of 1
        ##   .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        ##   .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        ##   .. .. .. .. .. ..@ i       : int [1:2286884] 70 166 178 326 363 410 412 492 494 495 ...
        ##   .. .. .. .. .. ..@ p       : int [1:2701] 0 781 2133 3264 4224 4746 5528 6311 7101 7634 ...
        ##   .. .. .. .. .. ..@ Dim     : int [1:2] 32738 2700
        ##   .. .. .. .. .. ..@ Dimnames:List of 2
        ##   .. .. .. .. .. .. ..$ : chr [1:32738] "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" ...
        ##   .. .. .. .. .. .. ..$ : chr [1:2700] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
        ##   .. .. .. .. .. ..@ x       : num [1:2286884] 1 1 2 1 1 1 1 41 1 1 ...
        ##   .. .. .. .. .. ..@ factors : list()
        ##   .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        ##   .. .. .. .. .. ..@ i       : int [1:2286884] 70 166 178 326 363 410 412 492 494 495 ...
        ##   .. .. .. .. .. ..@ p       : int [1:2701] 0 781 2133 3264 4224 4746 5528 6311 7101 7634 ...
        ##   .. .. .. .. .. ..@ Dim     : int [1:2] 32738 2700
        ##   .. .. .. .. .. ..@ Dimnames:List of 2
        ##   .. .. .. .. .. .. ..$ : chr [1:32738] "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" ...
        ##   .. .. .. .. .. .. ..$ : chr [1:2700] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
        ##   .. .. .. .. .. ..@ x       : num [1:2286884] 1 1 2 1 1 1 1 41 1 1 ...
        ##   .. .. .. .. .. ..@ factors : list()
        ##   .. .. .. ..@ scale.data   : num[0 , 0 ]
        ##   .. .. .. ..@ key          : chr "rna_"
        ##   .. .. .. ..@ assay.orig   : NULL
        ##   .. .. .. ..@ var.features : logi(0)
        ##   .. .. .. ..@ meta.features:'data.frame':   32738 obs. of  0 variables
        ##   .. .. .. ..@ misc         : list()
        ##   ..@ meta.data   :'data.frame': 2700 obs. of  3 variables:
        ##   .. ..$ orig.ident  : Factor w/ 1 level "PBMC analysis": 1 1 1 1 1 1 1 1 1 1 ...
        ##   .. ..$ nCount_RNA  : num [1:2700] 2421 4903 3149 2639 981 ...
        ##   .. ..$ nFeature_RNA: int [1:2700] 781 1352 1131 960 522 782 783 790 533 550 ...
        ##   ..@ active.assay: chr "RNA"
        ##   ..@ active.ident: Factor w/ 1 level "PBMC analysis": 1 1 1 1 1 1 1 1 1 1 ...
        ##   .. ..- attr(*, "names")= chr [1:2700] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
        ##   ..@ graphs      : list()
        ##   ..@ neighbors   : list()
        ##   ..@ reductions  : list()
        ##   ..@ images      : list()
        ##   ..@ project.name: chr "PBMC analysis"
        ##   ..@ misc        : list()
        ##   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        ##   .. ..$ : int [1:3] 4 1 0
        ##   ..@ commands    : list()
        ##   ..@ tools       : list()

=== "Dimension of expression matrices"

    ``` r
    dim(pbmc_small@assays$RNA@counts)
    ```

        ## [1] 32738  2700

    ``` r
    dim(pbmc_small@assays$RNA@data)
    ```

        ## [1] 32738  2700

=== "Cell Metadata"

    ``` r
    kable(head(pbmc_small@meta.data), "simple") #Preview of the cell metadata
    ```

    |                  | orig.ident    | nCount_RNA | nFeature_RNA |
    |------------------|:--------------|-----------:|-------------:|
    | AAACATACAACCAC-1 | PBMC analysis |       2421 |          781 |
    | AAACATTGAGCTAC-1 | PBMC analysis |       4903 |         1352 |
    | AAACATTGATCAGC-1 | PBMC analysis |       3149 |         1131 |
    | AAACCGTGCTTCCG-1 | PBMC analysis |       2639 |          960 |
    | AAACCGTGTATGCG-1 | PBMC analysis |        981 |          522 |
    | AAACGCACTGGTAC-1 | PBMC analysis |       2164 |          782 |

We can observe several slots via the `str` command:  

- `assays`: general slot that will include the different information of each
  study. They are composed of several things:
    - the starting expression matrix (`@counts`), usually raw counts or raw UMI
    - the one that will "undergo" all the modifications (filters,
      normalization, etc) (`@data`)
    - dataframe that will be created when scaling the data (`@scale_data`)
    - prefix used for each calculation that will use this assay (study)
      (`@key`)
    - vector of gene names that will be determined to have a variable
      expression (`@var.features`)
    - dataframe associated with genes with different metadata
      (`@meta.features`)

- `meta.data` : gathers all the information about the cells. At the beginning
  Seurat will calculate the size of the library (nCount_RNA, or the total
  number of UMI) and the number of detected genes (nFeatures_RNA) for each
  cell. If a dataframe is given in the `meta.data` parameter of the
  `CreateSeuratObject` function, its columns will be added after those
  calculated by Seurat.  

- `active.assay` : study used by default
- `active.ident` : default cell identity, here the name of the given project,
  also stored under the column `orig.ident` in the metadata

!!! note
    Navigation in the different slots is done via `@` or `$`. Each main slot is
    accessible via the `@`, *i.e.* `object@main slot` to go further in the
    slots tree, most often complex objects are accessible with a `@` (dgCMatrix,
    dataframe) and lists, vectors are accessible via `$`. If in doubt, you can
    refer to the result of the `str` command and use the character in front of
    each slot name.
