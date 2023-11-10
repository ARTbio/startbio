<!-- ## Dataframes and Matrices -->

Both `data.frame` and `matrix` are two-dimensional objects, they consist of rows and columns.
The main difference is that `matrix` can only store one class of data (either character or numeric),
while `data.frame` can store different classes of data (numeric, character and factor).


## Matrices

A matrix is a data structure used for mathematical computations, linear algebra operations, and storing homogeneous data. Learn more about matrices with chapter [5.3](https://rstudio-education.github.io/hopr/r-objects.html#matrices) of Garett’s book.

### Creation of Matrices

We can use `matrix()` function to create a matrix, or simply bind vectors as rows or columns in a matrix.

```r
matrix(1:6)
matrix(
  data = 1:6,
  nrow = 3, ncol = 2, # only need to specify one of these two parameters
  byrow = FALSE, # fill by column by default
  dimnames = list(paste0("row", 1:3), paste0("col", 1:2)) # can provide the rownames and colnames in a list
)

cbind(1:5, 6:10, 11:15) # bind vectors by columns
rbind(1:5, 6:10, 11:15) # bind vectors by rows
```

### Access and Modification of Matrices

To access the elements of a matrix, we use the operator `[]` and specify the index or the name of column(s) and/or row(s) separated by a comma:

```r
my_mat <- matrix(1:6, nrow = 3, dimnames = list(paste0("row", 1:3), paste0("col", 1:2)))

my_mat[1, ] # the 1st row, same as my_mat["row1", ]
my_mat[, 1] # the 1st column, same as my_mat[, "col"]

my_mat[c(1, 3), 2] # the 2nd elements of the 1st and the 3rd row, same as my_mat[c("row1", "row3"), "col2"]
```

To modify element(s) of a matrix, we can affect new value(s) to wanted position(s) using the index or colname/rowname.
When the provided new values have a different length than the original ones,
R will return an error, except you want to replace element(s) by a single new value.

```r
my_mat2 <- matrix(1:15, nrow = 3)
my_mat2[2, 1] <- 99
my_mat2[c(2, 3), c(4, 5)] <- 21:24 # replace by column
my_mat2[1, ] <- 2 # replace all values of the 1st row by 2
```

### Matrix Operation

We can do all kinds of matrix operation with R, for instance:

```r
X <- matrix(c(9, 2, -3, 2, 4, -2, -3, -2, 16), 3, byrow = TRUE)
Y <- matrix(0:8, ncol = 3)

X %*% Y # matrix multiplication
t(X) # matrix transpose
det(X) # calculate the determinant of a matrix
solve(X) # inverse X
diag(X) # matrix diagonals

colSums(X)
rowSums(X)
colMeans(X)
rowMeans(X)
```

### To Go Further

You may have heard of a sparse matrix, it is also a matrix but contains a lot of zeros.

For example, the expression matrix of single cell RNAseq data from `Seurat` object is stored using this special class `dgCMatrix` (see the scRNAseq tutorial [here](../scRNAseq_basics/import.md)).
We can create a sparse matrix using the R package <code>[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)</code>.

```r
# install.packages("Matrix")
library(Matrix)
sparse_mat <- sparseMatrix(
  i = c(1, 3:8), # position of rows
  j = c(2, 9, 6:10), # position of columns
  x = 7 * (1:7) # non-zero values to fill in the sparse matrix
)
sparse_mat
```


## Data Frames

A data framematrix is a two-dimensional data structure that organizes data into rows and columns, similar to a table or spreadsheet. It is particularly useful for working with heterogeneous datasets where different columns can have different data types.that stores data in a grid-like format consisting of rows and columns. It’s composed of vectors arranged in columns, see chapter [5.8](https://rstudio-education.github.io/hopr/r-objects.html#data-frames) of Garett’s book and sections [8.2.3 to 8.2.4](https://bookdown.org/ndphillips/YaRrr/creating-matrices-and-dataframes.html#data.frame) of Philips’ book for more details. 


### Tibble {#sec-tibble}
### data.table


### Manipulating Matrices and Dataframes

Please read with attention the sections [8.3 to 8.6](https://bookdown.org/ndphillips/YaRrr/matrix-and-dataframe-functions.html) of Philips’ book. 

**Exercises to manipulate matrices and dataframes**


