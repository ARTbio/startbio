<!-- ## Data Frames and Matrices -->

Both `data.frame` and `matrix` are two-dimensional objects, they consist of rows and columns.
The main difference is that `matrix` can only store one class of data (either character or numeric),
while `data.frame` can store different classes of data (numeric, character and factor).


## Matrices

A `matrix` is a data structure used for mathematical computations, linear algebra operations, and storing homogeneous data. Learn more about matrices with chapter [5.3](https://rstudio-education.github.io/hopr/r-objects.html#matrices) of Garett’s book.

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

If you use just one index without specifying any commas, you will get the Nth element of the matrix in column order.

```r
my_mat[3] # return the 3rd element
my_mat[5] # return the 5th element, so the element in row 2 and column 2
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

A `data.frame` is similar to a Excel spreadsheet.
It is particularly useful for working with heterogeneous datasets where different columns can have different data types.
It can be considered as a list of vectors of equal length, arranged in columns.
See chapter [5.8](https://rstudio-education.github.io/hopr/r-objects.html#data-frames) of Garett’s book and sections [8.2.3 to 8.2.4](https://bookdown.org/ndphillips/YaRrr/creating-matrices-and-dataframes.html#data.frame) of Philips’ book for more details. 

### Creation of Data Frames

Similar to how to create a matrix, we can bind named vectors through the function `data.frame()` to create a data frame with a mixture of numeric and character columns.

```r
my_df <- data.frame(
  "id" = 1:5,
  "age" = c(21, 25, 18, 35, 27),
  "sex" = c("female", "female", "male", "male", "male"),
  stringsAsFactors = FALSE # by default for R > 4.0
)
rownames(my_df) <- paste0("sample", 1:5) # data.frame can have unique rownames
```

By default, `data.frame` requires the columns are of equal length.
If it is not the case and when it is possible, `data.frame` will recycle the elements of the shorter vector(s).

```r
data.frame(x = 1:4, y = 1:2)
data.frame(x = 1:4, y = 1:3) # will return an error
```

You can use `rownames()` and `colnames()` to name/rename the names of row or column.
You can use the `as.data.frame()` function to convert a `vector`, a `list` or a `matrix` into a `data.frame`.

### Access and Modification of Data Frames

Similar to the way that we use for matrix, we can access the elements in a data.frame by using the index or the name of rows or columns.

```r
my_df[, c(2, 3)] # get the 2nd and the 3rd columns by column index
my_df[c(2, 3)] # the same as above but NOT suggested

my_df[, c("age", "sex")] # get the 2nd and the 3rd columns by colnames
my_df[c("age", "sex")] # the same as above but NOT suggested
```

We can inverse select the column by adding `-` before the column index or name:

```r
my_df[, -c(2, 3)] # get all columns except the 2nd and the 3rd
my_df[, -c("age", "sex")]
```

The same logique to access to rows:

```r
my_df[c("sample1", "sample5"), ]
my_df[c(1, 5), ]

my_df[-c("sample1", "sample5"), ]
my_df[-c(1, 5), ]
```

Particularly for columns, the `[[` and `$` operators can be used to select a single column and return the values in a vector.
The main difference is that `$` does not allow index, while `[[` allow both column name and index.

```r
my_df$age
my_df[["age"]]
my_df[[2]]
```

How to modify a `data.frame`?
We can use `$` or `cbind` to add a named vector of equal length as the other columns or a named vector length of 1 as a new column.

```r
my_df$new_col <- letters[1:5]
# or
cbind(my_df, "new_col" = letters[1:5])
my_df$new_col2 <- "cohort1"
```

To delete colummns in a `data.frame`, we can simply affect the wanted columns to `NULL`.

```r
my_df$new_col <- NULL
```

For other possible manipulations in `matrix` and `data.frame`, please refer to the sections [8.3 to 8.6](https://bookdown.org/ndphillips/YaRrr/matrix-and-dataframe-functions.html) of Philips’ book.

### To Go Further

Besides the traditional `data.frame` class, there are other "enhanced" version of `data.frame`, for example:

* the `tibble` ([tbl-df](https://search.r-project.org/CRAN/refmans/tibble/html/tbl_df-class.html)), which is the central data structure used by differents packages from <code>[tidyverse](https://www.tidyverse.org/packages/)</code>.
Read more about `tibble` [here](https://r4ds.had.co.nz/tibbles.html) of Hadley Wickham's book "R for Data Science".
* the <code>[data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)</code>, which was developped to manipulate large data structure in a fast and memory efficient way.


## Exercises

* Create an object of matrix with 3 rows and 4 columns, fill with numbers 1 to 12 by row,
name the rows with "r1", "r2", "r3" and the columns with "c1", "c2", "c3", "c4".

```r
my_mat <- matrix(
  1:12, nrow = 3, ncol = 4, byrow = TRUE,
  dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3", "c4"))
)
```

* Extract the 2nd row

```r
my_mat[2, ]
my_mat["r2", ]
```

* Extract the 2nd row but keep it in matrix format

```r
my_mat[2, , drop = FALSE]
my_mat["r2", , drop = FALSE]
```

* Extract the 2nd row using a logical vector

```r
my_mat[c(FALSE, TRUE, FALSE), ]
```

* What are the positions for the numbers that are multiples of 3?

```r
which((my_mat / 3) %in% seq(12/3))
```

* Add a column "c5" containing the values "a", "b", "c". What happens after this add?

```r
my_mat <- cbind(my_mat, "c5" = c("a", "b", "c"))
```

* Now delete the added column and convert the matrix to numeric mode.

```r
my_mat <- my_mat[, -5]
# or
my_mat <- my_mat[, !colnames(my_mat) %in% "c5"]

mode(my_mat) <- "numeric"
```

* Replace the element bigger than 10 by 99

```r
my_mat[my_mat > 10] <- 99
```

* Transforme the matrix to a `data.frame`

```r
my_df <- as.data.frame(my_mat)
```

* Use the rownames to create a new column "id".

```r
my_df$id <- rownames(my_mat)
```

* Which rows has duplicated values?

```r
apply(my_df, MARGIN = 1, FUN = function(i) length(unique(i)) < ncol(my_df))
```

* Create a new column "total" which calcul the sum of column "c1" to "c4" by row.

```r
my_df$total <- rowSums(my_df[, c("c1", "c2", "c3", "c4")])
```

* Change the column order to put the "id" in the first column.

```r
my_df <- my_df[, c("id", setdiff(colnames(my_df), "id"))]
```

* Remove the rownames.

```r
rownames(my_df) <- NULL
```

* Add a new row which contains the sum of each column (except the "id" column, put `NA` in the new row for this column).

```r
my_df <- rbind(my_df, c(NA, colSums(my_df[, -1])))
```

