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
##      [,1]
## [1,]    1
## [2,]    2
## [3,]    3
## [4,]    4
## [5,]    5
## [6,]    6

my_mat <- matrix(
  data = 1:6,
  nrow = 3, ncol = 2, # only need to specify one of these two parameters
  byrow = FALSE, # fill by column by default
  dimnames = list(paste0("row", 1:3), paste0("col", 1:2)) # can provide the rownames and colnames in a list
)
my_mat
##      col1 col2
## row1    1    4
## row2    2    5
## row3    3    6

cbind(1:5, 6:10, 11:15) # bind vectors by columns
##      [,1] [,2] [,3]
## [1,]    1    6   11
## [2,]    2    7   12
## [3,]    3    8   13
## [4,]    4    9   14
## [5,]    5   10   15

rbind(1:5, 6:10, 11:15) # bind vectors by rows
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    2    3    4    5
## [2,]    6    7    8    9   10
## [3,]   11   12   13   14   15
```

### Access and Modification of Matrices

To access the elements of a matrix, we use the operator `[]` and specify the index or the name of column(s) and/or row(s),
separated by a comma, for example with privously created `my_mat`:

```r
my_mat[1, ] # the 1st row, same as my_mat["row1", ]
## col1 col2 
##    1    4 

my_mat[, 1] # the 1st column, same as my_mat[, "col"]
## row1 row2 row3 
##    1    2    3

my_mat[c(1, 3), 2] # the 2nd elements of the 1st and the 3rd row, same as my_mat[c("row1", "row3"), "col2"]
## row1 row3 
##    4    6 
```

If you use just one index without specifying any commas, you will get the Nth element of the matrix in column order.

```r
my_mat[3] # return the 3rd element
## [1] 3

my_mat[5] # return the 5th element, so the element in row 2 and column 2
## [1] 5
```

!!! tip "About dimensionality..."
    Did you notice that the results are all vectors?
    To learn how to preserve the dimensionality, please check the section [4.2.5](https://adv-r.hadley.nz/subsetting.html#simplify-preserve) of Hadley Wickham's "Advanced R".


To modify element(s) of a matrix, we can affect new value(s) to wanted position(s) using the index or colname/rowname.
When the provided new values have a different length than the original ones,
R will return an error, except you want to replace element(s) by a single new value.

```r
my_mat2 <- matrix(1:15, nrow = 3)
my_mat2
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    4    7   10   13
## [2,]    2    5    8   11   14
## [3,]    3    6    9   12   15

my_mat2[2, 1] <- 99
my_mat2
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    4    7   10   13
## [2,]   99    5    8   11   14
## [3,]    3    6    9   12   15

my_mat2[c(2, 3), c(4, 5)] <- 21:24 # replace by column
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    4    7   10   13
## [2,]   99    5    8   21   23
## [3,]    3    6    9   22   24

my_mat2[1, ] <- 2 # replace all values of the 1st row by 2
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    2    2    2    2    2
## [2,]   99    5    8   21   23
## [3,]    3    6    9   22   24
```

### Matrix Operation

We can do all kinds of matrix operation with R, for instance:

```r
X <- matrix(c(9, 2, -3, 2, 4, -2, -3, -2, 16), 3, byrow = TRUE)
X
##      [,1] [,2] [,3]
## [1,]    9    2   -3
## [2,]    2    4   -2
## [3,]   -3   -2   16

colSums(X) # calculate the sum by column 
## [1]  8  4 11

rowSums(X) # calculate the sum by row
## [1]  8  4 11

colMeans(X) # calculate the average by column
## [1] 2.666667 1.333333 3.666667

rowMeans(X) # calculate the average by row
## [1] 2.666667 1.333333 3.666667

t(X) # matrix transpose
##      [,1] [,2] [,3]
## [1,]    9    2   -3
## [2,]    2    4   -2
## [3,]   -3   -2   16

det(X) # calculate the determinant of a matrix
## [1] 464

solve(X) # inverse X
##             [,1]        [,2]       [,3]
## [1,]  0.12931034 -0.05603448 0.01724138
## [2,] -0.05603448  0.29094828 0.02586207
## [3,]  0.01724138  0.02586207 0.06896552

diag(X) # matrix diagonals
## [1]  9  4 16

Y <- matrix(0:8, ncol = 3)
Y
##      [,1] [,2] [,3]
## [1,]    0    3    6
## [2,]    1    4    7
## [3,]    2    5    8

X %*% Y # matrix multiplication
##      [,1] [,2] [,3]
## [1,]   -4   20   44
## [2,]    0   12   24
## [3,]   30   63   96
```

### To Go Further

You may have heard of a sparse matrix, it is also a matrix but contains a lot of zeros (usually more than 2/3 of all values).

For example, the single cell RNAseq gives cell-level expression resolution,
and it is likely that only a small fraction of known genes are expressed in a single cell.
Therefore, single-cell RNAseq expression data are stored using a special class `dgCMatrix` developped for sparse matrix,
where only non-zero values are stored to save memory usage. 
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
## 8 x 10 sparse Matrix of class "dgCMatrix"
##                              
## [1,] . 7 . . .  .  .  .  .  .
## [2,] . . . . .  .  .  .  .  .
## [3,] . . . . .  .  .  . 14  .
## [4,] . . . . . 21  .  .  .  .
## [5,] . . . . .  . 28  .  .  .
## [6,] . . . . .  .  . 35  .  .
## [7,] . . . . .  .  .  . 42  .
## [8,] . . . . .  .  .  .  . 49
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
my_df
##   id age    sex
## 1  1  21 female
## 2  2  25 female
## 3  3  18   male
## 4  4  35   male
## 5  5  27   male

rownames(my_df) <- paste0("sample", 1:5) # name rows
my_df
##         id age    sex
## sample1  1  21 female
## sample2  2  25 female
## sample3  3  18   male
## sample4  4  35   male
## sample5  5  27   male
```
You can use `rownames()` and `colnames()` to name/rename the rows or columns.
But the rownames and colnames should be unique, which allow us to have acces to the exact wanted value(s). 

By default, `data.frame` requires the columns are of equal length.
If it is not the case and when it is possible, *i.e.* the length of the longest column is a multiple of the lenghth of shorter column(s) , `data.frame` will recycle the elements of the shorter column(s).

```r
data.frame(x = 1:4, y = 1:2)
##   x y
## 1 1 1
## 2 2 2
## 3 3 1
## 4 4 2

data.frame(x = 1:4, y = 1:3) # different length
## Error in data.frame(x = 1:4, y = 1:3) : 
##   arguments imply differing number of rows: 4, 3
```

You can use the `as.data.frame()` function to convert a `vector`, a `list` or a `matrix` into a `data.frame`.

### Access and Modification of Data Frames

Similar to the way that we use for matrix, we can access the elements in a data.frame by using the index or the name of rows or columns.

```r
my_df[, c(2, 3)] # get the 2nd and the 3rd columns by column index
##         age    sex
## sample1  21 female
## sample2  25 female
## sample3  18   male
## sample4  35   male
## sample5  27   male

## Alternative ways
my_df[c(2, 3)] # the same as above but NOT suggested
my_df[, c("age", "sex")] # get the 2nd and the 3rd columns by colnames
my_df[c("age", "sex")] # the same as above but NOT suggested
```

We can inverse select the column by adding `-` before the column index or name:

```r
my_df[, -c(2, 3)] # get all columns except the 2nd and the 3rd
## [1] 1 2 3 4 5
```

!!! question
    How to maintain dimensionality when subsetting results in selecting only one column?

The same logique to access to rows:

```r
my_df[c("sample1", "sample5"), ]
## or
my_df[c(1, 5), ]
##         id age    sex
## sample1  1  21 female
## sample5  5  27   male

my_df[-c(1, 5), ] # reverse selection
##         id age    sex
## sample2  2  25 female
## sample3  3  18   male
## sample4  4  35   male
```

Particularly for columns, the `[[` and `$` operators can be used to select a single column and return the values in a vector.
The main difference is that `$` does not allow index, while `[[` allow both column name and index.

```r
my_df$age
## [1] 21 25 18 35 27

## Alternative ways
my_df[["age"]]
my_df[[2]]
```

How to modify a `data.frame`?
We can use `$` or `cbind` to add a named vector of equal length as the other columns or a named vector length of 1 as a new column.

```r
my_df$new_col <- letters[1:nrow(my_df)]
## or
my_df <- cbind(my_df, "new_col" = letters[1:nrow(my_df)])
my_df
##         id age    sex new_col
## sample1  1  21 female       a
## sample2  2  25 female       b
## sample3  3  18   male       c
## sample4  4  35   male       d
## sample5  5  27   male       e

my_df$new_col2 <- "cohort1"
my_df
##         id age    sex new_col new_col2
## sample1  1  21 female       a  cohort1
## sample2  2  25 female       b  cohort1
## sample3  3  18   male       c  cohort1
## sample4  4  35   male       d  cohort1
## sample5  5  27   male       e  cohort1
```

How about joining two `data.frame`s to get a bigger one? `merge()` if your friend.
The merge is based on either a specific column or the `row.names`.

```r
my_df2 <- data.frame(
  id = 1:10,
  status = rep(c("case", "control"), each = 5)
)
my_df2
##    id  status
## 1   1    case
## 2   2    case
## 3   3    case
## 4   4    case
## 5   5    case
## 6   6 control
## 7   7 control
## 8   8 control
## 9   9 control
## 10 10 control

merge(
  x = my_df, y = my_df2,
  by = "id", # the name of the column to use for merging
  all.x = TRUE # the merging is based on the rows of the data frame provided in "x"
)
##   id age    sex new_col new_col2 status
## 1  1  21 female       a  cohort1   case
## 2  2  25 female       b  cohort1   case
## 3  3  18   male       c  cohort1   case
## 4  4  35   male       d  cohort1   case
## 5  5  27   male       e  cohort1   case
```

To delete colummns in a `data.frame`, we can simply affect the wanted columns to `NULL`.

```r
my_df$new_col <- NULL
my_df
##         id age    sex new_col2
## sample1  1  21 female  cohort1
## sample2  2  25 female  cohort1
## sample3  3  18   male  cohort1
## sample4  4  35   male  cohort1
## sample5  5  27   male  cohort1
```

For other possible manipulations in `matrix` and `data.frame`, please refer to the sections [8.3 to 8.6](https://bookdown.org/ndphillips/YaRrr/matrix-and-dataframe-functions.html) of Philips’ book.

---
