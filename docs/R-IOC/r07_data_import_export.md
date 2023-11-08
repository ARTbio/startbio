## Data Import and Export

We can deal with different formats of data with R, such as
the text files (.csv, .tsv, .txt), the Excel files (.xlx, .xlsx), and the R data file formats (.rds, .Rdata).
It's also possible to read or write other softwares-speficied formats, for example
"SAS", "SPSS" or "Minitab" files.


## Importing Data

In order to avoid error when you import the file,
it's suggested to name the data columns following the naming conventions (see [Best Practices](r04_bestpractices.md) section).

### Read Text Files (.csv/.tsv/.txt)

With the functions from the basic package "utils":

```r
my_file <- "path/to/my_file"

## read a comma separated file ".csv"
read.csv(
  file = my_file,
  header = TRUE,    # whether the file has a colnames in the 1st row
  sep = ",",        # what is the field separator
  quote = "\"",     # the character(s) for quotes
  dec = ".",        # the character for decimal points
  fill = TRUE,      # in case of rows with unequal length, whether to add empty field
  comment.char = "" # character used to indicated rows as comment lines, use empty string to turn off the interpretation of comment lines
)

## read a tab separated file ".tsv"
read.delim(
  file = my_file,
  header = TRUE,
  sep = "\t",       # what is the field separator
  quote = "\"",
  dec = ".",
  fill = TRUE,
  comment.char = ""
)

## read a text file ".txt"
read.table(
  file = my_file,
  header = TRUE,
  sep = "",          # what is the field separator
  quote = "\"'",     # the character(s) for quotes
  dec = ".",
  comment.char = "#"
)
```

There are other useful arguments that we didn't mentionned here,
such as `na.strings` (characters to be interpreted as `NA` values), `colClasse` (type of columns), *etc.*, 
please check the document with `?read.table`.

All these functions from `utils` will return a [data.frame](r06_df_matrices.md).

Apart from the `utils` packages, we can use other package for instance 
<code>[readr](https://cran.r-project.org/web/packages/readr/index.html)</code> and
<code>[vroom](https://cran.r-project.org/web/packages/vroom/index.html)</code> to achieve the same goal.


### Read Excel Files (.xlx, .xlsx)

Excel files are very often used to store clinical and biological experiment data.
We can use the <code>[readxl](https://cran.r-project.org/web/packages/readxl/index.html)</code>
to import them into R:

```r
library("readxl")
read_xls(
  path = my_file,
  sheet = NULL,                 # sheet to read, can be the sheet number or the sheet name
  range = NULL,                 # a cell range to read
  col_names = TRUE,             # use the 1st row as column names
  col_types = NULL,             # type of columns, use "NULL" to guess automatically the type of each column
  na = "",                      # characters to be interpreted as `NA` values
  trim_ws = TRUE,               # trim surrounding spaces
  skip = 0,                     # number of rows to skip before reading
  n_max = Inf,                  # maximum number of rows to read
  guess_max = min(1000, n_max), # maximum number of rows to use for guessing column type
  .name_repair = "unique"       # argument passed to tibble::as_tibble, default is to ensure unique and not empty column names
)

read_xlsx(
  path = my_file,
  sheet = NULL,
  range = NULL,
  col_names = TRUE,
  col_types = NULL,
  na = "",
  trim_ws = TRUE,
  skip = 0,
  n_max = Inf,
  guess_max = min(1000, n_max),
  progress = readxl_progress(),
  .name_repair = "unique"
)

```

All these functions will return a [tibble](https://tibble.tidyverse.org/reference/tibble-package.html).

!!! danger "Caution"
    The coloring of cells in Excel files CANNOT be handled.

!!! tip
    The merged cells can be handled by using the `openxlsx::read.xlsx` with specifying `fillMergedCells = TRUE`,
    the value in a merged cell is given to all cells within the merge.


### Read R Data Format (.rds, .Rdata)

### Other Software-specified Format (SAS, SPSS, Minitab, *etc.*)

The R package <code>[foreign](https://cran.r-project.org/web/packages/foreign/index.html)</code>
was developped to import these data. 

## Exporting Data
