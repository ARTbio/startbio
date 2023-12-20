<!-- ## Data Import and Export -->

We can deal with different formats of data with R, such as
the text files (.csv, .tsv, .txt), the Excel files (.xlx, .xlsx), and the R data file formats (.RDS, .RData).
It is also possible to read or write other program-speficied formats, for example
"SAS", "SPSS" or "Minitab" files.

There are R functions which allow users to download files from the Internet, for example the `download.file()` from the `utils` package or the similar function `curl_download()` from the <code>[curl](https://cran.r-project.org/web/packages/curl/index.html)</code> package.

```r
download.file(
  url = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000841/ScoringFiles/PGS000841.txt.gz",
    # example of the polygenic risk score file for the trait BMI
  destfile = "path/to/out_dir"
)

# install.packages("curl")
curl::curl_download(
  url = "http://url_to_wanted_file",
  destfile = "path/to/out_dir"
)
```

## Importing Data

In order to avoid error when you import the file,
it's suggested to name the data columns following the naming conventions (see [Best Practices](r04_bestpractices.md) section).

### Read Text Files (.csv/.tsv/.txt)

With the functions from the basic package "utils":

```r
## read a comma separated file ".csv"
my_csv <- read.csv(
  file = "path/to/my_file.csv",
  header = TRUE,    # whether the file has a colnames in the 1st row
  sep = ",",        # the character used in the file to separate columns
  quote = "\"",     # the character(s) for quotes
  dec = ".",        # the character for decimal points
  fill = TRUE,      # in case of rows with unequal length, whether to add empty field
  comment.char = "" # character used to indicated rows as comment lines, use empty string to turn off the interpretation of comment lines
)

## read a tab separated file ".tsv"
my_delim <- read.delim(
  file = "path/to/my_file.tsv",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  dec = ".",
  fill = TRUE,
  comment.char = ""
)

## read a text file ".txt"
### sometimes read.table does not work well, be careful of the parameter settings
my_table <- read.table(
  file = "path/to/my_file.txt",
  header = TRUE,
  sep = " ",
  quote = "\"'",     # the character(s) for quotes
  dec = ".",
  comment.char = "#"
)
```

!!! note
  In R, the backslash `\` is used to escape the character after it.
  As we use `""` to pass value to the argument `quote`,  meanwhile the `"` is the quoting character used in the file to read, we need to "protect" the quoting character by the backslash to let R know the `"` between the `""` is a real character to be evaluated. 


There are other useful arguments that we didn't mentionned here,
such as `na.strings` (characters to be interpreted as `NA` values), `colClasse` (type of columns), *etc.*, 
please check the document with `?read.table`.

All these functions from `utils` will return a [data.frame](r06_df_matrices.md).

Apart from the `utils` package, we can use other packages for instance 
<code>[readr](https://cran.r-project.org/web/packages/readr/index.html)</code> and
<code>[vroom](https://cran.r-project.org/web/packages/vroom/index.html)</code> to achieve the same goal.


### Read Excel Files (.xlx, .xlsx)

Excel files are very often used to store clinical and biological experiment data.
We can use the <code>[readxl](https://cran.r-project.org/web/packages/readxl/index.html)</code>
to import them into R:

```r
# install.packages("readxl")
library("readxl")
my_xls <- read_xls(
  path = "path/to/my_file.xls",
  sheet = 1,                    # sheet to read, can be the sheet number or the sheet name
  range = "B3:D87",             # a cell range to read from
  col_names = TRUE,             # use the 1st row as column names
  col_types = NULL,             # type of columns, use "NULL" to guess automatically the type of each column
  na = "",                      # characters to be interpreted as `NA` values
  trim_ws = TRUE,               # trim surrounding spaces
  skip = 0,                     # number of rows to skip before reading
  n_max = Inf,                  # maximum number of rows to read
  guess_max = min(1000, n_max), # maximum number of rows to use for guessing column type
  .name_repair = "unique"       # argument passed to tibble::as_tibble, default is to ensure unique and not empty column names
)

my_xlsx <- read_xlsx(
  path = "path/to/my_file.xlsx",
  sheet = NULL,
  range = NULL,
  col_names = TRUE,
  col_types = NULL,
  na = "",
  trim_ws = TRUE,
  skip = 0,
  n_max = Inf,
  guess_max = min(1000, n_max),
  .name_repair = "unique"
)

```

All these functions will return a [tibble](https://tibble.tidyverse.org/reference/tibble-package.html),
which is a more efficient version of data.frame created and used with the `tidyverse` packages (See "Tibble" section in the [tidyverse](r10_tidyverse.md) chapter).

!!! danger "Caution"
    The coloring of cells in Excel files CANNOT be handled.

!!! tip
    - `excel_sheets()` is useful to list all sheets without openning the file.
    - The merged cells can be handled by using the `openxlsx::read.xlsx()` with specifying `fillMergedCells = TRUE`,
    the value in a merged cell is given to all cells within the merge.
    

### Read R Data Format (.RDS, .RData)

`.RDS` is used to save a single R object and `.RData` is used to store multiple R objects.
We can use `readRDS` to import the `.RDS` file and `load` to open the `.RData` file.

```r
my_rds <- readRDS(file = "path/to/my_file.RDS")
load(file = "path/to/my_file.RData")
```

Please note that you don't need to assign the loaded RData to any object,
with `load`, you will import all objects stored in the RData file with their original name.
R will not show what objects were loaded into the working session,
if your environment has an object with the same name as one object from the RData file,
it will be overwritten by what you've loaded.

### Other Program-specified Format (SAS, SPSS, Minitab, *etc.*)

The R package <code>[foreign](https://cran.r-project.org/web/packages/foreign/index.html)</code>
was developped to import these data. For example:

```r
install.packages("foreign")
library("foreign")
my_sav <- read.spss(file = "path/to/my_file.sav")
my_xpt <- read.xport(file = "path/to/my_file.xpt")
my_mtp <- read.mtp(file = "path/to/my_file.mtp")
```

!!! tip
    More detailed information about the data import and export in R can be found in the chapter [R Files](https://rstudio-education.github.io/hopr/dataio.html) of the eBook [Hands-On Programming with R](https://rstudio-education.github.io/hopr/)


## Exporting Data

As we can import different formats of data into R, we can also export them from R and save it into various formats.
Based on what we've seen in the above section, the functions used to save data are usually named in the same way as the data importing functions, by changing "read" by "write" or "save".

### Write Text Files

```r
my_object <- data.frame(
  x = 1:3,
  y = letters[1:3]
)
write.csv(
  x = my_object,
  file = "path/to/my_file.csv",
  quote = TRUE,                 # whether to quote columns in the output file
  sep = ",",                    # the field separator
  eol = "\n",                   # the character to print at the end of the line
  na = "NA",                    # the character to mark missing value
  dec = ".",                    # the decimal point character
  row.names = FALSE,            # whether to write rownames in the output file
  col.names = TRUE              # whether to write colnames in the output file
)
write.table(x = my_object, file = "path/to/my_file.txt")
```

### Write Excel Files

Several packages are available to export data frame to Excel `.xlsx` format, for example 
<code>[writexl](https://cran.r-project.org/web/packages/writexl/index.html)</code> and
<code>[openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)</code>.
Here we illustrate with the `writexl` package:

```r
# install.packages("writexl")
library("writexl")

write_xlsx(
  x = my_object,
  path = "path/to/my_file.xlsx",
  col_names = TRUE
)
```

We can also store a list of data.frame into an Excel file with multiple sheet with the same function as follow: 

```r
list_df <- list(
  "df1" = data.frame(
    x = 1:3,
    y = letters[1:3]
  ),
  "df2" = data.frame(
    x = 4:6,
    y = letters[4:6]
  )
)
write_xlsx(
  x = list_df,
  path = "path/to/multi_sheets.xlsx",
  col_names = TRUE
)
```

### Write R Data Format

```r
a <- 1
b <- c(1, "abc", 5)
df <- data.frame(
  x = 1:3,
  y = letters[1:3]
)

saveRDS(a, file = "path/to/my_object.RDS")
save(a, b, df, file = "path/to/my_objects.RData")
```

---