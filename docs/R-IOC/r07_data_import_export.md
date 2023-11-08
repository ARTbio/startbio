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
  header = TRUE,    # whether the file has a colnames in the 1st row
  sep = "\t",       # what is the field separator
  quote = "\"",     # the character(s) for quotes
  dec = ".",        # the character for decimal points
  fill = TRUE,      # in case of rows with unequal length, whether to add empty field
  comment.char = "" # character used to indicated rows as comment lines, use empty string to turn off the interpretation of comment lines
)

## read a text file ".txt"
read.table(
  file = my_file,
  header = TRUE,     # whether the file has a colnames in the 1st row
  sep = "",          # what is the field separator
  quote = "\"'",     # the character(s) for quotes
  dec = ".",         # the character for decimal points
  comment.char = "#" # character used to indicated rows as comment lines, use empty string to turn off the interpretation of comment lines
)
```

There are other useful arguments that we didn't mentionned here,
such as `na.strings` (characters to be interpreted as `NA` values), `colClasse` (type of columns), *etc.*, 
please check the document with `?read.table`.

Apart from the `utils` packages, we can use other package such as 
<code>[readr](https://cran.r-project.org/web/packages/readr/index.html)</code> and
<code>[vroom](https://cran.r-project.org/web/packages/vroom/index.html)</code> to achieve the same goal.


### Read Excel Files (.xlx, .xlsx)




## Exporting Data


https://cran.r-project.org/web/packages/foreign/index.html