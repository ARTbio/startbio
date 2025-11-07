## What's the tidyverse?

The `tidyverse` is a set of R packages for data manipulation and visualisation. You can learn 
more on their [website](https://www.tidyverse.org/). It contains the following R packages : 

- `dplyr`: A Grammar of Data Manipulation
- `tidyr`: Tidy Messy Data
- `stringr`: Simple, Consistent Wrappers for Common String Operations
- `tibble`: Simple Data Frames
- `ggplot2`: Create Elegant Data Visualisations Using the Grammar of Graphics (coming [next](r09_viz_ggplot2.md))
- `readr`: Read Rectangular Text Data (seen previously [here](r07_data_import_export.md))
- `forcats`: Tools for Working with Categorical Variables (Factors)
- `purrr`: Functional Programming Tools

You can install and load all these packages with following commands : 

```r
install.packages("tidyverse")
library(tidyverse)
```

It's a collection that is really useful and powerful in data science. We are going to rely
on the book [R for Data Science (2e edition)](https://r4ds.hadley.nz/) (O'Reilly Book on tidyverse, "**R4DS**" in short).

!!! abstract "Take a break & Read"
    For an introduction of the tidyverse, please go read the 
    [introduction](https://r4ds.hadley.nz/intro) of R4DS.

### Tibbles, the "new" data.frame

Tidyverse packages are based on the manipulation of a new type of variable, the `"tibble"`.
It's a variant of `data.frame`. Don't worry, you can manipulate `tibble` the same way as
`data.frame` that you learn previously (brackets `[ ]`, `$`, and with basics functions such 
as : `colnames`, `rownames`, `str`, *etc*.). 

What does it look like? 

```r
# An example of a tibble from tidyverse
starwars
## # A tibble: 87 × 14
##    name   height  mass hair_color skin_color eye_color birth_year sex   gender homeworld species films
##    <chr>   <int> <dbl> <chr>      <chr>      <chr>          <dbl> <chr> <chr>  <chr>     <chr>   <lis>
##  1 Luke …    172    77 blond      fair       blue            19   male  mascu… Tatooine  Human   <chr>
##  2 C-3PO     167    75 NA         gold       yellow         112   none  mascu… Tatooine  Droid   <chr>
##  3 R2-D2      96    32 NA         white, bl… red             33   none  mascu… Naboo     Droid   <chr>
##  4 Darth…    202   136 none       white      yellow          41.9 male  mascu… Tatooine  Human   <chr>
##  5 Leia …    150    49 brown      light      brown           19   fema… femin… Alderaan  Human   <chr>
##  6 Owen …    178   120 brown, gr… light      blue            52   male  mascu… Tatooine  Human   <chr>
##  7 Beru …    165    75 brown      light      blue            47   fema… femin… Tatooine  Human   <chr>
##  8 R5-D4      97    32 NA         white, red red             NA   none  mascu… Tatooine  Droid   <chr>
##  9 Biggs…    183    84 black      light      brown           24   male  mascu… Tatooine  Human   <chr>
## 10 Obi-W…    182    77 auburn, w… fair       blue-gray       57   male  mascu… Stewjon   Human   <chr>
## # ℹ 77 more rows
## # ℹ 2 more variables: vehicles <list>, starships <list>
## # ℹ Use `print(n = ...)` to see more rows
```
Compared to a `data.frame` version of `starwars` : 

```r
# First rows of a converted tibble in data.frame
head(as.data.frame(starwars))
##             name height mass  hair_color  skin_color eye_color birth_year    sex    gender homeworld
## 1 Luke Skywalker    172   77       blond        fair      blue       19.0   male masculine  Tatooine
## 2          C-3PO    167   75        <NA>        gold    yellow      112.0   none masculine  Tatooine
## 3          R2-D2     96   32        <NA> white, blue       red       33.0   none masculine     Naboo
## 4    Darth Vader    202  136        none       white    yellow       41.9   male masculine  Tatooine
## 5    Leia Organa    150   49       brown       light     brown       19.0 female  feminine  Alderaan
## 6      Owen Lars    178  120 brown, grey       light      blue       52.0   male masculine  Tatooine
  species
## 1   Human
## 2   Droid
## 3   Droid
## 4   Human
## 5   Human
## 6   Human
##                                                                                                                                       films
## 1                                           The Empire Strikes Back, Revenge of the Sith, Return of the Jedi, A New Hope, The Force Awakens
## 2                    The Empire Strikes Back, Attack of the Clones, The Phantom Menace, Revenge of the Sith, Return of the Jedi, A New Hope
## 3 The Empire Strikes Back, Attack of the Clones, The Phantom Menace, Revenge of the Sith, Return of the Jedi, A New Hope, The Force Awakens
## 4                                                              The Empire Strikes Back, Revenge of the Sith, Return of the Jedi, A New Hope
## 5                                           The Empire Strikes Back, Revenge of the Sith, Return of the Jedi, A New Hope, The Force Awakens
## 6                                                                                     Attack of the Clones, Revenge of the Sith, A New Hope
##                             vehicles                starships
## 1 Snowspeeder, Imperial Speeder Bike X-wing, Imperial shuttle
## 2                                                            
## 3                                                            
## 4                                             TIE Advanced x1
## 5              Imperial Speeder Bike                         
## 6          
```

!!! abstract "Take a break & Read"
    For a detailed description of `tibbles`, please go read the [section 10](https://r4ds.had.co.nz/tibbles.html) 
    of "R4DS".


### Manipulating data

Some tidyverse packages such as `dplyr`, `tidyr` and `stringr` are used for manipulating 
your data. You can add/remove columns or rows, filter or manipulating strings and tables. 

!!! abstract "Take a break & Read"
    You can read the [section 3 Data Transform](https://r4ds.hadley.nz/data-transform) 
    of "R4DS" where you can see the basic commands for manipulating 
    `data.frames` and `tibbles` thanks to `dplyr`.

#### Pipe operator

In the previous section, you must have noticed some weird code: `|>`. It's called a 
*pipe operator*. Tidyverse used this weird grammar to improve code readability when
you combine multiple functions to accomplish a certain task. 

The following code lines are equivalent: 

```r
# with the pipe operator
starwars |> 
  filter(species == "Droid") |>
  head()

# without the pipe operator
head(filter(starwars, species == "Droid"))
```

But for those cases where multiple, possibly more complex, functions can be combined,
it does get confusing and difficult to read.
You can imagine "Russian dolls" using pipes when they are not stacked,
this is the best way to understand and see how many dolls we have and what they look like.

!!! abstract "Take a break & Read"
    In order to be an expert of pipe operator, please go read again [section 3.4](https://r4ds.hadley.nz/data-transform#sec-the-pipe) 
    but also [section 4.3](https://r4ds.hadley.nz/workflow-style#sec-pipes) of "R4DS" that show you the best practices.

#### Tidy data

Using tibbles is not enough to make full use of the tidyverse's functionalities, you need to
tidy your data. 
Tidy a data.frame or a tibble is a manipulation of the columns in order to obtain one column = 
one variable and one row = one observation.

For instance here is the difference between tidy and untidy data : 

```r
# Untidy data
untidy_df <- data.frame(Sample = paste0("Sample", 1:7),
                        T_cells = c(72, 0, 12, 11, 4, 10, 164), 
                        NK_cells = c(118, 24, 2, 0, 30, 4, 0),
                        Endothelial_cells = c(212, 49, 0, 29, 23, 4, 125)
                        )
untidy_df
##    Sample T_cells    NK_cells Endothelial_cells
## 1 Sample1      72         118               212
## 2 Sample2       0          24                49
## 3 Sample3      12           2                 0
## 4 Sample4      11           0                29
## 5 Sample5       4          30                23
## 6 Sample6      10           4                 4
## 7 Sample7     164           0               125
```

This format is often used when operating Excel tables, but it also has some inconveniences.
What are the numbers stands for? Potatoes? Okay, I may overstate it, but for complicated
tables it may be an issue and it makes it harder to manipulate untidy data. For example, if 
you need to visualize the number of cells for each sample but also for each cell type, it's 
not possible to do so easily in R. Instead we are going to favor this architecture: 

```r
# Tidy data
tidy_df <- untidy_df |> 
  pivot_longer(cols = contains("cells"),       # Tidy all columns that starts with "Sample"
               names_to = "Cell_types",        # Resume to a new column called "Sample"
               values_to = "Nbr_of_cells")     # Store the numeric value to a column called 
tidy_df
# A tibble: 21 × 3
##    Sample  Cell_types        Nbr_of_cells
##    <chr>   <chr>                    <dbl>
##  1 Sample1 T_cells                     72
##  2 Sample1 NK_cells                   118
##  3 Sample1 Endothelial_cells          212
##  4 Sample2 T_cells                      0
##  5 Sample2 NK_cells                    24
##  6 Sample2 Endothelial_cells           49
##  7 Sample3 T_cells                     12
##  8 Sample3 NK_cells                     2
##  9 Sample3 Endothelial_cells            0
## 10 Sample4 T_cells                     11
## # ℹ 11 more rows
## # ℹ Use `print(n = ...)` to see more rows
```

The R function `pivot_longer` was used to tidy the data.frame, because it's a tidyverse function,
the resulting value of the variable `tidy_df` is now a tibble. As you can see, we have less
columns and more rows but now each row describes one observation.

!!! abstract "Take a break & Read"
    To understand more about the power of tidy data, let's go read 
    [section 5](https://r4ds.hadley.nz/data-tidy) of "R4DS".

#### Transform data

Okay now your data is ready, you can use the pipe operator with your eyes closed, it's time to
take a closer look at `dplyr` and `stringr`. Thanks to these packages, you will be able to manipulate
and transform tables as you wish. And bonus! It will be useful when you want to visualise your data! 

!!! abstract "Take a break & Read"
    Please read carefully the [sections 12 to 19](https://r4ds.hadley.nz/transform) of "R4DS.

## Cheatsheets

You can retrieve [overviews](https://www.tidyverse.org/packages/) of all tidyverse packages,
you can also download their [cheatsheets](https://posit.co/resources/cheatsheets/), a small document that resumes all main functions and their utilisation for each package.