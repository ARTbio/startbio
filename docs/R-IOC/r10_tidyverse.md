## What's the tidyverse ?

The `tidyverse` is a set of R package for data visualisation and manipulation. You can learn 
more on their [website](https://www.tidyverse.org/). It contains the following R packages : 

- `dplyr` : A Grammar of Data Manipulation
- `tidyr` : Tidy Messy Data
- `stringr` : Simple, Consistent Wrappers for Common String Operations
- `tibble` : Simple Data Frames
- `ggplot2` : Create Elegant Data Visualisations Using the Grammar of Graphics
- `readr` : Read Rectangular Text Data
- `forcats` : Tools for Working with Categorical Variables (Factors)
- `purrr` : Functional Programming Tools

You can install and load all these packages with the following commands : 

```
install.packages("tidyverse")
library(tidyverse)
```
It's a collection that is really usefull and powerfull in data science. 

### Tibbles, the "new" data.frame

Tidyverse packages are based on the manipulation of a new type of variable, the `"tibble"`.
It's a variant of `data.frame`. Don't worry, you can manipulate `tibble` the same way as
`data.frame` that you learn previously (brackets `[ ]`, `$`, and with basics functions such as : `colnames`, `rownames`, `str`, *etc*...). 

How it looks ? 

```
# an example of a tibble from tidyverse
starwars
```

```
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

```
#First rows of a converted tibble in data.frame
head(as.data.frame(starwars))
```

```
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
    R for Data Science (O'Reilly Book tidyverse) [section 10](https://r4ds.had.co.nz/tibbles.html)


### Manipulating

R for Data Science Second Edition (O'Reilly Book tidyverse) Transform [sections 12 to 19](https://r4ds.hadley.nz/transform).