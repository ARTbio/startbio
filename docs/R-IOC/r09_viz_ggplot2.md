<!-- ## Visualization with ggplot2 -->

Visualization is used throughout data analysis, from controlling distribution to presenting final results.

With R, we can use the default plotting functions from the R package <code>graphics</code>
(`plot()`, `hist()`, `boxplot()`, *etc.*).
Read more about these functions in the [chapter 11 and 12](https://bookdown.org/ndphillips/YaRrr/plotting1.html) of Philips’ book.

In this tutorial, we will introduct the <code>[ggplot2](https://ggplot2.tidyverse.org/index.html)</code> package to make more flexible and beautiful plots.

## The Compositions of A ggplot

* data: what to visualize
* mapping: the properties of a graph ("aesthetics"), *e.g.*, the abscissa, the ordinate, the legend, the facets, *etc.*
* coordinates: interpretation of the "aesthetics" from `x` and `y` to define the position in the graph
* geometries: graphical interpretation of the "aesthetics" from `x` and `y`, *e.g.*, points, lines, or polygons
* statistics: calculation and transformation of data, *e.g.*, counting observations for a histogram
* scales: graphical translation of data, *e.g.*, associate colors to a variable, modify the presenting scales of axes
* facets: the grouping to be carried out
* theme: the style of a graph

## How to Build A ggplot

> All ggplot2 plots begin with a call to `ggplot()`, supplying default data and aesthethic mappings, specified by `aes()`. You then add layers, scales, coords and facets with `+`.

```r
library("ggplot2")
p0 <- ggplot(
  data = iris, mapping = aes(x = Sepal.Length, y = Petal.Length)
) # initiate a plot for "iris" dataset, display "Sepal.Length" on the abscissa and "Petal.Length" on the ordinate
p0

p1 <- p0 + geom_point(aes(color = Species)) # the data will be shown as dots in the graph
p1

p2 <- p1 + stat_smooth(method = "lm") # add a linear regression model line calculated based on x and y
p2

p3 <- p2 + scale_x_continuous(breaks = seq(4, 8, by = 0.5)) # change the breaks' position
p3

p4 <- p3 + facet_wrap(facets = vars(Species)) # show the graph by Species

p4 + theme_light() # use the "light" theme
```

Please check the [reference manuel](https://ggplot2.tidyverse.org/reference/index.html) of `ggplot2` for the documentation of all functions.

For more examples, please check 
- [chapter 5](https://egallic.fr/Enseignement/R/Book/graphiques.html) of a R course notes from the Aix-Marseille Université
- [chapter 2 and 3](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/making-beautiful-plots.html) of Brendan's book


**reserve for exercises**