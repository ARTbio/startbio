<!-- ## Visualization with ggplot2 -->

Visualization is used throughout data analysis, from controlling distribution to presenting final results.

With R, we can use the default plotting functions from the R package <code>graphics</code>
(`plot()`, `hist()`, `boxplot()`, *etc.*).
Read more about these functions in the [chapters 11 and 12](https://bookdown.org/ndphillips/YaRrr/plotting1.html) of Philips’ book.

In this tutorial, we will introduct the <code>[ggplot2](https://ggplot2.tidyverse.org/index.html)</code> package to make more flexible and beautiful plots.


## The Compositions of A ggplot

* **data**: what to visualize
* **mapping**: the properties of a graph ("aesthetics"), *e.g.*: the abscissa, the ordinate, the legend, the facets, *etc.*
* **coordinates**: interpretation of the "aesthetics" from `x` and `y` to define the position in the graph
* **geometries**: graphical interpretation of the "aesthetics" from `x` and `y`, *e.g.*: points, lines, or polygons
* **statistics**: calculation and transformation of data, *e.g.*: counting observations for a histogram
* **scales**: graphical translation of data, *e.g.*: associate colors to a variable, modify the presenting scales of axes
* **facets**: the grouping to be carried out
* **theme**: the style of a graph

## How to Build A ggplot

> All ggplot2 plots begin with a call to `ggplot()`, supplying default data and aesthethic mappings, specified by `aes()`. You then add layers, scales, coords and facets with `+`.

```r
library("ggplot2")

str(iris) # data structure of "iris" dataset

# initiate a plot for "iris" dataset, 
# display "Sepal.Length" on the abscissa and "Petal.Length" on the ordinate
p0 <- ggplot(
  data = iris,
  mapping = aes(x = Sepal.Length, y = Petal.Length)
)
p0

# the data will be shown as dots in the graph
p1 <- p0 + geom_point()
p1

# add a linear regression model line calculated based on x and y
p2 <- p1 + stat_smooth(method = "lm")
p2

# change the breaks' position
p3 <- p2 + scale_x_continuous(breaks = seq(4, 8, by = 0.5))
p3

# show the graph by Species
p4 <- p3 + facet_wrap(facets = vars(Species))

# use the "light" theme
p5 <- p4 + theme_light()
p5
```

## A Plot With More Detail?

```r
p_box <- ggplot( # init plot
  data = iris,
  mapping = aes(x = Species, y = Petal.Length)
) +
  geom_boxplot( # add a layer of boxplot
    mapping = aes(color = Species), # colored by species
    outlier.shape = NA # hide outlier points
  ) + 
  scale_color_viridis_d(begin = 0.2, end = 0.8) + # replace boxplot color by viridis palette
  geom_point( # add a layer of dots
    position = position_jitter(seed = 123) # use jitter position to avoid overlapping
    alpha = 0.5 # make the points transparent
  ) +
  stat_summary(# add summary of average value with specified form (a red point of shape 17 and size 2)
    fun = mean, shape = 17, geom = "point", size = 2, color = "red"
  ) +
  labs( # tweak labels
    x = NULL, # remove abscissa title
    y = "Petal Length (cm)", # change ordinate title
    title = "The distribution of iris' petal length" # add a title
  ) +
  theme_minimal() + # use the minimal theme
  theme( # extra tweaks on theme
    legend.position = "none", # hide legend
    axis.text.x = element_text(face = "italic", angle = 30) # show abscissa text at 30° angle with italic font face
  )
p_box
```

Please check the [reference manual](https://ggplot2.tidyverse.org/reference/index.html) of `ggplot2` for the documentation of all functions. For more examples, please check:

* [chapter 5](https://egallic.fr/Enseignement/R/Book/graphiques.html) of a R course notes from the Aix-Marseille Université
* [chapters 2 and 3](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/making-beautiful-plots.html) of Brendan's book

## Volcano Plot & Heatmap

The volcano plot and the heatmap are two widely used figure type to show biological research results.

Check the chapter [19.11](https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html) of Sarah's book
for a concrete example of how to build a Volcano plot for differentially expression analysis results.

Heatmap need a bit more data manipulation before draw it with ggplot2.
For instance, we want to visualize a set of 10 genes of 6 samples (3 control and 3 treated):

```r
## prepare a toy dataset
set.seed(123)
exp_mat_ctrl <- matrix(rexp(30, rate = 0.1), ncol = 3)
exp_mat_trt <- matrix(rexp(30, rate = 0.8), ncol = 3)
exp_mat <- cbind(exp_mat_ctrl, exp_mat_trt)
colnames(exp_mat) <- c(
  paste0("ctrl_", 1:ncol(exp_mat_ctrl)),
  paste0("trt_", 1:ncol(exp_mat_trt))
)
rownames(exp_mat) <- paste0("gene_", 1:nrow(exp_mat))

## transform the data into "long" format (tidydata)
exp_df <- as.data.frame(exp_mat)
exp_df$gene_name <- rownames(exp_df)
# install.packages("tidyr") # we need the 'gather' function from this package
exp_df_long <- tidyr::gather(
  exp_df,
  key = "sample", # new column name to store the sample ID
  value = "exp_value", # new column name to store the value of each sample
  -gene_name # the column to skip when gathering
)

## visualize the data
p_heatmap <- ggplot(exp_df_long, aes(x = sample, y = gene_name)) +
  geom_tile(aes(fill = exp_value)) + 
  scale_fill_gradient(high = "red", low = "blue")
p_heatmap
```

There is a built-in function in R `stats::heatmap()` to draw the graph directly.
But you can have more control on the figure (style, color, position, *etc.*) if you use ggplot2.

## Other Chart Types

Pleas check the [R graph gallery](https://r-graph-gallery.com/index.html) for more (complex, even dynamic) examples of different chart types.


## Export Graphs

`ggplot2` has an implemented function `ggsave()` to export the plots in a various formats (.png, .jpeg, .pdf, .svg, *etc.*),
by default it will save the last plotted graph if you don't specify.

```r
ggsave(
  plot = p5,
  filename = "path/to/my_plot.png",
  height = 6.3, width = 4.7, units = "in", dpi = 200
)
```

## Which Type of Plot?

* One variable: boxplot, histgram, pie chart, density plot
* Two quantitative variables: scatter plot (dots plot)
* Two qualitative variables: (nested) boxplot
* One quantitative and one qualitative: boxplot, violin plot

[The eBook of Claus](https://clauswilke.com/dataviz/) is interesting to have look for the general ideas of plot type to use and how to do a better visualization (not limited to ggplot2 figures).

And you can find a the cheat sheet for ggplot2 [here](https://rstudio.github.io/cheatsheets/data-visualization.pdf)

---