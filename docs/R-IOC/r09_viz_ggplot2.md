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

## Exercises

Let's play with the dataset `diamonds` provided in the `ggplot2` package,
it contains prices of more than 50,000 round cut diamonds, with 10 variables.
Use `?diamonds` to get the full description and `str(diamonds)` to have a glimpse of the data structure.

* Create a plot to visualize the `price` and the `carat`, colored by the quality of the `cut`

```r
library("ggplot2")
data("diamonds")
ggplot(
  data = diamonds,
  mapping = aes(x = carat, y = price)
) +
  geom_point(aes(color = cut))
```

* Change the shape and the size of the points

```r
ggplot(
  data = diamonds,
  mapping = aes(x = carat, y = price)
) +
  geom_point(aes(color = cut), shape = 2, size = 0.5) # possible shape `?pch`
```

* Create a histogram of `price` by the diamonds' `color`

```r
ggplot(data = diamonds) +
  geom_histogram(aes(x = price, color = color, fill = color))
```

* What happens if you add `position = "dodge"` in the `geom_histogram()` function ?

```r
ggplot(data = diamonds) +
  geom_histogram(
    aes(x = price, color = color, fill = color),
    position = "dodge"
  )
```

* Do the same figure but only for diamonds with prices higher than 10,000$.

```r
ggplot(data = subset(diamonds, price > 10000)) +
  geom_histogram(
    aes(x = price, color = color, fill = color),
    position = "dodge"
  )
```

* Draw a density plot of prices by group of `clarity`.

```r
ggplot(data = diamonds) +
  geom_density(aes(x = price, color = clarity))
```

* Visualize the diamonds' `carat` and width (`y`), colored by `clarity` and use `color` as facet.

```r
ggplot(data = diamonds) +
  geom_point(aes(x = carat, y = y, color = clarity)) +
  facet_wrap(facets = vars(color))
```

* Add a 2nd facet by using the `cut` and use free scales for both axes in the facets.

```r
ggplot(data = diamonds) +
  geom_point(aes(x = carat, y = y, color = clarity)) +
  facet_wrap(facets = vars(color, cut), scales = "free")
```

* What happens if you use `facet_grid()` (with appropriate arguments for facets, check `?facet_grid`) instead of `facet_wrap()`?

```r
ggplot(data = diamonds) +
  geom_point(aes(x = carat, y = y, color = clarity)) +
  facet_grid(rows = vars(color), cols = vars(cut), scales = "free")
```

### Bonus for Heatmap

* Use the previous built `p_heatmap`, try to add clustering tree on the figure.

Hints:
  - we need first have the dendrogram data
  - the R package {[ggdendro](https://andrie.github.io/ggdendro/)} can help you to draw the dendrogram data as a ggplot with `geom_segment()`
  - the R package {[patchwork](https://patchwork.data-imaginist.com)} is simple and useful to combine multiple ggplots
  (imagine we cut the plane on 4 parts:
  top-left for the sample-level dendrogram, top-right remains empty,
  bottom-left for the heatmap, bottom-right for the gene-level dendrogram
  )

```r
## compute dendrogram data
dendro_gene <- stats::as.dendrogram(
  stats::hclust(
    d = stats::dist(exp_mat, method = "euclidean"),
    method = "ward.D2"
  )
)
dendro_sample <- stats::as.dendrogram(
  stats::hclust(
    d = stats::dist(t(exp_mat), method = "euclidean"),
    method = "ward.D2"
  )
)

## build dendrogram ggplot
install.packages("ggdendro")
gg_dendro_gene <- ggplot() +
  geom_segment(
    data = ggdendro::segment(ggdendro::dendro_data(dendro_gene, type = "rectangle")),
    mapping = aes(x = y, y = x, xend = yend, yend = xend),
    linewidth = 0.5
  ) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5)))

gg_dendro_sample <- ggplot() +
  geom_segment(
    data = ggdendro::segment(ggdendro::dendro_data(dendro_sample, type = "rectangle")),
    mapping = aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5)))

## need to order samples/genes as how they were grouped in the dendrogram
p_heatmap <- ggplot(
  exp_df_long,
  aes(
    x = factor(sample, levels = colnames(exp_mat)[stats::order.dendrogram(dendro_sample)]),
    y = factor(gene_name, levels = rownames(exp_mat)[stats::order.dendrogram(dendro_gene)])
  )
) +
  geom_tile(aes(fill = exp_value)) + 
  scale_fill_gradient(high = "red", low = "blue") +
  labs(x = NULL, y = NULL)

## gather all ggplots
patchwork::wrap_plots(
  list(gg_dendro_sample, p_heatmap, gg_dendro_gene),
  design = "A#\nBC",
  guides = "collect",
  widths = c(2/3, 1/3),
  heights = c(1/3, 2/3)
)
```
