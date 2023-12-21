## Data visualization with ggplot2

We reach out almost the end of this IOC, where you've laid the groundwork with fundamental concepts.
Let's spice up your skill set with ggplot2!

`ggplot2` is the magic wand for creating cool charts and graphs,
it introduces you to the artistry of crafting meaningful and expressive visualizations.
It's not just about stats; it's about making your data speak visually. So, gear up to add that extra flair to your reports and impress your audience with data storytelling!

Unlock the power of data visualization by using the [ggplot2](r09_viz_ggplot2.md) page in the reference manual.


## Let's practice

For each week, you'll have a set of exercices that you must render in a R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script. Please note that
a Rscript has the extension `.R` but it's not supported by Google Form so you must add
the extension `.txt` so your filename will be : `NAME_week7_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself !**

#### Basic exercises with `diamonds`

Let's play with the dataset `diamonds` provided in the `ggplot2` package,
it contains prices of more than 50,000 round cut diamonds, with 10 variables.
Use `?diamonds` to get the full description and `str(diamonds)` to have a glimpse of the data structure.

- [x] Create a plot to visualize the `price` and the `carat`, colored by the quality of the `cut`.
- [x] Change the shape and the size of the points.
- [x] Create a histogram of `price` by the diamonds' `color`.
- [x] How to make the bars in histogram side by side?
- [x] Do the same figure but only for diamonds with prices higher than 10,000$.
- [x] Draw a density plot of prices by group of `clarity`.
- [x] Visualize the diamonds' `carat` and width (`y`), colored by `clarity` and use `color` as facet.
- [x] Add a 2nd facet for the `cut`, how to make the scales vary across both columns and rows?

#### Bonus for heatmap

- [x] Use the previously built `p_heatmap` from the [ggplot2 reference](r09_viz_ggplot2.md), try to add clustering tree on the figure.

Hints:
  - we need first have the dendrogram data
  - the R package {[ggdendro](https://andrie.github.io/ggdendro/)} can help you to draw the dendrogram data as a ggplot with `geom_segment()`
  - the R package {[patchwork](https://patchwork.data-imaginist.com)} is simple and useful to combine multiple ggplots
  (imagine we cut the plane on 4 parts:
  top-left for the sample-level dendrogram, top-right remains empty,
  bottom-left for the heatmap, bottom-right for the gene-level dendrogram
  )
 

Please be aware of the best practices for your Rscript, we will be attentive to them !

Now you can fill the following quizz : [Quizz of week 7]().


**Thank you for your attention and see you next week :clap: :clap: :clap:**