## Learning the two dimensional objects

Now you are familiar with the basics of R,
we will learn two more complexe data structures, the `matrix` and the `data.frame`.


### Matrices

A matrix is a fundamental two-dimensional data structure that organizes data into rows and columns.
Matrices are homogeneous, meaning they store elements of the same data type, making them efficient for mathematical operations.

Please check the [matrix](r06_df_matrices.md) part in the reference manual to learn more about it.

### Data Frames

The data frame is another essential two-dimensional data structure in R.
Unlike matrices, `data.frame` can represent heterogeneous data.

You can find how to create and manipulate it in the [data.frame](r06_df_matrices.md) part of the reference manual.


## Let's practice

For each week, you'll have a set of exercices that you must render in a R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script. Please note that
a Rscript has the extension `.R` but it's not supported by Google Form so you must add
the extension `.txt` so your filename will be : `NAME_week4_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself !**

- [x] Create an object of matrix named `my_mat` with 3 rows and 4 columns, fill with numbers 1 to 12 by row,
name the rows with "r1", "r2", "r3" and the columns with "c1", "c2", "c3", "c4".
- [x] Extract the 2nd row of `my_mat`.
- [x] Extract the 2nd row of `my_mat` but keep it in matrix format.
- [x] Extract the 2nd row of `my_mat` using a logical vector.
- [x] What are the positions for the numbers that are multiples of 3 in `my_mat`?
- [x] Based on `my_mat`, add a column "c5" containing the values "a", "b", "c". What happens after this add?
- [x] Now delete the added column of `my_mat` and convert the matrix to numeric mode.
- [x] Replace the element bigger than 10 by 99 in `my_mat`.
- [x] Transforme the matrix `my_mat` to a `data.frame` named `my_df`.
- [x] Use the rownames to create a new column "id" for `my_df`.
- [x] Which rows has duplicated values in `my_df`?
- [x] Create a new column named "total" in `my_df`, which calculates the sum of column "c1" to "c4" by row.
- [x] Change the column order to put the "id" in the first column in `my_df`.
- [x] Remove the rownames of `my_df`.
- [x] Add a new row in `my_df` which contains the sum of each column (except the "id" column, put `NA` in the new row for this column).

Please be aware of the best practices for your Rscript, we will be attentive to them !

Now you can fill the following quiz: [Quizz of week 4](https://forms.gle/9ge6VxjL9dmFapJK6).


**Thank you for your attention and see you next week :clap: :clap: :clap:**