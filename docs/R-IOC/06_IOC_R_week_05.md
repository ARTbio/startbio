## Improving your coding skills

You had an overview of the different main variable structures and how to manipulate them.
We propose to go further and learn how to improve your coding skills by using conditions
and manipulate several variables at the same time.

### Conditions

For a more advanced R scripting, you may want execute your code only when a specific 
requirement is fulfilled. And that's what conditions are for ! To take the mystery out
of functions `if`, `else` and `ifelse`, please go read carefully this [section](./r08_ifelse.md) 
of the reference manual.

### The `apply` family

One of the best practices to do in a code is to **never** duplicate code. Let's imagine that you
need to compute a mean for each column of a matrix, you may think about the "easiest" way is :

```r
#create a matrix
my_mat <- matrix(1:15, ncol = 3)

#compute mean for each column of my_mat
mean(my_mat[,1])
mean(my_mat[,3])
mean(my_mat[,3])
```

Fortunately `my_mat` is only composed of 3 columns, how about 300 ? Beyond the fact that it's
not optimized, it's tedious to do, you can insert errors when copy/pasting code and it makes
your script not really appealing to read and difficult to understand sometimes. We hope that
we succeed by helping you understand the importance of using the `apply` family functions. 
Take your time to read the following [section](./r11_apply.md), it's a big step ! `apply` 
is not known for being gentle with novices, but you'll see that once you understand the 
principle, it changes your life when you're coding !

## Let's Practice

For each week, you'll have a set of exercises that you must render in an R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script.
Please note that an Rscript has the extension `.R` but it's not supported by Google Form.
To avoid this inconvenience, you need to add the `.txt` extension to make your file named as: `NAME_week5_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

- [x] Create a matrix with several columns of numeric values (use `rnorm` for example) and use the apply function to calculate the max of each column.
- [x] Use the apply function to find the minimun value in each row of your matrix.
- [x] Write a function that takes a DNA sequence as input and checks if it contains any invalid characters (i.e., characters other than A, T, C, or G). If it does, print an error message, otherwise, print "Valid DNA sequence".
- [x] Write an update version of the function created in week 3 to also created false DNA sequences thanks to another parameter where when true it takes `nucleotides <- c("A", "T", "C", "G")`and otherwise `nucleotides <- c("A", "T", "C", "G", "X")`
- [x] By using apply and its subfonctions, create a list with 4 sequences where you select : 
    - [x] The length the same way as question 8 of week 3 
    - [x] The veracity randomly (no hard copy!) 
- [x] Test the validity of your sequences using apply and its subfunctions. 
- [x] Create a function called `which_season` that takes the month (integer) and returns the season
- [x] Create the variable `my_airquality` from available dataframe `airquality`.
- [x] Add the column `season` to `my_airquality` thanks to `which_season`
- [x] Compute the number of rows for each season
- [x] With `lapply`, create a list of numeric vectors (use `rnorm` for example) and calculate the sum of each vector.
- [x] For each element of this list, plot a simple histogram where you add a vertical line that represente the mean of the distribution. 
- [x] Create two numeric vectors of equal length and use mapply to calculate the element-wise product of the two vectors.
- [x] Create a custom function that takes to argument the day and month and determinate more accurately the season. Apply it for your first dataset.
- [x] Compare the result with those from `which_season`. If the result is equal, change the value in the column `season` to be in uppercase, otherwise don't change the value (or change to be in lowercase if it's already in uppercase).

Please be aware of the best practices for your Rscript, we will be attentive to them!

Now you can fill the following quiz: [Quiz of week 5](https://forms.gle/FRomLC2PCjYyydbZA).


**Thank you for your attention and see you next week :clap: :clap: :clap:**