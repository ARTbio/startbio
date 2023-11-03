A vector is a fundamental data structure in R and can be considered as a basic building
block for storing and manipulating data. In R, vectors are one-dimensional datasets
that consist of elements of the same data type. Unlike matrices and data frames,
which are two-dimensional, vectors represent a single sequence of values without any row
or column structure.

```
mycharacter_vector <- c("un", "deux", "trois", "quatre", "cinq", "six")
mynumeric_vector <- c(1, 2, 3, 4, 5, 6)
myboolean_vector <- c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
```

are character, numeric and logical vectors, respectively.

and you can watch their content just by typing and entering their names:
```
mycharacter_vector
[1] "un"     "deux"   "trois"  "quatre" "cinq"   "six"   
mynumeric_vector
[1] 1 2 3 4 5 6
myboolean_vector
[1]  TRUE FALSE FALSE  TRUE FALSE  TRUE
```

Learn more about vectors with the
[chapter 5.1](https://rstudio-education.github.io/hopr/r-objects.html#atomic-vectors)
of Garrett Grolemund’s Hands-on Programing with R’s book.
The R guide by Nathaniel D. Phillips shows how to create vectors its section
[5.2](https://bookdown.org/ndphillips/YaRrr/vectors.html) and how to manipulate them in
its sections [6 and 7](https://bookdown.org/ndphillips/YaRrr/vectorfunctions.html).

If you want to test your knowledge, don’t hesitate to do the exercises proposed by
Nathaniel D. Phillips where you can retrieve the solutions at sections
[18.2 to 18.4](https://bookdown.org/ndphillips/YaRrr/chapter-5-scalers-and-vectors.html). 


### Attributes

An attribute is a piece of information that you can attach to an atomic vector
(or any R object).
The attribute will not affect any of the values in the object, and it will not appear
when you display your object. You can think of an attribute as "metadata"; it is just a
convenient place to put information associated with an object.
R will normally ignore this metadata, but some R functions will check for specific
attributes. These functions may use the attributes to do special things with the data.

Attributes can be set and accessed individually with attr().

For instance


```
sd <- 1:10
attr(sd, "comment") <- "This is a bad name for a vector because sd is also a built-in function name in R"

sd
[1]  1  2  3  4  5  6  7  8  9 10

attr(sd, "comment")
[1] "This is a bad name for a vector because sd is also a built-in function name in R"
```

The three most important predefined attributes are:

#### 1. Names

A character vector giving each element a name. See the
[section on names](http://adv-r.had.co.nz/Data-structures.html#vector-names)
in the [Advanced R](http://adv-r.had.co.nz/) by Hadley Wickham. Names are undoubtedly
attributes that you will use extensively. Thus, let's give quickly to a vector a names attribute:


```
results <- c(1.2, 3.4, 8.02)

results
[1] 1.20 3.40 8.02

names(results) <- c("Rep1", "Rep2", "Rep3")

results
Rep1 Rep2 Rep2 
1.20 3.40 8.02
```

You can also give names to a vector when you create it:

```
other_results <- c(Rep1 = 5.6, Rep2 = 0.9, Rep3 = 5.7, Rep4 = 7.65)

other_results
Rep1 Rep2 Rep3 Rep4 
5.60 0.90 5.70 7.65
```

Finally, you can give names to a vector by creating a modified copy of a vector:


```
other_results <- setNames(other_results, c("replicat1", "replicat2", "replicat3", "replicat4"))

other_results
replicat1 replicat2 replicat3 replicat4 
     5.60      0.90      5.70      7.65
```

#### 2. Dimensions

Used to turn vectors into matrices and arrays. This will be covered later


#### 3. Class

Used to implement the S3 object system; this will be covered later too.

### Exercises on vectors
Exercises to manipulate vectors with operators and functions may be found
[here](https://www.learnbyexample.org/r-operators/)


### Factors

A factor is an R object that represents an R vector in a compact way.

To build the factor of a vector, R, behind the scene,  (i) determines the distinct values that are present in the vector, which are called _levels_ (ii) assigns a unique integer value to each level and (iii) stores for each value of the initial vector, the corresponding level, expressed with its assigned integer value.


#### Example-1

Let’s start from the vector

`pidigits <- c(3, 1, 4, 1, 5, 9, 3`)

and generate its factor pidigits_factor by


```
pidigits_factor <- factor(pidigits)
```


You may now compare the two objects:


```
> pidigits
[1] 3 1 4 1 5 9 3
> pidigits_factor
[1] 3 1 4 1 5 9 3
Levels: 1 3 4 5 9
```


At first glance, the two objects are similar, except that the `pidigits_factor` has now levels values stored along the initial vector, as a _metadata_. In addition, if you look at the structure of the pidigits_factor object, you will see that it has been built as explained above:


```
> str(pidigits_factor)
 Factor w/ 5 levels "1","3","4","5",..: 2 1 3 1 4 5 2
```


Thus, the value 2 has been assigned to the level "3" of Pi and comes at first in the initial vector, the value 1 has been assigned to the level "1" of Pi and comes at the second and fourth positions of the original vector, respectively, the value 3 has been assigned to the level "4" of Pi, etc… 


#### Example-2

You have 4 mice and you differentiate males and females in a vector


```
mice <- c("female", "male","male", "female")
```


And you build the factor object mice_factor with:


```
mice_factor <- factor(mice)
```


Then


```
> mice_factor
[1] female male   male   female
Levels: female male
```


And


```
> str(mice_factor)
 Factor w/ 2 levels "female","male": 1 2 2 1
```


As you see, the levels are "female" and "male", internally encoded by the value 1 and 2, respectively. Thus the internal sequence `1 2 2 1` returns the vector


```
[1] female male   male   female
```


To learn how to manipulate factors, go through sections [5.5.2](https://rstudio-education.github.io/hopr/r-objects.html#factors) of Garrett Grolemund’s guide and [R’s factor](https://www.geeksforgeeks.org/r-factors/) of GeeksforGeeks page. [https://hbctraining.github.io/Intro-to-R-flipped/lessons/02_introR-syntax-and-data-structures.html](https://hbctraining.github.io/Intro-to-R-flipped/lessons/02_introR-syntax-and-data-structures.html) 

Last, note that a factor is an object with two attributes, the attribute class and the attribute levels. Thus:


<table>
  <tr>
   <td><code>> class(mice_factor)</code>
<p>
<code>[1] "factor"</code>
<p>
<code>> levels(mice_factor)</code>
<p>
<code>[1] "female" "male"</code>
<p>
Or,
<p>
<code> </code>
<p>
<code>> attr(mice_factor, "class")</code>
<p>
<code>[1] "factor"</code>
<p>
<code>> attr(mice_factor, "levels")</code>
<p>
<code>[1] "female" "male"  </code>
   </td>
  </tr>
  <tr>
   <td>Note that the class attribute "factor" says R that mice_factor behave differently from regular vectors
   </td>
  </tr>
</table>


**Exercises to manipulate factors**

[To be done]


### Variable Manipulation

Now that you have mastered the differences between variable’s types and structures, you need to manipulate them.

Data manipulation 



*  O’Reilly Book 1 Hands-On Programming with R :
    * [https://rstudio-education.github.io/hopr/r-notation.html#selecting-values](https://rstudio-education.github.io/hopr/r-notation.html#selecting-values) chapter 6.1 manipulating vectors
    * [https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets](https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets) chapter 6.4 : subsetting lists with brackets and $
* R pirate Guide Nathaniel D. Phillips : 
    * [https://bookdown.org/ndphillips/YaRrr/vectorindexing.html](https://bookdown.org/ndphillips/YaRrr/vectorindexing.html) vectors and brackets chapter 7
    * [https://bookdown.org/ndphillips/YaRrr/slicing-dataframes.html](https://bookdown.org/ndphillips/YaRrr/slicing-dataframes.html) matrices and dataframes with brackets and basic functions, chapter 8.5 and 8.6

---
## Operators

In R programming, operators are symbols or special characters that perform specific operations on data. They allow you to manipulate values, perform calculations, compare values, and combine expressions. Here are some commonly used operators in R:



* Arithmetic Operators:
    * Addition: +
    * Subtraction: -
    * Multiplication: *
    * Division: /
    * Exponentiation: ^ or **
    * Modulo (remainder): %%
* Assignment Operators:
    * To assign value to a variable, use <- \
Note that although rarely used you can also use right assignments such as 5 -> myvar. You may also use = for variable assignment, but in order to follow the best practices in R programming, don’t do it and keep the = sign for argument assignment in functions.
    * To assign value to a function argument, use =
* Comparison Operators:
    * Equal to: ==
    * Not equal to: !=
    * Greater than: >
    * Less than: <
    * Greater than or equal to: >=
    * Less than or equal to: <=
* Logical Operators:
    * Logical AND: & or &&
    * Logical OR: | or ||
    * Logical NOT: !
* Membership Operators*:
    * %in%: Checks if an element is present in a vector or list.
    * 
* Miscellaneous Operators*:
    * Function call: ()
    * Indexing: []
    * Sequence generation: :
    * Access attributes: $

These operators can be used in combination with variables, literals, and expressions to perform a wide range of operations in R programming. By understanding and utilizing these operators effectively, you can manipulate and transform data, perform calculations, control program flow, and make comparisons in your R code.

You will find more examples on R operators [here](https://www.datamentor.io/r-programming/operator).

**Exercises to manipulate operators**


## R Functions

As soon as you start learning a programming language, you hear about functions. Functions are "self contained", named modules of code that accomplish a specific task. They usually take in some variable(s) of a specific data structure (simple scalar, vector, dataframe, list, etc.), process it (them), and return a result.


### Built-in functions

R has _many_ "built-in" functions and you already used some of them.

Thus, charToRaw(), typeof(), class(), sum(), mean(), max(), sin(), log2(), etc. are built-in functions provided by R as soon as you run the software.

A very useful function when you learn R (and later too !) is the help() function. Thus,

help("sum") will return documentation on the sum() function. Note that there is a specific shortcut for the help() function, the question mark. Thus ?sum returns the same as help("sum"). All functions coming from R and R packages have a help page that explains how to use the function. You may find extended explanations on the help function and help pages [here](https://rstudio-education.github.io/hopr/packages.html#getting-help-with-help-pages).


### Custom user functions

Importantly, you can also create, use and reuse your own functions.

A typical custom function has the following structure:


```
add_three <- function(x){
    y <- x + 3
   return(y)
}
```


The important elements here are:



* Function name (`add_three`): this is the name that you want to call your function. It should be something pretty short, easy to remember and as meaningful as possible. As when we create any variables or objects in R, we use the arrow <- to assign this name to our function… Yes, a function is an R object as any R objects, and you will see that a function may be used in some occasion as a variable !
* "function" and arguments (`function(x)`): we tell R that we want to create a function using `function()`. Within the parentheses, we can specify the number of arguments that we want our function to have. It doesn’t matter what we name our arguments within the parentheses (here,` x`), as long as we use the same names in the body of the function. If you want to have multiple arguments, it would look something like this: `function(arg1, arg2, arg3, ...)`. Later, when you put your function to use, you’ll have to specify values for the arguments.
* Curly brackets: `{` and `}` come after `function(argument)` and need to bracket the actual function code that you’re writing. Note that the best practice is to put the first bracket` `{ on the same line as the function assignment and the last bracket } alone on its own (last) line.
* Body of the function: this is the code in the function between the curly brackets that executes the task that you want. A new variable, y, is created to store the x + 3 value. Note that a good practice is to indent the body of the function for better readability.
* The return value (`return(y)`): Also inside the curly brackets, but usually at the end, this is the result that the function prints for you when it’s done running. Here the function returns the value of y (aka, x + 3). Note that the return statement is not mandatory, when this statement is missing, the function will return the content of the last variable manipul

_Tips: be careful when you create your own variables or functions, do not use reserved names that refer to R built-in functions like log, min, max, sum, mean, sd (for standard deviation), se (for standard error), and many others. If you do so, the function will be replaced by the value you assigned it for._

_If you are not sure, you can check if a name is already used by typing the name on your console. If there is a match, then pick another name ;)_

**Exercises on using help**: 



* calculate the log, log in base 2 and log in base 10 of the value 1.
* round the fraction 9/7 with 2 or 4 decimal numbers 


## Best Programming Practices

Like any other language, R has syntax conventions that it is not mandatory to follow in order to get the code working but make your code readable, primarily by yourself as well as by others.

On the one hand, don't be presumptuous: unless you're an alien, you won't be able to understand your R codes in 2 months or even next week if you do not comment on them, extensively stating your purposes and your algorithmic options.

On the other hand, always try to apply to the form of your code the same logic that governs its content! This is key to its readability. However, keep in mind that your sense of readability is rarely that of others. This is why it is necessary to comply with rules that are developed by a large community of developers and generally accepted by consensus.

For a nice summary of good practices in R coding, please see the section [10](https://www.r4epi.com/coding-best-practices.html) of R for epidemiology by Brad Cannell. You can also check the sections [2.4](https://bookdown.org/ndphillips/YaRrr/reading-and-writing-code.html) and [4.3](https://bookdown.org/ndphillips/YaRrr/a-brief-style-guide-commenting-and-spacing.html) of the Nathaniel D. Phillips‘s guide and the [post](http://adv-r.had.co.nz/Style.html) of Hadley Wickham which gives good and bad examples of R coding practices.

As a last word, we admit that the coding rules in R are less precisely described than for other languages (Python for example, not to name it). If how to code an R instruction seems ambiguous to you, look at what others are doing (StackOverflow is your friend) and choose the style of the majority!

**Exercice yourself with session2_part1 exercices**



* Create a variable called `surname` with the string _Marilyn Monroe_
* How would you name a variable containing the mean temperature ?** **
* Which variable names are correct?
    * 6h_day1 <- 
    * day1_6h
    * day16h
    * Day1 6h
    * 6h@day1
* Which assignments are correct?
    * `x <- 1,4`
    * `x <- 5`
    * `y <- x + 3`
    * y < -`
    * `a <- Marilyn`
    * `wishes <- "happy birthday"`
    * `val <- TRUE`
    * `TEMP <- `
    * `pi <- 3`
* What will be the value of `a` after this code? 

```

a <- 10

a + 10

a 

```



* Are those codes correct? 

``` 

a <- 5

b <- 17

ab <- a * b 

```


## Two dimensional objects


### Matrices

A matrix is a two-dimensional data structure that stores data in a grid-like format consisting of rows and columns. It is a data structure used for mathematical computations, linear algebra operations, and storing homogeneous data. Learn more about matrices with chapter [5.3](https://rstudio-education.github.io/hopr/r-objects.html#matrices) of Garett’s book.


### Data Frames

A data framematrix is a two-dimensional data structure that organizes data into rows and columns, similar to a table or spreadsheet. It is particularly useful for working with heterogeneous datasets where different columns can have different data types.that stores data in a grid-like format consisting of rows and columns. It’s composed of vectors arranged in columns, see chapter [5.8](https://rstudio-education.github.io/hopr/r-objects.html#data-frames) of Garett’s book and sections [8.2.3 to 8.2.4](https://bookdown.org/ndphillips/YaRrr/creating-matrices-and-dataframes.html#data.frame) of Philips’ book for more details. 


### Manipulating Matrices and Dataframes

Please read with attention the sections [8.3 to 8.6](https://bookdown.org/ndphillips/YaRrr/matrix-and-dataframe-functions.html) of Philips’ book. 

**Exercises to manipulate matrices and dataframes**


## Lists

Lists are variables that store heterogeneous information in one direction. chapter [5.7](https://rstudio-education.github.io/hopr/r-objects.html#lists) of Garett’s book

Manipulating list : 



* [https://data-flair.training/blogs/r-list-tutorial/](https://data-flair.training/blogs/r-list-tutorial/) blog post by 
