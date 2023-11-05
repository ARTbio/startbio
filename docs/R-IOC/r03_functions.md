As soon as you start learning a programming language, you hear about functions. Functions
are "self contained", named modules of code that accomplish a specific task. They usually
take in some variable(s) of a specific data structure (simple scalar, vector, dataframe,
list, etc.), process it (them), and return a result.

## Built-in functions

R has _many_ "built-in" functions and you already used some of them.

Thus,

- `charToRaw()`
- `str()`
- `typeof()`
- `class()`
- `sum()`
- `mean()`
- `max()`
- `sin()`
- `log2()`
- ...

are built-in functions provided by R as soon as you run the R environment.

A very useful function when you learn R (and later too !) is the `help()` function. Thus,
```
help("sum")
```
returns documentation on the sum() function. Note that the question mark acts as shortcut
to the help() function. Thus,
```
?sum
```
returns the same as `help("sum")`.

All functions coming from R and R packages have a help page that explains how to use them.

You will find extended explanations on the help function and help pages
[here](https://rstudio-education.github.io/hopr/packages.html#getting-help-with-help-pages).

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

- The function name (`add_three`)
    
    This is the name you will use to call the function. It should be something pretty short,
    easy to remember and as meaningful as possible.
    
    As when we create any variables or objects in R, we use the arrow `<-` to assign the
    name to the function. Indeed, a function is an R object as any R objects, and you will
    see that a function may be assigned in some occasion to a variable !
    
- the reserved token `function` and argument(s) (function`(x)`):
    
    We tell R that we want to create a function using `function()`.
    
    Within the parentheses, we can specify the number of arguments that we want our
    function to take in. It does not matter how we name the arguments (here,` x`), as
    long as we consistently use these argument names in the body of the function.
    If you need multiple arguments, it will look like this:
    ```
    mysuperfunc(arg1, arg2, arg3){
        ...
        ...
    }
    ```
    
    When you latter on call the function, you will have to specify values for the
    arguments. In this instance, this call will look like:
    ```
    myresult <- mysuperfunc("mytitle", 12.8, another_variable)
    ```
    where here the `mysupefunc` function take a string, a float number, and the content of
    the `another_variable` variable as `arg1`, `arg2` and `arg3`, respectively.
    
- The curly brackets: `{` and `}` come after `function(arguments)` and frame the
  actual function code.
    
    Note that the good practice is to put the first bracket `{` on the same line as the
    function assignment and the last bracket `}` alone on its own (last) line.
    
-  The body of the function is the code between the curly brackets.
    
    In the above example, `y`, is created to store the `x + 3` value.
    
    Note that the good practice is to indent the body of the function for readability.
    
    
- The `return` statement
    
    Here, `return(y)`)
    
    The return statement is usually at the end of the body function. This is the result
    that the function will return when executed, which can be assigned to a variable or
    directly used as an argument to another function.
    
    Here the function returns the value of y (aka, x + 3). Note that because the return
    statement is not mandatory (also we highly recommend to always use it), a function
    will return the last variable called in the function block.
    
    Here, we could have written our function as
    ```
    add_three <- function(x){
        y <- x + 3
    }
    ```
    without change, because y is the last called variable.
    
    But, let us say it again: for readability by others (and by yourself in a few weeks...)
    always put a return statement in your functions.


!!! warning ":warning:"
    be careful, when you create your own functions, to not using names reserved for R
    built-in functions such as `log`, `min`, `max`, `sum`, `mean`, `sd`
    (for standard deviation), `se` (for standard error), and many others.
    
    If you do so, you have to be aware that the new user function has precedence over the
    corresponding R built-in function and will be executed in place of it.
    If you are not sure, you can check if a name is already used by typing the name in your
    console.
    
    for instance:
    ```
    > sum
    function (..., na.rm = FALSE)  .Primitive("sum")
    > log
    function (x, base = exp(1))  .Primitive("log")
    ```
    If there is a match, then pick another name, unless you whish to mask the Built-in
    function on purpose !

**Exercises where using `help()` may help !**: 

- calculate the Ln, log in base 2 and log in base 10 of the value 1.
- round the fraction 9/7 with 2 or 4 decimal numbers


## Best Programming Practices

Like any other language, R has syntax conventions that it is not mandatory to follow in
order to get the code working but make your code readable, primarily by yourself as well
as by others.

On the one hand, don't be presumptuous: unless you're an alien, you won't be able to
understand your R codes in 2 months or even next week if you do not comment on them,
extensively stating your purposes and your algorithmic options.

On the other hand, always try to apply to the form of your code the same logic that
governs its content. This is key to its readability.

Keep also in mind that your sense of readability is rarely that of others. This is why it
is necessary to comply with rules that are developed by a large community of developers
and generally accepted by consensus.

For a nice summary of good practices in R coding, please see the section
[10](https://www.r4epi.com/coding-best-practices.html) of R for epidemiology by Brad Cannell.

You can also check the sections
[2.4](https://bookdown.org/ndphillips/YaRrr/reading-and-writing-code.html) and
[4.3](https://bookdown.org/ndphillips/YaRrr/a-brief-style-guide-commenting-and-spacing.html)
of the Nathaniel D. Phillips‘s guide and the [post](http://adv-r.had.co.nz/Style.html) of
Hadley Wickham which gives good and bad examples of R coding practices.

As a last word, we admit that the coding rules in R are less precisely described than for
other languages (Python for example, not to name it). If how to code an R instruction
seems ambiguous to you, look at what others are doing (StackOverflow is your friend)
and choose the style of the majority!

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
