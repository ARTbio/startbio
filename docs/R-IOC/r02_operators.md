In R programming, operators are symbols or special characters that perform specific
operations on data. They allow you to manipulate values, perform calculations,
compare values, and combine expressions. Here are some commonly used operators in R:

- Arithmetic Operators:
    - Addition: `+`
    - Subtraction: `-`
    - Multiplication: `*`
    - Division: `/`
    - Exponentiation: `^` or `**`
    - Modulo (remainder): `%%`
- Assignment Operators:
    - To assign value to a variable, use `<-`
      
        Note that although rarely used you can also use right assignments such as `5 -> myvar`.
        
        You may also use `=` for variable assignment, but don't do it !
        
        In order to follow the
        best practices in R programming, keep the `=` sign for argument assignment in
        functions.
        
    - To assign value to a function argument, use `=`

- Comparison Operators:
    - Equal to: `==`
    - Not equal to: `!=`
    - Greater than: `>`
    - Less than: `<`
    - Greater than or equal to: `>=`
    - Less than or equal to: `<=`
- Logical Operators:
    - Logical AND: `&` or `&&`
    - Logical OR: `|` or `||`
    - Logical NOT: `!` (for instance, the logical "different" is encoded with `!=`)

- Membership Operators:
    - `%in%`: Checks if an element is present in a vector or list.

- Miscellaneous Operators-:
    - Function call: `()`
    - Indexing: `[]`
    - Sequence generation: `:` (for instance `5:8` return `[1] 5 6 7 8`)
    - Access attributes: `$`

These operators can be used in combination with variables, literals, and expressions to
perform a wide range of operations in R programming.

By understanding and utilizing these operators effectively, you can manipulate and
transform data, perform calculations, control program flow, and make comparisons in your
R code.

You will find more examples on R operators
[here](https://www.datamentor.io/r-programming/operator).

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
