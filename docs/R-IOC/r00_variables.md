### Introduction

A variable is a named (`x`, `mydata`, `results`) placeholder to store data (`12`, `"my boss"`,
`[1]  4.59 12.40 20.21 28.02 35.83 43.64 51.45 59.26 67.07 74.88`).
Variables are at the core of any programming language. Instead of directly manipulating
the data, you manipulate variables, which you can think of as abstractions of the data.
It should be noted that sometimes, R programmers use the term "object" in place of "variable".
No worries, we are still talking about variables!

A variable has a type (`integer`, `character`, *etc.*) and may have different structures
(`scalar`, `vector`, `dataframe`, *etc.*)

For an introduction to variables in R, read the section
[1.10 - Variables](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#variables){:target="_blank"}
and [1.11 - Vectors](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#vectors){:target="_blank"}
from Brendan R. E.


### Type of variable

The notion of type of a variable is pretty intuitive.
You are indeed familiar with most of the variable types.
The type of a variable may be:

1. Numeric. E.g.: `pi`, `1.0`, `2.7` or `2`
* Integer. E.g.: `1`, `2`, `45893`. Note that in R, a variable assigned to an integer value,
has the type "numeric" by default.
If you want to give it the type integer (which use less space in memory), you have to do
it actively by typing `x <- 7L` or `x <- as.integer(7)`
* Complex. E.g.: `2 + 3i`
* Character. E.g.: `"a"`, `"X"`, `"ARTbio"`, `"I had a dream"`. Note that characters are declared
as character using double quotes.
* Logical. Also called boolean. Takes only two possible value: `TRUE` and `FALSE`
* Raw. Store any piece of information as raw bytes, using the function `charToRaw()`.
For instance, `A <- charToRaw("ARTbio")`

While the type of a variable is most often obvious, you can check it out using the function `typeof()`:

```r
myvariable <- "ARTbio"
typeof(myvariable)
```
returns
```
[1] "character"
```

For a quick review on the types of the variables, see
[R Data Types](https://www.programiz.com/r/data-types){:target="_blank"}.

### Structure of variable

You may be less familiar with the notion of structure of an R variable.
Take your time here because it is important to understand that the data you place in
a variable may have different structures and how differences of structures will determine
what you can and cannot do with a variable.

The main structures of variables in R are:

1. Vectors
2. Factors
3. Matrix
4. Data Frame
5. Lists

These structures are detailed in the
[HBC training], section
["R Syntax and Data Structures"](https://hbctraining.github.io/Intro-to-R-flipped/lessons/02_introR-syntax-and-data-structures.html){:target="_blank"}.

Note that we will come back extensively to the data frames and lists later on.

Now that you have dug into the R variables you may also read the sections
[5.1](https://rstudio-education.github.io/hopr/r-objects.html){:target="_blank"} to
[5.8](https://rstudio-education.github.io/hopr/r-objects.html#data-frames){:target="_blank"}
of "Hands-On Programming with R" by Garrett Grolemund. This will recapitulate most of the
notions introduced in the "Variables" section and help you to reinforce your
comprehension of these notions.


### Scalars

The scalar represents a single value and is the simplest form of data that is
manipulated and stored in R. It can be of any existing R type
(numeric, logical, character,...). It's important to note that even though a scalar
represents a single value, R treats scalars as vectors of length 1, which means many
vector operations can be performed on scalars as well. Most of the assignments
that you have seen or manipulated so far are scalars. You can read a bit more about it
in ["YaRrr! The Pirate’s Guide to R"](https://bookdown.org/ndphillips/YaRrr/scalars.html){:target="_blank"}, by Nathaniel D. Phillips.

---
