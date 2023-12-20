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

`pidigits <- c(3, 1, 4, 1, 5, 9, 3)`

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
