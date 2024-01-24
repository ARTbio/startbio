List are a variable structure that can be related to atomic vectors. The main difference is lists
can store heterogeneous information. 

For instance, a vector can contain only one data type, either `num`, `chr`, `int`, *etc*...
A list can store variables with different structures and different types. 

You can read a perfectly well explained introduction of list from the chapter [5.7](https://rstudio-education.github.io/hopr/r-objects.html#lists) of Garett’s book.


## Manipulating list

### How to navigate within a list? 

To manipulate lists, you can use the square brackets `[ ]` but not quite the same way that 
you saw for vectors. 

Let's go with our list : 

```
mylist <- list(100:130, "R", list(TRUE, FALSE))
mylist
## [[1]]
## [1] 100 101 102 103 104 105 106 107 108 109 110 111 112
## [14] 113 114 115 116 117 118 119 120 121 122 123 124 125
## [27] 126 127 128 129 130
## 
## [[2]]
## [1] "R"
##
## [[3]]
## [[3]][[1]]
## [1] TRUE
##
## [[3]][[2]]
## [1] FALSE
```

It's a list composed of 3 elements : 

- a numeric vector
- a string character
- a list itself composed of two elements : 
    - the logical value `TRUE`
    - the logical value `FALSE`

As you might seen, when we visualize `mylist` it have different square brackets (`[ ]`) in front
of each row print in the R console. It's the first lead to be able to navigate within a list.

So if we want to retrieve the second element of a vector, we'll go `myvec[2]` and for list it's :

```
mylist[[2]]

## [1] "R"
```

And if we want to retrieve the `TRUE` value, we'll write : 

```
mylist[[3]][[1]]

## logi TRUE
```

**But why are we using a double square brackets?**

Because if we use a simple square brackets you only filter the list, *i.e.* by using simple
square brackets you'll retrieve a *smaller* list instead of the element. 

```
mylist[2]

## [[1]]
## [1] "R"

str(mylist[2]) #structure when using simple square brackets

## List of 1
##  $ : chr "R"

str(mylist[[2]]) #structure when using double square brackets

## chr "R"
```

In a nutshell : 

- `[ ]` : to filter a list
- `[[ ]]` : to retrieve an element of a list

### Usefull small functions 

If a list has "one dimension" just like vectors, you can also use all functions that
manipulate variable with one dimension. Here is the most usefull :

- `length()` : to know how many elements are in the list
-  `names()` : name each element of a list, improve the manipulation


```
length(mylist)

## [1] 3
```

#### Naming elements

To name the different elements of a list, you can use the function `names()` the same way as
you do for vectors : 

```
names(mylist) <- c("a_vector", "a_string", "a_list")

```

Or directly at the creation of the list : 

```
mylist <- list(a_vector = 100:130, 
               a_string = "R", 
               a_list = list(TRUE, FALSE))
```

**What difference does it make?**

```
mylist

## $a_vector
##  [1] 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123
## [25] 124 125 126 127 128 129 130
##
## $a_string
## [1] "R"
##
## $a_list
## $a_list[[1]]
## [1] TRUE
##
## $a_list[[2]]
## [1] FALSE
```

Now each element of the list `mylist` has a name, and you can manipulate the element not with
their position in the list but based on their name. For example : 

```
mylist$a_vector
# or
mylist[["a_vector"]]

## [1] 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123
## [25] 124 125 126 127 128 129 130
```
