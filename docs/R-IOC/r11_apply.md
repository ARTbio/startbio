
Functions such as `apply` and its derivatives (`lapply`, `mapply`, etc...) are 
very important in R that allow to run code more efficiently. If you are familiar
to loop in programmation, you know that it enables to run the same code repetitively.
But loops are not suitable for R, this is why we often use these functions.

## Apply

The main function is `apply` and it *applies* a function for each row and/or columns
of a two dimensional object.

!!! abstract "Take a break & Read"
    In order to learn more about `apply` please go read carefully the section 1-2
    of the [Chapter 4](https://ademos.people.uic.edu/Chapter4.html) of Erin Sovansky Winter


## Other Apply functions

If you want to apply a function on other object than a two dimensional variable,
you may be interested in `lapply` and `mapply` for example. It performs a function
for each element of a vector or a list. 

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


## for each element of mylist compute the number of item
lapply(mylist, length)
## [[1]]
## [1] 31
## 
## [[2]]
## [1] 1
##
## [[3]]
## [1] 2
```

The main difference between `lapply` and `mapply` is there is only one thing that 
differs for each lap in `lapply` (the element vector or list). Whereas for `mapply`,
you can define a different set of parameters for each lap. `lapply` returns a list
and `mapply` a vector. 

```
mylist2 <- list(test = letters[1:3], test2 = letters[4:6])

#Concatenate each vector with an underscore
lapply(X = mylist2,
       FUN = paste,
       collapse = "_")

## $test
## [1] "a_b_c"
##
## $test2
## [1] "d_e_f"

#Concatenate each vector with an underscore for the first element and a dash for the second
mapply(FUN = paste,
       mylist2,
       collapse = c("_", "-"))

##    test   test2 
## "a_b_c" "d-e-f" 
```

!!! abstract "Take a break & Read"
    In order to learn more about oher functions please go read carefully the section 3-6
    of the [Chapter 4](https://ademos.people.uic.edu/Chapter4.html) of Erin Sovansky Winter