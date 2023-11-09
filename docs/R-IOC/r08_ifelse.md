## Conditions and If statements

The use of conditions in R is very important. It allows us to use logical variables 
and "if" statements to manipulate variables or even to execute R code R code only 
in certain modalities.

!!! abstract "Take a break & Read"
    For an introduction in the differents logical conditions and if statements please look at [W^3schools If/Else parts](https://www.w3schools.com/r/r_if_else.asp). Thanks to this reference you can also test the differents notions directly in this website. 

### Comparison operators

You previously saw this kind of operators but they are really at the base of conditions and If statements. You can see a more detailled and well described documentation by reading again the 
[7.2](https://bookdown.org/ndphillips/YaRrr/logical-indexing.html) section of R pirate Guide
by Nathaniel D. Phillips.

### If Statements

Now that you masteurize the comparison operators we can look forward into the conditional
statements. 
You need to know 3 R functions : `if`, `else` and `ifelse`. 


The structure of `if` statement will always be the same :  

```
if(cond){
    code to execute if `cond` returns `TRUE`
} else {
    code to execute if `cond` returns `FALSE`
}
```

But there is a way to use a statement within a R command "one line" thanks to the R function `ifelse` :

```
my_vec <- 1:10
my_test <- ifelse(test = my_vec < 5, #the value in my_vec is inferior to 5
                  yes = my_vec,      #if `test = TRUE` it returns the value
                  no = 5)            #if `test = TRUE` it returns `5`

my_test
## [1] 1 2 3 4 5 5 5 5 5 5
```

Thanks to this command, you can really imagine manipulate more complex variables such as
dataframe or matrices. 

```

```