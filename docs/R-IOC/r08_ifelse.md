## Conditions and If statements

The use of conditions in R is very important. A logical variable is a variable that can only take two values : `TRUE` and `FALSE`. A if statement indicates that the block code between `{ and }` will only be executed if the condition between parenthesis is fulfilled (or `TRUE`, to say it in short).


!!! abstract "Take a break & Read"
    For an introduction in the differents logical conditions and if statements please look at [W3schools If/Else parts](https://www.w3schools.com/r/r_if_else.asp). Thanks to this reference you can also test the differents notions directly in this website. 

### Comparison operators

You previously saw this kind of operators but they are really at the base of conditions and If statements. You can see a more detailled and well described documentation by reading again the 
[7.2](https://bookdown.org/ndphillips/YaRrr/logical-indexing.html) section of R pirate Guide
by Nathaniel D. Phillips.

### If Statements

Now that you masteurize the comparison operators we can look forward into the conditional
statements. 
You need to know 2 statements : `if` and `else` and 1 R function : `ifelse`. 


The structure of `if` statement will always be the same :  

```
if(cond){
    code to execute if `cond` returns `TRUE`
} else {
    code to execute if `cond` returns `FALSE`
}
```

!!! Warning
    Note that you can use `if` without `else` but never in the opposite way.
    ```
    if(cond){
        code to execute if `cond` returns `TRUE`
    } 
    ```

If you have more than two conditions, you will need to combine if statements. But if you are
familiar with python you may know the existence of `elif`, unfortunately it doesn't exist in
R. You'll write instead : 

```
my_val <- -5
if(my_val > 0){
    print("my value is a positive value.")
} else if(my_val < 0){
    print("my value is a negative value.")
} else {
    print("my value is a null value.")
}

## [1] "my value is a negative value."
```

!!! note
    In this case, R will test if our value is negative only if the test `my_val > 0` 
    returns `FALSE`. 

These statements are only working when your condition return a single value (`TRUE` or `FALSE`). 
But sometimes you'll need to test the same statement multiple times, so there is a way to use a
statement within a R command "one line" thanks to the R function `ifelse` :

```
my_vec <- 1:10
my_test <- ifelse(test = my_vec < 5, #the value in my_vec is inferior to 5
                  yes = my_vec,      #if `test = TRUE` it returns the value
                  no = 5)            #if `test = TRUE` it returns `5`

my_test
## [1] 1 2 3 4 5 5 5 5 5 5
```

The R function `ifelse` make it easier to write code when manipulating values inside a 
variable. Thanks to this, you can really imagine manipulate more complex variables such 
as dataframe or matrices. 

```
my_df <- data.frame(gene = c("A1BG", "EPC1", "MTMR7", "SLC20A2", "ZZZ3"),
                    mean_exp = c(1, 4, 0, 10, 3))
my_df
##      gene mean_exp
## 1    A1BG        1
## 2    EPC1        4
## 3   MTMR7        0
## 4 SLC20A2       10
## 5    ZZZ3        3

#create a new column based on another
my_df$is_expressed <- ifelse(my_df$mean_exp == 0,
                             "not_detected",
                             ifelse(my_df$mean_exp < 5,
                             "low detection",
                             "high_detection"))

my_df
##      gene mean_exp   is_expressed
## 1    A1BG        1  low detection
## 2    EPC1        4  low detection
## 3   MTMR7        0   not_detected
## 4 SLC20A2       10 high_detection
## 5    ZZZ3        3  low detection
```