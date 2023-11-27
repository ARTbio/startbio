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

* How to check if there are same elements from `v1` (`v1 <- c(1, 2, 3, 4, 5)`) and `v2` (`v2 <- c(8, 3, 7, 9)`)

```r
v1 <- c(1, 2, 3, 4, 5)
v2 <- c(8, 3, 7, 9)
v1 == v2 # warning of not equal length!
v1 %in% v2
```

* Do all elements from `v1` are greater than 3?

```r
v1 > 3
```

* Does any element from `v1` is greater than 8 AND any element from `v2` is greater than 8?

```r
any(v1 > 8) & any(v2 > 8)
```

* Try `c(TRUE, FALSE) & TRUE`, `c(TRUE, FALSE) & c(TRUE, FALSE)`, `c(TRUE, FALSE) && TRUE` and
`c(TRUE, FALSE) && c(TRUE, FALSE)` in the R terminal,
can you tell how to use properly `&` and `&&`? 

```r
c(TRUE, FALSE) & TRUE # element wise evaluation
c(TRUE, FALSE) & c(TRUE, FALSE)

c(TRUE, FALSE) && TRUE # error, `&&` is used to evaluate condition of length 1 
c(TRUE, FALSE) && c(TRUE, FALSE) # same error
TRUE && FALSE
```

---