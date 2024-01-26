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

- Miscellaneous Operators:
    - Function call: `()`
    - Indexing: `[]`
    - Sequence generation: `:` (for instance, `5:8` returns `[1] 5 6 7 8`)
    - Access attributes: `$`

These operators can be used in combination with variables, literals, and expressions to
perform a wide range of operations in R programming.

By understanding and utilizing these operators effectively, you can manipulate and
transform data, perform calculations, control program flow, and make comparisons in your
R code.

You will find more examples on R operators
[here](https://www.datamentor.io/r-programming/operator).

---