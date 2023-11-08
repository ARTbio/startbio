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
