## Learning the basics of R

You are more familiar with the different tools for this IOC, especially RStudio. 
It's good because now we can go to the heart of the matter, R !

### Variables

We'll start slowly but surely by learning what's a variable. You need to understand
what's a variable, how to create one and manipulate them. It's a key concept for R. 

You'll need to go read carefully the [variable](./r00_variables.md) page in the 
reference manual.

### Functions

The other key concept of R is the use of functions. An R function is a set of R command
that achieve a specific task. 

To learn more about them, how to use and create a function, go to the [function](./r03_functions.md) page in the reference manual. 

### Best Practices

Last but not least for this week, the best practices of R programming. Some people may
disagree but for us in ARTbio, it's a very important topic. R is a language and like
all languages, there are rules to correctly write it.
Those rules are called [best practices](./r04_bestpractices.md) and you can learn more on their page in the reference manual.

### Bonus: Automatic Reporting (R Markdown/Quarto)

Do you want to simplify your life by generating automatic reports?
As a bonus, we'll introduce you a simple and practical way -- using R Markdown (or Quarto)!

R Markdown and Quarto provides a streamlined and efficient way to generate dynamic and reproducible reports in data analysis and research.
They are both markup languages that **integrate code, text and output in a single document**.
In addition, you can choose different output formats (Word, Powerpoint, PDF, HTML, *etc.*) to write reports, presentations or even articles.
With the ability to incorporate R code directly into the document, these tools ensuring that reports can be easily updated with new data or changes in analysis.

Please check the [automatic reporting](./r11_auto_reporting.md) page in the reference manual for more details.

## Let's Practice

For each week, you'll have a set of exercises that you must render in an R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script.
Please note that an Rscript has the extension `.R` but it's not supported by Google Form.
To avoid this inconvenience, you need to add the `.txt` extension to make your file named as: `NAME_week1_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

- [x] Create a variable called `my_var` that contain your favorite color.
- [x] Create a variable called `surname` with the string _Marilyn Monroe_.
- [x] Create a variable with the number 9.
- [x] What's its type?
- [x] Change it to character.
- [x] Calculate the Ln, log in base 2 and log in base 10 of the value 1.
- [x] Round the fraction 9/7 with 2 and then 4 decimal numbers.
- [x] Create a function that takes a value and substract the number 4.
- [x] Test your function for the values : 12, 5.6 and 0.

Please be aware of the best practices for your Rscript, we will be attentive to them!

Now you can fill the following quiz: [Quiz of week 1](https://forms.gle/Y6enoxKSH5Nfa14w9).


**Thank you for your attention and see you next week :clap: :clap: :clap:**

### To go further

You need more practice? You can test your R with the amazing Pirate's Guide to R of
Nathaniel D. Phillips :

- [Exercices for basic R](https://bookdown.org/ndphillips/YaRrr/test-your-r-might.html)
- [Exercices for custom function](https://bookdown.org/ndphillips/YaRrr/test-your-r-might-6.html)

Solutions are available in [Chapter 18](https://bookdown.org/ndphillips/YaRrr/solutions.html).