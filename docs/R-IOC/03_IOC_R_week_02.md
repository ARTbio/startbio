## Learning the basics of R (2)

Based on what you learned last week, we will continue playing with R,
integrating a few new basic concepts and trying to import or export data in R.

### Vectors

Let's start with the simplest structure in R -- the vectors!
It is the foundation stone of other data structure in R.

Please use the [vector](r01_vectors.md) page in the reference manual to learn about its characteristics and manipulation.

### Operators

R has different groups of operators, for arthmetic operations, for logical operations, for assignment, *etc.*
What exactly are the operators and when to use them?
You can find answers on the [operator](r02_operators.md) page in the reference manual.

### Data Import & Export

With the help of previous basic concepts, we know you are now keen to use R to manipulate your data.
Wait a moment! How to import your data into R? You want to import data from other softwares into R?
The [data import and export](r07_data_import_export.md) page in the reference manual is your friend.
And you will of course find how to export data from R on this page.


## Let's Practice

For each week, you'll have a set of exercises that you must render in an R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script.
Please note that an Rscript has the extension `.R` but it's not supported by Google Form.
To avoid this inconvenience, you need to add the `.txt` extension to make your file named as: `NAME_week2_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){ style="width:75px"} **Do it yourself!**

- [x] Create factor of exam grades "A", "B", "C", "D". What is the current reference level?
- [x] Now set the grade "B" as the reference level.
- [x] The grade "D" is no longer used in exam grades, please delete it from the vector and drop this unused level.
- [x] How to check if there are same elements in `v1` (`v1 <- c(1, 2, 3, 4, 5)`) and `v2` (`v2 <- c(8, 3, 7, 9)`)
- [x] Are all elements in `v1` greater than 3?
- [x] Is any element in `v1` greater than 8 AND is any element in `v2` greater than 8?
- [x] Try `c(TRUE, FALSE) & TRUE`, `c(TRUE, FALSE) & c(TRUE, FALSE)`, `c(TRUE, FALSE) && TRUE` and
`c(TRUE, FALSE) && c(TRUE, FALSE)` in the R terminal, can you tell how to use properly `&` and `&&`? 
- [x] Download the data in any folder of your choice using this url:
https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000841/ScoringFiles/PGS000841.txt.gz
- [x] Read (*i.e.*, import into R) the downloaded file and observe what you got.
- [x] How many lines of comment (also called metadata) should we skip to get the data?
- [x] Re read the file again with appropriate parameters of `read.delim()`.
- [x] Save the readed table in `.csv` format and in Excel `.xlsx` format.
- [x] The comment lines are sometime useful, in this example we can get the information of the downloaded polygenic score (PGS). Try to read only the comment lines in R and transforme it into a `data.frame`.
- [x] Save the PGS information table in an `.RDS`.
- [x] Save both PGS score table and the information table in an `.RData`.
- [x] Save both PGS score table and the information table in a single Excel `.xlsx` file.
- [x] Read the cells A8 to C10 of the first sheet of the previous saved Excel file.

Please be aware of the best practices for your Rscript, we will be attentive to them!

Now you can fill the following quiz: [Quiz of week 2](https://forms.gle/GX4gkqkARvns1mrU6).


**Thank you for your attention and see you next week :clap: :clap: :clap:**