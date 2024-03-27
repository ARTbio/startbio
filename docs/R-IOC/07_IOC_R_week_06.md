## A whole new world

You already learn how to manipulate dataframe but now we are going a step further.
The `tidyverse` is a set of packages for more complex data manipulation with a nicer layout
and sometimes even more intuitive. On your marks, get set, go to the section [tidyverse](./r10_tidyverse.md) !

## Let's Practice

For each week, you'll have a set of exercises that you must render in an automatic rapport. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your RMD/QMD.
Please note that an Rscript has the extension `.Rmd`/`.qmd` but it's not supported by Google Form.
To avoid this inconvenience, you need to add the `.txt` extension to make your file named as: `NAME_week6_script.Rmd.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

- [x] 1. Create a tibble `my_phones` from the available data.frame `WorldPhones`. Beware of the rownames ! We don't want to lose them
- [x] 2. Make it tidy and write it in `my_phones`
- [x] 3. Filter the tibble to retrieve only the European (don't write it in `my_phones`).
- [x] 4. Select the tibble to retrieve only the column that contains region names (don't write it in `my_phones`).
- [x] 5. Replace "." by an underscore in region name 
- [x] 6. Replace truncated region names (like `Amer`) by the full continent name using stringr's pattern matching functions.
- [x] 7. Group the data based on the region
- [x] 8. Compute the mean of number of telephones per region
- [x] 9. Add a column with a normalized number of telephone per region (Reminder : `norm_val = val - mean(val) / sd(val)`)
- [x] 10. Resume the information to check if the mean equal 0 and the sd equal 1. What do you get ?
- [x] 11. Do the same for the last 4 questions (from Q7 to 10) but by grouping on the year you can change the created column names
- [x] 12. Resume the information to retrieve the Year with the most phones for each region

Please be aware of the best practices for your Rmarkdown or quarto document, we will be 
attentive to them!

Now you can fill the following quiz: [Quiz of week 6](https://forms.gle/aNvCqoyqZbajxtZQ9).


**Thank you for your attention and see you next week :clap: :clap: :clap:**