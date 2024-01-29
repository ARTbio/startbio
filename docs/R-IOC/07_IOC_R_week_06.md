## A whole new world

You already learn how to manipulate dataframe but now we are going a step further.
The `tidyverse` is a set of packages for more complex data manipulation with a nicer layout
and sometimes even more intuitive. On your marks, get set, go to the section [tidyverse](./r10_tidyverse.md) !

## Let's Practice

For each week, you'll have a set of exercises that you must render in an R script. 
After that you need to complete the following google form to answer some MCQ (Multiple
Choice Questions) where the final question is to deposit your R script.
Please note that an Rscript has the extension `.R` but it's not supported by Google Form.
To avoid this inconvenience, you need to add the `.txt` extension to make your file named as: `NAME_week6_script.R.txt`. 

![](images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

- [x] Create a tibble `my_phones` from the available data.frame `WorldPhones`. Beware of the rownames ! We don't want to lose them
- [x] Make it tidy and write it in `my_phones`
- [x] Filter the tibble to retrieve only the European (don't write it in `myphones`).
- [x] Select the tibble to retrieve only the region (don't write it in `myphones`).
- [x] Replace "." by underscore in region name
- [x] Replace truncated region names (`Amer`) by the full continent name using stringr's pattern matching functions.
- [x] Group the data based on the region
- [x] Compute the mean of number of telephones per region
- [x] Add a column with a normalized number of telephone per region (Reminder : `norm_val = val - mean(val) / sd(val)`)
- [x] Resume the information to check if the mean equal 0 and the sd equal 1. What do you get ?
- [x] Do the same for the last 4 questions but by grouping on the year you can change the created column names
- [x] Resume the information to retrieve the Year with the most phones for each region

Please be aware of the best practices for your Rscript, we will be attentive to them!

Now you can fill the following quiz: [Quiz of week 6](https://forms.gle/aNvCqoyqZbajxtZQ9).


**Thank you for your attention and see you next week :clap: :clap: :clap:**