## Séance 5 : _Advanced manipulation of Data Objects_



* Filtration, sélection et structuration de données d’intérêt
* _dplyr_
* Tutorials : 
    * Brendan R. E. Ansell : Introduction to R - tidyverse
        * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/manipulating-data-with-dplyr.html](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/manipulating-data-with-dplyr.html) : chapter 5 manipulating data with dplyr

Séance 5  : _Structures de contrôle_



* Condition if/else
    * https://www.w3schools.com/r/r_if_else.asp
* Exercices

Séance 6 et 7 : _Les fonctions_



* Création de nouvelles fonctions
* Fonctions avancées : _apply_ and co
* Tests et exceptions
* Exercices
* Tutorials
    * Notes de cours de R - Version 2, Ewen Gallic  : 
        * Apply & co : [https://egallic.fr/Enseignement/R/Book/boucles.html#boucles-vectorisation-apply](https://egallic.fr/Enseignement/R/Book/boucles.html#boucles-vectorisation-apply)
    * R for Data Science (O'Reilly Book tidyverse) 1er : 
        * Custom function : [https://r4ds.had.co.nz/functions.html](https://r4ds.had.co.nz/functions.html) chapter 19

Séance 8 et 9 : _Les graphiques en R_



* concepts de base et fonctions génériques
* _ggplot2_
* Tutorials : 
    * Brendan R. E. Ansell : Introduction to R - tidyverse
        * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/making-beautiful-plots.html](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/making-beautiful-plots.html) chapter 4.5
    * Notes de cours de R - Version 2, Ewen Gallic (français) : 
        * [https://egallic.fr/Enseignement/R/Book/graphiques.html](https://egallic.fr/Enseignement/R/Book/graphiques.html) : chapter 5
    * Introduction to R - tidyverse Brendan R. E. Ansell 
        * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/what-the-factor.html](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/what-the-factor.html) : rappel factor + factor in plots  → section 8
* Exercices

Séance 10 : Récapitulatif général du savoir faire pour application à l'analyse SingleCell RNA-Seq


        Data structure - Factor : [https://r4ds.hadley.nz/factors.html](https://r4ds.hadley.nz/factors.html) : R for Data Science (O'Reilly Book tidyverse) 2e : Chapter 17



    * Tutorials : 
        * O’Reilley Book 1 Hands-On Programming with R
            * [2 The Very Basics | Hands-On Programming with R](https://rstudio-education.github.io/hopr/basics.html#objects)
            * [https://rstudio-education.github.io/hopr/r-objects.html#atomic-vectors](https://rstudio-education.github.io/hopr/r-objects.html#atomic-vectors)
        * Introduction à R (Français)
            * [https://egallic.fr/Enseignement/R/Book/introduction.html](https://egallic.fr/Enseignement/R/Book/introduction.html) (Chapter 1)
            * 
        * Brendan R. E. Ansell : Introduction to R - tidyverse
            * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#variables](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#variables)
            * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#vectors](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#vectors)
        * Guide R Joseph Larmarange (Français)
            * [https://larmarange.github.io/guide-R/bases/vecteurs.html](https://larmarange.github.io/guide-R/bases/vecteurs.html) + utile pour les fonctions de bases aussi
            * 


## Session 3

Best practices : 



* [https://www.r4epi.com/coding-best-practices.html](https://www.r4epi.com/coding-best-practices.html) R for epidemiology - Brad Cannell, chapter 10
* Advanced R by Hadley Wickham [http://adv-r.had.co.nz/Style.html](http://adv-r.had.co.nz/Style.html)

Séances 2 et 3 : _Les objets R_



* Conventions & bonnes pratiques
    * Revenir sur le nommage des variables
    * Syntaxe virgule, espaces…
    * Pas de `?fonction` dans le script
    * Les commentaires
* Présentation et manipulation des différents types d’objets : vecteurs, matrices, data frames et listes
* Object R
    * Vecteur
        * Référencement des éléments
        * names
        * length
    * Factor
        * Levels
        * str
    * Matrices 
        * Dimensions [ , ]
    * Dataframes 
        * Référencement des éléments de la dataframe
        * colnames/rownames
        * Dim
        * 
    * List
        * [[ ]], $  
        * Names
        * unlist
    * Tutorials : 
        * O’Reilly Book 1 Hands-On Programming with R : 
            * [https://rstudio-education.github.io/hopr/r-objects.html](https://rstudio-education.github.io/hopr/r-objects.html) (Chapter 5)
            * [https://rstudio-education.github.io/hopr/r-notation.html#selecting-values](https://rstudio-education.github.io/hopr/r-notation.html#selecting-values) chapter 6.1 manipulating vectors
            * [https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets](https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets) chapter 6.4 : subsetting lists with brackets and $
        * Brendan R. E. Ansell : Introduction to R - tidyverse 
            * [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#data-frames](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#data-frames) : dataframes
            * Factors but it may be too early for these explanations : [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/what-the-factor.html](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/what-the-factor.html)
        * Notes de cours de R - Version 2, Ewen Gallic : 
            * [https://egallic.fr/Enseignement/R/Book/donn%C3%A9es.html#structures-de-base](https://egallic.fr/Enseignement/R/Book/donn%C3%A9es.html#structures-de-base)
            * [https://egallic.fr/Enseignement/R/Book/donn%C3%A9es.html#manipulation-des-donn%C3%A9es](https://egallic.fr/Enseignement/R/Book/donn%C3%A9es.html#manipulation-des-donn%C3%A9es)
        * R for Data Science (O'Reilly Book tidyverse) 2e : 
            * Factors : [https://r4ds.hadley.nz/factors.html](https://r4ds.hadley.nz/factors.html)
* Exercices

Data structure - Data-frames : [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#data-frames](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#data-frames) : Brendan R. E. Ansell : Introduction to R - tidyverse - Chapter 1.14

[https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets](https://rstudio-education.github.io/hopr/r-notation.html#dollar-signs-and-double-brackets) 



* Fonctions de bases
    * Son fonctionnement
        * Usage
        * Paramètre
        * Détails
        * Exemple
    * ?, help
    * typeof()
    * Vecteurs, c()
    * Head, length, str, rep, seq, table
    * Mean, min, max, paste, 
    * is.[...] 
    * Tutorials : 
        * O’Reilley Book 1 Hands-On Programming with R
            * [https://rstudio-education.github.io/hopr/basics.html#functions](https://rstudio-education.github.io/hopr/basics.html#functions) : chapter 2.3 + 2.5 (chapter 2.4 seems too early, it’s about custom functions)
        * Brendan R. E. Ansell : Introduction to R - tidyverse : [https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#help](https://bookdown.org/ansellbr/WEHI_tidyR_course_book/welcome-to-r.html#help) : How to use the help section (chapter 1.9) very short
        *  R for Data Science (O'Reilly Book tidyverse) 2e : 
            * [https://r4ds.hadley.nz/workflow-basics.html](https://r4ds.hadley.nz/workflow-basics.html) (chapter 3.1-3.3 variable assignation and small intro & chapter 3.4 small intro and examples about functions)
            * 
    * Exercices
* Quizz général sur les notions de la session 1

idées questions



1. Which variable names are incorrect : 


## Packages

[[https://hbctraining.github.io/Intro-to-R-flipped/lessons/04_introR_packages.html](https://hbctraining.github.io/Intro-to-R-flipped/lessons/04_introR_packages.html)]

R works in a model of **package** **“add-ons”** where packages are collections of functions, data sets, and other resources that extend the functionality of the base R system.

As R is used by many communities (biologists, economists, statisticians, meteorologists…), several **R packages **were designed to address specific tasks or provide specialized tools for the different domains or areas of interest. 

Therefore, many (most) of R packages are available in the official **[CRAN package repository](https://cran.r-project.org/web/packages/available_packages_by_date.html)** but we will also use packages from** [Bioconductor](https://bioconductor.org/) **(open source software for Bioinformatics). 

When you start an R session, default packages are automatically loaded allowing you to use R as a calculator and access basic mathematical functions. Some of the key default packages include:



* **base**: It provides basic data structures (vectors, matrices, data frames, lists), operators for calculations (+, -, *, /), control structures (if-else, loops), and functions for mathematical operations, such as log(), sin(), cos(), and many more.
* **stats**: This package contains functions for statistical analysis and modeling. It includes methods for probability distributions, hypothesis testing, regression, analysis of variance, and other statistical procedures. Some commonly used functions from this package are mean(), sd() or t.test().
* **graphics**: The graphics package provides functions for creating plots and visualizations. It includes functions for basic plots (scatter plots, bar plots, histograms), customization options, and annotation features. The functions plot(), hist() or boxplot() are examples of functions available in this package.
* **utils**: This package contains utility functions that assist in data manipulation, file handling, and other miscellaneous tasks. Functions like head(), tail(), str(), and install.packages() are part of this package._ _

To be able to use the functionalities of a R package, you need to load the package in your current R environment. Here is a small video to show you how to install, load R packages.


### How to install an R package
