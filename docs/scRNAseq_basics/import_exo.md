![](images/yourmove.png){width="500"}

Now it's your time to shine ! We are going to put into pratice what we 
have just seen. By using a more complicated use case, you are going to 
reproduce the whole scRNAseq analysis with Seurat.  

# Dataset test

The dataset for this analysis will be single cell RNAseq from chick embryos 
from [de Lima *et al*](https://doi.org/10.1038/s41467-021-24157-x). We are 
going to analyse only single cells of E10 stage of embryos. You can download 
the dataset at the GEO accession [GSE166981](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE166981).

??? tip "Do you need help to find the data ?"
    !!! quote "Click to reveal spoiler: <span class="spoiler">Hello ! I'm a spoiler !</span>"

    First tip : <span class="spoiler">Look to supplementary file....</span>  
    Second tip: <span class="spoiler">To see what's inside the tar archive you can click on `(custom)`</span>  
    Still no clue ? : <span class="spoiler">You need the three files with *"E10"* in the filename. You need to check the three files and then click on download. You will have the archive `GSE166981_RAW.tar` in your downloads.</span>  
    
Once you download the data, all you have to do is import it onto the Rstudio server.

## Biomart is your friend

Don't forget to import biomaRt in order to help you annotate your genes.
Make sure to tweak parameters to fit this new dataset. 

!!! question "Trouble to use biomaRt ?"
    Here is some links to help you with biomaRt :

    - [Vignette of the R package](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html)
    - [Short presentation](https://docs.google.com/presentation/d/1ck41d_0a6bMEreTfeeES67RExEc2pp_OXQvqbJ3ZdhU/edit?usp=sharing) of the use of the R package, comparing it with the Ensembl interface

# Quizz time !

Before moving on, please take this quiz :
[Quizz n°1](https://docs.google.com/forms/d/e/1FAIpQLSdEXJyXY-2pF6Lmu5YhJr1wQvBhcDAkFoV26TVGrinqFzENDQ/viewform?usp=pp_url)

To complete this quizz you'll need to :

1. Retrieve the chick dataset
2. Import data in Rstudio
3. Import data in your global environment
4. Create a Seurat Object
5. Create an annotation table of chick genes using `biomaRt`. 