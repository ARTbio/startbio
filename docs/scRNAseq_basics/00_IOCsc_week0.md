## First Steps with Seurat

Please go read the following pages : 

- [Introduction to Single-Cell RNAseq analysis](introduction.md)  
- [Initialization of R analysis](intro_seurat.md)  
- [Import data and intialization of Seurat object](import.md)  

---

![](../R-IOC/images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

Now it's your time to shine ! We are going to put into pratice what we 
have just seen. By using a more complicated use case, you are going to 
reproduce the whole scRNAseq analysis with Seurat.  

# Dataset test

The dataset for this analysis will be single cell RNAseq from zebrafish embryos 
from [Metikala *et al*](https://doi.org/10.1371/journal.pone.0254024). You can download 
the dataset at the GEO accession [GSE152982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152982).

??? question "Do you need help to find the data ?"
    ??? tip "First tip : "
        Look to supplementary file....
        ??? tip "Second tip : "
            To see what's inside the tar archive you can click on `(custom)`
    
Once you download the data, all you have to do is import it onto the Rstudio server.

## Biomart is your friend

Don't forget to import biomaRt in order to help you annotate your genes.
Make sure to tweak parameters to fit this new dataset. 

!!! question "Trouble to use biomaRt ?"
    Here is some links to help you with biomaRt :

    - [Vignette of the R package](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html)
    - [Short presentation](https://docs.google.com/presentation/d/1ck41d_0a6bMEreTfeeES67RExEc2pp_OXQvqbJ3ZdhU/edit?usp=sharing) of the use of the R package, comparing it with the Ensembl interface


# Render your RMD/QMD

To complete this week you'll need to :

- [x] 1. Retrieve the zebrafish dataset
- [x] 2. Import data in Rstudio
- [x] 3. Import data in your global environment
- [x] 4. Create a Seurat Object
- [x] 5. Create an annotation table of zebrafish genes using `biomaRt`. 

Add your RMD/QMD in your Trello card.

**Thank you for your attention and see you next week :clap: :clap: :clap:**

----