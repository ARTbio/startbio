## Data Preprocessing

The preprocessing is the most important part of a single cell analysis because you
can skew your result if you filter too much **or too little** and you must really
understand what's going on these steps.

Please go read the following pages to learn more about it : [Preprocessing.](preprocessing.md).

The preprocesing is composed of : 

- Filtering of low quality barcodes
- Barcode Normalization
- Selection of most variable features

---

![](../R-IOC/images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

Now it's your time to shine ! We are going to put into pratice what we 
have just seen. By using a more complicated use case, you are going to 
reproduce the whole scRNAseq analysis with Seurat.  

# Dataset test

The dataset for this analysis will be single cell RNAseq from zebrafish embryos 
from [Metikala *et al*](https://doi.org/10.1371/journal.pone.0254024). You need
to continue your RMD/QMD from last week. 

# Render your RMD/QMD

To complete this week you'll need to :

- [x] 1. Filtering the low quality barcodes **and explain each cutoff**
- [x] 2. Normalize the data
- [x] 3. Identify the most variable genes

!!! warning "IMPORTANT"
    Please note you **must** be explicative in your cutoff choices and detailled
    each step of your thoughts. 
    In general, try to explain in your own words, each step of your analysis !

Add your RMD/QMD in your trello card.

**Thank you for your attention and see you next week :clap: :clap: :clap:**

----