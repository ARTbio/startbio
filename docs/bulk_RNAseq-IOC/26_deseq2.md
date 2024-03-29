# Analysis of differential gene expression in PRJNA630433 using `DESeq2`

## DESeq2 Analysis

To begin, navigate to the history `PRJNA630433 FeatureCounts Counting on HISAT2 bam
alignments` and copy the three dataset collections of counts generated by FeatureCounts:
`Dc FeatureCounts counts`, `Mo FeatureCounts counts` and `Oc FeatureCounts counts` into a
new history that you will name `PRJNA630433 DESeq2 analysis`

Then, search for `DESeq2` in the tool search bar

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `DESeq2` settings"
    - how
        
        --> Select datasets per levels
    - 1: Factor
        
        --> Tissue
    - 1: Factor level
        
        Note that there will be three factor levels in this analysis: Dc, Mo and Oc.
        
        --> Oc
        
    - Counts file(s)
        
        --> select the data collection icon, then `15: Oc FeatureCounts counts`
    - 2: Factor level
        
        --> Mo
        
    - Counts file(s)
        
        --> select the data collection icon, then `10: Mo FeatureCounts counts`
    - 3: Factor level (you must click on :heavy_plus_sign: `Insert Factor level`)
        
        --> Dc
        
    - Counts file(s)
        
        --> select the data collection icon, then `5: Mo FeatureCounts counts`
    - (Optional) provide a tabular file with additional batch factors to include in the model.
        
        --> Leave to `Nothing selected`
    - Files have header?
        
        --> Yes
    - Choice of Input data
        
        --> Count data
    - Advanced options
        
        --> No, leave folded
    - Output options
        
        --> Unfold and check `Output all levels vs all levels of primary factor (use when
        you have >2 levels for primary factor)` in addition to the already checked
        `Generate plots for visualizing the analysis results`
        
        --> Leave `Alpha value for MA-plot` to 0,1: note that this option is used for
        plots and does not impact DESeq2 results
    - `Run Tool`

??? warning "Note on the order of Factors levels in the DESeq2 html form"
    As specified in the help section of the DESeq2 html form, the order of the Factors
    levels matters ! See why in that section.
    
    In a nutshell, the Factor level you put as last in the form, will be taken as the
    reference Factor level.
    
    Thus in our use case, the condition `Mo` will serve as reference condition for
    differential gene expression in the DESeq2 analysis.

## Inspect DESeq2 plots

There is a lot of information here which we will discuss online or in live

## Add a missing header to DESeq2 tabular outputs

If you have a look to three datasets in the collection `DESeq2 result files on data 4,
data 3, and others`, you'll see that a header indicating what is the content of the 7
columns is missing. This lack of header is inconfortable when you are not very familiar
with DE analyses.

Indeed, this header should be
```
GeneID	Base_mean	log2FC	StdErr	Wald-Stats	P-value	P-adj
```

Fortunately, there is a nice tool in Galaxy to quickly add this header.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Add Header` settings"
    - List of Column headers (comma delimited, e.g. C1,C2,...) 
        
        --> `GeneID,Base_mean,log2FC,StdErr,Wald-Stats,P-value,P-adj`
    - Data File (tab-delimted)
        
        --> Select the data collection icon, then `DESeq2 result files on data 4, data 3,
        and others`
    - `Run Tool`

:warning: Rename the new collection `DESeq2 Results Tables`

---
