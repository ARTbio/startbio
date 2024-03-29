# Analysis of differential gene expression in PRJNA630433 using `limma`

## limma Analysis

To begin, navigate to the history `PRJNA630433 FeatureCounts Counting on HISAT2 bam
alignments` and copy the three dataset collections of counts generated by FeatureCounts:
`Dc FeatureCounts counts`, `Mo FeatureCounts counts` and `Oc FeatureCounts counts` into a
new history that you will name `PRJNA630433 limma analysis`

Then, search for `limma` in the tool search bar

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `limma` settings"
    - Differential Expression Method
        
        --> `limma-voom`
    - Apply voom with sample quality weights?
        
        --> `No`
    - Count Files or Matrix?
        
        --> Separate Count Files
    - 1: Factor/Name
        
        --> Tissue
    - 1: Factor/1: Group
        
        Note that there will be three Groups (ie factor levels) in this analysis: Dc, Mo and Oc.
        
        --> Oc
        
    - Counts file(s)
        
        --> select the data collection icon, then `15: Oc FeatureCounts counts`
    - 2: Factor/2: Group
        
        --> Mo
        
    - Counts file(s)
        
        --> select the data collection icon, then `10: Mo FeatureCounts counts`
    - 3: Factor level (you must click on :heavy_plus_sign: `Insert Group`)
        
        --> Dc
        
    - Counts file(s)
        
        --> select the data collection icon, then `5: Mo FeatureCounts counts`
    - Use Gene Annotations?
        
        --> `No`
     - Input Contrast information from file?
       
       --> `No`
    - 1: Constrast
        
        --> `Mo-Dc`
    - 2: Constrast (click :heavy_plus_sign: `Insert Contrast`)
        
        --> `Oc-Dc`
    - 3: Constrast (click :heavy_plus_sign: `Insert Contrast`)
        
        --> `Oc-Mo`
    - Filter Low Counts
        
        --> No, leave folded
    - Output options
        
        --> Leave folded
    - Advanced options
        
        --> Put `P-Value Adjusted Threshold` to 0.1 (to be consistent with DESeq settings)
        
        --> Leave other advanced options unchanged
    - `Run Tool`

??? warning "Note on the order of Factors levels (Groups) in the limma html form"
    In contrast to DESeq2, the order of the Factors levels (Groups) does not matter with
    the limma approach.
    
    This is because here you specify manually the comparison formulas. Yet, in these
    formula, the order of the levels matters !
    
    Thus when we specify `Mo-Dc` this implies specifically that we consider the Dc as the
    reference level: we "subtract" the test level `Mo` from the reference level `Dc`

## Inspect limma plots

There is a lot of information here which we will discuss online or in live. You should
also compare these plots side by side with the plots generated by edgeR or DESeq,
especially edgeR since the format of plot reporting is very similar between limma and
edgeR.

## Here, no need for adding header to limma tabular outputs !
However, note that the loom headers are not exactly the same as the edgeR or DESeq2 headers.

---