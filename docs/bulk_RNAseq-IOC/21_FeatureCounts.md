## ![](images/tool_small.png){width="30" align="absbottom"} Use of `FeatureCounts` tool on `PRJNA630433` datasets

Before using `FeatureCounts` ensure that you have ready:

- [x] A bam dataset, or a collection of bam file.
- [x] A GTF file corresponding to your reference genome
- [x] The knowledge of you library design (strandness, single or paired-ends and orientation of reads)

and

- [x] Copy the appropriate datasets (or collections) in a new history which you will name
  `PRJNA630433 FeatureCounts Counting on HISAT2 bam alignments`
  In the PRJNA630433 use case, this corresponds to 3 collections (Dc, Mo and Oc HISAT2 
  alignments), as well as the GTF `Mus_musculus.GRCm38.102.chr.gtf`, all present in your
  history `HISAT alignments`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `FeatureCounts` settings"
    - Alignment file
        
        --> Click on the collection icon and select `Dc HISAT2 alignments (BAM)`
    - Specify strand information
        
        --> Standed (Reverse)
    - Gene annotation file
        --> A GFF/GTF file in your history
        
        --> Mus_musculus.GRCm38.102.chr.gtf
        
    - GFF feature type filter
        
        --> exon    - GFF gene identifier
        
        --> gene_id_
    - On feature level
        
        --> No (keep default). If you select "yes", the counting will be done at the exon
            level, since your `GFF feature type filter` is `exon`
    - Output format *
        
        --> Gene-ID "\t" read-count (MultiQC/DESeq2/edgeR/limma-voom compatible)
    - Create gene-length file
        
        --> Yes
    - Does the input have read pairs?
        
        --> No, single-end
    - Advanced options
        
        --> Leave folded, no advanced options !
    - `Run Tool`

### Repeat the exact same operation twice for the collections Mo and Oc HISAT2 alignments

:bulb: use the rerun functionality !

## ![](images/tool_small.png){width="30" align="absbottom"} `MultiQC`

The MultiQC tool can use to nicely summarise the FeatureCounts countings

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `MultiQC` settings"
    - `1: Results`/ `Which tool was used generate logs?`
        
        --> FeatureCounts
    - `1: Results`/ `Output of FeatureCounts`
        
        --> Click on the collection icon, then select the three collections generated by
        featureCounts and suffixed with `: Summary`
    - click `Execute`(or `Run Tool` in the latest Galaxy version)

![](images/oeil.png){width="25" align="absbottom"} examine the results by clicking the eye
icon of the generated collection `MultiQC... ...others:Webpage`
---
