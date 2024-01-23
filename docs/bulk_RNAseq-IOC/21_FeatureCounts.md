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
        
        --> exon
    - GFF gene identifier
        
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
    - `Execute`

### Repeat the exact same operation twice for the collections Mo and Oc HISAT2 alignments

:bulb: use the rerun functionality !

### Rename your collections
For each FeatureCounts run, 3 collection are generated whose name is generic, but
self-explanatory.

:warning: As usual, take time to rename your collections such these names reflects their
specific content. In this case, this can be a prefix like `Dc FeatureCounts counts`,
`Dc FeatureCounts summary` and Dc `Dc FeatureCounts Feature Length`

# `MultiQC`

----
![](images/tool_small.png)

  1. In your history `HISAT2` or `STAR`
  2. Select the `MultiQC` tool with the following parameters:
      1. `1: Results`/ `Which tool was used generate logs?`: **STAR** `or` **HISAT2** (depending on your analysis track)
      2. `STAR or HISAT output`: shift-click select all files with the extension **.log**
      3. Click on the `+ Insert Result` button
      4. `2: Results`/ `Which tool was used generate logs?`: **featureCounts**
      5. `Output of FeatureCounts`: shift-click select all files with the extension **: Summary**
  3. `Execute`

----
![](images/oeil.png) Examine the results and ![](images/lamp.png)




