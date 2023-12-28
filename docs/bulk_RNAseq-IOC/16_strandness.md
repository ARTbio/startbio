![](images/galaxylogo.png)
# Estimation of the strandness

In practice, with Illumina RNA-seq protocols you will most likely deal with either:

  - Unstranded RNAseq data
  
  - Stranded RNA-seq data produced with - kits and dUTP tagging (ISR)

This information should be provided with your FASTQ files, ask your sequencing facility!

If not, try to find it on the site where you downloaded the data or
in the corresponding publication.

Another option is to estimate these parameters with a tool called `Infer Experiment` from
the `RSeQC` tool suite. This tool takes the output of your mappings (BAM files), selects
a subsample of your reads and compares their genome coordinates and strands with those of
the reference gene model (from an annotation file).

Based on the strand of the genes, it can gauge whether sequencing is strand-specific,
and if so, how reads are stranded.

## Use of `Infer Experiment` tool

----
![](images/tool_small.png)

### `Convert GTF to BED12` tool to convert the GTF file to BED

1. Go to your history `STAR` or `HISAT2`
2. Select the tool `Convert GTF to BED12`
    1. `GTF File to convert`: Drosophila_melanogaster.BDGP6.95.gtf
3. `Execute`

----
![](images/tool_small.png)

### `Infer Experiment` tool to determine the library strandness

1. In the same history `STAR` or `HISAT2`
2. Select the tool `Infer Experiment`
    1. `Input .bam file`: mapped.bam files (outputs of RNA STAR or HISAT2 tools)
    2. `Reference gene model` : BED12 file (output of Convert GTF to BED12 tool)
    3. `Number of reads sampled from SAM/BAM file (default = 200000)`: 200000
3. `Execute`

----
![](images/tool_small.png)

###  Summarize results with `MultiQC` tool

1. Select the tool `MultiQC`
2. Which tool was used generate logs?
    - `RSeQC`
3.  RSeQC output (Type of RSeQC output?)
    - infer_experiment
4. Select the 7 datasets of type `Infer Experiment on ...`
3. `Execute`

----
![](images/lamp.png)

Infer Experiment tool generates one file with information on:

- Paired-end or single-end library
- Fraction of reads failed to determine
- 2 lines:
    - For single-end
        - Fraction of reads explained by “++,–” (SF in previous figure)
        - Fraction of reads explained by “+-,-+” (SR in previous figure)
    - For paired-end
        - Fraction of reads explained by “1++,1–,2+-,2-+” (SF in previous figure)
        - Fraction of reads explained by “1+-,1-+,2++,2–” (SR in previous figure)
    
If the two “Fraction of reads explained by” numbers are close to each other (*i.e.* a mix of SF and SR),
we conclude that the library is not a strand-specific dataset (U in previous figure).

As it is sometimes quite difficult to find out which settings correspond to those of
other programs, the following table might be helpful to identify the library type:

|Library type              |Infer Experiment | TopHat           | HISAT            | htseq-count |featureCounts |
|--------------------------|-----------------|------------------|------------------|-------------|--------------|
|Paired-End (PE) - SF      |1++,1–,2+-,2-+   |FR Second Strand  |Second Strand F/FR|yes          |Forward (1)   |
|PE - SR                   |1+-,1-+,2++,2–   |FR First Strand   |First Strand R/RF |reverse      |Reverse (2)   |
|Single-End (SE) - SF      |+,–              |FR Second Strand  |Second Strand F/FR|yes          |Forward (1)   |
|SE - SR                   |+-,-+            |FR First Strand   |First Strand R/RF |reverse      |Reverse (2)   |
|PE, SE - U                |undecided        |FR Unstranded     |default           |no           |Unstranded (0)|
