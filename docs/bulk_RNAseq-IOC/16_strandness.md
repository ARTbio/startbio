![](images/galaxylogo.png)
# Estimation of the strandness

In practice, with Illumina RNA-seq protocols you will most likely deal with either:

  - Unstranded RNAseq data
  - Stranded RNA-seq data produced with - kits and dUTP tagging (ISR)

This information, here called "the strandness of the libraries" should be provided by your
sequencing platform along with your FASTQ files. If you cannot find the information, ask
for it, it is always better than guessing ! If you are working on published data, the
strandness of the libraries can often be deduced from the kit reference for library preparation.

In the absence of strandness information, it is still possible to make a (very) good guess
using a tool called `Infer Experiment` from the `RSeQC` tool suite.

This tool takes the output of your mappings (BAM files), selects
a subsample of your reads and compares their genome coordinates and strands with those of
the reference gene model (from an annotation file).

Based on the strand of the genes, it can gauge whether sequencing is strand-specific,
and if so, how reads are stranded.

!!! info "a paradoxical situation"
    At this point you can notice a sort of paradox: splice-aware aligner need to know the
    strandness of the library, but if you do not know this strandness, you must use the
    "infer experiment" tool which itself will analyze a bam alignment...
    
    If so, you are right !
    
    However:
    
    1. It is not necessary to use a splice-aware aligner to generate this first set of
    alignment. A "simple" bowtie or bwa alignment will be enough and less resource-demanding.
    2. You can even use HISAT2 and specify "unstranded" for the library design. Even if it
    is not true, this is not an issue since the bam file will be re-analyzed by the tool
    `infer experiment`.
    3. You can anyway use the STAR aligner to generate this fist BAM file, since STAR does
    not require specifying whether the library is stranded or not.
    4. Finally, all your libraries have most likely the same design, unless they have not been
    generated at the same time. Therefore, analyzing a **single** sequencing dataset should
    be enough !

## Preliminary read alignment from a single sequencing dataset

To keep it sample, we are going to do this preliminary alignment using the BWA aligner.
Although BWA is not commonly used in transcriptome analysis, it will be the opportunity
for you to use it once.

- [x] First of all, create a new history and rename it `PRJNA630433 Strandness analysis`
- [x] Using the data library. Follow the main menu Shared Data --> Data Libraries -->
  IOC_bulk_RNAseq --> PRJNA630433 --> FASTQ files, and check out the  `SRR11688218` dataset.
- [x] Select the `Export to History`tab and `as Datasets`
- [x] On the next panel that pops up, select the newly created history `PRJNA630433
  Strandness analysis` (which should be already selected). Click the `Import` button.
- [x] If you are fast enough, click on the green pop up area which will return you to the
  working history where the dataset is now imported and visible. Otherwise, click on the
  house icon which will get you to the same history.

![](images/tool_small.png)

- [x] Now, select the `Map with BWA-MEM` tool and check that his version is `0.7.17.2`. If it
  is not, change it using the version icon (3 stacked cubes) at the top-right of the tool
  form.
!!! info "![](images/tool_small.png){width="25" align="absbottom"} Map with BWA-MEM tool settings"
    - Will you select a reference genome from your history or use a built-in index?
        
        --> Use a built-in genome index
    - Using reference genome
        
        --> GRCm38
    - Single or Paired-end reads
        
        --> single
    - Select fastq dataset
        
        --> SRR11688218
    - Leave other paramaters as is
    - Click the `Execute` button
    
    The tool should run about 2-3 mins before returning a green bam dataset, which you can
    examine with the eye icon
     
- [x] Before using the `Infer Experiment` tool, you will need the GTF annotation file
  "Mus_musculus.GRCm38.102.chr.gtf" from the data library `IOC_bulk_RNAseq`. Import it as
  you already did with the `SRR11688218` dataset.


## ![](images/tool_small.png){width="30" align="absbottom"} Use of `Infer Experiment` tool

Unfortunalty, the `Infer Experiment` tool does not with annotations in GTF format, but rather
in another format: the BED12 format. No worries, there is a converter tool for this.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Convert GTF to BED12` settings"
    - GTF File to convert
        
        --> Mus_musculus.GRCm38.102.chr.gtf
    - Advanced options
        
        --> Use default options
    - Click `Execute` button
    
    This will return a bed12 dataset "Convert GTF to BED12 on data 2: BED12" to be used In
    the next step

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Infer Experiment` settings"
    - Input BAM file
        
        --> Map with BWA-MEM on data 1 (mapped reads in BAM format)
    - Reference gene model
        
        --> "Convert GTF to BED12 on data 2: BED12" (the output of Convert GTF to BED12 tool)
    - Number of reads sampled
        --> 200000
        
    - Minimum mapping quality
        
        --> 30
    - `Execute`

The output of the `Infer Experiment` tool is the following text file:
```
This is SingleEnd Data
Fraction of reads failed to determine: 0.0299
Fraction of reads explained by "++,--": 0.0035
Fraction of reads explained by "+-,-+": 0.9666
```

In most cases, when a read in mapped to the (+) strand of the genome, the parental gene is
located on the (-) strand of the genome (thus "+-") or conversely when a read in mapped to
the (-) strand of the genome, the parental gene is located on the (+) strand of the genome
(thus "-+").

In conclusion, in this study, the library is **stranded** and reads correspond to the reverse
of transcribed RNA.


![](images/tool_small.png)

###  Summarize results with `MultiQC` tool

For a marginal benefit, you can also use the `MultiQC` tool
!!! info "![](images/tool_small.png){width="25" align="absbottom"} `MultiQC` settings"
    - Which tool was used generate logs?
        
        --> RSeQC
    - RSeQC output (Type of RSeQC output?)
        
        --> infer_experiment
    - Select dataset `Infer Experiment on ...`
    - `Execute`

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
