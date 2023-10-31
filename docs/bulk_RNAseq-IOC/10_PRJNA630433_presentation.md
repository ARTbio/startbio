# Primary myeloid cell proteomics and transcriptomics: importance of β-tubulin isotypes for osteoclast function
David Guérit, Pauline Marie, Anne Morel, Justine Maurin, Christel Verollet, Brigitte
Raynaud-Messina, Serge Urbach and Anne Blangy 

J Cell Sci (2020) 133 (10): jcs239772.

[link to the article](https://doi.org/10.1242/jcs.239772)

In this study, the authors analysed RNA expression in three different myeloid cell types
of mouse:

- Osteoclasts (OCs, 4 replicates)
- Monocyte-derived immature dendritic cells (DCs, 4 replicates)
- Bone marrow macrophages (MOs, 4 replicates) 

The raw data, 50nt single reads Illumina, were deposited to the Small Read Archive of the
EBI, under the study accession ID
[PRJNA630433](https://www.ebi.ac.uk/ena/browser/view/PRJNA630433).

We are going to use these data as a common case of RNAseq analysis.

You are expected to conduct a canonical analysis of the PRJNA630433 data, _in addition_ to
the analysis of your data.

??? info "Material and methods section of the article (this can help much)"
    RNAseq analyses
    
    The 12 RNA samples (2 μg) were processed and analyzed in parallel by Fasteris SA
    (Switzerland), according to the HiSeq Service Stranded Standard Protocol’
    (https://support.illumina.com/sequencing/sequencing_instruments/ hiseq-3000.html).
    
    The stranded mRNA libraries were sequenced by HiSeq 4000 Illumina technology,
    generating single reads of 1×50 bp. Adapter sequences were removed from the obtained
    1×50 bp reads and adapter trimmed reads were used for further analysis.
    
    About 30 million raw reads were obtained per sample (from 26,717,590 to 36,916,924),
    with around 99% of the reads mapping on reference mouse genome GRCm38.
    Multiple mapping percentages ranged between 23.62 and 32.43% according to sample
    (Fig. S1C).
    
    Sequence mapping (Mus musculus genome GRCm38, from iGenome downloaded on the
    2017-07-13), normalization and estimation of transcript abundances (FKPM) were
    performed using the Tuxedo suite of short read mapping tools (Bowtie v2.0.5, Tophat
    v2.0.6, Samtools 1.2 and Cufflinks v2.1.1).
    
    Differential expression analysis was performed with DESeq2 R package from Bioconductor
    v2.13. For each comparison by pairs, the mean of the normalized counts obtained for
    the four replicates within each group of samples was calculated as well as the log2
    fold change. The p, adjusted for multiple testing with the Benjamini-Hochberg
    procedure, was used to qualify fold changes as significant (padj<0.05).

12 fastq files can be retrieved from the SRA PRJNA630433 study.
Each fastq file constitutes a separate biological replicate of the corresponding condition
indicated in the following table.

| Condition                |replicate| id. in EBI SRA|
|--------------------------|---------|---------------|
|DC                        |1        |SRR11688218    |
|DC                        |2        |SRR11688221    |
|DC                        |3        |SRR11688224    |
|DC                        |4        |SRR11688228    |
|MO                        |1        |SRR11688219    |
|MO                        |2        |SRR11688222    |
|MO                        |3        |SRR11688225    |
|MO                        |4        |SRR11688227    |
|OC                        |1        |SRR11688220    |
|OC                        |2        |SRR11688223    |
|OC                        |3        |SRR11688226    |
|OC                        |4        |SRR11688229    |

We well also need a GTF annotation file for the Mus musculus genome, version mm38/GRCm38.
Interestingly, the link given in the mat and met section of the article for the GTF file
==does not work anymore== (the Illumina
[iGenome site](https://support.illumina.com/sequencing/sequencing_software/igenome.html)).
This well illustrate the need to reference thoroughly the external
data you are using in your scientific reports.

In this particular case, we will go to a more reliable source to retrieve the GTF for the
GRCm38 version of the Mus musculus genome: the
[Ensembl database](https://www.ensembl.org/Mus_musculus/Info/Index)

In the next section, we will use various ways to upload all these data (fastq and GTF) in
your Galaxy account.
----
