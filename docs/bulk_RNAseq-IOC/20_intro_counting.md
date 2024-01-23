# Counting strategy ![](images/lamp.png)

![](images/readcounts.png)

## Count the number of reads per annotated gene
To compare the expression of single genes between different conditions
(e.g. with or without *Pasilla* depletion), an essential first step is to quantify
the number of reads per gene.

![](images/gene_counting.png)

From the image above, we can compute:

### Number of reads per exons

|Gene   |Exon   |Number of reads |
|-------|-------|----------------|
|gene1  |exon1  |        3       |
|gene1  |exon2  |        2       |
|gene2  |exon1  |        3       |
|gene2  |exon2  |        4       |
|gene2  |exon3  |        3       |

- The gene1 has 4 reads, not 5 (gene1 - exon1 + gene1 - exon2) because of the splicing of the last read.
- The gene2 has 6 reads (3 spliced reads)

## Counting tools
Two main tools could be used for that: HTSeq-count
([Anders et al, Bioinformatics, 2015](https://academic.oup.com/bioinformatics/article/31/2/166/2366196))
or featureCounts ([Liao et al, Bioinformatics, 2014](https://academic.oup.com/bioinformatics/article/31/2/166/2366196)).

FeatureCounts is considerably faster and requires far less computational resources.

HTSeq-counts was originally developed by Simon Anders (the developer of DESeq2 and DEXSeq).
Thus, there is certainly a guaranty of quality. However, it is less easy to use and requires
several additional steps, in particular to carefully control the GTF/GFF files taken as
an input by HTSeq-counts.

!!! note "Three important points"
    - [x] Here again (in the counting task), the genome annotation file, either a GTF or a GFF3, is
    central, regardless of the chosen counting software.
    - [x] The library strandness is absolutly required for an accurate gene read counting.
    - [x] If you are working with paired-end sequencing datasets, you will not count reads.
      Instead, you will count ==fragments==. If you count reads in this situation, you will
      overestimate gene expressions !
---
