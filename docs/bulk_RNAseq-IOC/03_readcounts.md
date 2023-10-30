<center>
![](images/lamp.png)

**==In Reference-based RNAseq analysis, the read counts are treated as proxies to RNA steady state levels==**
</center>

![](images/readcounts.png)

## This idea implies 3 steps

### 1. Map reads to a reference genome with aligners such as:

- [x] TopHat2
- [x] HiSat2
- [x] STAR

→ These aligners are “splice aware”

→ They generate a **BAM Alignment file**

### 2. Use a software tool and genome annotations to count the reads spanning genes or transcripts

- [x] Input files for the counting tools are mostly the **BAM Alignment file** generated at step `1.`
- [x] Mostly used counting tools:
    - FeatureCounts
    - HTseqCounts
    - Some aligners (STAR, Salmon, ...) are embedding their own counting procedure
- [x] Genome annotations
    
    They are an **essential** element of the whole analysis: *Keep in mind* that you
    only count what you told the counting tool to count!
    
    Most often genome annotations are provided in a GTF format. GFF3 or BED12 may also be
    used by some counting tools.
    
    Given the importance of annotations in Reference-based Expression analysis it is highly
    recommended that you take some time to well understand the structure of a GTF file, its
    links with the genome version as well as its own version. In particular annotation versions
    evolve more rapidly that the genome version to which they are linked. There are several
    years between different releases of genome assembly, whereas different versions of the
    genome annotation may be separated by only several months (curation, new gene discovery,
    etc.)

### 3. Use a tools to detect ==statistically significant== changes of read counts between conditions

Most used statistics tools (R packages):

- [x] DESeq2
- [x] EdgeR
