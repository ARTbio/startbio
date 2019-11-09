## ![](images/lamp.png) Reference-base Expression analysis: the key idea

![](images/readcounts.png)

### Map reads to a reference genome with aligners

- TopHat
- TopHat2
- HiSat
- HiSat2
- STAR

→ These aligners are “splice aware”

→ They generate a **BAM Alignment file**

### Use a read counting software and annotation information (GTF, GFF3, BED, …) to count the read spanning a gene / transcript

The input file for this counting software is the **BAM Alignment file**

# Read counts are *proxies* to RNA steady state levels
