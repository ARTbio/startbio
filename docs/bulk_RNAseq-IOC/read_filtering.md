 ![](images/lamp.png)
# Focus on quality control & “filtering” in RNAseq analysis

## It is tempting to **filter** the data to get “good counts”

- low quality alignments
- PCR duplicates

## But..

- Why low quality reads should be skipped if they were aligned ? Is the implicit hypothesis
"low quality read are miss-mapped" a likely hypothesis ?

- When we remove PCR duplicates (exact same sequence and exact same location), are we sure
that we remove *PCR duplicates* ? What are the metrics that support the implicit hypothesis that
read with same sequence & same location are PCR duplicates ?

Reflect of miRNA sequencing...

