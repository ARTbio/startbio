 ![](images/lamp.png)

It is tempting to **filter** the data to get “good reads” and **discard** "bad reads"
(false reads ??), including: 

- Reads with low quality alignments
- Reads suspected to be PCR duplicates

**HOWEVER**

  ==Discarding reads is a radical decision because it changes the counts==

- [ ] Why low quality reads should be skipped if they were aligned ?
- [ ] Is the implicit hypothesis "low quality read are miss-mapped" a likely hypothesis ?


- [ ] When we remove PCR duplicates (exact same sequence and exact same location), are we sure
      that we remove *PCR duplicates* ? What are the metrics that support the implicit hypothesis that
      read with same sequence & same location are PCR duplicates ?


- [ ] You may also reflect on the case of miRNA sequencing...

??? tip "Spoiler alert"
    **For the questions mentionned above, we are pretty reluctant to discard reads for
    quality reasons, when it comes to bulk RNAseq analysis.**

