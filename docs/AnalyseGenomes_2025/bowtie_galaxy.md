---
hide:
  - toc
---
## Import data

- Rename the `Unnamed history` to `Bowtie` using the pencil icon
- Go to `Upload Data` (to the left bar) and select `Paste/Fetch Data`
- Paste the following content
```
https://ftp.flybase.net/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz
https://psilo.sorbonne-universite.fr/index.php/s/HYLtfo9d2eD3Q2A/download/GRH-103_R1.fastq.gz
```
- And click the `start` button

- Check the imported datasets in the history bar
- Check the content of the imported datasets by clicking the eye icon in each dataset

## Install required packages
==Required packages (`bowtie` and `samtools`) are already installed in your Galaxy server==

## Clip fastq reads from their sequence adapter and output clipped sequences in a fasta format
- type "clip adapter" in the search toolbar box
- select the `Clip adapter` Galaxy toolbar
- Fill the tool form as following, indicating which file to clip, the min and max sizes of the
  reads you wish to keep in the processed dataset, that you want a fasta output, do no want
  N in the retrieved clipped reads, and that the adapter in the dataset is the Illumina
  TruSeq adapter.
!!! info ":wrench: Clip adapter parameters"
    - **Source file**: `2: GRH-103_R1.fastq.gz`
    - **min size**: `18`
    - **max size**: `36`
    - **Select output format**: `fasta`
    - **Accept reads containing N?**: `reject`
    - **Source**: `Use a built-in adapter (select from the list below)`
    - **Select Adapter to clip**: `Illumina TruSeq TGGAATTCTCGGGTGCCAAGTGGAAT`

![clip tool](images/clip.png){width="500"}

- Click the `Execute` icon

Check the result in the history:

- how many clipped sequences ? --> click on the dataset to deploy it
- which format ?
- How do the sequences look like ? --> click on the eye icon

## Prepare dmel_r6.54 bowtie index

==No need to prepare the bowtie index, the next tool will do it for us on the fly==

## Align the clipped fasta reads to dmel.r6.54 using `bowtie`

- In the search toolbar box, type `bowtie`
- Select the tool `sR_bowtie for small RNA short reads`
!!! info ":wrench: sR_bowtie for small RNA short reads parameters"
    - **Input fasta or fastq file: reads clipped from their adapter**: `Clipped GRH-103_R1.fastq.gz-then-fasta` 
    -  **What kind of matching do you want to do?**: `Match on DNA as fast as possible, ...`
    - **Number of mismatches allowed**: `0`
    - **Will you select a reference genome from your history or use a built-in index?**: `Use one from the history`
    - **Select a fasta file, to serve as index reference**: `dmel-all-chromosome-r6.54.fasta`
    -  **Select output format**: `bam`
    - **additional fasta output**: `both aligned and unaligned`

Examine the output datasets (`Bowtie Output`, `Matched reads` and `Unmatched reads`)


## Convert SAM file to BAM file and sort the alignments by chromosome positions

==This is automatically done by Galaxy==