---
hide:
  - toc
---
## Import data

We first create a working directory for our bowtie alignment and import the required input
data in it:
```
mkdir ~/bowtie_work && cd ~/bowtie_work
```
```
wget https://ftp.flybase.net/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz \
     https://psilo.sorbonne-universite.fr/index.php/s/HYLtfo9d2eD3Q2A/download/GRH-103_R1.fastq.gz
```
Check the imported files using:
```
ll
```

We also need to uncompress the `.gz` files
```
gunzip *.gz
```
you can check the result by
```
ll -rt
```

## Install required packages
We will need the `bowtie` and `samtools` programs:
```
apt update && apt install -y bowtie samtools
```

## Clip fastq reads from their sequence adapter and output clipped sequences in a fasta format
```
cat GRH-103_R1.fastq | \
perl -ne 'if (/^([GATC]{18,})TGGAATT/){$count++; print ">$count\n"; print "$1\n"}' \
> clipped_GRH-103.fa
```
Check the result with
```
grep -c ">" clipped_GRH-103.fa
```
and
```
wc -l clipped_GRH-103.fa
```

## Prepare dmel_r6.54 bowtie index
The following command line is masked. Before unmasking it, you can try to find the
appropriate command line using the `man` command or the `--help` argument
??? question "Bowtie indexing command line"
    ```
    time bowtie-build --threads 7 dmel-all-chromosome-r6.54.fasta dmel.r6.54
    ```
    Note the `time` here is to indicate the time consumed to index the genome, it is optional.

## Align the clipped fasta reads to dmel.r6.54 using `bowtie`
```
time bowtie dmel.r6.54 -f clipped_GRH-103.fa \
                       -v 0 \
                       -k 1 \
                       -p 7 \
                       --al dmel_matched_GRH-103.fa \
                       --un unmatched_GRH-103.fa \
                       -S \
                       > GRH-103.sam
```
## Convert SAM file to BAM file and sort the alignments by chromosome positions
```
samtools view -Sb -@ 7 GRH-103.sam | samtools sort -@ 4 -o GRH-103.bam
```
Check the result using
```
samtools view GRH-103.bam | more
```
