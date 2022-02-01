## [Bowtie](http://bowtie-bio.sourceforge.net/) align reads on _indexed_ genomes

![](images/bowtie.png){width="600"}

## Drosophila genome index

This task has already been performed for you. Here is the code used by the "root" user.

:warning: **Please do not execute the following code, we do not have enough time**

```
root@instance-1:~# mkdir dmel && chmod 777 dmel
root@instance-1:~# cd /home/dmel/
root@instance-1:~# wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.49_FB2013_01/fasta/dmel-all-chromosome-r5.49.fasta.gz
root@instance-1:~# chmod 444 dmel-all-chromosome-r5.49.fasta.gz
root@instance-1:~# gunzip dmel-all-chromosome-r5.49.fasta.gz

root@instance-1:~# apt install bowtie
root@instance-1:~# apt install samtools

root@instance-1:~# bowtie-build dmel-all-chromosome-r5.49.fasta Dmel_r5.49
# This last step took about 3:21 min and returned the following informations:
```

??? info "output of the `bowtie-build dmel-all-chromosome-r5.49.fasta Dmel_r5.49` command"
    ```
    Returning from Ebwt constructor
    Headers:
    len: 162367812
    bwtLen: 162367813
    sz: 40591953
    bwtSz: 40591954
    lineRate: 6
    linesPerSide: 1
    offRate: 5
    offMask: 0xffffffe0
    isaRate: -1
    isaMask: 0xffffffff
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 5073995
    offsSz: 20295980
    isaLen: 0
    isaSz: 0
    lineSz: 64
    sideSz: 64
    sideBwtSz: 56
    sideBwtLen: 224
    numSidePairs: 362429
    numSides: 724858
    numLines: 724858
    ebwtTotLen: 46390912
    ebwtTotSz: 46390912
    reverse: 0
    Total time for backward call to driver() for mirror index: 00:03:21
    ```

:warning:
Thus, for every user of the server, there is a Dmel_r5.49 genome bowtie index available at
`/home/dmel/Dmel_r5.49`

## Bowtie alignement of the clipped reads

```
bowtie /home/dmel/Dmel_r5.49 -f clipped_GKG13.fasta \
-v 1 -k 1 -p 4 \
--al droso_matched_GKG-13.fa \
--un unmatched_GKG13.fa \
-S > GKG13_bowtie_output.sam
```
??? info "The bowtie alignment command explained"
    - `bowtie dmel/Dmel_r5.49 -f clipped_GKG13.fasta` # tells bowtie where is the index and the input clipped_GKG13.fasta
    - `-v 1 -k 1 -p 4` # These are bowtie options
    - `--al droso_matched_GKG-13.fa` # aligned reads will be in the droso_matched_GKG-13.fa file
    - `--un unmatched_GKG13.fa` # Unaligned reads will be in the unmatched_GKG13.fa file
    - `-S > GKG13_bowtie_output.sam` # tells bowtie to return an alignement file in SAM format (-S) GKG13_bowtie_output.sam

### Bowtie Outputs

```
ls -laht
```

returns a number of files, among which the bowtie alignment outputs:

- GKG13_bowtie_output.sam
- droso_matched_GKG-13.fa
- unmatched_GKG13.fa

Please, look at the content of these files:

```
less GKG13_bowtie_output.sam
```
```
less droso_matched_GKG-13.fa
```
```
less unmatched_GKG13.fa
```

The last 2 steps we need to perform are

1. Compressing the SAM alignment file to a BAM alignment file
```
samtools view -Sb -@ 4 GKG13_bowtie_output.sam > GKG13_bowtie_output.bam # bam compression
```
2. Sorting the alignement **by position on the chromosomes**
```
samtools sort -@ 4 GKG13_bowtie_output.bam -o GKG13_bowtie_output.sort.bam # bam sorting
```

