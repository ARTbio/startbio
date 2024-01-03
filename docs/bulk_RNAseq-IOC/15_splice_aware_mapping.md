## Splice-Aware Aligners

For RNAseq analysis, it is common to speak of "Splice-Aware" aligners.

This is in particular
mandatory if you work with an eukaryote organism where mature messenger RNAs are made of
joint exons coming from genome regions separated by introns. Indeed, in this situation,
mRNA derived sequencing reads maybe split between distant genomic regions and distance
between two paired reads may be much higher than expected.

Actually, splice-aware aligners are just BWA-base aligners wrapped in additional code to take
into accounts split or distant pair alignments.

Importantly, if you are working with a model organism with available genome annotations,
splice-aware aligners will heavily rely on these annotations. Therefore, splice-Aware
aligners will most of the time work with GTF (or GFF3) input files, in addition to the
fastq files and the genome reference index.

However if your working organism is not a model organism, splice-aware aligners are still
useful, since the will reconstruct de novo the exon-exon junctions identified in the
sequencing reads. Indeed they have often been used to discover new mRNA isoforms !

![](images/splice_aware_alignment.png)

## software

Historically, the first popular splice-aware aligner has been TopHat and TopHat2,
based on bowtie and bowtie2 aligners, respectively.

Nowadays, the two popular splice-aware aligners are

- HISAT2 (based on bowtie2)
- STAR (with its own aligner implementation).
  Note that in the case of STAR, you have the possibility to build index already incorporating
  GTF informations. It is also possible to provide GTF information at the runtime of the
  STAR alignment.
