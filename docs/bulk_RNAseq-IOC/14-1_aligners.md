## Aligners softwares

The main alignment softwares are currently:

- BWA
- Bowtie
- STAR

They are all based on the Burrows-Wheeler Algorithm.
This implies to build a genome index in which the genome is recoded using the BWA, ensuring
very fast read alignments.

BWA-based aligner are CPU- and IO-demanding. In contrast they usually are not demanding in
RAM (with maybe the exception of STAR, for index building)

Aligners take FASTQ (FASTQ.gz) filesas well as a genome reference index
appropriately built as inputs.

They return BAM files which are compressed SAM files (Simple Alignment/Map).

The SAM format is really at the heart of RNAseq analyses, because it contains ==all== the
information needed to profile gene expressions from sequencing datasets.

==**Therefore, we highly recommend** to take a few hours to look at all the details of the SAM
format==, which can be found in the [GitHub repository](https://github.com/samtools/hts-specs).
You can start with [Sequence Alignment/Map format specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf),
and also have a closer look at
[Sequence Alignment/Map optional fields specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)

## Pseudo-aligners

Other aligners rather operate using a pseudo-alignment mode based on graphs of k-mers.

These include [Kallisto](https://cyverse-leptin-rna-seq-lesson-dev.readthedocs-hosted.com/en/latest/section-8.html)
and [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)