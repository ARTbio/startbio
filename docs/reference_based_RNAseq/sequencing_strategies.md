# ![](images/lamp.png) Inserts and sequencing strategies

![](images/sequencing_strategies.png)

You can retrieve three different informations :

1. The relative orientation of reads :
    - `I` : Inwards
    - `M` : Matching
    - `O` : Outwards
2. The strandedness of the library :
    - `S` : Stranded
    - `U` : Unstranded
3. The strand origin of reads :
    - `F` : read 1 (or single-end read) comes from the forward strand
    - `R` : read 1 (or single-end read) comes from the reverse strand


## in practice, with Illumina paired-end RNAseq protocols you will either deal with:

### Unstranded RNAseq data
**IU type from above. Also called fr-unstranded in TopHat/Cufflinks nomenclature**

### Stranded RNAseq data produced with Illumina TrueSeq RNAseq kits
**ISR type from above or fr-firststrand in TopHat/Cufflinks nomenclature**

