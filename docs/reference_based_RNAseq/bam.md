# Inspection of BAM files

----
![](images/oeil.png)

Click on the small eye icon of a Bam dataset (generated either with `RNA STAR` or `HISAT2`)
![](images/view-dataset.png)

The header contains the chromosome specifications (their name and length) and other informations
such as the software that generation the Bam file and the command line to run the software.

A BAM file (or a SAM file, the non compressed version) consists of:

A header section (the lines starting with @) containing metadata, in particular the chromosome names and lengths (lines starting with the @SQ symbol)
An alignment section consisting of a table with 11 mandatory fields, as well as a variable number of optional fields:

|Col |	Field |Type	  |      Brief Description               |
|----|--------|-------|--------------------------------------|
|1	 | QNAME  |String | Query template NAME                  |
|2	 | FLAG	  |Integer| bitwise FLAG                         |
|3	 | RNAME  |String |References sequence NAME              |
|4	 | POS    |Integer|1-based leftmost mapping POSition     |
|5	 | MAPQ   |Integer|MAPping Quality                       |
|6	 | CIGAR  |String |CIGAR String                          |
|7	 | RNEXT  |String |Ref. name of the mate/next read       |
|8	 | PNEXT  |Integer|Position of the mate/next read        |
|9	 | TLEN   |Integer|observed Template LENgth              |
|10	 | SEQ    |String |segment SEQuence                      | 
|11	 | QUAL   |String |ASCII of Phred-scaled base QUALity+33 |

----
