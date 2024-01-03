## Reference Genomes

- Fasta format

- Assembly version, generally, associated to a number and a date of assembly

- A same assembly may be provided by various organisation (Genome Resource Consortium, Ensembl, NCBI, UCSC, etc)
    
    This will be the same DNA sequence but formats may differ:
    
    - by the name of the chromosomes (chr1, 1, NC_000001.11, ...)
    - by the presence (or the absence) of unmapped contigs and haplotypes

## exemple 1: human genome

- GRCh37/hg19 - juil 2007
- GRCh38/hg38 - déc 2011
- GRCh39/hg39 - juin 2020 (repeat ++)

This various versions (or "releases") may in addition contain

- chromosomal regions "Aplotypes" (HLA, HBV inserts, etc…)
- unmapped contigs (regions which are significant assembly of reads, but are not assigned to a specific chromosome)

## exemple 2: mouse genome

<table>
  <tr>
   <td><strong>Release name</strong>
   </td>
   <td><strong>Date of release</strong>
   </td>
   <td><strong>Equivalent UCSC version</strong>
   </td>
  </tr>
  <tr>
   <td>GRCm39
   </td>
   <td>June 2020
   </td>
   <td>mm39
   </td>
  </tr>
  <tr>
   <td>GRCm38
   </td>
   <td>Dec 2011
   </td>
   <td>mm10
   </td>
  </tr>
  <tr>
   <td>NCBI Build 37
   </td>
   <td>Jul 2007
   </td>
   <td>mm9
   </td>
  </tr>
  <tr>
   <td>NCBI Build 36
   </td>
   <td>Feb 2006
   </td>
   <td>mm8
   </td>
  </tr>
  <tr>
   <td>NCBI Build 35
   </td>
   <td>Aug 2005
   </td>
   <td>mm7
   </td>
  </tr>
  <tr>
   <td>NCBI Build 34
   </td>
   <td>Mar 2005
   </td>
   <td>mm6
   </td>
  </tr>
</table>

## Annotations

It is important to note that annotations of genomes (GTF, GFF, etc.) although generally
equivalent, are strictly linked to their genome version because they refere to the DNA
sequences using the format of the release. This is why a GTF annotation file downloaded
from Ensembl is not interchangeable with a GTF annotation file from the UCSC or from another
organisation.

Moreover, since genome annotations may be considered as genome metadata (data on data), it is
normal and expected that genome annotation versions are different from genome versions and
that they are released at a faster pace.
