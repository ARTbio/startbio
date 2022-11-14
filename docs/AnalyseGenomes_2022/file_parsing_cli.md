---
hide:
  - toc
---

### Initial Format (EMBL flat file)
```
ID   INE1    standard; DNA; INV; 611 BP.
XX
AC   U66884;
XX
DR   FLYBASE; FBte0000312; Dmel\INE-1.
XX
SY   synonym: mini-me
SY   synonym: DINE
SY   synonym: narep1
SY   synonym: Dr. D
XX
FT   source          U66884:4880..15490
XX
CC   This is presumably a dead element.
CC   Derived from U66884 (e1371475) (Rel. 52, Last updated, Version 6).
CC   Michael Ashburner, 28-Sep-2001.
CC   Any changes to original sequence record are annotated in an FT line.
XX
SQ   Sequence 611 BP; 193 A; 123 C; 93 G; 202 T; 0 other;
     TATACCCGTT ACTAGATTCG TTGAAATGAA TGTAACAGGC AGAAGGAAGC GTCTTAGACC        60
     ATATATAGTA TATACATACA TGTATATTCT TGATCAGGAT CAATAGCCGA GTCGATCTTG       120
     CCATATCCGT CTGTCCGTAT GAACGTCGAG ATCTCAGGAA CTATAAAAGC TAGAAGGTTT       180
     AGATTCAGCA TACAGAGACA AAGACGCAAG TAGCCATGCC CACTCTAACG TCCACAAACA       240
     GCGCAAAACT ATCACGCCCA CACTTTTGAA AAATGTGTTG TTCTTTTCAC ATTCTGATTA       300
     GTCTTTTACA TTTCTATCGA TTTCCAAAAA AAAACTTTTT GCCAACGCCC TAAAACCGCC       360
     CAAAACTCCG ACACCCACAT TTGTAAAAAA TTGTTGGGAA TTTTTTTCAT AAATTTATTA       420
     GTTTATTATT TATTATAAAT TTAAGTTTAT ATCGATTTGC CGACAACATA TTTTAATTTT       480
     TTTTCTCATT TTATCTTTTA TCTATCGATA TCCCAGAAAA ATTGTGCAAT TTCGCATTCA       540
     CACTAGCTGA GTAACGGGTA TCTGATAGTC GGGAAACTCG ACTATAGCAT TCTCTCTTTT       600
     TGAAATTGCG G                                                            611
//
```
### Target Format
```
>INE1
TATACCCGTTACTAGATTCGTTGAAATGAATGTAACAGGCAGAAGGAAGCGTCTTAGACC
ATATATAGTATATACATACATGTATATTCTTGATCAGGATCAATAGCCGAGTCGATCTTG
CCATATCCGTCTGTCCGTATGAACGTCGAGATCTCAGGAACTATAAAAGCTAGAAGGTTT
AGATTCAGCATACAGAGACAAAGACGCAAGTAGCCATGCCCACTCTAACGTCCACAAACA
GCGCAAAACTATCACGCCCACACTTTTGAAAAATGTGTTGTTCTTTTCACATTCTGATTA
GTCTTTTACATTTCTATCGATTTCCAAAAAAAAACTTTTTGCCAACGCCCTAAAACCGCC
CAAAACTCCGACACCCACATTTGTAAAAAATTGTTGGGAATTTTTTTCATAAATTTATTA
GTTTATTATTTATTATAAATTTAAGTTTATATCGATTTGCCGACAACATATTTTAATTTT
TTTTCTCATTTTATCTTTTATCTATCGATATCCCAGAAAAATTGTGCAATTTCGCATTCA
CACTAGCTGAGTAACGGGTATCTGATAGTCGGGAAACTCGACTATAGCATTCTCTCTTTT
TGAAATTGCGG
```
## Import the dataset
Create a dedicated screen session:
```
screen -S parsing
```
Then:
```
mkdir ~/file_parsing && \
cd ~/file_parsing && \
wget https://raw.githubusercontent.com/bergmanlab/drosophila-transposons/9b28cdbe9d2b3ef895df37f8495b33104677e516/releases/transposon_sequence_set_v9.5.embl.txt
```
## Reformat the file using sequencial command lines:

```
grep -P "(^ID)|(^ +[GATCNgatcn ]+\d+)" transposon_sequence_set_v9.5.embl.txt > transposon_sequence_set_v9.5.fa
```

```
sed -i.bak -E "s/^ID   />/" transposon_sequence_set_v9.5.fa
```

```
sed -i.bak2 -E "s/(>[^ ]+) .+/\1/g" transposon_sequence_set_v9.5.fa
```

```
sed -i.bak3 -E "s/([GATCNgatcn]+) /\1/g" transposon_sequence_set_v9.5.fa
```

```
sed -i.bak4 -r "s/^ +//g" transposon_sequence_set_v9.5.fa
```

```
sed -i.bak5 -r "s/ +[0-9]+//g" transposon_sequence_set_v9.5.fa
```
## Check the conversion

- Download the file reference for the conversion (ie, a file that we know is correctly converted...)
```
wget https://raw.githubusercontent.com/ARTbio/AnalyseGenome/main/Exercises/transposon_sequence_set_v9.5.fa
```
- check the content, what do you see ?
```
ll -tr
```
- compute the difference between your conversion and the reference conversion
```
diff transposon_sequence_set_v9.5.fa transposon_sequence_set_v9.5.fa.1
```

