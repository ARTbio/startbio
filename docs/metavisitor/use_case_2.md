# Histories for Use Cases 2-1, 2-2

Now that are more familiar with manipulations in Galaxy with the Use Cases 1-1 to 1-4, we will
describe the other Use Case analyses more concisely. If you experience lack of skills in basic
Galaxy operations (tool usage, copy of datasets, etc), do not hesitate to go back and
examine the previous chapters step by step.

---
#### Discovery of novel viruses

Here we are going to use Metavisitor workflows to discover new viruses infecting a
laboratory colony of *Anopheles coluzzii* mosquitoes.

- Workflow for Use Case 2-1:

    Takes reads from EBI SRA ERP012577 (small RNAs), assembles contigs, blastx them against
    **vir2**, selects contigs hitting *Dicistroviridae* proteins, re-assembles selected
    contigs to align them to the *Drosophila* C virus (DCV) genome and integrates them to its sequence

- Worflow for Use Case 2-2:

    Takes reads from EBI SRA ERS977505 (mRNA), assembles contigs, and blastx them against
    **vir2**
---

## Input data for Use Cases 2-1 and 2-2

As for the previous Use Case 1, the first step is to collect all input data in an history
that we will name `Input data for Use Cases 2-1 and 2-2`

 - Create a new history
 - Rename this history `Input data for Use Cases 2-1 and 2-2`
 - For the small RNA sequence datasets (ERP012577) in this study, we are going to use
 another tool to upload to the Galaxy Metavisitor server: the `EBI SRA ENA SRA`tool which
 in the "Get data" section of the left tool bar.
    - Click on this tool and enter ERP012577 in the search field that shows up in the
    European Nucleotide Archive web page, and search. Click on the `ERP012577` link.
    In the column "Submitted files (galaxy)" of the table, click on the first "fastq file 1".
    This action should send you back to your Galaxy page automatically and you see the
    fastq dataset loading (yellow dataset in the history bar).
    - Repeat the exact same operation, for the three other "fastq file 1".
    - At final you should have uploaded four fastq datasets corresponding to the sequencing
    runs "post_infected_rep1.fastq", "post_infected_rep2.fastq", "post_non-infected_rep1.fastq"
    and "post_non-infected_rep2.fastq"
    - Select a dataset to make sure its datatype is `fastqsanger.gz`. If not click on the
    crayon button of any one of these four datasets and select the `Datatypes` tab and set
    it to `fastqsanger.gz`.
    - Repeat this operation for the other 3 datasets.
 - Create a dataset collection as [previously explained](https://artbio.github.io/metavisitor/use_case_1/#history-for-remapping-in-use-cases-1-123)
 (step 9) and name it `Small RNA reads ERP012577`
 - For the RNA sequence datasets (ERS977505) that will be used in Use Case 2-2, use again
 the `EBI SRA ENA SRA`tool which in the "Get data" section of the left tool bar.
    - Click on this tool and enter ERS977505 in the search field that shows up in the
    European Nucleotide Archive web page, and search. Click on the `ERS977505` link
    (Sample 1 result found). In the column "Submitted files (galaxy)" of the table, click
    on the first "fastq file 1". This action should send you back to your Galaxy page
    automatically and you see the fastq dataset loading (yellow dataset in the history bar).
    - Repeat the exact same operation for the other "fastq file 1" and the two other
    "fastq file 2"
    - In the end you should have uploaded four additional fastq datasets corresponding to
    the sequencing runs "IP-isoT-1_AGTCAA_L001_R_1.fastq", "IP-isoT-1_AGTCAA_L001_R_2.fastq",
    "IP-isoT-2_ATGTCA_L002_R_1.fastq" and "IP-isoT-2_ATGTCA_L002_R_2.fastq"
 - Create a dataset collection as explained in the previous chapter and name it
 `long read RNAseq datasets`. Note that we are not handling the files as paired-ends, thus
 use the `Build dataset list` command and not the `Build list of dataset pairs` command.
 - Use the `Retrieve FASTA from NCBI`, paste `phix174[title]` in the "Query to NCBI in
 entrez format" field and select `Nucleotide` for the NCBI database. This will upload 154
 fasta sequences from phix174.
 - Use the wheel icon at the top of the history bar to copy `nucleotide vir2 blast database`,
 `protein vir2 blast database` and `P. berghei` **from** the history `References` **to**
 the current history `Input data for Use Cases 2-1 and 2-2`. If you don't remember well
 how to copy datasets between histories, you may read again the explanation
 [here](https://artbio.github.io/metavisitor/use_cases_input_data/#history-with-input-data-for-use-cases-1-1-1-2-1-3-and-1-4)
 (step 3)

**_Your are now ready for generating Uses Cases 2-1 and 2-2_**

## History for Use Case 2-1

- Stay in the current history `Input data for Use Cases 2-1 and 2-2` !
- In the `Workflow` menu, select the workflow `Metavisitor: Workflow for Use Case 2-1` and
directly select `Run` (you may also look at the workflow using the `edit` option).
- Be careful at selecting `Small RNA reads ERP012577` in step 1 (Input Dataset Collection).
- **Be careful** in selecting `P. berghei` in **step 2**.
- **Be careful** in **step 3** select :
```
Retrieve FASTA from NCBI (Nucleotide) with queryString 'phix174[title]'
```

- In step 4, the option `protein vir2 blast database` is forced, because the workflow is expecting of protein blast database in this step and only one dataset with this datatype is available in the history
- Click the `Send results to a new history` checkbox and rename the history to "History for Use Case 2-1".
- Run Workflow !

You may follow the link to the new history when the workflow has started.

## History for Use Case 2-2

- If you are not already in, go back to the history `Input data for Use Cases 2-1 and 2-2`
- In the `Workflow` menu, select the workflow `Metavisitor: Workflow for Use Case 2-2` and
directly select `Run` (you may also look at the workflow using the `edit` option)
- Be careful at selecting `long read RNAseq datasets` in step 1 (Input Dataset Collection)
- In step 2, the option `protein vir2 blast database` is forced, because the workflow is
expecting of protein blast database in this step and only one dataset with this datatype
is available in the history
- Click the `Send results to a new history` checkbox and rename the history to "History
for Use Case 2-2".
- Run Workflow.

## Re-mapping of the small RNA reads (ERP012577) to the AnCV genome (KU169878).
The previous Workflow for Use Case 2-2 allowed to assemble a large contig of 8919 nt which
significantly matched structural and non-structural polyproteins of Drosophila C Virus and
Cricket Paralysis Virus in blastx alignments (see the dataset
`blastx Filter sequences by length on data 17 vs 'protein BLAST database from data 2'`
of the history). This large contig corresponds to the genome of a new Anopheles C Virus
deposited to the NCBI nucleotide database under accession number KU169878 (see the
[companion Metavisitor article](http://dx.doi.org/10.1101/048983) and
[Carissimo et al](http://dx.doi.org/10.1371/journal.pone.0153881)).

Here, we are going to perform manually a few steps, before using another workflow in the
history 2-2 to remap the ERP012577 small RNA reads to the AnCV genome.

- Look at the `blast analysis, by subjects` dataset and copy the name of the 8919 nt
contig that aligned to DCV and CrPV sequences. It is noteworthy that the names may vary
from one Oases run to another because the Oases algorithm is not totally deterministic.
In the [companion Metavisitor article](http://dx.doi.org/10.1101/048983), this name was
Locus_69_Transcript_1/1_Confidence_0.000_Length_8919.
- Copy this name, find the tool `Pick Fasta sequences with header satisfying a
query string` in the Galaxy tool bar, and paste the name in the field `Select sequences
with this string in their header` of the tool form. Select the dataset `Oases_optimiser on
data 21: Denovo assembled transcripts` as a source file, and run the tool.
- Now, we are going to change the header of the previously extracted fasta sequences using
the tool `Regex Find And Replace`.
    - Select the previous dataset `Concatenated datasets` as input dataset for this tool.
    Click on `+ Insert Check`. Use `>.*Confidence(.*)_Length_8919` (or the equivalent you
    extracted) as *Find Regex* and `>Anopheles_C_Virus|KU169878_confidence\1` as *Replacement*.

- Copy the dataset collection `Small RNA reads ERP012577` from the history `Input data for
Use Cases 2-1 and 2-2` into the *current* history `Use Case 2-2`. You may have to refresh
the history bar to see this collection and the attached datasets popping up.

We are now ready to run the workflow.

----

- In the workflow menu, pick up the workflow
`Metavisitor: Workflow for remapping in Use Cases 2-1,2` and select the `run` option.
- In the workflow form, ensure that `Small RNA reads ERP012577` are selected in step 1 and
`Regex Find And Replace on data 28` is selected in step 2 (this should be the case if you
followed the instructions).
- This time, *do not* check the box `Send results to a new history` and directly click the
`Run workflow`button.

This workflow will provide you with a graphical view of ERP012577 small RNA mapping to the
AnCV genome.
