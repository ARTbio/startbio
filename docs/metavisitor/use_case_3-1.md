Now that you are familiar with manipulations in Galaxy with the Use Cases 1-1 to 1-4 described in detail in the previous chapters, we will describe the other Use Case analyses more concisely. If you experience lack of skills in basic Galaxy operations (tool usage, copy of datasets, etc), do not hesitate to go back and examine the [previous chapters](use_cases_input_data) step by step.

---
#### Virus detection in human RNAseq libraries

In the Use Cases 3-X we'll use Metavisitor to detect viruses in RNA sequencing dataset of human patients from 3 different studies.

In Use Case 3-1 we use Metavisitor to detect and assemble HIV genomes from patients Innate lymphoid cells sequencing data EBI SRP068722.

---

## Input data for Use Case 3-1

As for the previous Use Cases 1 and 2, the first step is to collect all the input data in a history that we will name `Input data for Use Case 3-1`.

1. Create a new history
    - Rename this history `Input data for Use Case 3-1`
    - We are going to upload 40 datasets form the EBI ENA SRP068722 :

        Go to the upload files menu and select `Paste/Fetch data`. Copy-Paste the following table (excluding the headers):

        SRR id | Patient id|
        -----------|-----------------|
        SRR3111582 | patient 0450-318|
        SRR3111583 | patient 0450-318|
        SRR3111584 | patient 0450-318|
        SRR3111585 | patient 0450-318|
        SRR3111586 | patient 0450-318|
        SRR3111587 | patient 0450-318|
        SRR3111588 | patient 0387-272|
        SRR3111589 | patient 0387-272|
        SRR3111590 | patient 0387-272|
        SRR3111591 | patient 0387-272|
        SRR3111592 | patient 0387-272|
        SRR3111593 | patient 0387-272|
        SRR3111594 | patient 0629-453|
        SRR3111595 | patient 0629-453|
        SRR3111596 | patient 0629-453|
        SRR3111597 | patient 0629-453|
        SRR3111598 | patient 0629-453|
        SRR3111599 | patient 0629-453|
        SRR3111600 | patient 0444-312|
        SRR3111601 | patient 0444-312|
        SRR3111602 | patient 0444-312|
        SRR3111603 | patient 0444-312|
        SRR3111604 | patient 0500-355neg|
        SRR3111605 | patient 0500-355neg|
        SRR3111606 | patient 0292-xxxneg|
        SRR3111607 | patient 0292-xxxneg|
        SRR3111608 | patient 0394-274|
        SRR3111609 | patient 0394-274|
        SRR3111610 | patient 0218-162neg|
        SRR3111611 | patient 0218-162neg|
        SRR3111612 | patient 0311-217HIVneg|
        SRR3111613 | patient 0311-217HIVneg|
        SRR3111614 | patient 0440-307neg|
        SRR3111616 | patient 0440-307neg|
        SRR3111617 | patient 0518-370neg|
        SRR3111618 | patient 0518-370neg|
        SRR3111619 | patient 0560-420neg|
        SRR3111620 | patient 0560-420neg|
        SRR3111621 | patient 0575-419neg|
        SRR3111622 | patient 0575-419neg|

    - Click the `Start` button Name and rename the dataset "Use-Case_3-1_information".
    - Use the tool `Cut columns from table`. In the "Cut columns field" write `c1` and make sure you select "Use-Case_3-1_information" file in the "From" field before executing. Rename the output "Use-Case_3-1_accessions".
    - Use the tool `Download and Extract Reads in FASTA/Q format from NCBI SRA`, select `List of SRA accession, one per line`from `select input type` and "Use-Case_3-1_accessions" in sra accession list. Click the `Execute` button.
    - When the tool is finished running you should have 2 new dataset collections in your history, one of them is empty. Delete the empty collection and verify that you have 40 pairs of datasets in the second collection.
    - If you are missing some sequences you'll have to re-do the steps above with only the missing identifiers. Once done, merge the collections using the tool `Merge Collections`.
    - Use the `Concatenate multiple datasets tail-to-head` tool and select "Paired collection" as type of data. Set the paired collection as input and select "Concatenate pairs of datasets" as type of concatenation. Execute the tool.
    - Rename the outputed collection to `SRP068722` and delete the previous one by clicking the `X` button and selecting "Permanently Delete Datasets".

2. Copy the `vir2 nucleotide BLAST database` from the `References` history to the current history `Input data for Use Case 3-1`.
3. Now we still have to associate sequencing dataset coming from the same patient. We are going to use the tool `Tag elements from file` to add the patient information as metadata.
    - Click on the `Tag elements from file` tool and select the collection "SRP068722" in "Input Collection" and "Use-Case_3-1_information" in "Tag collection elements according to this file". Execute the tool. Rename the new dataset collection `SRP068722_with_patient_information`.
    - Select the `Apply Rule to Collection` and set "SRP068722_with_patient_information" as "Input Collection". Click on the "Edit" button at the right of the form.
        - Click the "Column" button and select `Add Column from Metadata` from the list.
        - In the "From" list select "Tags". Then click the "Apply" button.
        - Click the "Rules" button and select `Add / Modify Colmn Definitions` from the list.
        - Click the "Add Definition" button and select the `List identifier(s)` from the list.
        - In the "Select a column" list select "B" then click on `... Assign Another Column` and select "A". Click the "Apply" button.
        - Click the "Save" and execute the tool.
        - Select the `Concatenate multiple datasets tail-to head` tool. In "What type of data do you wish to concatenate?" select "Nested collection". In "Input nested collection" select "SRP068722_with_patient_information (re-organized)". Execute the tool.
        - Rename the resulting collection "patient collection".
4. We are done. You can now permanently delete "SRP068722_with_patient_information",  "SRP068722_with_patient_information (re-organized)" and "SRP068722". This will save you some disk space.

## History for Use Case 3-1
1. Stay in the history `Input data for Use Case 3-1`
- pick the workflow `Metavisitor: Workflow for Use Case 3-1` in the workflows menu, and select the `run` option.
- For Step 1 (Fever Patient Sequences collection), select `patient collection` (this should be already selected).
- For Step 2, select the `nucleotide vir2 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `History for Use Case 3-1`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started may take time to show up.

---

## Results

The results for this use case differ whether you use [Metavisitor](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0168397) or Metavisitor2.
There are failed (red), empty datasets in this history. These datasets correspond to patients who didn't have any sequence matchng the viral database. You will notice that only `patient 0629-453` has contigs matching HIV sequences. However, this patient is a false positive. In support to this conclusion, you can:

- Copy the name of the contig
- Select `Pick Fasta sequences`
- Paste the contig name in the "Select sequences with this string in their header" section
- Select dataset `85: Oases viral contigs` and run the tool
- In a new browser tab go to the [Blastn web page](https://blast.ncbi.nlm.nih.gov/Blast.cgi?LINK_LOC=blasthome&PAGE_TYPE=BlastSearch&PROGRAM=blastn)
- Copy paste the contig sequence in the query section and mae sure the `Nucleotide collection nr` is selected as database before running blast

The sequence does not match viruses but cloning vectors.
