#### Use Case 3-2's aim

In this Use Case, Metavisitor is used to search for the presence of viruses and identify them in RNA sequencing data of serums of children suffering from fevers of unknown origins. To compare Metavisitor's results to [Yozwiak *et al*](https://dx.doi.org/10.1371%2Fjournal.pntd.0001485).

---
## Input data for Use Case 3-2

As for the previous Use Cases 1, 2 and 3-1, the first step is to collect all input data in an history that we will name `Input data for Use Case 3-2`.

- Create a new history
- Rename this history `Input data for Use Case 3-2`
- Using the tool `Extract reads in FASTQ/A format from NCBI SRA`, we are going to upload 43 paired end datasets. Indeed, these 43 datasets correspond to 86 fastq paired-ended sequence files. In addition, some datasets derive from the same patient; in those cases we will merge those datasets using the tool `Concatenate multiple datasets tail-to-head` and delete and purge the original datasets.
- Open the "Upload datasets" menu. Click on the `Paste/Fetch data` button, name the file "Use-Case_3-2_SRR_information" and copy-paste the following text:

    SRR id | Patient id|
    -------|------------|
    SRR453487 | patient 566|
    SRR453437 | patient 438|
    SRR453443 | patient 401|
    SRR453458 | patient 401|
    SRR453430 | patient 382|
    SRR453491 | patient 377|
    SRR453499 | patient 375|
    SRR453484 | patient 350|
    SRR453464 | patient 349|
    SRR453506 | patient 345|
    SRR453417 | patient 344|
    SRR453490 | patient 335|
    SRR453478 | patient 331|
    SRR453465 | patient 330|
    SRR453480 | patient 330|
    SRR453489 | patient 329|
    SRR453505 | patient 329|
    SRR453498 | patient 322|
    SRR453446 | patient 321|
    SRR453427 | patient 315|
    SRR453440 | patient 315|
    SRR453438 | patient 282|
    SRR453450 | patient 275|
    SRR453460 | patient 274|
    SRR453485 | patient 270|
    SRR453448 | patient 266|
    SRR453424 | patient 263|
    SRR453457 | patient 263|
    SRR453510 | patient 193|
    SRR453456 | patient 187|
    SRR453425 | patient 186|
    SRR453469 | patient 186|
    SRR453481 | patient 183|
    SRR453531 | patient 180|
    SRR453474 | patient 179|
    SRR453509 | patient 171|
    SRR453451 | patient 168|
    SRR453495 | patient 161|
    SRR453504 | patient 161|
    SRR453500 | patient 159|
    SRR453493 | patient 156|
    SRR453444 | patient 131|
    SRR453426 | patient 78|

- Click the `Start` button.
- Select the `Cut columns from a table` tool and set the "Cut columns" parameter to "c1" and select "Use-Case_3-2_SRR_information" as input file in the "From" list. Execute the tool and rename the new dataset collection "Use_Case_3-2_accessions".
- Select the tool `Download and Extract Reads in FASTA/Q format from NCBI SRA`and select "List of SRA accession, one per line" in "select input type" and "Use_Case_3-2_accessions" in "sra accession list". Execute the tool. The data downloading step might take 40 minutes to 1h. Delete "Single-end data (fastq-dump)".
- Select the `Concatenate multiple datasets tail-to-head` tool and set "Paired collection" in "What type of data do you wish to concatenate?" and "Pair-end data (fastq-dump)" as "Input paired collection to concatenate". In "What type of concatenation do you wish to perform?" select "Concatenate pairs of datasets (outputs an unpaired collection of datasets)". Execute the tool.
- Select `Tag elements from file` tool and set "Concatenation by pairs" as "Input Collection" and "Use-Case_3-2_SRR_information" as "Tag collection elements according to this file". Execute the tool.
- Select `Apply Rule to Collection` tool and set "data 1, data 144, and others (Tagged)" as "Input Collection" and click the "Edit" button.
    - Click the "Column" button and select "Add Column from Metadata" from the list.
    - Select "Tags" from the "For" list and click the "Apply" button.
    - Click the "Rules" button and select "Add / Modify Column Definitions".
    - Click "Add Definitions" button and select "List identifier(s)" from the list.
    - Select "B" from the "Select a column" list.
    - Click "*... Assign Another Column*" and select "A" from the "Select column" list. Click the "Apply" button and the "Save" button. Execute the tool.
- Select the `Conatenate multiple datasets tail-to-head` tool ans set "Nested collection" in "What type of data do you wish to concatenate?" and select the "... (re-organized)" dataset collection in "Input nested collection". Click the "Execute" button.
- Rename the output collection as "Tractable Patient Datasets".
- Copy the `vir2 nucleotide BLAST database` from the `References` history to the current history `Input data for Use Case 3-2`.

## History for Use Case 3-2
1. Stay in the history `Input data for Use Case 3-2`
- pick the workflow `Metavisitor: Workflow for Use Case 3-2` in the workflows menu, and select the `run` option.
- For Step 1 (Fever Patient Sequences collection), select `Tractable Patient Datasets` (this should be already selected).
- For Step 2, select the `nucleotide vir2 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `History for Use Case 3-2`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up; you can even have a "504 Gateway Time-out" warning. This is not a serious issue: just go in your `User` -> `Saved history` menu, you will see you `History for Use Case 3-2` running and you will be able to access it.

As a last note, the workflow for Use Case 3-2 may take a long time. Be patient.
