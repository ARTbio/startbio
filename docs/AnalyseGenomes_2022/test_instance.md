Here is the procedure for rapidly testing whether your instance is correctly deployed.
The main issues that you can encounter with a Galaxy server are improperly installed tool
dependencies. In galaxy, most of this dependencies (codes, packages, modules that tools
need to access for run) are installed using the conda packages & environments manager.

In order to test a significant sample of these tools and dependencies, the strategy here is
to import a test history in Galaxy and a workflow which will take inputs from this history.
The workflow is then run and must produce only "green" datasets in order to get the test
validating.

In addition, if the test does not pass, information in the "red" dataset will be very useful
to fix the issues experienced in your instance.

### 1. Import the test history

In the history panel (Menu `User` --> `histories`), click the ++"Import from file"++ button
at the top-right corner of the panel.

Paste the url of the test history archive (.tar.gz) in the field (checkbox `Export URL from
another Galaxy instance` checked)
```
https://storage.googleapis.com/analyse-genome-coupon-1/%20Test-History-sampleRNAseq.tar.gz
```

Wait for a dozen of secondes

### 2. Go to the `your histories` link.

There, you will see a new history named `imported from archive: sample_RNAseq`.

### 3. Go to that history `imported from archive: sample_RNAseq`

### 4. from the menu `Copy Datasets`, copy the `dmel-all-r6.18.gtf` from your history `References`
or whatever you named it).

### 5. Go to the menu ++"Workflows"++ of Galaxy. Here, you will notice a workflow `Analyse RNAseq`
that the Ansible playbook `GalaxyKickStart` has preloaded for you !

At the right side of this workflow name, there is an arrow to trigger the workflow execution.
Trigger it !

### 6. In the workflow form,

!!! note "Fill the form of :wrench: **Workflow: Analyse RNAseq**"
    - **Send results to a new history**: `Yes`
    - **History name**: Analyse RNAseq Test 1
    - **WT Collection**: `8: sample WT`
    - **SF Collection**: `7: sample SF1`
    - **dmel GTF**: `9: dmel-all-r6.18.gtf` (It must be if you correctly copied the dataset
      from the reference history)
    - **Click the ++"Run Workflow"++ button**
    
### 7. The workflow you take few minutes to run. You can follow the operation in the new history
which was created.

### 8. Test Results

If all dataset are green at the end of the workflow run: The test is ok.
If the workflow stops with some red datasets, look carefully at the error and bug icons of
these datasets: they contain useful information which you can escalate to your trainers for
help.

!!! info "Correctly reporting an error during the training"
    Remember the best, cleanest way to report an error is to [raise an issue](https://github.com/ARTbio/Run-Galaxy/issues/new/choose)
    in the GitHub repository ARTbio/Run-Galaxy, with **a maximum of detail** such as
    
    - Detailed description of the issue
    - logs and alert messages which you can copy between two lines that will contain only three
      back ticks: 
      
      
      **\`\`\`**
      
      Put your code/log here
      
      **\`\`\`**
    
    - Screen shots !


