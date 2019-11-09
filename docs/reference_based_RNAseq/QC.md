![](images/galaxylogo.png)

# Quality Control

## FastQC tool to analyse the fastq (or fastq.gz) datasets

----
![](images/tool_small.png)

  1. Create a new history and name it `Quality Control`
  
  2. Copy again all fastq.gz files from the data library into this history. You should
  have 11 datasets in your history
  
  3. Select the `fastqc` tool.
  
  4. In the `Short read data from your current history` menu, select the `multiple datasets` button. ![](images/multiple-datasets.png)
  
  5. Shift-Click to select all 11 datasets
  
  6. Click `Execute`
  ----
   
  ![](images/oeil.png)
  
  - Look at the results of `FastQC`: These are the datasets named `FastQC on data xx: Webpage`
  ----
  
## MultiQC to aggregate and have a general view of sequence qualities in the project

----
  ![](images/tool_small.png)
  
  1. Select the `MultiQC`tool (you can use the search bar).
  
  2. `Which tool was used generate logs?` : Select `FastQC`
  
  3. `Type of FastQC output?` : Select `Raw data`
  
  4. `FastQC output` Cmd-Click (discontinuous, multiple selection) the *11* files named
  `FastQC on xx: RawData`
  
  5. Click `Execute`
  ----
  
  ![](images/oeil.png)
  
  Look at the result of `MultiQC`, dataset named `MultiQC on ...: Webpage`
  
  - Pay attention to the General Statistics that indicate the read sizes.
  - Pay attention to the `Sequence Quality Histograms`. What can you say about the
  quality of the samples ?
  - Have a look to the `Adapter Content` section.
  ----
