# Data

The original data is available at NCBI Gene Expression Omnibus (GEO)
under accession number [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508).
It is also mirrored at the EBI Small Read Archive under the accession number
[SRP001537](https://www.ebi.ac.uk/ena/data/view/SRP001537)

The data was generated through deep Sequencing of mRNA from the Drosophila melanogaster
S2-DRSC cells that have been RNAi depleted of mRNAs encoding RNA binding proteins.

In the tutorial, we are going to focus on 7 datasets generated to study the effect of the
*Pasilla* gene inactivation by RNAi knock-down.

- 4 untreated samples: GSM461176, GSM461177, GSM461178, GSM461182
- 3 treated samples (Pasilla gene depleted by RNAi): GSM461179, GSM461180, GSM461181 

Each sample constitutes a separate biological replicate of the corresponding condition
(treated or untreated).

Two of the treated and two of the untreated samples are from a paired-end sequencing assay,
while the remaining samples are from a single-end sequencing experiment. Thus the following
table will be (very) useful in our analysis since each of the 7 datasets are designated
with (i) its original ID in GEO (or EBI SRA) (ii)  its condition (untreated or treated)
and (iii) the sequencing technology used (single read or paired-end).

| id. in GEO               | id. in EBI SRA           |
|--------------------------|--------------------------|
|GSM461176_untreat_single  |SRR031709_untreat_single  |
|GSM461177_untreat_paired  |SRR031714_untreat_paired  |
|GSM461178_untreat_paired  |SRR031716_untreat_paired  |
|GSM461179_treat_single    |SRR031718_treat_single    |
|GSM461180_treat_paired    |SRR031724_treat_paired    |
|GSM461181_treat_paired    |SRR031726_treat_paired    |
|GSM461182_untreat_single  |SRR031728_untreat_single  |

----
![](images/galaxylogo.png)

## Data upload

We will take benefit of this mandatory stage, to review various possibilities to upload
datasets in Galaxy. Specifically, we will review two options for uploading the gtf annotations
for the Drosophila genome dm6 in a Galaxy history. We will also have a look to a third option
that allows specifically to directly transfer FASTQ sequence files from the EBI SRA to a Galaxy history.

Transfers of Big Files take time, especially when the Internet connection speed is moderate to low...
To avoid consuming too much time on this task, you will have the possibility to import the full set
of the 11 FASTQ files in one of your histories, from a data library that has been pre-set in your Galaxy
server for this training session.

### Uploading data from your local computer

----
![](images/tool_small.png)

1. Download from the Ensembl database the sample [Drosophila_melanogaster.BDGP6.95.gtf.gz](ftp://ftp.ensembl.org/pub/release-95/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.95.gtf.gz) to your computer.
2. Upload this local file Drosophila_melanogaster.BDGP6.95.gtf.gz to your Galaxy history using the upload/Download
Galaxy interface that pops up if you click the upload icone  ![](images/upload_button.png)
----

### Importing data via links is more efficient and reliable !

----
[](images/tool_small.png)

The previous strategy is not efficient. Indeed, we can directly transfert the Drosophila_melanogaster.BDGP6.95.gtf.gz
from its primary location in the Ensembl database server to your Galaxy History !

----

1. Copy its URL below 
    
```
ftp://ftp.ensembl.org/pub/release-95/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.95.gtf.gz
```
    
2. and paste it in the `Paste/Fetch data` tab of the Galaxy upload interface.

3. In addition, select `gtf` in the `Type` menu.

4. Press the start button.

----

### Importing data via the `EBI SRA ENA SRA`

Finally there is a tool to specifically fetch fastq sequence file from the EBI SRA to Galaxy

The sample `GSM461178/SRR031716` was sequenced using a `paired-end` strategy (both ends of fragments
in the library are sequenced, giving rise to 2 read files, a forward read fastq file and a reverse
read fastq file).

We are going to download the fastq.gz files directly from
the EBI SRA using the tool `EBI SRA ENA SRA` in the `Get data` tool submenu.

----
 ![](images/tool_small.png)
----
  
  1. Click on the tool `EBI SRA ENA SRA` (you can select it rapidly using the search bar)
  
  2. In the search box of the EBI SRA website, enter `SRR031716`
  
  3. Two categories of results are retrieved, Experiment and Run.
    What we want to get are the files from the sequencing runs. Thus, click the
    SRR031716 link in the Run section (1 results found).
  
  4. Click on "File 1" in the `FASTQ files (Galaxy)` Column.
    You will be switched back to the Galaxy interface, and the download of the
    SRR031716_1.fastq.gz file will start immediately as a yellow dataset in the history right panel.
    
    Without waiting for the complete download of SRR031716_1.fastq.gz, you can repeat
    the previous steps 1, 2, 3 and 4. Just Click on `File 2` instead of `File 1` in step 4

Then, to save time, stop the tools (by clicking the small cross) and go to the next section.

----
    
### Importing data from data libraries

For collaborative work, Galaxy offers data libraries, where datasets can be stored and
available to one or multiple users.

This is what we are going to use to import rapidly all the input data you need for this RNAseq
analysis.

All datasets have been preloaded in the data library named `RNAseq`.

To access this library and import its content in your histories:

  ----
  ![](images/tool_small.png)
  
  1. Click the menu `Données partagées` (`Shared data`) and select the submenu
  `Bibliothèque de Données` (`Data libraries`).
  
  2. Navigate to the data library `RNAseq`
  
  3. Select all datasets
  
  4. Click the `To History` button and select `as Datasets`
  
  5. In the pop up window, `or create new` and type `Input data` to transfer the datasets
  in a new history with this name.
  
  6. Click on the green box to navigate to this new history (or click on the main menu `analyse data`)
  and start using these datasets.
----
