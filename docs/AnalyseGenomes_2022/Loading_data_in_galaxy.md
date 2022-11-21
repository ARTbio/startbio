For the course "Analyse des Génomes", we need three types of datasets

- The reference sequences that will be used to align sequencing reads (full genome, miRNA, transposons, etc.)
- libraries of sequencing reads from small RNAs (for analysis of piRNAs)
- Librairies of sequencing reads from mRNA (for Gene differential expression analysis)

All these data have been deposited in the storage server
[Psilo](https://psilo.sorbonne-universite.fr/index.php/s/yHSoKGZKMeJkeXa) at
Sorbonne-Université.

### Get data "by URL"
As these data _**are available through a URL**_ (Universal Resource Location) we will use
as before the menu `Paste/Fetch Data` of the `Upload Data` menu.

??? warning "There are other methods to upload data in Galaxy !"
    * You can transfer data from your local machine (the one where your keyboard is plugged !)
      to Galaxy
    * You can upload data to your Galaxy FTP account and then transfer these data from your
    Galaxy FTP directory to one of your Galaxy histories.

    
#### 1. Upload of reference files as a batch of multiple URLs :heavy_plus_sign: Programmatic file naming

As you have already uploaded single files using their url, we are going to use a more
powerful procedure which is appropriate when uploading numerous files.

Before all, create a new history by clicking the :heavy_plus_sign: icon in the history header
![](images/history_header.png){ width="300"} and immediately renaming the new history as
`References`.

- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- This time, Click the `Rule-based`tab !
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:

!!! info ":candy: URLs of references (genome and RNA classes)"
    ```
    https://psilo.sorbonne-universite.fr/index.php/s/7YqGeFTxTgxtafy/download/dmel-all-chromosome-r6.18.fasta	dmel-r6.18
    https://psilo.sorbonne-universite.fr/index.php/s/w3soG64bRcZd5CJ/download/Dmel_all-transposon_merge.fasta	transposons
    https://psilo.sorbonne-universite.fr/index.php/s/2Y7BfGNZQokMDfT/download/Dmel_piRNA_clusters.fasta	piRNA_clusters
    https://psilo.sorbonne-universite.fr/index.php/s/s9reeBC79Pw6jfN/download/dmel-all-gene-r6.18.fasta	genes
    https://psilo.sorbonne-universite.fr/index.php/s/JndyRYqeWE9d8eb/download/dmel-all-intron-r6.18.fasta	introns
    https://psilo.sorbonne-universite.fr/index.php/s/aZAPqW9PJ2ncnXD/download/dmel-all-miRNA-r6.18.fasta	miRNAs
    https://psilo.sorbonne-universite.fr/index.php/s/dKKHKqHWLj45QWx/download/dmel-all-miscRNA-r6.18.fasta	miscRNAs
    https://psilo.sorbonne-universite.fr/index.php/s/cDy6iqnzkomkPyd/download/dmel-all-ncRNA-r6.18.fasta	ncRNA
    https://psilo.sorbonne-universite.fr/index.php/s/riKojdtpwB9xrKK/download/dmel-all-r6.18.gtf	dmel-all-r6.18.gtf
    https://psilo.sorbonne-universite.fr/index.php/s/F46CrxzzZxcrESa/download/dmel-all-transcript-r6.18.fasta	transcripts
    https://psilo.sorbonne-universite.fr/index.php/s/JM6aQkDMKF6YZRE/download/dmel-all-tRNA-r6.18.fasta	tRNAs
    https://psilo.sorbonne-universite.fr/index.php/s/A8dB5qPW3KgxNey/download/PlacW.fasta	PlacW
    ```

- Click the `Build` button
- In the `Build Rules ...` pannel that opens, click the ![](images/plus_rules.png){ width="80"}
and choose `Add/Modify Column Definitions`
- Click a first time on `Add Definition` and Select `URL`. Leave the URL column to `A`
- Click a second time on `Add Definition`, select `Name` and choose the column `B` for `Name`
- Now, click the `Apply` button
- And to finish the job, click on the dark-blue button `Upload`

<center>:tada:	:confetti_ball:	:balloon:</center>

#### 2. Upload of small RNA sequencing datasets :heavy_plus_sign: Programmatic dataset naming.

- Create a new history using the :heavy_plus_sign: icon of the history menu, and rename it
  `Small RNA sequence datasets`
- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- Click the `Rule-based`tab as we just did with the reference datasets
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:
!!! info ":ice_cream: small RNAseq datasets"
    ```
    https://psilo.sorbonne-universite.fr/index.php/s/HYLtfo9d2eD3Q2A/download/GRH-103_R1.fastq.gz	GRH-103
    https://psilo.sorbonne-universite.fr/index.php/s/C3o48iyyaeYw9gk/download/GRH-105_R1.fastq.gz	GRH-105
    ```
??? note "supplementary smallRNAseq datasets"
    You only need GRH-103 and GRH-105 datasets for the training.
    However, you can also experiment with additional datasets if you wish, to be downloaded
    as before using the following URL/name table:
    ```
    https://psilo.sorbonne-universite.fr/index.php/s/xgFNxA6fb7GgaLL/download/GRH-104_R1.fastq.gz	GRH-104
    https://psilo.sorbonne-universite.fr/index.php/s/Cp5ToXyj3X5jZxT/download/GRH-106_R1.fastq.gz	GRH-106
    ```

- Click the `Build` button
- In the `Build Rules ...` pannel that opened, click the ![](images/plus_rules.png){ width="80"}
and choose `Add/Modify Column Definitions`
- Click a first time on `Add Definition` and Select `URL`. Leave the URL column to `A`
- Click a second time on `Add Definition`, select `Name` and choose the column `B` for `Name`
- Now, click the `Apply` button
- select the Type "fastqsanger.gz" at the bottom of the panel

    ![](images/type_fastqsanger_gz.png){ width="200"}

- To finish the job, click on the dark-blue button `Upload`
<center>:tada:	:confetti_ball:	:balloon: :tada:	:confetti_ball:	:balloon:</center>

#### 3. RNAseq datasets (for gene differential expression analysis)

- Create a new history in Galaxy and rename it `RNA sequence datasets`
- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- Click the `Rule-based`tab as we just did with the reference datasets
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:
!!! info ":doughnut: RNAseq datasets"
    ```
    https://psilo.sorbonne-universite.fr/index.php/s/iaQxyLJzAXknCYp/download/SF1_R1.fastq.gz	SF1
    https://psilo.sorbonne-universite.fr/index.php/s/7HTaFnwfdQPCkL3/download/SF2_R1.fastq.gz	SF2
    https://psilo.sorbonne-universite.fr/index.php/s/ig6yZgEg9FPgRiw/download/SF3_R1.fastq.gz	SF3
    https://psilo.sorbonne-universite.fr/index.php/s/gyDgrq8S98C6KyL/download/WT1_R1.fastq.gz	WT1
    https://psilo.sorbonne-universite.fr/index.php/s/J9o73xmQWk43P6A/download/WT2_R1.fastq.gz	WT2
    https://psilo.sorbonne-universite.fr/index.php/s/YFnT3KzPNLdH3N4/download/WT3_R1.fastq.gz	WT3
    ```

- Click the `Build` button
- In the `Build Rules ...` pannel that opened, click the ![](images/plus_rules.png){ width="80"}
and choose `Add/Modify Column Definitions`
- Click a first time on `Add Definition` and Select `URL`. Leave the URL column to `A`
- Click a second time on `Add Definition`, select `Name` and choose the column `B` for `Name`
- Click the `Apply` button
- select the Type "fastqsanger.gz" at the bottom of the panel

    ![](images/type_fastqsanger_gz.png){ width="200"}

- And to finish the job, click on the dark-blue button `Upload`

<center>:tada:	:confetti_ball:	:balloon: :tada:	:confetti_ball:	:balloon:
:tada:	:confetti_ball:	:balloon: :tada:	:confetti_ball:	:balloon:</center>


#### 4. Uncompress datasets [_Section 4 should be optionnal_]

At this stage, we have uploaded small RNA and RNA sequencing datasets as `fastqsanger.gz`.
To simplify the subsequent analyzes we are going to uncompress all these datasets, whose
datatype will therefore become `fastqsanger`.

##### Procedure for a single dataset

  1. Go to your `small RNA input datasets` history (or whatever you named it).
  2. Click on the pencil icon ![](images/pencil.png){ width="70"} of the first dataset.
  3. Click on the tab `Convert` ![](images/convert.png){ width="100"}, _**NOT**_ on the tab `datatype`
  ![](images/datatypes.png){ width="100"}.
  
    ??? warning "Why 'Convert file' is different from 'Change Datatype' ?"
        - Let's imagine a Galaxy dataset whose name is `Hamlet`
        - the _content_ of this dataset is:
        ```
        To be, or not to be, that is the question:
        ```
        - Would you agree that the `datatype` of this dataset is ==`english`==? I think so.
        - Let's put it all together in the form of:
        ```
        @name: Hamlet
        @datatype: english
        @content:
        To be, or not to be, that is the question:
        ```
        
        **Now, what if you change the `Datatype` of this dataset from `english` to `french`
        using the `edit attribute` panel? ==This -->==**
        ```
        @name: Hamlet
        @datatype: french
        @content:
        To be, or not to be, that is the question:
        ```
        **This does not seem correct ! Do you aggree ?**
        
        **If you `Convert` instead this dataset from `english` to `french`, you will have
        ==This -->==**
        ```
        @name: Hamlet
        @datatype: french
        @content:
        Être ou ne pas être, telle est la question
        ```
        **It is looking better, isn't it ?**
        
        **In contrast, if your starting dataset was as this:**
        ```
        @name: Hamlet
        @datatype: english
        @content:
        Être ou ne pas être, telle est la question
        ```
        **There, you would "just" change the Datatype of the dataset from `english` to `french` and
        get**:
        ```
        @name: Hamlet
        @datatype: french
        @content:
        Être ou ne pas être, telle est la question
        ```
  
  4. Select `Convert compressed file to uncompressed`
  5. Click on ![](images/convert_datatype.png){ width="120"}
  
  ==**-->**== A new dataset is created. During the decompression job, its name looks like
  `5: Convert compressed file to uncompressed. on data 1`. But when the job finishes, the
  name of the dataset changes to more self-explanatory: `5: GRH-103 uncompressed`.

##### Repeat the same procedure for every small RNAseq dataset.

##### Repeat the same procedure for every RNAseq dataset.
==_Naturally, you can launch as many jobs as you need in the same time_==

##### When all datasets are decompressed

- Delete the compressed datasets (by clicking on the cross icon of datasets).
- Rename the uncompressed datasets by removing the `uncompressed` suffix.
- Purge the deleted datasets. This is done by clicking the wheel icon of the **top**
history menu, and selecting `Purge Deleted Datasets` in the **Datasets Actions** section.
    
    ![](images/wheel.png){ width="250" }
    
    - :warning: If you do not perform this last action, the deleted datasets remain on your
      instance disk !

#### 5. Dataset collections :milky_way: :alien:

If we have enough time, we are going to organize our various datasets using an additional
structure layer: the **Galaxy Collection**.

A Galaxy Collection is a container object which is very convenient to treat together multiple
equivalent datasets, such as a list of sequencing dataset, of text labels, of fasta sequences,
etc.

For those of you who are a bit familiar with Python language, a Galaxy Collection is actually
just a dictionary, whose `keys` are the names of the datasets in the collection (in Galaxy
these names are referred to as `element identifiers`), and `values` are the paths to the
corresponding datasets. Well, a dictionary as I said :stuck_out_tongue_winking_eye:

##### A. Making a collection of the small RNA sequence datasets.

For clarity, we are going first to _copy_ the small RNA sequence dataset from their initial
history to a **new** history.

- Go to your small RNAseq sequence datasets.
- Click on the wheel icon of the history top menu
    
    ![](images/wheel.png){ width="200" }

- Select the submenu `Copy Datasets` in the section `Dataset Actions`
- In the pop-up panel, `Source History:`, check-in the 4 small RNA sequencing datasets
- In the same pop-up panel, `Destination History:`, field `New history named`, write
  ```
  small RNAs in collection
  ```
- Click the `Copy History Items` button.
- Still on the same pop-up panel, at the top in a green area, you have now a :link: to the
  new history that was created and where the datasets were copied. Click on that link !
    
    ??? info "When you copy datasets in that way..."
        The new datasets actually do not take any space on your disk. New symbolic links to
        the actual files are only created.

- Now, that your are in the new history, click on the checkbox icon in the top area of the
  history.
    
    ![](images/history_checkbox.png){ width="250" }
    
- Check-in the 4 small RNA datasets
- In the menu `All 4 selected` (also in the top area of the history), select
  `Build Dataset List`
  
  ![build list](images/build_list.png){ width="250" }
  
- In the pop-up panel, just write a meaningful name in the field `Name`, something like
  ```
  Small RNA collection
  ```
- Reorganize the datasets order by dragging the datasets in the list. For instance, you may
  prefer here the order GRH-103 --> GRH-104 --> GRH-105 --> GRH-106, from top to botom.
  In this specific use case, you can also use the `alphabetic sorting` icon.
- Press the button `Create Collection`

??? question "What do you see when you click on name of the new dataset collection?"
    You see the content of the collection, with datasets identified with names called
    `element_identifiers.
    
    Click on  the `<< History` link, to come back to the normal history view.


??? question "what do you see if you click the `crossed eye` icon at the right corner ![](images/hidden.png){ width="150"} ? "
    You see the actual dataset contained in the Collection. If you click on `unhide` for
    each of these datasets, you will actually see both the container collection and the contained
    datasets !

##### B. Making 2 collections RNA sequence datasets.
For RNAseq datasets, collections are also very convenient. However, it is even better to
anticipate the type of analysis that you are going to perform. Indeed, you are going to
compare 3 "test" (mutant, treated, whatever...) datasets with 3 control datasets.

Therefore, we are going to organise the RNAseq datasets as 2 collections: a collection `WT`
and a collection `SF`.

- Go back to your RNAseq input datasets history
- As before, _copy_ the 6 RNAseq dataset to a new history which you will name
  `RNAseq dataset Collections`
- This time, create first a collection by only checking the three datasets `WT1`, `WT2`
  and `WT3`, which you will name:
  ```
  WT
  ```
- Create also a second collection by only checking the three datasets `SF1`, `SF2`
  and `SF3`, which you will name:
  ```
  SF
  ```
