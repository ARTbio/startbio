For the training, we need three types of datasets

- The reference sequences that will be used to align sequencing reads (full genome, miRNA, transposons, etc.)
- libraries of sequencing reads from small RNAs (for analysis of piRNAs)
- Librairies of sequencing reads from mRNA (for Gene differential expression analysis)

All these data have been deposited in 2 differents repositories. A first one is a so-called
S3 Amazon bucket. The second one is a
[Nextcloud server](https://usegalaxy.sorbonne-universite.fr/nextcloud) located at
Sorbonne-Université. You may get your input data from one or the other repositories.

### Get data "by URL"
We are going to focus on one method to upload data in galaxy, which is applicable when these
data _**are available through a URL**_ (Universal Resource Location).

??? info "The other methods to upload data in Galaxy are:"
    * transfering data from your local machine (the one that is running your web browser)
      to Galaxy
    * uploading data to your Galaxy FTP account and then transfering these data from your
    Galaxy FTP directory to one of your Galaxy histories. We are not going to use them in
    this training, and invite you to look at one of the "Galaxy tours" available
    in the menu `Help` :arrow_forward: `Interactive tours`

#### 1. Single URL, simple trial.

- Click the `Upload Data` button at the top-left corner of the Galaxy interface:

![](images/galaxy_upload_button.png){: style="width:200px"}

- Stay with the regular tab and click the `Paste/Fetch data` button

![](images/regular_upload.png){: style="width:600px"}

- Paste the following url in the open text field,
```
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=PlacW.fasta
```
- Paste `PlacW.fasta` in the name text field (instead of `New File`)
- Finally, press the dark-blue `Start` button.

==-->== a dataset should appear soon in your current history and turn green when the upload is
complete.
    
#### 2. Upload of reference files as a batch of multiple URLs :heavy_plus_sign: Programmatic file naming

Delete the previously uploaded dataset, we are going to re-upload it in a batch.

- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- This time, Click the `Rule-based`tab !
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:
```
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-r6.18.gtf	dmel-all-r6.18.gtf
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-miscRNA-r6.18.fasta	miscRNA
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=PlacW.fasta	PlacW
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-ncRNA-r6.18.fasta	ncRNA
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-miRNA-r6.18.fasta	miRNA
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-intron-r6.18.fasta	introns
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-gene-r6.18.fasta	genes
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=Dmel_piRNA_clusters.fasta	piRNA_clusters
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=Dmel_all-transposon_merge.fasta	all-transposons
https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/B433xtdmdQqdFYd/download?path=%2F&files=dmel-all-chromosome-r6.18.fasta	dmel-all-chromosome-r6.18
```

??? info "To balance the load of the data servers, a half of the trainees may also use the S3 Amazon bucket"
    ```
    https://mydeepseqbucket.s3.amazonaws.com/References/PlacW.fasta	PlacW
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-ncRNA-r6.18.fasta	ncRNA
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-miscRNA-r6.18.fasta	miscRNA
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-miRNA-r6.18.fasta	miRNA
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-intron-r6.18.fasta	introns
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-gene-r6.18.fasta	genes
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-chromosome-r6.18.fasta	dmel-r6.18
    https://mydeepseqbucket.s3.amazonaws.com/References/Dmel_piRNA_clusters.fasta	piRNA_clusters
    https://mydeepseqbucket.s3.amazonaws.com/References/Dmel_all-transposon_merge.fasta	transposons
    https://mydeepseqbucket.s3.amazonaws.com/References/dmel-all-r6.18.gtf	dmel-all-r6.18.gtf
    ```
- Click the `Build` button
- In the `Build Rules ...` pannel that opened, click the ![](images/plus_rules.png){ width="80"}
and choose `Add/Modify Column Definitions`
- Click a first time on `Add Definition` and Select `URL`. Leave the URL column to `A`
- Click a second time on `Add Definition`, select `Name` and choose the column `B` for `Name`
- Now, click the `Apply` button
- And to finish the job, click on the dark-blue button `Upload`
- After the upload is complete, rename the history "References"

<center>:tada:	:confetti_ball:	:balloon:</center>

#### 3. Upload of small RNA sequencing datasets :heavy_plus_sign: Programmatic dataset naming.

Before all, create a new history by clicking the **+** icon in the history header
![](images/history_header.png){ width="300"} and immediately renaming the new history as
**"Small RNA sequence datasets"**.

- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- Click the `Rule-based`tab as we just did with the reference datasets
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:
??? info "from the Nextcloud server"
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=GRH-103_R1.fastq.gz	GRH-103
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=GRH-104_R1.fastq.gz	GRH-104
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=GRH-105_R1.fastq.gz	GRH-105
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=GRH-106_R1.fastq.gz	GRH-106
    ```
Or
??? info "from the S3 Amazon bucket"
    ```
    https://mydeepseqbucket.s3.amazonaws.com/smRNAseq/GRH-103_R1.fastq.gz	GRH-103
    https://mydeepseqbucket.s3.amazonaws.com/smRNAseq/GRH-104_R1.fastq.gz	GRH-104
    https://mydeepseqbucket.s3.amazonaws.com/smRNAseq/GRH-105_R1.fastq.gz	GRH-105
    https://mydeepseqbucket.s3.amazonaws.com/smRNAseq/GRH-106_R1.fastq.gz	GRH-106
    ```
- Click the `Build` button
- In the `Build Rules ...` pannel that opened, click the ![](images/plus_rules.png){ width="80"}
and choose `Add/Modify Column Definitions`
- Click a first time on `Add Definition` and Select `URL`. Leave the URL column to `A`
- Click a second time on `Add Definition`, select `Name` and choose the column `B` for `Name`
- Now, click the `Apply` button
- select the Type "fastqsanger.gz" at the bottom of the panel

    ![](images/type_fastqsanger_gz.png){ width="200"}

- And to finish the job, click on the dark-blue button `Upload`
<center>:tada:	:confetti_ball:	:balloon: :tada:	:confetti_ball:	:balloon:</center>

#### 4. RNAseq datasets (for gene differential expression analysis)

- Create a new history in Galaxy and rename it `RNA sequence datasets`
- Click the `Upload Data` button at the top-left corner of the Galaxy interface.
- Click the `Rule-based`tab as we just did with the reference datasets
- Leave **Upload data as** `Datasets` and **Load tabular data from** `Pasted Table`
- In the text field `Tabular source data to extract collection files and metadata from`,
paste the following Tabular source data:
??? info "from the Nextcloud server"
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=WT1_R1.fastq.gz	WT1
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=WT2_R1.fastq.gz	WT2
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=WT3_R1.fastq.gz	WT3
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=SF1_R1.fastq.gz	SF1
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=SF2_R1.fastq.gz	SF2
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/LqKb3Qmy8m9RXtk/download?path=%2F&files=SF3_R1.fastq.gz	SF3
    ```
Or
??? info "from the S3 Amazon bucket"
    ```
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/WT1_R1.fastq.gz	WT1
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/WT2_R1.fastq.gz	WT2
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/WT3_R1.fastq.gz	WT3
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/SF1_R1.fastq.gz	SF1
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/SF2_R1.fastq.gz	SF2
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq/SF3_R1.fastq.gz	SF3
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


#### 5. Uncompress datasets

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
        @datatype: english
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

#### 6. Dataset collections :milky_way: :alien:

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
- In the menu `Pour toute la sélection` (also in the top area of the history), select
  `Build Dataset List`
- In the pop-up panel, just write a meaningful name in the field `Name`, something like
  ```
  Small RNA collection
  ```
- Press the button `Create Collection`

??? question "What do you see when you click on name of the new dataset collection? (please not the :heavy_multiplication_x:...)"
    You see the content of the collection, with datasets identified with names called
    `element_identifiers.
    
    Click on the `recycling` icon ![](images/recycle.png){ width="20"}, or, the `< back to
    the Small RNA Collection` link, to come back to the normal history view.


??? question "what do you see if you click the `hidden` hyperlink at the top right corner ![](images/hidden.png){ width="150"} ? "
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

#### This is the end of this training session, you deserve :coffee: or :beer: ! 



