## Prepare your input data history

!!! note "Step-1: Transfer datasets from the Data library to a new history"
    1. Rename your `Unnamed history` to
    ```
    Input Dataset and collections
    ```
    2. Go to menu `Shared Data`--> `Data Libraries` (`Données Partagées` --> `Bibliothèque de Données`)
       ![](images/import_data.png){: style="width:450px"}
    
    3. Choose `Mouse Genetics` library
    4. Select the 4 fastq files (A_R1.fastq, A_R2.fastq, B_R1.fastq and B_R2.fastq)
    5. Select the `To History` tab --> `as datasets`
       ![](images/import_library_datasets.png){: style="width:600px"}
    
    6. Select your freshly renamed `Input Dataset and collections` in the `select history` menu
    7. Click `Import` button
    8. After the import, navigate directly to this history by clicking the `green warning`

!!! note "Step-2: Create data collections in you `Input Datasets and Collections` history"
    
    1. Prepare two `collections` from your raw input datasets.
        1. Toggle the "checkbox" mode by clicking the small checkbox icon at the top of the history bar
            ![](images/checkbox_mode.png){: style="width:250px"}
        2. Select the 2 A fastq files OR the 2 B fasq files (not all 4 files, choose as you feel it!)
        3. Select `Build List of Dataset Pairs` from the tab `Pour toute la sélection`
            ![](images/build_paired-dataset-collection.png){: style="width:250px"}
        4. in the pop up window, replace `_1` by `_R1` and `_2` by `_R2`
        5. Click the `Pair these datasets` tab
            ![](images/pair_datasets.png)
        6. Name your new "paired dataset" collection with a single element `A_fastq` (or `B_fastq`if you chose the B fastq file at the previous step) and click on `Create list`
        7. Back to your history, that is still in "checkbox" mode, select the 4 fastq files, and repeat the operation to produce this time a collection of 2 paired-sequences element, which you will name this time `patient sequences`

!!! note "Step-3: Send dataset collections in a new history"
    1. Select the `Copy datasets`in the history "wheel" menu 
        ![](images/copy_datasets.png){: style="width:250px"}
    - Select the first collection with a single element (A or B) that you first prepared
    - in the `destination history` area, fill the `New history named` field with
      ```
      Single sequence dataset analysis` and click the `Copy History Items
      ```
    - Click the "Copy" button
    - Click the link that shows up to navigate directely to this new history !
