??? example "Step-1: Transfer datasets from the Data library to a new history"
    1. Rename your `Unnamed history` to
    ```
    Input Datasets and Collections
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

??? example "Step-2: Create a data collections in your `Input Datasets and Collections` history"
    
    1. Toggle the "checkbox" mode by clicking the small checkbox icon at the top of the history bar
            ![](images/checkbox_mode.png){: style="width:250px"}
    2. Select the 4 fastq files
    3. Select `Build List of Dataset Pairs` from the tab `Pour toute la sélection`
        ![](images/build_paired-dataset-collection.png){: style="width:350px"}
    4. in the pop up window, replace `_1` by `_R1` and `_2` by `_R2`
    5. Click the `Pair these datasets` tab
        
        ![](images/pair_datasets.png)
        
    6. Name your new "paired dataset" collection 
       ```
       Patient data
       ```
       and click on `Create list`
    7. Back to your history, you can now untoggle the "checkbox" mode

??? example "Step-3: Send dataset collections in a new history"
    1. Select the `Copy datasets`in the history "wheel" menu 
        
        ![](images/copy_datasets.png){: style="width:250px"}
        
    - Select the collection `Patient data` which you prepared
    - in the `destination history` area, fill the `New history named` field with
      ```
      Translocation analysis
      ```
    - Click the `Copy History Items`
    - Click the link that shows up to navigate directely to this new history !
