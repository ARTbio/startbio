## Upload input data to you newly deployed Galaxy instance
For the training, we need three types of datasets

- The reference sequences that will be used to align sequencing reads (full genome, miRNA, transposons, etc.)
- libraries of sequencing reads from small RNAs (for analysis of piRNAs)
- Librairies of sequencing reads from mRNA (for Gene differential expression analysis)

- Click the main menu `User` --> `Saved Histories`
- Press the top right button (above history list) `Import from file`
- copy this url :
```
https://galaxy.pasteur.fr/history/export_archive?id=4c5da5ad7355ff42
```
    
- repeat the same operation with: 
```
https://galaxy.pasteur.fr/history/export_archive?id=eb4c1d5564c9f78c
```
and
```
https://galaxy.pasteur.fr/history/export_archive?id=69a1b70d1c4a6bdb
```

??? bug "In case of emergency"
    
    ### NEXTCLOUD
    
    Small RNAseq
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/wC4DrxHN3gLtx4c/download
    ```
    Reférences
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/Sri6HKiwCn4RbSq/download
    ```
    RNAseq
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/4FWrWxZf72KDNji/download
    ```
    
    ### Amazon S3
    RNAseq
    ```
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq
    ```
    small RNAseq
    ```
    https://mydeepseqbucket.s3.amazonaws.com/smallrnaseqs
    ```
    references
    ```
    https://mydeepseqbucket.s3.amazonaws.com/references
    ```
    
