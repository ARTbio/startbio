## Add the URL for GDC downloads to the UUIDs of files in the metadata table 

The first column in the metadata file `gdc_sample_sheet.yyyy-mm-dd.tsv` has `File ID` as
a header. `FILE ID` actually stands for UUID, which in turn means Unique Universal Identifier.
This UUID allows to find files in the GDC repository file system, provided that the correct
URL is formed with the UUID.

Thus, to download the file whose `File Name` (second column) is `2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz`
it suffices to form the following URL
```
https://api.gdc.cancer.gov/data/0d3ddab0-4c6e-4f96-88dc-3c9813d9c292
```
by adding the prefix string `https://api.gdc.cancer.gov/data/`
to the UUID of the file `2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz`, namely
`0d3ddab0-4c6e-4f96-88dc-3c9813d9c292`

You can already test this URL by copying and pasting it in your web browser. This should trigger
the download of a file named `2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz`

If you are able to use a linux terminal in your computer, you can also use the same URL
in the following command line:

```
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/0d3ddab0-4c6e-4f96-88dc-3c9813d9c292'
```

This will also trigger locally the download of the `2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz`
file.

In addition, you can modify this command line by concatenating several UUIDs separated by commas,
resulting in the download of the several corresponding files. Thus, this command:

```
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/0d3ddab0-4c6e-4f96-88dc-3c9813d9c292,cc802aa8-e299-40c3-9539-a854ef950ff6,85310996-75c9-4b5b-a16e-406cde0a3373,0dc459fe-ac56-42af-9c99-3fe3e14bc156'
```
triggers the download of the files:
```
2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz
6f475a4d-53ea-4110-b80f-d519289d41f9.FPKM.txt.gz
7d07636d-e4cb-4642-83ee-c04b21d15227.FPKM.txt.gz
a405cab4-28d0-4cf3-b3c0-2007539002f3.FPKM.txt.gz
```
In a gzipped archive `gdc_download_20200219_152427.247343.tar.gz`

So far, so good. However we agree that we should not repeat this operation for > 5 files and
we need a more automated procedure.

### SED in rescue

First of all, we are going to replace all UUIDs in the first column of the `gdc_sample_sheet.yyyy-mm-dd.tsv`
metadata file, by the corresponding URL for downloads.

This can be done using regular expressions. If your text editor has regex replace function,
it is fine using it. Here we show how to do it in a terminal with the `sed` program, which will
work in linux and MacOSX environments:

```
sed -i '.bak' -E  's#^([^F])#https://api.gdc.cancer.gov/data/\1#' /home/chris/DOWNLOADS/gdc_sample_sheet.2020-02-19.tsv
```
where

`-i` option triggers the backup of the file in case of problems in the search/replace procedure

`-E` is for using modern regular expressions

`s#` is the substitute fonction  (s) of sed

`#^([^F])#` is the search expression between 2 #. It searches for a motif
at the beginning of a line ^ that **is not** an F (allowing to skip the first header line)
and capture it in the \1 variable

`#https://api.gdc.cancer.gov/data/\1#` is the replacement expression, \1 being the content
of the variable \1 being captured by the parenthesis in the search expression.

and

`/home/chris/DOWNLOADS/gdc_sample_sheet.2020-02-19.tsv` is the path to the file to be edited
by sed.

This sed command should generate de following transformed `gdc_sample_sheet.2020-02-19.tsv`

```
File ID	File Name	Data Category	Data Type	Project ID	Case ID	Sample ID	Sample Type
https://api.gdc.cancer.gov/data/0d3ddab0-4c6e-4f96-88dc-3c9813d9c292	2a3814fb-e9b1-404c-95da-db348ecfa14a.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PAEFGT	TARGET-20-PAEFGT-03A	Primary Blood Derived Cancer - Peripheral Blood
https://api.gdc.cancer.gov/data/cc802aa8-e299-40c3-9539-a854ef950ff6	6f475a4d-53ea-4110-b80f-d519289d41f9.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PASMHY	TARGET-20-PASMHY-03A	Primary Blood Derived Cancer - Peripheral Blood
https://api.gdc.cancer.gov/data/85310996-75c9-4b5b-a16e-406cde0a3373	7d07636d-e4cb-4642-83ee-c04b21d15227.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PASWAJ	TARGET-20-PASWAJ-04A	Recurrent Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/0dc459fe-ac56-42af-9c99-3fe3e14bc156	a405cab4-28d0-4cf3-b3c0-2007539002f3.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PAMVKZ	TARGET-20-PAMVKZ-09A	Primary Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/e8a0e021-96cd-4ced-b881-da6b31351ac3	a9fd9d1c-348f-41a9-9da2-a72981ef1ae5.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PANFMG	TARGET-20-PANFMG-09A	Primary Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/5b28c25d-8775-403d-b286-dc88fa3a2d17	bd29f662-4f53-43ce-ad0f-ebcb62546109.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PAMYGX	TARGET-20-PAMYGX-09A	Primary Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/23735039-36e9-4029-bb57-e2ccab9598ed	f49764c5-6b18-435a-90de-64f114a4ce89.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PASVVS	TARGET-20-PASVVS-04A	Recurrent Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/30a6e42d-81de-4837-abb4-84a3687f60de	2f601196-beaf-4a1e-91ee-72ef1de17f98.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PANHYK	TARGET-20-PANHYK-09A	Primary Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/6694923c-6ac6-4ba8-bddc-b35f438924b9	50b79e4d-0e67-4d79-914f-924e0ea920d9.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PAPWYK	TARGET-20-PAPWYK-09A	Primary Blood Derived Cancer - Bone Marrow
https://api.gdc.cancer.gov/data/a0f31602-503b-4964-b68c-1ef996536a74	b5514b40-5b8b-4faf-a202-49438c386032.FPKM.txt.gz	Transcriptome Profiling	Gene Expression Quantification	TARGET-AML	TARGET-20-PARGVC	TARGET-20-PARGVC-03A	Primary Blood Derived Cancer - Peripheral Blood
```

the original file being backed up in gdc_sample_sheet.2020-02-19.tsv.bak


Et Voilà... We can now use this metada file to upload all expression data file we have selected in Galaxy.



    

