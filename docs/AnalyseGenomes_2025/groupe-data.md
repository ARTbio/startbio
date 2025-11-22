Galaxy has numerous tools to analyse tables and variables they contains.
Here, we are going use the tool `Group data by a column and perform aggregate
operation on other columns. (Galaxy Version 2.1.4)` to rapidly extract information
from a complex GTF file.

## The input data

The input data is a GTF file annotation of the dmel 6.59 *D. melanogaster* genome which contains
only annotations for exons on chromosomes X, 2, 3 and 4.
This file is available on the GitHub repository `ARTbio/AnalyseGenome` at this URL:
```
https://github.com/ARTbio/AnalyseGenome/raw/refs/heads/main/Exercises/dmel.gtf.gz
```

1. Create a new history `Grouping exercise`, open the `Upload` menu (top-left corner, just below the Galaxy logo).
    
    Click the :pencil:`Paste/Fetch data` button, and paste the above URL in the central field of the panel.
    Instead of the pre-filled `New File`, type `dmel.gtf`, and click the `Start` button.
    
    - [ ] Once the upload is complete, you will notice that the dataset's data type is `gtf`, not `gtf.gz`. This is because Galaxy uncompresses most uploaded data on the fly. Notable exceptions to this default behavior are `fastq.gz` and `fasta.gz` files. Since these files are often very large, they are kept in their compressed form to save disk space on the server and can be processed directly.
    
    - [ ] Have a look to the content of the file. Standard GTF ! 9 columns, the last column contains 4 types of annotations:
        `gene_id`, `gene_symbol`, `transcript_id` and `transcript_symbol`.
 
2. The first thing we are going to do is simplifying the 9th column, keeping only the `gene_id` information.
    - [x] Select the  :wrench:`Column Regex Find And Replace` in the left tools bar.
    - [x] Select the `column 9` for the dataset `1: dmel.gtf`
    - [x] Click the :heavy_plus_sign:`Insert Check`
    - [x] In the **Find Regex** enter the regular expression
          ```
          ;.+
          ```
    - [x] In the Replacement field, do not enter anything
    - [x] Click the `Run`Button.
    - [x] In the resulting dataset `2: Column Regex Find And Replace on data 1`, check the column 9 (**Attributes**).
          What do you see ?

3. **First Grouping operation.**
     
     Now, we are going to group the data in the 9th column (we collapse lines by unique gene_id), while counting the "events"
     on column 3 (exon), and randomly picking a value on column 1 (the chromosome name, which is expected to be the same for
     all exons of a given gene)
     
     - [x] Select the  :wrench:`Group data by a column and perform aggregate operation on other columns` in the left tools bar.
     - [x] The **Select data** should already be `2: Column Regex Find And Replace on data 1`
     - [x] In the **Group by column** menu, select the 9th column.
     - [x] Click once on the :heavy_plus_sign:`Insert Operation`
     - [ ] In the **Type** menu, select `Count`
     - [ ] In the **On column** menu, select `Column: 3`
     - [x] Click one more time the :heavy_plus_sign:`Insert Operation`
     - [ ] In the **Type** menu, select `Randomly pick` (bottom of the menu)
     - [ ] In the **On column** menu, select `Column: 1`
     - [x] Click the `Run`Button !
     - [x] What do you see in the generated dataset ? --> It is a 3-column table. The first column contains unique
           gene_ids. The second column contains the computed number of exons for this gene. The third column contains
           the chromosome name for this gene.

4. **Second Grouping operation.**

    Finally, we are going the group again the data. This time, we group on the third column (chromosome names) while counting
    the variable (gene_id) in the 1 column.

     - [x] Select the  :wrench:`Group data by a column and perform aggregate operation on other columns` in the left tools bar.
     - [x] The **Select data** should already be `3: Group on data 2`
     - [x] In the **Group by column** menu, select the 3rd column.
     - [x] Click once on the :heavy_plus_sign:`Insert Operation`
     - [ ] In the **Type** menu, select `Count`
     - [ ] In the **On column** menu, select `Column: 1`
     - [x] Click the `Run`Button !
     - [x] What do you see in the generated dataset ? --> It is a 2-column table. The first column contains unique
           chromosome arms. The second column contains the computed number of genes for the chromosome arm.

           **:sunglasses:** ?

    




