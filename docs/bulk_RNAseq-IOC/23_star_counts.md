## ![](images/tool_small.png){width="30" align="absbottom"} Using `RNA STAR` for both alignment and read counting

We have already used the STAR aligner. But, for the sake of simplicity, we did not used its
integrated fonction which allows to counts reads after alignments, still using the appropriate
GTF input file.

This is what we are going to do in this section.

At first, navigate to the history `STAR Alignments` which we previously generated in the
section [STAR alignments](../18_star/#star-alignments).

From this history, copy (using the menu `copy datasets` item in the wheel history menu)

- [x] The three fastq.gz collections `5: Dc`, `10: Mo`, and `15: Oc`
- [x] and the GTF file `Mus_musculus.GRCm38.102.chr.gtf`

in a new history that you will name `STAR alignments AND counting`.

Navigate to this new history and run `RNA STAR` with the following settings


!!! info "![](images/tool_small.png){width="25" align="absbottom"} RNA STAR settings"
    - Single-end or paired-end reads
        
        --> Single-end
    - RNA-Seq FASTQ/FASTA file
        
        --> select the collection icon and then the collection `5: Dc`
    - Custom or built-in reference genome
        
        --> Use a built-in index
    - Reference genome with or without an annotation
        
        --> use genome reference without builtin gene-model but provide a gtf
    - Select reference genome
        
        --> GRCm38_w/o_GTF
    - Gene model (gff3,gtf) file for splice junctions
        
        --> Mus_musculus.GRCm38.102.chr.gtf
    - In `Output filter criteria`, **Exclude the following records from the BAM output**
        
        --> check Select all
    - Per gene/transcript output
        
        --> This time, select `Per gene read counts (GeneCounts)`
    - `Output filter criteria`, **Exclude the following records from the BAM output**
        
        --> check `Select all`

    The tool will run during several minutes, generating four new dataset collections, whose
    name is self-explanatory. Take benefit of the run time, to rename at least 3 of these
    collections with more meaningful names:
    
    - `RNA STAR on collection 5: log` --> `Dc STAR log`
    - `RNA STAR on collection 5: mapped.bam` --> `Dc RNA STAR mapped.bam`
    - `RNA STAR on collection 5: reads per gene` --> `Dc nbre of reads per gene (STAR)`
    
    :warning: Reminder: we understand it is a bit borring to rename datasets but these
    renaming operations are essential to the readibility of your histories.

### Re-run the RNA STAR tool for the collections:
- [x] `10: Mo`
- [x] `15: Oc`

This time, each run of `RNA STAR` generate a 5th dataset collection named
`RNA STAR on collection X: reads per gene`.

Rename these collections `Dc STAR counts`, `Mo STAR counts` and `Oc STAR counts`,
respectively.

## ![](images/tool_small.png){width="30" align="absbottom"} Mapping statistics with MultiQC tool

You can re-run MultiQC on the 3 RNA STAR log collection but note that we already permormed
this operation in the history `STAR alignments` with the section [18_star](../18_star/#mapping-statistics-with-multiqc-tool)

??? info "![](images/tool_small.png){width="25" align="absbottom"} MultiQC settings"
    - 1: Results
    - Which tool was used generate logs?
        
        --> STAR
    - Click "Insert STAR output"
    - Type of STAR output?
        
        --> Log
    - STAR log output
        
        --> Click first the collection icon ![](images/library_icon.png){width="75" align="absbottom"}
        
        --> Select the 3 collections `Dc`, `Mo` and `Oc RNA STAR log`, holding down the
        ++command++ key
    - Leave the other settings as is
    - Press `Execute` !

This is the occasion to use the `window manager`which you can trigger by clicking this
icon ![](images/window_manager_icon.png){width="250" align="absbottom"} (becomes yellow
when activated).

- [x] Click first on the eye of the collection `MultiQC on ... and others: Webpage` in the history
`STAR alignments AND counting`.
- [x] The web report opens in a floatting window in the center of the screen.
- Switch to the history `HISAT Alignments` using the history switch menu at the top of the history:
  
  ![](images/switch_history_icon.png){width="300" align="absbottom"}
  
- [x] Click on the eye of the collection `MultiQC on ... and others: Webpage` in the history
  `HISAT Alignments`.
- [x] You can now compare the results from both aligners, sided by side in the center of the
  screen.

## ![](images/tool_small.png){width="30" align="absbottom"} Adapt the format of STAR counts collections

