# Manipulation of DESeq2 data for visualisation and comparisons

Now we would like to extract the most differentially expressed genes in the various
conditions, and then visualize them using an heatmap of the normalized counts for each
sample.

We will proceed in several steps:

- [x] For each package, extract the normalized counts of genes for each sample (all three
  packages, DESeq2, edgeR and limma, provide this functionality.
- [x] For each package, extract the most differentially expressed genes at a given log2FC
  threshold (let's say 2, corresponding to a 4x or 1/4x fold time expression), and at a
  given p-adjusted value (let's say p-adj < 0.01). We will keep these gene lists apart to
  build latter a venn diagram for comparison of the three tools.
- [x] Plot heatmaps of normalized counts

## Extract the most differentially expressed genes (PRJNA630433 / DESeq2)

Basically, we navigate in the DESeq history of the PRJNA630433 use-case and we repeat a
DESeq2 run, asking in addition for a **rLog-Normalized** counts output.

??? info "![](images/tool_small.png){width="25" align="absbottom"} `DESeq2` settings"
    Basically, the same as before, except that we ask for a Normalized counts file
    
    - how
        
        --> Select datasets per levels
    - 1: Factor
        
        --> Tissue
    - 1: Factor level
        
        Note that there will be three factor levels in this analysis: Dc, Mo and Oc.
        
        --> Oc
        
    - Counts file(s)
        
        --> select the data collection icon, then `15: Oc FeatureCounts counts`
    - 2: Factor level
        
        --> Mo
        
    - Counts file(s)
        
        --> select the data collection icon, then `10: Mo FeatureCounts counts`
    - 3: Factor level (you must click on :heavy_plus_sign: `Insert Factor level`)
        
        --> Dc
        
    - Counts file(s)
        
        --> select the data collection icon, then `5: Mo FeatureCounts counts`
    - (Optional) provide a tabular file with additional batch factors to include in the model.
        
        --> Leave to `Nothing selected`
    - Files have header?
        
        --> Yes
    - Choice of Input data
        
        --> Count data
    - Advanced options
        
        --> No, leave folded
    - Output options
        
        --> ==This time, check the `Output rLog normalized table` box !==
        
        --> Unfold and check `Output all levels vs all levels of primary factor (use when
        you have >2 levels for primary factor)` in addition to the already checked
        `Generate plots for visualizing the analysis results`
        
        --> Leave `Alpha value for MA-plot` to 0,1: note that this option is used for
        plots and does not impact DESeq2 results
    - `Run Tool`

:warning: This time you can trash the DESeq2 plots and result files which we have already
generated.

:warning: Keep this output for latter, will use it for a clustered heatmap

## Generate top lists of DE genes

We will do that with the help of the tool `Filter data on any column using simple
expressions`. We will also use 3 other tools `Compute on rows`, `Column Regex Find And
Replace` and `Filter data on any column using simple expressions`

### Select genes with |log2FC > 2| and p-adj < 0.01 with ![](images/tool_small.png){width="30" align="absbottom"}`Filter data on any column using simple expressions`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Filter data on any column...` settings"
    - Filter
        
        --> DESeq2 Results Tables
    - With following condition
        
        --> abs(c3) > 2 and c7 < 0.01
    - Number of header lines to skip
        
        --> `1` (these tables have an added header !)
    - Click the `Run Tool` button

:warning: Rename the "filter on..." collection to `Top gene lists`

### Compute a boolean value by row

This is to determine whether genes in the lists are up or down-regulated

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Compute on rows` settings"
    - Input file
        
        --> `Top gene lists` (:warning: collection !)
    - Input has a header line with column names?
        
        --> `Yes`
    - 1: Expressions
    - Add expression
        
        --> `c3 > 0`
    - Mode of the operation
        
        --> `Append`
    - The new column name
        
        --> `Regulation`
    - Avoid scientific notation in any newly computed columns
        
        --> `No`
    - Click the `Run Tool` button

:warning: Look at the effect of evaluating the expression `c3 > 0` in the new column
`expression` in the output datasets.

### Transform `True` and `False` values to `up` and `down`, respectively

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Column Regex Find And Replace` settings"
    - Select cells from
        
        --> `Compute on collection 36 (or so)`
    - using column
        
        --> `8`
    - Check
        
        --> click the button :heavy_plus_sign:`Insert Check`
    - Find Regex
        
        --> `False`
    - Replacement
        
        --> `down`
    - Check
        
        --> click another time the button :heavy_plus_sign:`Insert Check`
    - Find Regex
        
        --> `True`
    - Replacement
        
        --> `up`
    - Click the `Run Tool` button

:warning: rename the collection `Column Regex Find And Replace on collection 40` with
`top gene lists - oriented`

### Split the lists in `up` and `down` regulated lists

This will be performed through 2 successive runs of the
![](images/tool_small.png){width="25" align="absbottom"} tool `Select lines that match an
expression`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Select lines that match an expression` settings"
    - Select lines from
        
        --> `top gene lists - oriented`
    - that
        
        --> `matching`
    - the pattern
        
        --> `\tup` (a tabulation immediately followed by the string *up*)
    - Keep header line
        
        --> `Yes`
    - Click the `Run Tool` button

:warning: Immediately rename the collection `Select on collection...` to `top up-regulated
gene lists`

Redo exactly the same operation with a single change in the setting of the
![](images/tool_small.png){width="25" align="absbottom"} tool `Select lines that match an
expression`

??? info "![](images/tool_small.png){width="25" align="absbottom"} `Select lines that match an expression` settings"
    - Select lines from
        
        --> `top gene lists - oriented`
    - that
        
        --> `matching`
    - the pattern
        
        --> `\tdown` (a tabulation immediately followed by the string *down*)
    - Keep header line
        
        --> `Yes`
    - Click the `Run Tool` button

:warning: Rename the collection `Select on collection...` to `top down-regulated
gene lists`

:warning: keep the last three generated collections for later comparison with edgeR and
limma tools

## Plotting an heatmap of the most significantly de-regulated genes

For this, we are going to collect and gather all significantly de-regulated genes in any of the
3 conditions, and to intersect (join operation) this list with the rLog normalized count
table precedently generated.

### Use ![](images/tool_small.png){width="25" align="absbottom"}`Advanced cut` to select the list of deregulated genes in all three comparisons

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `advanced cut` settings"
    - File to cut
        
        --> `Top gene lists` (this is a collection)
    - Operation
        
        --> `Keep`
    - Delimited by
        
        --> `Tab`
    - Cut by
        
        --> `fields`
    - List of Fields
        
        --> `Column 1`
    - First line is a header line
    - Click the `Run Tool` button

:warning: Rename this collection of single column datasets `top genes names`
### Next we concatenate the three datasets of the previous collection in a single dataset

We do that using the ![](images/tool_small.png){width="25" align="absbottom"}
`Concatenate multiple datasets tail-to-head while specifying how` tool

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Concatenate multiple datasets tail-to-head while specifying how` settings"
    - What type of data do you wish to concatenate?
        
        --> `Single datasets`
    - Concatenate Datasets
        
        --> :warning: Click on the collection icon and select `top genes names`
    - Include dataset names?
        
        --> `No`
    - Number of lines to skip at the beginning of each concatenation:
        
        --> `1`
    - Click the `Run Tool` button

:warning: Rename the return single dataset as `Pooled top genes`

### Next we extract *Uniques* gene names from the `Pooled top genes` dataset

You probably agree that the same gene may be deregulated in the three pair-wise comparisons
which we have performed with DESeq2.

Thus we need to eliminate the redundancy, using the tool
![](images/tool_small.png){width="25" align="absbottom"}`Unique
occurrences of each record`.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Unique occurrences of each record` settings"
    - File to scan for unique values
        
        --> `Pooled top genes`
    - Ignore differences in case when comparing
        
        --> `No`
    - Column only contains numeric values
        
        --> `No`
    - Advanced Options
        
        --> Leave as `Hide Advanced Options`
    - Click the `Run Tool` button

### Add a header the list of unique gene names associated we significant DE in any of the comparisons

We do this with the tools `Add Header`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Add Header` settings"
    - List of Column headers (comma delimited, e.g. C1,C2,...)
        
        --> `All_DE_genes`
    - Data File (tab-delimted)
        
        --> `Unique on data 1xx...`
    - Click the `Run Tool` button

:warning: Rename the generated dataset `All_DE_genes`

### Intersection (join operation) between the list of unique gene name associated with DE and the rLog-Normalized counts file.

This is the moment when we are going to use the `rLog-Normalized counts file on data...`
and intersect it (join operation) with the list of DE genes in all three condition.

To do this, we are going to use the tool
![](images/tool_small.png){width="25" align="absbottom"}`Join two files`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Join two files` settings"
    - 1st file
        
        --> `rLog-Normalized counts file on data...`
    - Column to use from 1st file
        
        --> `1`
    - 2nd File
        
        --> `All_DE_genes`
    - Column to use from 2nd file
        
        --> `1`
    - Output lines appearing in
        
        --> `Both 1st and 2nd files`
    - First line is a header line
        
        --> `Yes`
    - Ignore case
        
        --> `No`
    - Value to put in unpaired (empty) fields
        
        --> `NA`
    - Click the `Run Tool` button

:warning: Rename the output dataset `rLog-Normalized counts of DE genes`

### Plot a heatmap of the rLog-Normalized counts of DE genes in all three conditions

We do this using the ![](images/tool_small.png){width="25" align="absbottom"}`Plot
heatmap with high number of rows` tool

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Plot heatmap with high number of rows` settings"
    - Input should have column headers - these will be the columns that are plotted
        
        --> `rLog-Normalized counts of DE genes`
    - Data transformation
        
        --> `Plot the data as it is`
    - Enable data clustering
        
        --> `Yes`
        
    - Clustering columns and rows
        
        --> `Cluster rows and not columns`
    - Distance method
        
        --> `Euclidean`
        
    - Clustering method
        
        --> `Complete`
    - Labeling columns and rows
        
        --> `Label columns and not rows`
        
    - Coloring groups
        
        --> `Blue to white to red`
    - Data scaling
        
        --> `Scale my data by row`
    - tweak plot height
        
        --> `35`
    - tweak row label size
        
        --> `1`
    - tweak line height
        
        --> `24`
    - `Run Tool`
