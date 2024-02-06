# Manipulation of limma data for visualisation and comparisons

Now we would like to extract the most differentially expressed genes in the various
conditions, and then visualize them using an heatmap of the normalized counts for each
sample.

## Extract the most differentially expressed genes (PRJNA630433 / limma)

Basically, we navigate in the limma history of the PRJNA630433 use-case and we **repeat a
limma run**, asking in addition for a file containing the normalised counts, these are in
**log2 counts per million (logCPM)**.

Note that, as edgeR, limma returns **log2 counts per million (logCPM)**.

??? info "![](images/tool_small.png){width="25" align="absbottom"} `limma` settings"
    - Differential Expression Method
        
        --> `limma-voom`
    - Apply voom with sample quality weights?
        
        --> `No`
    - Count Files or Matrix?
        
        --> Separate Count Files
    - 1: Factor/Name
        
        --> Tissue
    - 1: Factor/1: Group
        
        Note that there will be three Groups (ie factor levels) in this analysis: Dc, Mo and Oc.
        
        --> Oc
        
    - Counts file(s)
        
        --> select the data collection icon, then `15: Oc FeatureCounts counts`
    - 2: Factor/2: Group
        
        --> Mo
        
    - Counts file(s)
        
        --> select the data collection icon, then `10: Mo FeatureCounts counts`
    - 3: Factor level (you must click on :heavy_plus_sign: `Insert Group`)
        
        --> Dc
        
    - Counts file(s)
        
        --> select the data collection icon, then `5: Mo FeatureCounts counts`
    - Use Gene Annotations?
        
        --> `No`
     - Input Contrast information from file?
       
       --> `No`
    - 1: Constrast
        
        --> `Mo-Dc`
    - 2: Constrast (click :heavy_plus_sign: `Insert Contrast`)
        
        --> `Oc-Dc`
    - 3: Constrast (click :heavy_plus_sign: `Insert Contrast`)
        
        --> `Oc-Mo`
    - Filter Low Counts
        
        --> No, leave folded
    - Output options
        
        --> Unfold and select `Output Normalised Counts Table?: Yes`
    - Advanced options
        
        --> Put `P-Value Adjusted Threshold` to 0.1 (to be consistent with DESeq settings)
        
        --> Leave other advanced options unchanged
    - `Run Tool`

:warning: Note that limma is nicer than edgeR (in Galaxy) and return a single extra dataset
 `limma on data ... and others: Normalised counts`. Thus, no need here to do collection
 manipulations.

Beside, as with edgeR, a `limma-voom_normcounts.tsv` also show up as a html link in the dataset `limma on
data 4, data 3, and others: Report`, that download directly to your local computer if you
click it.

:warning: Keep the dataset `limma on ... and others: Normalised counts` for latter, we will use
it for the clustered heatmap.

## Generate top lists of limma DE genes

### Select genes with |log2FC > 2| and p-adj < 0.01 with ![](images/tool_small.png){width="30" align="absbottom"}`Filter data on any column using simple expressions`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Filter data on any column...` settings"
    - Filter
        
        --> limma on data ... others: DE tables (:warning: this is a collection)
    - With following condition
        
        --> `abs(c2) > 2 and c6 < 0.01` :warning: this expression is different from the one
        used for DESeq2 tables because the column structure is different. Actually this is
        the same expression that the one used for edgeR.
    - Number of header lines to skip
        
        --> `1` (these tables have an added header !)
    - Click the `Run Tool` button

:warning: Rename the "filter on..." collection to `limma Top gene lists`

### Compute a boolean value by row

This is to determine whether genes in the lists are up or down-regulated

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Compute on rows` settings"
    - Input file
        
        --> `limma Top gene lists` (:warning: collection !)
    - Input has a header line with column names?
        
        --> `Yes`
    - 1: Expressions
    - Add expression
        
        --> `c2 > 0`
        tables
    - Mode of the operation
        
        --> `Append`
    - The new column name
        
        --> `Regulation`
    - Avoid scientific notation in any newly computed columns
        
        --> `No`
    - Click the `Run Tool` button

:warning: Look at the effect of evaluating the expression `c2 > 0` in the new column
`expression` in the output datasets.

### Transform `True` and `False` values to `up` and `down`, respectively

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Column Regex Find And Replace` settings"
    - Select cells from
        
        --> `Compute on collection 32 (or so)`
    - using column
        
        --> `8` (for limma there are now 8 columns...)
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

:warning: rename the collection `Column Regex Find And Replace on collection 44` with
`limma top gene lists - oriented`

### Split the lists in `up` and `down` regulated lists

This will be performed through 2 successive runs of the
![](images/tool_small.png){width="25" align="absbottom"} tool `Select lines that match an
expression`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Select lines that match an expression` settings"
    - Select lines from
        
        --> `limma top gene lists - oriented`
    - that
        
        --> `matching`
    - the pattern
        
        --> `\tup` (a tabulation immediately followed by the string *up*)
    - Keep header line
        
        --> `Yes`
    - Click the `Run Tool` button

:warning: Immediately rename the collection `Select on collection...` to `limma top up-regulated
gene lists`

Redo exactly the same operation with a single change in the setting of the
![](images/tool_small.png){width="25" align="absbottom"} tool `Select lines that match an
expression`

??? info "![](images/tool_small.png){width="25" align="absbottom"} `Select lines that match an expression` settings"
    - Select lines from
        
        --> `limma top gene lists - oriented`
    - that
        
        --> `matching`
    - the pattern
        
        --> `\tdown` (a tabulation immediately followed by the string *down*)
    - Keep header line
        
        --> `Yes`
    - Click the `Run Tool` button

:warning: Rename the collection `Select on collection...` to `limma top down-regulated
gene lists`

:warning: keep the last three generated collections for later comparison with edgeR and
DESeq2 tools

## Plotting an heatmap of the most significantly de-regulated genes

For this, we are going to collect and gather all significantly de-regulated genes in any of the
3 conditions, and to intersect (join operation) this list with the rLog normalized count
table precedently generated.

### Use ![](images/tool_small.png){width="25" align="absbottom"}`Advanced cut` to select the list of deregulated genes in all three comparisons

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `advanced cut` settings"
    - File to cut
        
        --> `limma Top gene lists` (this is a collection)
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

:warning: Rename this collection of single column datasets `limma top genes names`

### Next we concatenate the three datasets of the previous collection in a single dataset

We do that using the ![](images/tool_small.png){width="25" align="absbottom"}
`Concatenate multiple datasets tail-to-head while specifying how` tool

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Concatenate multiple datasets tail-to-head while specifying how` settings"
    - What type of data do you wish to concatenate?
        
        --> `Single datasets`
    - Concatenate Datasets
        
        --> :warning: Click on the collection icon and select `limma top genes names`
    - Include dataset names?
        
        --> `No`
    - Number of lines to skip at the beginning of each concatenation:
        
        --> `1`
    - Click the `Run Tool` button

:warning: Rename the return single dataset as `limma Pooled top genes`

### Next we extract *Uniques* gene names from the `Pooled top genes` dataset

You probably agree that the same gene may be deregulated in the three pair-wise comparisons
which we have performed with DESeq2.

Thus we need to eliminate the redundancy, using the tool
![](images/tool_small.png){width="25" align="absbottom"}`Unique
occurrences of each record`.

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Unique occurrences of each record` settings"
    - File to scan for unique values
        
        --> `limma Pooled top genes`
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
        
        --> `limma_All_DE_genes`
    - Data File (tab-delimted)
        
        --> `Unique on data 7x...`
    - Click the `Run Tool` button

:warning: Rename the generated dataset `limma_All_DE_genes`

### Intersection (join operation) between the list of unique gene name associated with DE and the rLog-Normalized counts file.

This is the moment when we are going to use the dataset `limma on ... and others: Normalised counts`
 and intersect it (join operation) with the list of DE genes in all three condition.

To do this, we are going to use the tool
![](images/tool_small.png){width="25" align="absbottom"}`Join two files`

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Join two files` settings"
    - 1st file
        
        --> select `limma on ... and others: Normalised counts`
    - Column to use from 1st file
        
        --> `1`
    - 2nd File
        
        --> `limma_All_DE_genes`
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

:warning: Rename the single-element output collection `limma Log2CPM Normalized counts of DE genes`

### Plot a heatmap of the Log2CPM Normalized counts of limma DE genes in all three conditions

We do this using the ![](images/tool_small.png){width="25" align="absbottom"}`Plot
heatmap with high number of rows` tool

!!! info "![](images/tool_small.png){width="25" align="absbottom"} `Plot heatmap with high number of rows` settings"
    - Input should have column headers - these will be the columns that are plotted
        
        --> Click the collection icon and select `limma Log2CPM Normalized counts of DE genes`
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
