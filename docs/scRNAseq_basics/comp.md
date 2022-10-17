# Comparing two populations

One last thing we might ask is what makes two populations distinct.
This could be two conditions or even two clusters. This is possible
with the `FindMarkers` function which works almost like `FindAllMarkers`
since it calls the former. Many of the parameters are equivalent,
however `FindMarkers` allows you to perform a differential expression
analysis between two populations which are defined with the `ident.1`
and `ident.2` parameters. By default, these are cell identities present
in the `active.ident` but we can select another variable contained in
the metadata using the `group.by` parameter. So if we want to study the
impact of gender in the cells of the platelet then we would run:

`FindMarkers(objectName, ident.1 = "female", ident.2 = "male", group.by = "sex", subset.ident = "Platelet")`

(if we had a "sex" column in the metadata slot). If we want to study
the impact of sex in all cells then we leave the `subset.ident` parameter
as default (i.e. `NULL`).

Here we will test the difference between our NK and CD8+ T clusters which
were very difficult to differentiate. The results are similar to
`FindAllMarkers` with the difference that there is no `gene` column
because as a differential analysis it will not be possible to have a
gene overexpressed in both NK and CD8+ T cells.

``` r
NK_CD8_diff_markers <- FindMarkers(pbmc_small,
                           ident.1 = "NK",
                           ident.2 = "CD8+ T")

## Merge markers results with biomart annotation
NK_CD8_diff_markers_annotated <- merge(x = NK_CD8_diff_markers,  #First df to merge
                                       y = annotated_hg19,       #Second df to merge
                                       by.x = 0,                 #Column name of first df used for matching lines, 0 for rownames
                                       by.y = "ensembl_gene_id", #Column name of second df used for matching lines
                                       all.x = TRUE)             #Keep all lines from first df even if there is no match with second df

## Filter dataset based on Fold change and p-value adjusted
NK_CD8_diff_markers_annotated_signif <- subset(NK_CD8_diff_markers_annotated,
                                               p_val_adj < 0.05 &
                                                 abs(avg_log2FC) >= 0.25)       #Filter dataframe based on p_val_adj column

## Sorting results by average log2(Fold Change)
NK_CD8_diff_markers_annotated_signif <- NK_CD8_diff_markers_annotated_signif %>%                 #Rearrange df with dplyr package
  arrange(desc(avg_log2FC))                  #Sort lines by descending the column avg_log2FC and by group

## Most DE gene marker for each cluster
kable(NK_CD8_diff_markers_annotated_signif[(c(1:3, (nrow(NK_CD8_diff_markers_annotated_signif)-2):nrow(NK_CD8_diff_markers_annotated_signif))),])
```

|     | Row.names       | p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | external_gene_name | description                                                                                                  | gene_biotype   | chromosome_name |
|:----|:----------------|------:|-----------:|------:|------:|----------:|:-------------------|:-------------------------------------------------------------------------------------------------------------|:---------------|:----------------|
| 1   | ENSG00000115523 |     0 |   3.183558 | 0.968 | 0.276 |         0 | GNLY               | granulysin \[Source:HGNC Symbol;Acc:4414\]                                                                   | protein_coding | 2               |
| 2   | ENSG00000100453 |     0 |   2.617544 | 0.955 | 0.256 |         0 | GZMB               | granzyme B (granzyme 2, cytotoxic T-lymphocyte-associated serine esterase 1) \[Source:HGNC Symbol;Acc:4709\] | protein_coding | 14              |
| 3   | ENSG00000011600 |     0 |   1.955812 | 0.904 | 0.162 |         0 | TYROBP             | TYRO protein tyrosine kinase binding protein \[Source:HGNC Symbol;Acc:12449\]                                | protein_coding | 19              |
| 216 | ENSG00000227507 |     0 |  -1.197050 | 0.263 | 0.621 |         0 | LTB                | lymphotoxin beta (TNF superfamily, member 3) \[Source:HGNC Symbol;Acc:6711\]                                 | protein_coding | 6               |
| 217 | ENSG00000113088 |     0 |  -1.388411 | 0.115 | 0.582 |         0 | GZMK               | granzyme K (granzyme 3; tryptase II) \[Source:HGNC Symbol;Acc:4711\]                                         | protein_coding | 5               |
| 218 | ENSG00000167286 |     0 |  -1.749184 | 0.083 | 0.885 |         0 | CD3D               | CD3d molecule, delta (CD3-TCR complex) \[Source:HGNC Symbol;Acc:1673\]                                       | protein_coding | 11              |

Here are the results for the three most over-expressed genes in NK cells
(`avg_log2FC` positive) and the 3 most over-expressed genes in CD8+ T
cells (`avg_log2FC` negative).

You can totally use the different methods of gene cluster analysis
on this one.
