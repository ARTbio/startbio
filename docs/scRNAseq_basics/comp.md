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

|     | Row.names       | p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | external_gene_name | description                                                                             | gene_biotype   | chromosome_name |
|:--|:-----|--:|----:|--:|--:|----:|:------|:---------------------------|:-----|:-----|
| 1   | ENSG00000163453 |     0 |   7.673621 | 0.551 | 0.006 |  0.00e+00 | IGFBP7             | insulin-like growth factor binding protein 7 \[Source:HGNC Symbol;Acc:5476\]            | protein_coding | 4               |
| 2   | ENSG00000135077 |     0 |   5.656534 | 0.314 | 0.009 |  0.00e+00 | HAVCR2             | hepatitis A virus cellular receptor 2 \[Source:HGNC Symbol;Acc:18437\]                  | protein_coding | 5               |
| 3   | ENSG00000175294 |     0 |   5.513915 | 0.109 | 0.000 |  2.03e-05 | CATSPER1           | cation channel, sperm associated 1 \[Source:HGNC Symbol;Acc:17116\]                     | protein_coding | 11              |
| 224 | ENSG00000167286 |     0 |  -3.971611 | 0.083 | 0.885 |  0.00e+00 | CD3D               | CD3d molecule, delta (CD3-TCR complex) \[Source:HGNC Symbol;Acc:1673\]                  | protein_coding | 11              |
| 225 | ENSG00000137078 |     0 |  -4.379914 | 0.006 | 0.232 |  1.01e-05 | SIT1               | signaling threshold regulating transmembrane adaptor 1 \[Source:HGNC Symbol;Acc:17710\] | protein_coding | 9               |
| 226 | ENSG00000172116 |     0 |  -5.098008 | 0.019 | 0.368 |  0.00e+00 | CD8B               | CD8b molecule \[Source:HGNC Symbol;Acc:1707\]                                           | protein_coding | 2               |

Here are the results for the three most over-expressed genes in NK cells
(`avg_log2FC` positive) and the 3 most over-expressed genes in CD8+ T
cells (`avg_log2FC` negative).

You can totally use the different methods of gene cluster analysis
on this one.
