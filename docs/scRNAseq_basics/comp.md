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
paged_table(NK_CD8_diff_markers_annotated_signif[(c(1:3, (nrow(NK_CD8_diff_markers_annotated_signif)-2):nrow(NK_CD8_diff_markers_annotated_signif))),])
```

Here are the results for the three most over-expressed genes in NK cells
(`avg_log2FC` positive) and the 3 most over-expressed genes in CD8+ T 
cells (`avg_log2FC` negative). 

You can totally use the different methods of gene cluster analysis 
on this one. 
