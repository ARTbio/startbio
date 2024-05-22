# Differential expression analyses

We identified the different clusters but we don't know to which cell types they correspond
to. To be able to characterize them, we often use differential expression analyses. Please,
go read carefully [introduction](./intro_markers.md) and [marker annotation](./marker_annot.md).

---

![](../R-IOC/images/toolbox-do-it-yourself.png){: style="width:75px"} **Do it yourself!**

# Render your RMD/QMD

To complete this week you'll need to :

- [x] 1. Perfom the differential expression analyses by taking all deregulated genes (up **and**
         down), with `min.pct = 0.1` and `logfc.threshold = 0`.
- [x] 2. Annotate the result with biomaRt
- [x] 3. In a new dataframe, retrieve only significant result, all genes with a p-value adjusted
         inferior to 0.05.
- [x] 4. How many genes are significantly deregulated by cluster ? 
- [x] 5. Group the dataframe by cluster
- [x] 6. Retrieve the 5 most upregulated genes and 5 most downregulated genes by cluster
- [x] 7. Plot a heatmap of these genes
- [x] 8. Visualise the most upregulated gene for each cluster in a violin plot, describe the plot

!!! warning "IMPORTANT"
    Please note you **must** be explicative in your cutoff choices and
    provide detailed explanations for each step of your thought process 
    In general, try to explain in your own words for each step of your analysis!

Add your RMD/QMD in your Trello card. You can also add the html version of your rapport.


**Thank you for your attention and see you next week :clap: :clap: :clap:**

----