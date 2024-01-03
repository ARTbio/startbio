![](images/lamp.png)

# Analysis of functional enrichment among the differentially expressed genes

We have extracted genes that are differentially expressed in treated (Pasilla gene-depleted)
samples compared to untreated samples. We would like to know if there are categories of
genes that are enriched among the differentially expressed genes.

Gene Ontology (GO) analysis is widely used to reduce complexity and highlight biological
processes in genome-wide expression studies.

However, standard methods give biased results on RNA-seq data due to over-detection
of differential expression for long and highly-expressed transcripts.

The goseq tool provides methods for performing GO analysis of RNA-seq data,
taking length bias into account. The methods and software used by goseq are equally
applicable to other category based tests of RNA-seq data, such as KEGG pathway analysis.