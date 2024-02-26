## Analyzing GO Enrichment from DEGs

Gene Ontology (GO) enrichment analysis is used to
identify `biological processes`, `cellular components`, and `molecular functions` that are
significantly over-represented (or under-represented) in a set of genes compared to a
background list. This is particularly valuable when analyzing differentially expressed
genes (DEGs) identified from RNA-seq or microarray experiments.

### Individual Gene Analysis (IGA)

- [x] **Concept**
  
  This approach tests each GO term individually for enrichment within the DEG list.

- [x] **Methods:**
    * **Hypergeometric test:** Calculates the probability of observing
the number of DEGs in a specific GO term by chance.
    * **Fisher's exact test:** Similar to
the hypergeometric test but suitable for smaller datasets.
  
- [x] **Limitations:**
    
    * Ignores the hierarchical structure of GO, potentially missing related terms.
    * Susceptible to multiple testing issues, requiring correction methods like Bonferroni
    adjustment.

- [x] **Advanced Considerations:**
    
    * **Multiple Testing Correction:** As mentioned, IGA is susceptible to multiple testing
    issues. Here are some commonly used correction methods:
        * *Bonferroni adjustment:* A conservative approach that controls the family-wise
        error rate (FWER) but can be overly stringent.
        * *Benjamini-Hochberg (BH) procedure:* Controls the false discovery rate (FDR) and
        is less conservative than Bonferroni.
        * *False discovery rate (q-value):* Provides a measure of significance adjusted
        for multiple testing.
    * **Gene Ontology Consortium (GOC) recommendations:** The GOC recommends using a
    combination of statistical significance (p-value) and fold change thresholds to
    identify relevant enriched terms, acknowledging the limitations of p-values alone.

### Gene Set Analysis (GSA)

- [x] **Concept:**
  
  Considers the entire set of DEGs and their relationships within the GO hierarchy.

- [x] **Methods:**
    
    * **Pathway analysis tools:** Tools like Enrichr, clusterProfiler, and GSEA analyze
    pre-defined gene sets like KEGG pathways and analyze enrichment within DEGs.
    * **GO-based GSA methods:** 
        * **Rank-based approaches:** Assign a rank to each gene based on its differential
        expression and analyze enrichment within ranked gene sets. (e.g., GSEA)
        * **Permutation-based approaches:** Randomly shuffle gene labels and recalculate
        enrichment scores to assess statistical significance. (e.g., fgsea)
        - Tools like GOseq, fgsea, and piano utilize various
        statistical models to account for the hierarchical structure of GO and identify
        enriched functional categories.

- [x] **Advantages of using GSA:**
    * Incorporates information about gene relationships within the GO hierarchy, leading
    to more biologically relevant insights.
    * Reduces the burden of multiple testing compared to individual GO term analysis.

### Advanced Methods for Deeper Exploration

- [x] **Cluster enrichment analysis:** Tools like CeaGO group related GO terms based on
semantic similarity and analyze enrichment within these clusters. This approach can reveal
broader functional themes beyond individual terms.
- [x] **Network analysis:** Integrating protein-protein interaction data with GO
annotations allows identifying functionally connected subnetworks enriched in DEGs. This
provides a network-based understanding of the underlying biological processes.

### Choosing the right method

The choice of method depends on factors like:

- [x] **Size of the DEG list:** For smaller lists, IGA might be sufficient, while larger lists
benefit from GSA approaches.
- [x] **Research question:** If interested in specific GO terms,
IGA might be suitable. For broader functional insights, GSA is preferred.

### Additional considerations

- [x] **Over-detection bias** standard methods give biased results on RNA-seq data due to
over-detection of differential expression for long and highly-expressed transcripts. The
==**goseq**== tool provides methods for performing GO analysis of RNA-seq data, taking length
bias into account. The methods and software used by goseq are equally applicable to other
category based tests of RNA-seq data, such as KEGG pathway analysis.
- [x] **Background gene list:** Choosing a relevant background list representing the genes not
differentially expressed is crucial for accurate enrichment analysis.
- [x] **Multiple testing correction:** Apply appropriate correction methods to account for
testing multiple GO terms simultaneously.
- [x] **Visualization:** Utilize graphical representations like bar
charts or heatmaps to visualize enriched GO terms and their significance levels.
---