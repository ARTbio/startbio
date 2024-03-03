# EGSEA

EGSEA, an acronym for Ensemble of Gene Set Enrichment Analyses, is a Bioconductor package
that utilizes the analysis results of eleven prominent GSE algorithms from the literature
to calculate collective significance scores for gene sets. These methods are currently:
ora, globaltest, plage, safe, zscore, gage, ssgsea, roast, fry, padog, camera, gsva. The
ora, gage, camera and gsva methods depend on a competitive null hypothesis while the
remaining seven methods are based on a self-contained hypothesis. EGSEA’s gene set
database, the EGSEAdata Bioconductor package, contains around 25,000 gene sets from 16
collections from MSigDB, KEGG and GeneSetDB. ==Supported organisms are human, mouse and rat,
however MSigDB is only available for human and mouse==.

An example
[EGSEA workflow](https://www.bioconductor.org/help/workflows/EGSEA123/) is available at
the Bioconductor workflows website.