# References

Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd,
Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of
Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21. doi:
10.1016/j.cell.2019.05.031. Epub 2019 Jun 6. PMID: 31178118; PMCID:
PMC6687398.

Neo Christopher Chung, John D. Storey, Statistical significance of
variables driving systematic variation in high-dimensional data,
Bioinformatics, Volume 31, Issue 4, 15 February 2015, Pages 545–554,
<https://doi.org/10.1093/bioinformatics/btu674>

Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, +7, Sayan Mukherjee,
Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L.
Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov. Gene set
enrichment analysis: A knowledge-based approach for interpreting
genome-wide expression profiles. PNAS, Volume 102, Number 43, Pages
15545-15550, (2005) <https://doi.org/10.1073/pnas.0506580102>

Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1α-responsive genes
involved in oxidative phosphorylation are coordinately downregulated in
human diabetes. Nat Genet 34, 267–273 (2003).
<https://doi.org/10.1038/ng1180>

Arthur Liberzon, Aravind Subramanian, Reid Pinchback, Helga
Thorvaldsdóttir, Pablo Tamayo, Jill P. Mesirov, Molecular signatures
database (MSigDB) 3.0, Bioinformatics, Volume 27, Issue 12, 15 June
2011, Pages 1739–1740, <https://doi.org/10.1093/bioinformatics/btr260>

The Gene Ontology Consortium, The Gene Ontology resource: enriching a
GOld mine, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021,
Pages D325–D334, <https://doi.org/10.1093/nar/gkaa1113>

Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the
unification of biology. Nat Genet 25, 25–29 (2000).
<https://doi.org/10.1038/75556>

Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and
Genomes, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000,
Pages 27–30, <https://doi.org/10.1093/nar/28.1.27>

Kanehisa, M. Toward understanding the origin and evolution of cellular
organisms. Protein Science. 2019; 28: 1947– 1951.
<https://doi.org/10.1002/pro.3715>

Minoru Kanehisa, Miho Furumichi, Yoko Sato, Mari Ishiguro-Watanabe, Mao
Tanabe, KEGG: integrating viruses and cellular organisms, Nucleic Acids
Research, Volume 49, Issue D1, 8 January 2021, Pages D545–D551,
<https://doi.org/10.1093/nar/gkaa970>

Xinxin Zhang, Yujia Lan, Jinyuan Xu, Fei Quan, Erjie Zhao, Chunyu Deng,
Tao Luo, Liwen Xu, Gaoming Liao, Min Yan, Yanyan Ping, Feng Li, Aiai
Shi, Jing Bai, Tingting Zhao, Xia Li, Yun Xiao, CellMarker: a manually
curated resource of cell markers in human and mouse, Nucleic Acids
Research, Volume 47, Issue D1, 08 January 2019, Pages D721–D728,
<https://doi.org/10.1093/nar/gky900>

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/rstudio-server_conda/conda/envs/rstudio-server_4.3.1/lib/libopenblasp-r0.3.24.so;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] presto_1.0.0           data.table_1.15.4      Rcpp_1.0.12           
    ##  [4] vroom_1.6.5            msigdbr_7.5.1          rmarkdown_2.26        
    ##  [7] knitr_1.46             enrichplot_1.22.0      org.Hs.eg.db_3.18.0   
    ## [10] AnnotationDbi_1.64.1   IRanges_2.36.0         S4Vectors_0.40.2      
    ## [13] Biobase_2.62.0         BiocGenerics_0.48.1    clusterProfiler_4.10.0
    ## [16] clustree_0.5.1         magrittr_2.0.3         dplyr_1.1.4           
    ## [19] plyr_1.8.9             biomaRt_2.58.0         RColorBrewer_1.1-3    
    ## [22] gridExtra_2.3          ggraph_2.1.0           ggplot2_3.5.1         
    ## [25] Seurat_5.0.3           SeuratObject_5.0.1     sp_2.1-3              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.4                matrixStats_1.3.0       spatstat.sparse_3.0-3  
    ##   [4] bitops_1.0-7            HDO.db_0.99.1           httr_1.4.7             
    ##   [7] backports_1.4.1         tools_4.3.1             sctransform_0.4.1      
    ##  [10] utf8_1.2.4              R6_2.5.1                lazyeval_0.2.2         
    ##  [13] uwot_0.2.2              withr_3.0.0             prettyunits_1.2.0      
    ##  [16] progressr_0.14.0        cli_3.6.2               spatstat.explore_3.2-7 
    ##  [19] fastDummies_1.7.3       scatterpie_0.2.2        labeling_0.4.3         
    ##  [22] spatstat.data_3.0-4     ggridges_0.5.6          pbapply_1.7-2          
    ##  [25] yulab.utils_0.1.4       gson_0.1.0              DOSE_3.28.2            
    ##  [28] R.utils_2.12.3          parallelly_1.37.1       limma_3.56.2           
    ##  [31] rstudioapi_0.16.0       RSQLite_2.3.6           generics_0.1.3         
    ##  [34] gridGraphics_0.5-1      ica_1.0-3               spatstat.random_3.2-3  
    ##  [37] GO.db_3.18.0            Matrix_1.6-5            ggbeeswarm_0.7.2       
    ##  [40] fansi_1.0.6             abind_1.4-5             R.methodsS3_1.8.2      
    ##  [43] lifecycle_1.0.4         yaml_2.3.8              qvalue_2.34.0          
    ##  [46] BiocFileCache_2.10.1    Rtsne_0.17              grid_4.3.1             
    ##  [49] blob_1.2.4              promises_1.3.0          crayon_1.5.2           
    ##  [52] miniUI_0.1.1.1          lattice_0.22-6          cowplot_1.1.3          
    ##  [55] KEGGREST_1.42.0         pillar_1.9.0            fgsea_1.28.0           
    ##  [58] future.apply_1.11.2     codetools_0.2-20        fastmatch_1.1-4        
    ##  [61] leiden_0.4.3.1          glue_1.7.0              ggfun_0.1.4            
    ##  [64] vctrs_0.6.5             png_0.1-8               treeio_1.26.0          
    ##  [67] spam_2.10-0             gtable_0.3.5            cachem_1.0.8           
    ##  [70] xfun_0.43               mime_0.12               tidygraph_1.2.3        
    ##  [73] survival_3.5-8          fitdistrplus_1.1-11     ROCR_1.0-11            
    ##  [76] nlme_3.1-164            ggtree_3.10.0           bit64_4.0.5            
    ##  [79] progress_1.2.3          filelock_1.0.3          RcppAnnoy_0.0.22       
    ##  [82] GenomeInfoDb_1.38.5     irlba_2.3.5.1           vipor_0.4.7            
    ##  [85] KernSmooth_2.23-22      colorspace_2.1-0        DBI_1.2.2              
    ##  [88] ggrastr_1.0.2           tidyselect_1.2.1        bit_4.0.5              
    ##  [91] compiler_4.3.1          curl_5.0.2              xml2_1.3.6             
    ##  [94] plotly_4.10.4           shadowtext_0.1.3        checkmate_2.3.0        
    ##  [97] scales_1.3.0            lmtest_0.9-40           rappdirs_0.3.3         
    ## [100] stringr_1.5.1           digest_0.6.35           goftest_1.2-3          
    ## [103] spatstat.utils_3.0-4    XVector_0.42.0          htmltools_0.5.8.1      
    ## [106] pkgconfig_2.0.3         highr_0.10              dbplyr_2.5.0           
    ## [109] fastmap_1.1.1           rlang_1.1.3             htmlwidgets_1.6.4      
    ## [112] shiny_1.8.1.1           farver_2.1.1            zoo_1.8-12             
    ## [115] jsonlite_1.8.8          BiocParallel_1.36.0     GOSemSim_2.28.1        
    ## [118] R.oo_1.26.0             RCurl_1.98-1.14         GenomeInfoDbData_1.2.11
    ## [121] ggplotify_0.1.2         dotCall64_1.1-1         patchwork_1.2.0        
    ## [124] munsell_0.5.1           ape_5.8                 babelgene_22.9         
    ## [127] viridis_0.6.5           reticulate_1.35.0       stringi_1.8.3          
    ## [130] zlibbioc_1.48.0         MASS_7.3-60.0.1         parallel_4.3.1         
    ## [133] listenv_0.9.1           ggrepel_0.9.5           deldir_2.0-4           
    ## [136] Biostrings_2.70.1       graphlayouts_1.0.1      splines_4.3.1          
    ## [139] tensor_1.5              hms_1.1.3               igraph_2.0.3           
    ## [142] spatstat.geom_3.2-9     RcppHNSW_0.6.0          reshape2_1.4.4         
    ## [145] XML_3.99-0.16.1         evaluate_0.23           tzdb_0.4.0             
    ## [148] tweenr_2.0.2            httpuv_1.6.15           RANN_2.6.1             
    ## [151] tidyr_1.3.1             purrr_1.0.2             polyclip_1.10-6        
    ## [154] future_1.33.2           scattermore_1.2         ggforce_0.4.1          
    ## [157] xtable_1.8-4            RSpectra_0.16-1         tidytree_0.4.6         
    ## [160] later_1.3.2             viridisLite_0.4.2       tibble_3.2.1           
    ## [163] aplot_0.2.2             memoise_2.0.1           beeswarm_0.4.0         
    ## [166] cluster_2.1.6           globals_0.16.3