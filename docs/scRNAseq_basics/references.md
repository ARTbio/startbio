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

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.4 LTS
    ##
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/rstudio-server_conda/conda/envs/rstudio-server_4.1.0/lib/libopenblasp-r0.3.15.so
    ##
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
    ##
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets
    ## [8] methods   base
    ##
    ## other attached packages:
    ##  [1] vroom_1.5.7           msigdbr_7.5.1         rmarkdown_2.14
    ##  [4] knitr_1.39            enrichplot_1.12.3     org.Hs.eg.db_3.13.0
    ##  [7] AnnotationDbi_1.54.1  IRanges_2.26.0        S4Vectors_0.30.2
    ## [10] Biobase_2.52.0        BiocGenerics_0.38.0   clusterProfiler_4.0.5
    ## [13] clustree_0.4.4        magrittr_2.0.3        dplyr_1.0.9
    ## [16] plyr_1.8.7            biomaRt_2.48.3        RColorBrewer_1.1-3
    ## [19] gridExtra_2.3         ggraph_2.0.5          ggplot2_3.3.6
    ## [22] sp_1.4-7              SeuratObject_4.1.0    Seurat_4.1.1
    ##
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2             reticulate_1.22        tidyselect_1.1.2
    ##   [4] RSQLite_2.2.14         htmlwidgets_1.5.4      grid_4.1.0
    ##   [7] BiocParallel_1.26.2    Rtsne_0.16             scatterpie_0.1.7
    ##  [10] munsell_0.5.0          codetools_0.2-18       ica_1.0-2
    ##  [13] future_1.25.0          miniUI_0.1.1.1         withr_2.5.0
    ##  [16] spatstat.random_2.2-0  colorspace_2.0-3       GOSemSim_2.18.1
    ##  [19] progressr_0.10.1       filelock_1.0.2         highr_0.9
    ##  [22] rstudioapi_0.13        ROCR_1.0-11            tensor_1.5
    ##  [25] DOSE_3.18.3            listenv_0.8.0          labeling_0.4.2
    ##  [28] GenomeInfoDbData_1.2.6 polyclip_1.10-0        bit64_4.0.5
    ##  [31] farver_2.1.0           downloader_0.4         treeio_1.16.2
    ##  [34] parallelly_1.32.0      vctrs_0.4.1            generics_0.1.2
    ##  [37] xfun_0.31              BiocFileCache_2.0.0    R6_2.5.1
    ##  [40] GenomeInfoDb_1.28.4    graphlayouts_0.8.0     gridGraphics_0.5-1
    ##  [43] bitops_1.0-7           spatstat.utils_2.3-1   cachem_1.0.6
    ##  [46] fgsea_1.18.0           assertthat_0.2.1       promises_1.2.0.1
    ##  [49] scales_1.2.0           rgeos_0.5-9            gtable_0.3.0
    ##  [52] globals_0.15.0         goftest_1.2-3          tidygraph_1.2.1
    ##  [55] rlang_1.0.2            splines_4.1.0          lazyeval_0.2.2
    ##  [58] checkmate_2.1.0        spatstat.geom_2.4-0    yaml_2.3.5
    ##  [61] reshape2_1.4.4         abind_1.4-5            backports_1.4.1
    ##  [64] httpuv_1.6.5           qvalue_2.24.0          tools_4.1.0
    ##  [67] ggplotify_0.1.0        ellipsis_0.3.2         spatstat.core_2.4-4
    ##  [70] ggridges_0.5.3         Rcpp_1.0.8.3           progress_1.2.2
    ##  [73] zlibbioc_1.38.0        purrr_0.3.4            RCurl_1.98-1.6
    ##  [76] prettyunits_1.1.1      rpart_4.1.16           deldir_1.0-6
    ##  [79] pbapply_1.5-0          viridis_0.6.2          cowplot_1.1.1
    ##  [82] zoo_1.8-10             ggrepel_0.9.1          cluster_2.1.3
    ##  [85] RSpectra_0.16-1        data.table_1.14.2      scattermore_0.8
    ##  [88] DO.db_2.9              lmtest_0.9-40          RANN_2.6.1
    ##  [91] fitdistrplus_1.1-8     matrixStats_0.62.0     hms_1.1.1
    ##  [94] patchwork_1.1.1        mime_0.12              evaluate_0.15
    ##  [97] xtable_1.8-4           XML_3.99-0.9           compiler_4.1.0
    ## [100] tibble_3.1.7           shadowtext_0.1.2       KernSmooth_2.23-20
    ## [103] crayon_1.5.1           htmltools_0.5.2        tzdb_0.3.0
    ## [106] ggfun_0.0.6            mgcv_1.8-40            later_1.3.0
    ## [109] aplot_0.1.4            tidyr_1.2.0            DBI_1.1.2
    ## [112] tweenr_1.0.2           dbplyr_2.1.1           MASS_7.3-57
    ## [115] rappdirs_0.3.3         babelgene_22.3         Matrix_1.4-1
    ## [118] cli_3.3.0              igraph_1.3.1           pkgconfig_2.0.3
    ## [121] plotly_4.10.0          spatstat.sparse_2.1-1  xml2_1.3.3
    ## [124] ggtree_3.0.4           XVector_0.32.0         yulab.utils_0.0.4
    ## [127] stringr_1.4.0          digest_0.6.29          sctransform_0.3.3
    ## [130] RcppAnnoy_0.0.19       spatstat.data_2.2-0    Biostrings_2.60.2
    ## [133] leiden_0.4.2           fastmatch_1.1-3        tidytree_0.3.9
    ## [136] uwot_0.1.11            curl_4.3.2             shiny_1.7.1
    ## [139] lifecycle_1.0.1        nlme_3.1-157           jsonlite_1.8.0
    ## [142] limma_3.48.3           viridisLite_0.4.0      fansi_1.0.3
    ## [145] pillar_1.7.0           lattice_0.20-45        KEGGREST_1.32.0
    ## [148] fastmap_1.1.0          httr_1.4.3             survival_3.3-1
    ## [151] GO.db_3.13.0           glue_1.6.2             png_0.1-7
    ## [154] bit_4.0.4              ggforce_0.3.3          stringi_1.7.6
    ## [157] blob_1.2.3             memoise_2.0.1          ape_5.6-2
    ## [160] irlba_2.3.5            future.apply_1.9.0
