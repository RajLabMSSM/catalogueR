# Colocalization

``` r

library(catalogueR)
```

## Introduction

Here we aim to robustly test whether the genetic signals underlying two
datasets (e.g. GWAS vs. eQTL in the same locus) are indeed the same,
using a methodology called colocalization.

Specifically, we use [`coloc`](https://github.com/chr1swallace/coloc),
which infers the probability that each SNP is causal in a given locus in
each of the datasets. It then tests the hypothesis that those signals
show substantially similar association distributions.

## List datasets

We have previously queried *eQTL Catalogue* using several Parkinson’s
disease GWAS loci (with
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md)).

Let’s first gather those files saved as TSV files.

``` r

gwas.qtl_paths <- catalogueR::eQTLcatalogue_example_queries(fnames = c(
  "BST1__Alasoo_2018.macrophage_IFNg+Salmonella.tsv.gz",
  "BST1__Alasoo_2018.macrophage_naive.tsv.gz",
  "BST1__Alasoo_2018.macrophage_Salmonella.tsv.gz"
))
```

## Run coloc

Since its original release,
[`coloc`](https://github.com/chr1swallace/coloc) has been updated so
that it can now model multiple causal variants within a given dataset
(see the
[paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008720)),
whereas previously it could assume only one causal variant. Thus, it may
be able to better estimate the colocalization probability between two
datasets.

``` r

coloc_QTLs <- catalogueR::COLOC_run(
  gwas.qtl_paths = gwas.qtl_paths,
  top_snp_only = TRUE,
  split_by_group = FALSE,
  method = "abf")
```

## Plot results

First, let’s plot only the results with \>80% colocalization
probability. This is a colocalization threshold commonly used in the
field.

``` r

coloc_plot <- catalogueR::COLOC_heatmap(
  coloc_QTLs = coloc_QTLs,
  coloc_thresh = .8)
```

If there are no colocalizations at \>=80% probability, you can lower the
threshold to explore the data (though not advisable for real analyses).

``` r

coloc_plot <- catalogueR::COLOC_heatmap(
  coloc_QTLs = coloc_QTLs,
  coloc_thresh = 0)
```

**Too many results:**

In other cases, you may have too many colocalized results to show in a
plot all at once. As a general rule, the more tests you run the higher
your probability threshold should be, since Bayesian methods do not
calculate P-values (and thus cannot use multiple-testing correction
methods like FDR). For example, if you ran many thousands of coloc tests
across many GWAS-eQTL pairs, it is advisable to raise your PP threshold
to 0.90 or even 0.99. This also has the benefit of narrowing down your
(sometimes many!) results so that you may focus on only the strongest
colocalizations.

## Session info

``` r

utils::sessionInfo()
```

    ## R Under development (unstable) (2026-03-12 r89607)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] catalogueR_2.0.1 BiocStyle_2.39.0
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
    ##   [3] magrittr_2.0.4              GenomicFeatures_1.63.1     
    ##   [5] farver_2.1.2                rmarkdown_2.30             
    ##   [7] fs_1.6.7                    BiocIO_1.21.0              
    ##   [9] ragg_1.5.1                  vctrs_0.7.1                
    ##  [11] memoise_2.0.1               Rsamtools_2.27.1           
    ##  [13] RCurl_1.98-1.17             base64enc_0.1-6            
    ##  [15] mixsqp_0.3-54               htmltools_0.5.9            
    ##  [17] S4Arrays_1.11.1             curl_7.0.0                 
    ##  [19] downloadR_1.0.0             SparseArray_1.11.11        
    ##  [21] sass_0.4.10                 bslib_0.10.0               
    ##  [23] htmlwidgets_1.6.4           basilisk_1.23.0            
    ##  [25] desc_1.4.3                  plyr_1.8.9                 
    ##  [27] cachem_1.1.0                GenomicAlignments_1.47.0   
    ##  [29] lifecycle_1.0.5             pkgconfig_2.0.3            
    ##  [31] Matrix_1.7-4                R6_2.6.1                   
    ##  [33] fastmap_1.2.0               MatrixGenerics_1.23.0      
    ##  [35] digest_0.6.39               reshape_0.8.10             
    ##  [37] patchwork_1.3.2             AnnotationDbi_1.73.0       
    ##  [39] S4Vectors_0.49.0            irlba_2.3.7                
    ##  [41] aws.signature_0.6.0         textshaping_1.0.5          
    ##  [43] GenomicRanges_1.63.1        RSQLite_2.4.6              
    ##  [45] filelock_1.0.3              httr_1.4.8                 
    ##  [47] abind_1.4-8                 compiler_4.6.0             
    ##  [49] bit64_4.6.0-1               S7_0.2.1                   
    ##  [51] BiocParallel_1.45.0         viridis_0.6.5              
    ##  [53] DBI_1.3.0                   R.utils_2.13.0             
    ##  [55] DelayedArray_0.37.0         rjson_0.2.23               
    ##  [57] piggyback_0.1.5             tools_4.6.0                
    ##  [59] otel_0.2.0                  zip_2.3.3                  
    ##  [61] R.oo_1.27.1                 glue_1.8.0                 
    ##  [63] restfulr_0.0.16             grid_4.6.0                 
    ##  [65] generics_0.1.4              BSgenome_1.79.1            
    ##  [67] gtable_0.3.6                tzdb_0.5.0                 
    ##  [69] susieR_0.14.2               R.methodsS3_1.8.2          
    ##  [71] tidyr_1.3.2                 data.table_1.18.2.1        
    ##  [73] hms_1.1.4                   xml2_1.5.2                 
    ##  [75] XVector_0.51.0              BiocGenerics_0.57.0        
    ##  [77] echoconda_1.0.0             pillar_1.11.1              
    ##  [79] stringr_1.6.0               dplyr_1.2.0                
    ##  [81] lattice_0.22-9              rtracklayer_1.71.3         
    ##  [83] bit_4.6.0                   tidyselect_1.2.1           
    ##  [85] Biostrings_2.79.5           coloc_5.2.3                
    ##  [87] knitr_1.51                  gridExtra_2.3              
    ##  [89] bookdown_0.46               IRanges_2.45.0             
    ##  [91] Seqinfo_1.1.0               SummarizedExperiment_1.41.1
    ##  [93] stats4_4.6.0                xfun_0.56                  
    ##  [95] Biobase_2.71.0              matrixStats_1.5.0          
    ##  [97] DT_0.34.0                   stringi_1.8.7              
    ##  [99] UCSC.utils_1.7.1            yaml_2.3.12                
    ## [101] echodata_1.0.0              evaluate_1.0.5             
    ## [103] codetools_0.2-20            cigarillo_1.1.0            
    ## [105] tibble_3.3.1                BiocManager_1.30.27        
    ## [107] cli_3.6.5                   echotabix_1.0.1            
    ## [109] reticulate_1.45.0           systemfonts_1.3.2          
    ## [111] jquerylib_0.1.4             Rcpp_1.1.1                 
    ## [113] GenomeInfoDb_1.47.2         dir.expiry_1.19.0          
    ## [115] png_0.1-8                   XML_3.99-0.22              
    ## [117] parallel_4.6.0              pkgdown_2.2.0              
    ## [119] ggplot2_4.0.2               readr_2.2.0                
    ## [121] blob_1.3.0                  basilisk.utils_1.23.1      
    ## [123] aws.s3_0.3.22               bitops_1.0-9               
    ## [125] viridisLite_0.4.3           VariantAnnotation_1.57.1   
    ## [127] scales_1.4.0                openxlsx_4.2.8.1           
    ## [129] purrr_1.2.1                 crayon_1.5.3               
    ## [131] rlang_1.1.7                 KEGGREST_1.51.1
