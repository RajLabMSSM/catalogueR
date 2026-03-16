# catalogueR

``` r

library(catalogueR)
```

## Introduction

*eQTL Catalogue* includes a large number of standardised QTL datasets
(110 datasets from 20 studies as of 7/4/2020). It actually contains more
than just eQTL data. For each dataset, the following kinds of QTLs can
be queried:

- **gene expression QTL**: `quant_method="ge"` (*default*) or
  `quant_method="microarray"`, depending on the dataset. `catalogueR`
  will automatically select whichever option is available.
- **exon expression QTL**: *under construction* `quant_method="ge"`
- **transcript usage QTL**: *under construction* `quant_method="tx"`
- **promoter, splice junction and 3’ end usage QTL**: *under
  construction* `quant_method="txrev"`

## Metadata

You can view a metadata table for all current datasets:

``` r

data("meta")
print(
  paste("There are", formatC(nrow(meta), big.mark = ","),
        "datasets in eQTL Catalogue.")
)
```

    ## [1] "There are 758 datasets in eQTL Catalogue."

``` r

knitr::kable(meta[1:10, ])
```

| study_id | dataset_id | study_label | sample_group | tissue_id | tissue_label | condition_label | sample_size | quant_method | ftp_path | ftp_cs_path | ftp_lbf_path | Tissue_group | System |
|:---|:---|:---|:---|:---|:---|:---|---:|:---|:---|:---|:---|:---|:---|
| QTS000001 | QTD000001 | Alasoo_2018 | macrophage_naive | CL_0000235 | macrophage | naive | 84 | ge | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000001/QTD000001.all.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000001/QTD000001.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000001/QTD000001.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000002 | Alasoo_2018 | macrophage_naive | CL_0000235 | macrophage | naive | 84 | exon | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000002/QTD000002.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000002/QTD000002.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000002/QTD000002.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000003 | Alasoo_2018 | macrophage_naive | CL_0000235 | macrophage | naive | 84 | tx | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000003/QTD000003.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000003/QTD000003.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000003/QTD000003.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000004 | Alasoo_2018 | macrophage_naive | CL_0000235 | macrophage | naive | 84 | txrev | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000004/QTD000004.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000004/QTD000004.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000004/QTD000004.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000005 | Alasoo_2018 | macrophage_naive | CL_0000235 | macrophage | naive | 84 | leafcutter | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000005/QTD000005.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000005/QTD000005.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000005/QTD000005.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000006 | Alasoo_2018 | macrophage_IFNg | CL_0000235 | macrophage | IFNg_18h | 84 | ge | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000006/QTD000006.all.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000006/QTD000006.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000006/QTD000006.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000007 | Alasoo_2018 | macrophage_IFNg | CL_0000235 | macrophage | IFNg_18h | 84 | exon | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000007/QTD000007.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000007/QTD000007.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000007/QTD000007.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000008 | Alasoo_2018 | macrophage_IFNg | CL_0000235 | macrophage | IFNg_18h | 84 | tx | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000008/QTD000008.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000008/QTD000008.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000008/QTD000008.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000009 | Alasoo_2018 | macrophage_IFNg | CL_0000235 | macrophage | IFNg_18h | 84 | txrev | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000009/QTD000009.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000009/QTD000009.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000009/QTD000009.lbf_variable.txt.gz> | macrophage | Blood |
| QTS000001 | QTD000010 | Alasoo_2018 | macrophage_IFNg | CL_0000235 | macrophage | IFNg_18h | 84 | leafcutter | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000010/QTD000010.cc.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000010/QTD000010.credible_sets.tsv.gz> | <ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000010/QTD000010.lbf_variable.txt.gz> | macrophage | Blood |

You can search through the metadata for datasets with certain keywords
(case-insensitive substrings across multiple columns).

``` r

qtl_datasets <- catalogueR::eQTLcatalogue_search_metadata(
  qtl_search = c("Alasoo_2018", "monocyte"))

print(qtl_datasets)
```

    ##  [1] "QTD000001" "QTD000002" "QTD000003" "QTD000004" "QTD000005" "QTD000006"
    ##  [7] "QTD000007" "QTD000008" "QTD000009" "QTD000010" "QTD000011" "QTD000012"
    ## [13] "QTD000013" "QTD000014" "QTD000015" "QTD000016" "QTD000017" "QTD000018"
    ## [19] "QTD000019" "QTD000020" "QTD000021" "QTD000022" "QTD000023" "QTD000024"
    ## [25] "QTD000025" "QTD000069" "QTD000081" "QTD000082" "QTD000083" "QTD000084"
    ## [31] "QTD000409" "QTD000410" "QTD000411" "QTD000412" "QTD000413" "QTD000414"
    ## [37] "QTD000415" "QTD000416" "QTD000417" "QTD000418" "QTD000419" "QTD000420"
    ## [43] "QTD000421" "QTD000422" "QTD000423" "QTD000424" "QTD000425" "QTD000426"
    ## [49] "QTD000427" "QTD000428" "QTD000429" "QTD000430" "QTD000431" "QTD000432"
    ## [55] "QTD000433" "QTD000499" "QTD000500" "QTD000501" "QTD000502" "QTD000503"
    ## [61] "QTD000504" "QTD000505" "QTD000506" "QTD000507" "QTD000508" "QTD000594"
    ## [67] "QTD000595" "QTD000596" "QTD000752" "QTD000753" "QTD000754" "QTD000755"
    ## [73] "QTD000756" "QTD000757" "QTD000758"

## Approach 1: Query with summary stats

Supply one or more paths to GWAS summary stats files (one per locus) and
automatically download any eQTL data within that range.

The files can be any of these formats, either gzip-compressed (`.gz`) or
uncompressed: `.csv`, `.tsv`, `space-separated`.

The summary stats files must have the following column names (order does
not matter):

- `SNP` (rsid for each SNP)
- `CHR` (chromosome; with or without the “chr” prefix)
- `POS` (basepair position)
- … (optional extra columns)

``` r

sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
gwas.qtl_paths <- catalogueR::eQTLcatalogue_query(
  sumstats_paths = sumstats_paths["BST1"],
  qtl_search = c("Alasoo_2018.macrophage_naive"),
  split_files = TRUE)
```

Because you selected the argument `split_files=TRUE`, the query results
have been distributed across multiple files and saved to disk. This is
useful when you do not want to load one massive `data.table` every time
you want to look at specific result subsets. In this scenario,
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md)
returns a list of those file paths (instead of the `data.table` itself).

To merge these files back together and import them into R, use:

``` r

GWAS.QTL <- catalogueR::merge_files(file_paths = gwas.qtl_paths)
```

## Approach 2: Query with coordinates

You can also make queries to *eQTL Catalogue* by manually specifying the
coordinates of the region you want to extract, as well as the
`unique_label` of the QTL dataset (see `data("meta")` for IDs).

``` r

GWAS.QTL_manual <- catalogueR::eQTLcatalogue_fetch(
  unique_label = "Alasoo_2018.macrophage_IFNg",
  chrom = 8,
  bp_lower = 21527069 - 500,
  bp_upper = 21527069 + 500)
```

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
