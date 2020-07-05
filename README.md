# catalogueR  

Rapid querying, colocalization, and plotting of summary stats from the eQTL Catalogue.  

## Intro    

The functions in **catalogueR** are partly derived from the following
[eQTL Catalogue tutorial](http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html).  
Additional eQTL Catalogue Resources:  
- [GitHub](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)  
- [In-depth API documentation](https://www.ebi.ac.uk/eqtl/api-docs/)  
- [FTP server](ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv) 

*NOTE*: [The ALT allele is always the effect allele in eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/).  

<hr>  

## Getting started  


### Installation

```R
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
devtools::install_github("RajLabMSSM/catalogueR")
```


### Command line dependencies  

- [tabix](http://www.htslib.org/doc/tabix.html) (Install via [conda](https://anaconda.org/bioconda/tabix) or from [source](http://www.htslib.org/download/))  


### Intro example  

Supply one or more paths to [GWAS] summary stats files (one per locus) and automatically download any eQTL data within that range. The files can be any of these formats, either *gzip*-compressed (`.gz`) or uncompressed: `.csv`, `.tsv`, `space-separated`  

<br>

The summary stats files must have the following column names (order doesn't matter):
  - `SNP` (rsid for each SNP)
  - `CHR` (chromosome; with or without the "chr" prefix is fine)
  - `POS` (basepair position)
  - ... (optional extra columns)

```R 
sumstats_paths <- example_sumstats_paths()

gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths = sumstats_paths,  
                                       qtl_search = c("myeloid","Alasoo_2018"),
                                       output_dir = "./catalogueR_queries", 
                                       split_files = T,  
                                       merge_with_gwas = T,
                                       force_new_subset = T,
                                       nThread=4)
GWAS.QTL <- gather_files(file_paths = gwas.qtl_paths)
# Interactive datatable of results 
## WARNING: Don't use this function on large datatables, might cause freezing.
createDT(head(GWAS.QTL))
``` 
  
<hr>  
 
## Author  

Brian M. Schilder  
Raj Lab  
Icahn School of Medicine at Mount Sinai  
New York, New York, USA  
