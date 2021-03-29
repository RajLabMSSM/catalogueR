# catalogueR  

Rapid querying, colocalization, and plotting of summary stats from the *eQTL Catalogue*. 

*eQTL Catalogue* currently contains 110 QTL datasets (full, genome-wide, standardized summary statistics and metadata) 
from 20 different studies (including GTEx V8), across many tissues/cell types/conditions (*updated: 7/5/20*).  

## Intro  

The functions in **catalogueR** are partly derived from the 
[*eQTL Catalogue* tutorial](http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html).  
Additional *eQTL Catalogue* resources:  
- [GitHub repository](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)  
- [In-depth API documentation](https://www.ebi.ac.uk/eqtl/api-docs/)  
- [Preprint](https://www.biorxiv.org/content/10.1101/2020.01.29.924266v1)  
- FTP server: ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv  

*NOTE*: [The ALT allele is always the effect allele in eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/).  

<hr>  

<br>  

## Documentation  

## [Website](https://rajlabmssm.github.io/catalogueR/)  
## [PDF](https://rajlabmssm.github.io/catalogueR/catalogueR_0.1.0.pdf)  

<br>

## Getting started  


### Installation 

### [1] *Optional step: conda*   

To ensure all dependendencies are installed and don't conflict with each other, 
you can create a *conda* environment by downloading [this *yaml*` file](https://github.com/RajLabMSSM/catalogueR/blob/master/inst/conda/catalogueR.yml) and 
entering the following in command line:  
```
conda env create -n <path_to_yml>/catalogueR.yml
```

### [2] Install **catalogueR** in R  

#### Tabix 

If you don't use the *conda* env, you will also need to make sure [tabix](http://www.htslib.org/doc/tabix.html) is installed.

**NOTE**: Different distributions of tabix may work with varying degrees of consistency. I recommend you try them in the following order of preference:

1. [brew](https://formulae.brew.sh/formula/htslib) (for Mac users).
2. [source](http://www.htslib.org/download/) 
3. [bioconda](https://anaconda.org/bioconda/tabix) 
  - As of Feb 2021, I noticed that the bioconda distribution of tabix will randomly fail at queries. It appears this hasn't been updated in over 2 years, so I do NOT recommend you use it.
  
#### data.table 

[*data.table*](https://github.com/Rdatatable/data.table) is a great R package for reading/writing/transforming tables far faster than base R methods. Unfortunately, in *data.table 1.14* they made some changes that cause errors reading/writing with compressed files (*.gz*) without additional setup steps on your computer. Therefore, I recommend you use versions 1.13 (the current CRAN binary distribution, NOT the source distribution, as of Feb 2021).  


```R
if(!"remotes" %in% installed.packages()){install.packages("remotes")}
remotes::install_github("RajLabMSSM/catalogueR")
```



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
[Raj Lab](http://rajlab.org)  
Icahn School of Medicine at Mount Sinai  
New York, New York, USA  
