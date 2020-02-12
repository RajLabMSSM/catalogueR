# catalogueR

## Intro    

The following extends and build upon the APIs provided by [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) (with which I am not affiliated):  
- [GitHub source code](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)  
- **FTP Server**: *ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv*  
- [In-depth API documentation](https://www.ebi.ac.uk/eqtl/api-docs/)  

A list of all current tabix-indexed QTL datasets is provided [here](https://github.com/RajLabMSSM/catalogueR/blob/master/resources/eQTLcatalogue_tabix_ftp_paths.tsv) (or [here]() for the original source).  


<hr>  

## Getting started  

### Clone this repo  
`cd <your_preferred_path>`  
`git clone https://github.com/RajLabMSSM/catalogueR.git`  
`cd catalogueR`

### Install required software

#### Command line 
- [tabix](http://www.htslib.org/doc/tabix.html) (Install via [conda](https://anaconda.org/bioconda/tabix) or from [source](http://www.htslib.org/download/))  

#### R  
- dplyr  
- ggplot2  
- readr  
- stringr  
- httr  
- jsonlite  
- tidyverse  
- coloc  
- biomaRt
- wiggleplotr  
- GenomicRanges  
- AnnotationDbi   
- data.table  
  
### Import *catalogueR* functions  

`source("./functions/catalogueR.R")`  

<hr>  

## Tutorial

## Import functions


### [Approach 1]  

Supply one or more paths to [GWAS] summary stats files (one per locus) and automatically download any eQTL data within that range. The files can be any of these formats, either *gzip*-compressed (`.gz`) or uncompressed: `.csv`, `.tsv`, `space-separated`  
<br>
The summary stats file must have the following column names (order doesn't matter). You can have as many additional columns as you want:  
  - `SNP` (rsid for each SNP)
  - `CHR` (chromosome)
  - `POS` (basepair position)

```
# Returns both the gwas_data you supplied and the queried QTL summary stats, 
## all merged into one data.table.  

gwas.qtl <- catalogueR.run(# Any number of summary stats files
                           sumstats_paths =
                           c("./example_data/Nalls23andMe_2019/BIN3_Nalls23andMe_2019_subset.tsv.gz",
                             "./example_data/Nalls23andMe_2019/BST1_Nalls23andMe_2019_subset.tsv.gz",
                             "./example_data/Nalls23andMe_2019/SNCA_Nalls23andMe_2019_subset.tsv.gz"),
                           
                           # Name of each locus (a vector of the same length as the sumstats_paths)  
                           # If none are supplied (NULL), names will be assigned based on the chromosomal coordinates.
                           loci_names=c("BIN3","BST1","SNCA"),
                           
                           # The file you want the merged gwas/qtl results to be saved in.
                           ## The file will be tab-delimited, but you can choose to leave it uncompressed by removing the '.gz' at the end of the path here.
                           output_path="./example_data/Nalls23andMe_2019/eQTL_Catalogue.tsv.gz",
                           
                           # This function will automatically search for any datasets that match your criterion.
                           ## For example, if you search "Alasoo_2018", it will query the datasets "Alasoo_2018.macrophage_naive", "Alasoo_2018.macrophage_IFNg",Alasoo_2018.macrophage_Salmonella", and "Alasoo_2018.macrophage_IFNg+Salmonella").
                           ## You can be more specific about which datasets you want to 
                           qtl_datasets=c("ROSMAP","Alasoo_2018"),
                           
                           # Tabix is about ~17x faster (default; =T) than the REST API (=F).
                           use_tabix=T,
                           
                           # Use multiple CPU cores to speed up your queries
                           nThread=4, 
                           
                           # Multi-thread across QTL datasets (good when you're querying lots of QTL datasets.)
                           multithread_qtl=T,
                           
                           # Multi-thread across loci (good when you have lots of gwas loci)
                           multithread_loci=F)
```


### [Approach 2]

Download a subset of QTL summary stats directly by specifying the coordinates you want to extract:  
 
```
# Returns a data.table with the QTL summary stats subset.  

gwas.qtl <- catalogueR.fetch(# Unique QTL id (<study>.<qtl_group>)
							 unique_id="Alasoo_2018.macrophage_IFNg",

                              # You can specify the QTL quantification method you want to use.
                              ## (options: "ge","exon", "tx","txrev","microarray")
                              ## NOTE: Not all methods are available for all datasets. 
                              ## So if for example you select "ge" for a microarray dataset,   
                              ## this function will default to a method that is available for that dataset (i.e. "microarray")
                              quant_method="ge",
                              
                              # Set to false since you're not using the gwas_data to infer coordinates.
                              infer_region=F, 
                              
                              # The number of threads to use when reading in the QTL data subset
                              nThread = 4,
                               
                             # Tabix is about ~17x faster (default; =T) than the REST API (=F).
                              use_tabix = T,
                              
                             # Chromosome to query
                             chrom = 8,
                             
                             # Minimum basepair position of the query window
                             bp_lower=21527069
                             
                             # Maximum basepair position of the query window
                             bp_upper=23525543)
```

<hr>

## Author  

Brian M. Schilder  
Raj Lab  
Icahn School of Medicine at Mount Sinai  
New York, New York, USA  
