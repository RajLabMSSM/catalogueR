<img src='https://github.com/RajLabMSSM/catalogueR/raw/dev/inst/hex/hex.png' height='300'><br><br>
[![](https://img.shields.io/badge/devel%20version-0.1.1-black.svg)](https://github.com/RajLabMSSM/catalogueR)
[![R build
status](https://github.com/RajLabMSSM/catalogueR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RajLabMSSM/catalogueR/actions)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/catalogueR.svg)](https://github.com/RajLabMSSM/catalogueR/commits/master)
[![](https://codecov.io/gh/RajLabMSSM/catalogueR/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/catalogueR)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
README updated: <i>Dec-03-2021</i>
</h5>

## `catalogueR`: Rapid querying, colocalization, and plotting of summary stats from the [*eQTL Catalogue*](https://www.ebi.ac.uk/eqtl/). *eQTL Catalogue* currently contains &gt;100 datasets from 20 different studies (including GTEx), across many tissues/cell types/conditions.

This R package is part of the *echoverse* suite that supports
[`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR): an automated
genomic fine-mapping pipeline.

If you use `catalogueR`, please cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021) echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline, *Bioinformatics*; btab658,
> <https://doi.org/10.1093/bioinformatics/btab658>

and

> N Kerimov, JD Hayhurst, K Peikova, et al. A compendium of uniformly
> processed human gene expression and splicing quantitative trait loci.
> Nat Genet 53, 1290–1299 (2021).
> <https://doi.org/10.1038/s41588-021-00924-w>

## Installation

### \[1\] *Optional step: conda*

To ensure all dependencies are installed and don’t conflict with each
other, you can create a *conda* environment by downloading [this *yaml*
file](https://github.com/RajLabMSSM/catalogueR/blob/master/inst/conda/catalogueR.yml)
and entering the following in command line:

    conda env create -n <path_to_yml>/catalogueR.yml

### \[2\] Install `catalogueR` in R

#### Tabix

If you don’t use the *conda* env, you will also need to make sure
[tabix](http://www.htslib.org/doc/tabix.html) is installed.

**NOTE**: Different distributions of *tabix* may work with varying
degrees of consistency. I recommend you try them in the following order
of preference:

1.  [brew](https://formulae.brew.sh/formula/htslib) (for Mac users).
2.  [source](http://www.htslib.org/download/)
3.  [bioconda](https://anaconda.org/bioconda/tabix)

-   As of Feb 2021, I noticed that the bioconda distribution of *tabix*
    will randomly fail at queries. It appears this hasn’t been updated
    in over 2 years, so I do NOT recommend you use it.

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("RajLabMSSM/catalogueR")
library(catalogueR)
```

## Documentation

### [Website](https://rajlabmssm.github.io/catalogueR)

### [Getting started](https://rajlabmssm.github.io/catalogueR/articles/catalogueR)

## Notes

-   The functions in `catalogueR` are partly derived from the [*eQTL
    Catalogue*
    tutorial](http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html).  

-   [The ALT allele is always the effect allele in eQTL
    Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/).

-   Additional *eQTL Catalogue* resources:

    -   [GitHub
        repository](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)  
    -   [In-depth API
        documentation](https://www.ebi.ac.uk/eqtl/api-docs/)
    -   FTP server: *<ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv>*

<hr>

## Contact

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian
M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department
of Neuroscience, Icahn School of Medicine at Mount Sinai</a>
