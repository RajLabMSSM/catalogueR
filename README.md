<img src='https://github.com/RajLabMSSM/catalogueR/raw/dev/inst/hex/hex.png' height='300'><br><br>
[![](https://img.shields.io/badge/devel%20version-0.1.1-black.svg)](https://github.com/RajLabMSSM/catalogueR)
[![R build
status](https://github.com/RajLabMSSM/catalogueR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RajLabMSSM/catalogueR/actions)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/catalogueR.svg)](https://github.com/RajLabMSSM/catalogueR/commits/master)
[![](https://app.codecov.io/gh/RajLabMSSM/catalogueR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/RajLabMSSM/catalogueR)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
README updated: <i>Sep-12-2022</i>
</h5>

## `catalogueR`: Rapid querying, colocalization, and plotting of summary stats from the [*eQTL Catalogue*](https://www.ebi.ac.uk/eqtl/). *eQTL Catalogue* currently contains \>100 datasets from 20 different studies (including GTEx), across many tissues/cell types/conditions.

This R package is part of the *echoverse* suite that supports
[`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR): an automated
genomic fine-mapping pipeline.

If you use `catalogueR`, please cite:

and

> N Kerimov, JD Hayhurst, K Peikova, et al. A compendium of uniformly
> processed human gene expression and splicing quantitative trait loci.
> Nat Genet 53, 1290–1299 (2021).
> <https://doi.org/10.1038/s41588-021-00924-w>

## Installation

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
