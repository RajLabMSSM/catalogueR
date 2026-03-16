![](https://github.com/RajLabMSSM/catalogueR/raw/master/inst/hex/hex.png "Hex sticker for catalogueR")\
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btab658-blue.svg)](https://doi.org/10.1093/bioinformatics/btab658)
[![](https://img.shields.io/badge/devel%20version-2.0.1-black.svg)](https://github.com/RajLabMSSM/catalogueR)
[![](https://img.shields.io/github/languages/code-size/RajLabMSSM/catalogueR.svg)](https://github.com/RajLabMSSM/catalogueR)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/catalogueR.svg)](https://github.com/RajLabMSSM/catalogueR/commits/master)\
[![R build
status](https://github.com/RajLabMSSM/catalogueR/workflows/rworkflows/badge.svg)](https://github.com/RajLabMSSM/catalogueR/actions)
[![](https://codecov.io/gh/RajLabMSSM/catalogueR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/RajLabMSSM/catalogueR)\
[![](https://codecov.io/gh/RajLabMSSM/catalogueR/branch/master/graphs/icicle.svg "Codecov icicle graph")](https://app.codecov.io/gh/RajLabMSSM/catalogueR/tree/master)\

#### \
Authors: *Brian Schilder, Jack Humphrey, Towfique Raj*\

## `catalogueR`: Rapid querying, colocalization, and plotting of summary stats from the eQTL Catalogue <https://www.ebi.ac.uk/eqtl/>. eQTL Catalogue currently contains \>100 datasets from 20 different studies (including GTEx), across many tissues/cell types/conditions.

This R package is part of the *echoverse* suite that supports
[`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR): an automated
genomic fine-mapping pipeline.

If you use `catalogueR`, please cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021). echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline. Bioinformatics, btab658.
> <https://doi.org/10.1093/bioinformatics/btab658>

and

> N Kerimov, JD Hayhurst, K Peikova, et al. A compendium of uniformly
> processed human gene expression and splicing quantitative trait loci.
> Nat Genet 53, 1290–1299 (2021).
> <https://doi.org/10.1038/s41588-021-00924-w>

## Installation

``` r

if(!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("RajLabMSSM/catalogueR")
library(catalogueR)
```

## Documentation

### [Website](https://rajlabmssm.github.io/catalogueR)

### [Getting started](https://rajlabmssm.github.io/catalogueR/articles/catalogueR)

## Notes

- The functions in `catalogueR` are partly derived from the [*eQTL
  Catalogue*
  tutorial](http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.md).\

- [The ALT allele is always the effect allele in eQTL
  Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/).

- Additional *eQTL Catalogue* resources:

  - [GitHub
    repository](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)\
  - [In-depth API documentation](https://www.ebi.ac.uk/eqtl/api-docs/)
  - FTP server: *<ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv>*

------------------------------------------------------------------------

## Contact

[Brian M. Schilder, Bioinformatician
II](https://bschilder.github.io/BMSchilder/)\
[Raj Lab](https://rajlab.org)\
[Department of Neuroscience, Icahn School of Medicine at Mount
Sinai](https://icahn.mssm.edu/about/departments-offices/neuroscience)
