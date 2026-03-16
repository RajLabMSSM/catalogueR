# Rapid querying, colocalization, and plotting of summary stats from the eQTL Catalogue.

The functions in **catalogueR** are partly derived from the following
[eQTL Catalogue
tutorial](http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.md).
Additional eQTL Catalogue Resources:

- [GitHub](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)

- [In-depth API documentation](https://www.ebi.ac.uk/eqtl/api-docs/)

- [eQTL Catalogue data portal](https://www.ebi.ac.uk/eqtl/)

## Details

**Notes on parallelization**: There's multiple levels to parallelize on.
You can only choose one level at a time:

- multithread_qtl=TRUE:

  Across QTL datasets.

- multithread_loci=TRUE:

  Across loci.

- multithread_tabix=TRUE:

  Within tabix files.

You can also get a speedup by using tabix instead of the rest API Test:
For 3 loci, and X QTL datasets:

- RESTful API:

  7.5 minutes.

- Tabix:

  27 seconds (\*clear winner! ~17x speedup). That said, if you're only
  query a small number of specific SNPs (rather than a large range), the
  RESTful API can sometimes be faster.

## See also

Useful links:

- <https://github.com/RajLabMSSM/catalogueR>

- Report bugs at <https://github.com/RajLabMSSM/catalogueR/issues>

## Author

**Maintainer**: Brian Schilder <brian_schilder@alumni.brown.edu>
([ORCID](https://orcid.org/0000-0001-5949-2191))

Authors:

- Jack Humphrey <Jack.Humphrey@mssm.edu>
  ([ORCID](https://orcid.org/0000-0002-6274-6620))

- Towfique Raj <towfique.raj@mssm.edu>
  ([ORCID](https://orcid.org/0000-0002-9355-5704))
