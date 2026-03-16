# Changelog

## catalogueR 2.0.1

### Bug fixes

- Add `skip_on_cran()` and `skip_if_offline()` to `test-COLOC_run` to
  prevent hard failures on CI when piggyback downloads are unavailable.
- Fix hardcoded `$POS` in `test-eQTLcatalogue_fetch` — now detects
  position column dynamically (POS or BP).
- Fix `construct_locus_name()` to detect position column dynamically.
- Add `"aut"` role to maintainer in DESCRIPTION (required by
  `--as-cran`).

## catalogueR 2.0.0

### New features

- Remove 15 deprecated function wrappers that have been deprecated since
  v1.0.x.

### Bug fixes

- Fix `dataset_id` column detection in eQTL Catalogue queries.
- Fix data.table recycling warnings.
- Update R requirement from `>= 3.6` to `>= 4.1`.
- Remove unused imports: echoconda, stringr, utils, rlang.
- Rd: convert `\itemize` to `\describe` for named lists, `citEntry` to
  `bibentry`, bare doi URLs to `\doi{}`, title case.

## catalogueR 1.0.3

### Bug fixes

- Local R CMD check fixes and compatibility updates.

## catalogueR 1.0.2

### Bug fixes

- Update to v3 of *eQTL Catalogue* API
  (<https://github.com/eQTL-Catalogue/eQTL-SumStats>) which specifies
  chromosomes as an arg instead of part of the API path.
- Change the base URL from *http* to *https*.
- Update
  [`catalogueR::meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)
  object.\`
- Standardize `rworkflows.yml` with canonical template: enable Docker on
  `ghcr.io`, set `write-all` permissions, use `GITHUB_TOKEN`, add
  `devel` branch trigger.
- Update R requirement from `>= 3.6.0` to `>= 4.1`.
- Remove tracked `.DS_Store` files.

## catalogueR 1.0.1

### New features

- Implement `rworkflows`.

## catalogueR 1.0.0

### New features

- Offload functions to echoverse deps:
  - `get_data()` –\> `echodata`
  - `get_os()` –\> `echodata`
  - `load_rdata()` –\> `downloadR`
  - `example_sumstats_paths()` –\> `echpdata::get_Nalls2019_loci()`
  - `find_consensus_SNPs` –\>
    [`echodata::find_consensus_snps()`](https://rdrr.io/pkg/echodata/man/find_consensus_snps.html)
  - `locus_plot` –\> `echoplot`
- Rename functions to fit *echoverse* conventions:
  - Remove “.” from any function names.
  - All eQTL Catalogue functions start with `eQTLcatalogue_`.
  - All COLOC functions start with `COLOC_`.
  - Deprecated all functions properly to let users know.
- Change arg in all functions:
  - `gwas_data` –\> `query_granges`
  - `genome_build` –\> `query_genome`
- Condensed multiple arguments into
  [`echotabix::construct_query()`](https://rdrr.io/pkg/echotabix/man/construct_query.html)
  - `gwas.data`, `chrom`, `bp_min`, `bp_max`
- New exported functions:
  - `COLOC_heatmap`
  - `COLOC_corplot`

### Bug fixes

- Updated function/arg names from *echoverse* deps.
- Updated `meta` data to include latest paths and datasets.
- `fix_ftp`: Replacing “ftp:” prefix with “http:” seems to drastically
  improve the ability of functions to query the files.

## catalogueR 0.1.1

### New features

- Added a `NEWS.md` file to track changes to the package.
- Split functions into separate files.
- Made README into .Rmd file.
- Changed license to GPL-3 and removed LICENSE file.
- Removed old `data.table` version requirement (issues with compressed
  files have since been resolved).
- Updated .Rbuildignore and .gitignore
- Added GitHub Actions via
  `biocthis::use_bioc_github_action(pkgdown_covr_branch = c("main","master"))`.\
- Replaced `T`/`F` with `TRUE`/`FALSE` in accordance with best coding
  practices.
- Replaced NAMESPACE file with Roxygen-generated version.
- Replaced internal functions with *echoverse* packages:
  - `echotabix`
  - `echoconda`
- Replaced built-in datasets with `echodata` datasets:
  - `BST1`
  - `LRRK2`
  - `MEX3C`
- Styled code with:
  `styler::style_pkg(transformers = styler::tidyverse_style(indent_by = 4))`
- Replace `printer` with `messager`.
- Created [*dev*
  branch](https://github.com/RajLabMSSM/catalogueR/tree/dev).
- Add GTEX_v8 metadata to `eQTL_Catalogue.list_datasets`.
- Example coloc results now accessed via `get_coloc_QTLs`:
  - Moved `coloc_QTLs_full` to Releases as it was \>12Mb.
  - Moved `coloc_QTLs` to Releases as it was \>1.5Mb.
- Moved eQTL Catalogue query results to Release via
  `example_eQTL_Catalogue_query_paths`.
- Replace liftover functions with
  [`echotabix::liftover`](https://rdrr.io/pkg/echotabix/man/liftover.html)
  (removed unused deps).
- Replaced `rsids_from_coords` with
  [`echodata::coords_to_rsids`](https://rdrr.io/pkg/echodata/man/coords_to_rsids.html)
  (removed unused deps).
- Made “REST” default `method` for all relevant functions.
- Deprecated:
  - `gather_files` in favor of more descriptive `merge_files`
  - `gather_colocalized_data` in favor of more descriptive
    `merge_colocalized_data`

### Bug fixes

- Reduced *NAMESPACE* by properly exporting functions.
- Removed `library` calls within functions.
- Added `requireNamespace` where suggests are used.
- Changed `eQTL_Catalogue.query` so that the default is
  `use_tabix=FALSE` due to instability of using `tabix` with the EBI
  server. See [here](https://github.com/RajLabMSSM/catalogueR/issues/5)
  for details.
- `run_coloc` automatically infers (effective) sample size and checks
  that MAF is present for all SNPs (after removing MAF\>=1 or \<=0 or
  ==NA).
- Removed example query results
  (e.g. `MEX3C__Alasoo_2018.macrophage_IFNg`) from built-in data as
  these were already included in *exdata*.

## catalogueR 0.1.0

First version of `catalogueR` created as a sister package of
[`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR).
