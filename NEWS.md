# catalogueR 0.1.1

## New features

* Added a `NEWS.md` file to track changes to the package.
* Split functions into separate files.
* Made README into .Rmd file.
* Changed license to GPL-3 and removed LICENSE file. 
* Removed old `data.table` version requirement
(issues with compressed files have since been resolved).
* Updated .Rbuildignore and .gitignore
* Added GitHub Actions via 
`biocthis::use_bioc_github_action(pkgdown_covr_branch = c("main","master"))`.  
* Replaced `T`/`F` with `TRUE`/`FALSE` in accordance with best coding practices.
* Replaced NAMESPACE file with Roxygen-generated version. 
* Replaced internal functions with *echoverse* packages:
  - `echotabix`
  - `echoconda` 
* Replaced built-in datasets with `echodata` datasets:
  - `BST1`
  - `LRRK2` 
  - `MEX3C`
* Styled code with: 
`styler::style_pkg(transformers = styler::tidyverse_style(indent_by = 4))` 
* Replace `printer` with `messager`. 
* Created [*dev* branch](https://github.com/RajLabMSSM/catalogueR/tree/dev). 
* Add GTEX_v8 metadata to `eQTL_Catalogue.list_datasets`. 
* Example coloc results now accessed via `get_coloc_QTLs`: 
  - Moved `coloc_QTLs_full` to Releases as it was >12Mb. 
  - Moved `coloc_QTLs` to Releases as it was >1.5Mb. 
* Moved eQTL Catalogue query results to Release via
`example_eQTL_Catalogue_query_paths`. 
* Replace liftover functions with `echotabix::liftover` 
(removed unused deps). 
* Replaced `rsids_from_coords` with `echodata::coords_to_rsids` 
(removed unused deps). 
* Made "REST" default `method` for all relevant functions.
* Deprecated:
  - `gather_files` in favor of more descriptive `merge_files`
  - `gather_colocalized_data` in favor of more descriptive
  `merge_colocalized_data`

## Bug fixes

* Reduced *NAMESPACE* by properly exporting functions. 
* Removed `library` calls within functions. 
* Added `requireNamespace` where suggests are used.
* Changed `eQTL_Catalogue.query` so that the default is `use_tabix=FALSE` due to 
instability of using `tabix` with the EBI server. See [here](https://github.com/RajLabMSSM/catalogueR/issues/5) for details. 
* `run_coloc` automatically infers (effective) sample size and checks that MAF
is present for all SNPs (after removing MAF>=1 or <=0 or ==NA).
* Removed example query results (e.g. `MEX3C__Alasoo_2018.macrophage_IFNg`) 
from built-in data as these were already included in *exdata*.


# catalogueR 0.1.0

First version of `catalogueR` created as a sister package of [`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR). 