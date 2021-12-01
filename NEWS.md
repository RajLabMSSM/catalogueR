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
* Replaced `T`/`F` with `FALSE`/`TRUE` in accordance with best coding practices.
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


## Bug fixes

* Reduced NAMESPACE by properly exporting functions. 
* Removed `library` calls within functions. 
* Added `requireNamespace` where suggests are used.
* Changed `eQTL_Catalogue.query` so that the default is `use_tabix=FALSE` due to 
instability of using `tabix` with the EBI server. See [here](https://github.com/RajLabMSSM/catalogueR/issues/5) for details.



# catalogueR 0.1.0

First version of `catalogueR` created as a sister package of [`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR). 