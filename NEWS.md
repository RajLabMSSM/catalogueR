# catalogueR 0.1.1

## New features

* Added a `NEWS.md` file to track changes to the package.
* Split functions into separate files.
* Made README into .Rmd file.
* Changed license to GPL-3 and removed LICENSE file. 
* Removed old `data.table` version requirement
(issues with compressed files have since been resolved).
* Update .Rbuildignore and .gitignore
* Add GitHub Actions via 
`biocthis::use_bioc_github_action(pkgdown_covr_branch = c("main","master"))`.  
* Replace `T`/`F` with `FALSE`/`TRUE` in accordance with best coding practices.
* Replace NAMESPACE file with Roxygen-generated version. 
* Replace internal functions with *echoverse* packages:
  - `echotabix`
  - `echoconda`
* Styled code with: 
`styler::style_pkg(transformers = styler::tidyverse_style(indent_by = 4))` 
* Replace `printer` with `messager`.


## Bug fixes

* Reduce NAMESPACE by properly exporting functions. 
* Remove `library` calls within functions. 
* Add `requireNamespace` where suggests are used.



# catalogueR 0.1.0

First version of `catalogueR` created as a sister package of [`echolocatoR`](https://github.com/RajLabMSSM/echolocatoR). 