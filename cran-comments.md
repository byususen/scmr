## Test environments

* Local: native R 4.3.3 on Ubuntu 24.04.3, x86_64 Linux.
* `R CMD check --no-manual scmr_0.3.0.tar.gz`: Status OK.
* Tests: 136 passed, 0 failed, 0 warnings, 0 skipped.
* GitHub Actions: workflow retained; not executed during local validation.

## Notes

The package contains synthetic-data simulation utilities. The full article
study is opt-in; package checking runs small integration cases. See
`dev/validation.md` in the repository for details. Other platforms and the PDF
manual build have not been checked locally. This is a development refactor,
not a claim of CRAN submission or acceptance.
