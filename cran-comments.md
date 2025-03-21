## Resubmission changes

First and foremost, I truly appreciate the time and effort you have taken to
review my package, and I am grateful for your valuable feedback. I sincerely
apologize for not being able to address your concerns due to the following
reasons:

1. Example Execution (\dontrun{} → \donttest{})

This package implements the SAHMI algorithm using polars and includes a dataset
(4.9M), which makes it impractical to add additional datasets. The algorithm
requires at least three FASTQ files, and the absence of such data makes it
difficult to provide meaningful runnable examples.

Additionally, the package depends on external command-line tools such as kraken2
and seqkit, as well as packages from `Additional_repositories`. Since CRAN does
not support install dependencies from `Additional_repositories`, it is
impossible to provide fully executable examples. 

Given these constraints, we believe the use of `\dontrun{}` is appropriate to
prevent execution failures due to missing data and dependencies.

2. Small Executable Examples
The package follows the workflow of the original algorithm, where each step
corresponds to a function. While we acknowledge CRAN’s preference for runnable
examples, the lack of small, self-contained input data and external dependencies
makes this infeasible.

3. Avoiding `install.packages()` in Functions, Examples, and Vignettes The
`install.packages()` calls in: `R/import-standalone-pkg.R`
`R/import-standalone-polars.R`

They were included to ensure a user-friendly experience, as `polars` is not yet
available on `CRAN`. Since there are no runnable examples, we believe this
approach remains acceptable, as it simplifies installation for users unfamiliar
with setting up `polars`. We hope `polars` will be available on CRAN in the
future, as it is widely used by many users.

Once again, I sincerely appreciate your time, patience, and valuable feedback. I
deeply regret any inconvenience caused and will continue working to align the
package with CRAN’s best practices.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
