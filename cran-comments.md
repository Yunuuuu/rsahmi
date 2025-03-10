## Resubmission changes 
This resubmission addresses all the points raised by the reviewer:  

1. Expanded `DESCRIPTION` field
   - Added more details on package functionality and implemented methods.
   - Included relevant references in the required format.

2. Added `\value{}` tags  
   - Included `\value{}` sections in all relevant `.Rd` files (`blsd.Rd`, `extractor.Rd`, `prep_dataset.Rd`, `taxa_counts.Rd`).
   - The `install_polars` function has been removed, as all other functions now prompt users to install `polars` if needed.
   - Clearly described function outputs, including their structure and meaning.

3. Improved examples in `.Rd` files  
   - Added examples for exported functions.
   - Since all functions require large sequencing data, examples are wrapped in `\dontrun{}` to prevent execution during checks.

4. Ensured compliance with CRAN file writing policies  
   Removed any default file paths in functions that write data (including `extract_kraken_output()` and `prep_dataset()`)

5. Non-standard file/directory found at top level: 'rsahmi.Rmd'
   remote the the file

Thank you for your time and consideration.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
