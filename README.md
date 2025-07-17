Single-cell Analysis of Host-Microbiome Interactions
================

> **⚠️ This package has been superseded by
> [`scmire`](https://github.com/WangLabCSU/scmire).**  
> While `rsahmi` remains available for archival purposes, future
> development and improvements will be focused on `scmire`, which offers
> an improved design for single-cell microbiome reconstruction and
> quantification.

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/rsahmi)](https://CRAN.R-project.org/package=rsahmi)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://cranlogs.r-pkg.org/badges/rsahmi)](https://cran.r-project.org/package=rsahmi)
[![DOI](https://zenodo.org/badge/605180683.svg)](https://doi.org/10.5281/zenodo.15076675)
<!-- badges: end -->

A computational resource designed to accurately detect microbial nucleic
acids while filtering out contaminants and false-positive taxonomic
assignments from standard transcriptomic sequencing of mammalian tissues
(See [SAHMI](https://github.com/sjdlabgroup/SAHMI)). This implementation
leverages the `rust` for fast and systematic microbial signal recovery
and denoising from host tissue genomic sequencing.

## Installation

You can install `rsahmi` from `CRAN` using:

``` r
# install.packages("pak")
pak::pak("rsahmi")
```

Alternatively, install the development version from
[r-universe](https://yunuuuu.r-universe.dev/rsahmi) with:

``` r
pak::repo_add("https://yunuuuu.r-universe.dev")
pak::pak("rsahmi")
```

or from [GitHub](https://github.com/Yunuuuu/rsahmi) with:

``` r
pak::pak("Yunuuuu/rsahmi")
```

You must also install
[kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual).

## Cite

If you use this package, please cite:

- Song, Yuxuan PhD; Peng, Yun PhD; Qin, Caipeng PhD; Jiang, Shan PhD;
  Lin, Jiaxing PhD; Lai, Shicong MD; Wu, Jilin PhD; Ding, Mengting PhD;
  Du, Yiqing PhD*; Yu, Luping MD*; Xu, Tao MD\*. Antibiotic use
  attenuates response to immune checkpoint blockade in urothelial
  carcinoma via inhibiting CD74-MIF/COPA: revealing cross-talk between
  anti-bacterial immunity and ant-itumor immunity. International Journal
  of Surgery 111(1):p 972-987, January 2025. \| DOI:
  10.1097/JS9.0000000000001901.
- Ghaddar B, Blaser MJ, De S. Denoising sparse microbial signals from
  single-cell sequencing of mammalian host tissues. Nat Comput Sci. 2023
  Sep;3(9):741-747. doi: 10.1038/s43588-023-00507-1. Epub 2023 Sep 18.
  PMID: 37946872; PMCID: PMC10634611.

## sessionInfo

``` r
sessionInfo()
#> R version 4.4.2 (2024-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.1 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libmkl_rt.so;  LAPACK version 3.8.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: Asia/Shanghai
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.4.2    fastmap_1.2.0     cli_3.6.5         tools_4.4.2      
#>  [5] htmltools_0.5.8.1 yaml_2.3.10       rmarkdown_2.29    knitr_1.50       
#>  [9] xfun_0.52         digest_0.6.37     rlang_1.1.6       evaluate_1.0.3
```
