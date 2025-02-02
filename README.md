Single-cell Analysis of Host-Microbiome Interactions
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Installation

Here, we use `pak` to install package, you can also use `remotes`

``` r
if (!requireNamespace("pak")) {
    install.packages("pak",
        repos = sprintf(
            "https://r-lib.github.io/p/pak/devel/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch
        )
    )
}
```

`rsami` depends on [polars](https://rpolars.github.io/index.html), you
can install it by:

``` r
pak::repo_add("https://rpolars.r-universe.dev")
pak::pkg_install("polars")
```

You can also install the development version of `polars` from
[Github](https://github.com/pola-rs/r-polars)

In this way, you must install `rustc` to compile `polars`

You can install `rustc` with `rustup`:

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Then install the `polars` package:

``` r
pak::pkg_install("pola-rs/r-polars")
```

You can install the development version of `rsahmi` from
[GitHub](https://github.com/Yunuuuu/rsahmi):

``` r
pak::pkg_install("Yunuuuu/rsahmi")
```

You must also install [seqkit](https://bioinf.shenwei.me/seqkit/) and
[kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual).

## Cite

If you use this package, please cite:

- Ghaddar B, Biswas A, Harris C, Omary MB, Carpizo DR, Blaser MJ, De S.
  Tumor microbiome links cellular programs and immunity in pancreatic
  cancer. Cancer Cell. 2022 Oct 10;40(10):1240-1253.e5. doi:
  10.1016/j.ccell.2022.09.009. PMID: 36220074; PMCID: PMC9556978.
- Ghaddar B, Blaser MJ, De S. Denoising sparse microbial signals from
  single-cell sequencing of mammalian host tissues. Nat Comput Sci. 2023
  Sep;3(9):741-747. doi: 10.1038/s43588-023-00507-1. Epub 2023 Sep 18.
  PMID: 37946872; PMCID: PMC10634611.
- Song, Yuxuan PhD; Peng, Yun PhD; Qin, Caipeng PhD; Jiang, Shan PhD;
  Lin, Jiaxing PhD; Lai, Shicong MD; Wu, Jilin PhD; Ding, Mengting PhD;
  Du, Yiqing PhD*; Yu, Luping MD*; Xu, Tao MD\*. Antibiotic use
  attenuates response to immune checkpoint blockade in urothelial
  carcinoma via inhibiting CD74-MIF/COPA: revealing cross-talk between
  anti-bacterial immunity and ant-itumor immunity. International Journal
  of Surgery 111(1):p 972-987, January 2025. \| DOI:
  10.1097/JS9.0000000000001901
