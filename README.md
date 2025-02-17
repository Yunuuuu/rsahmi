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

`rsahmi` depends on [polars](https://rpolars.github.io/index.html), you
can install it by:

``` r
if (Sys.info()["sysname"] == "Linux") {
    pak::pkg_install(
        "https://github.com/pola-rs/r-polars/releases/latest/download/polars__x86_64-pc-linux-gnu.gz"
    )
} else if (Sys.info()["sysname"] == "Darwin") {
    pak::pkg_install(
        "https://github.com/pola-rs/r-polars/releases/latest/download/polars__x86_64-apple-darwin20.tgz"
    )
} else {
    stop("Unsupported operating system")
}
```

You can also install the development version of `polars` from
[Github](https://github.com/pola-rs/r-polars) or
[R-universe](https://community.r-multiverse.org/polars):

In this way, you must install `rustc` to compile `polars`, you can
install `rustc` with `rustup`:

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Then install the `polars` package:

``` r
pak::pkg_install("pola-rs/r-polars")
```

Or

``` r
pak::repo_add("https://community.r-multiverse.org")
pak::pkg_install("polars")
```

Then you can install `rsahmi`:

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
