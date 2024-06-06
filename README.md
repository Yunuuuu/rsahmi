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
