Single-cell Analysis of Host-Microbiome Interactions
================

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
leverages the `polars` package for fast and systematic microbial signal
recovery and denoising from host tissue genomic sequencing.

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

You must also install [seqkit](https://bioinf.shenwei.me/seqkit/) and
[kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual).

Additionally, `rsahmi` relies on
[polars](https://rpolars.github.io/index.html). If you call any
functions that require `polars`, you will be prompted to install the
package.

## Workflow

``` r
library(ggplot2)
```

Prepare the output directory and input data

``` r
odir <- "rsahmi"
# For 10x Genomic data, `fq1` only contain barcode and umi, but the official
# didn't give any information for this. In this way, I prefer using `umi-tools`
# to transform the `umi` into fq2 and then run `rsahmi` with only fq2.
fq1 <-
    fq2 <- # can be `NULL`
    kraken_db <- # specify the kraken database
    if (dir.exists(odir)) dir.create(odir, recursive = TRUE)
```

### Taxonomic sequence classifier

We must run `kraken2` with `--use-names`, and `--report-minimizer-data`,
`classified_out` is not necessary in `rsahmi` package, but should be
specified if you want to use official `SAHMI` script.

``` r
blit::kraken2(
    fq1 = fq1,
    fq2 = fq2,
    classified_out = "classified.fq",
    blit::arg("--threads", 10L, format = "%d"),
    blit::arg("--db", kraken_db),
    "--use-names", "--report-minimizer-data",
    # Note: you should change `sample` into your sample name
    odir = file.path(odir, sample)
) |> blit::cmd_run()
```

### Extract microbiome kraken output

``` r
# Note: you should change `sample` into your sample name
#
# this'll extract all taxids specified in `taxon`
# which default to: c("d__Bacteria", "d__Fungi", "d__Viruses")
taxids <- rsahmi::extract_taxids(file.path(odir, sample, "kraken_report.txt"))
# then we extract the kraken2 output for these taxon.
rsahmi::extract_kraken_output(
    file.path(odir, sample, "kraken_output.txt"),
    taxids = taxids,
    dir = file.path(odir, sample)
)
```

### Extract microbiome kraken reads

``` r
# Note: you should change `sample` into your sample name
rsahmi::extract_kraken_reads(
    kraken_out = file.path(odir, sample, "kraken_microbiome_output.txt"),
    reads = c(fq1, fq2),
    odir = file.path(odir, sample),
    threads = 10L,
    # try to change `seqkit` argument into your seqkit path
    # If `NULL`, the internal will detect it in your `PATH` environment variable
    seqkit = NULL
)
```

**Following workflow will need all samples (multiple), you should run
above command for your all samples, then proceed to the subsequent
steps.**

### Prepare datasets for signal denoising and taxa counting

``` r
# paths: the output directory of all samples (from above steps)
paths <- fs::dir_ls(odir, type = "directory")
# this is the time-consuming step, so the internal will save the results in `odir`
sahmi_datasets <- lapply(paths, function(dir) {
    rsahmi::prep_dataset(
        fa1 = file.path(dir, "kraken_microbiome_reads.fa"),
        # if you have paired sequence, please also specify `fa2`,
        # !!! Also pay attention to the file name of `fa1` (add suffix `_1`) if
        #     you use paired reads.
        # fa2 = file.path(dir, "kraken_microbiome_reads_2.fa"),
        kraken_report = file.path(dir, "kraken_report.txt"),
        kraken_out = file.path(dir, "kraken_microbiome_output.txt"),
        odir = dir
    )
})
# If you have run above command, you can use following code
# to read the datasets
# sahmi_datasets <- lapply(paths, rsahmi::read_dataset)
```

### Identify contaminants

``` r
library(polars) # rsahmi dependents on `polars`
truly_microbe <- rsahmi::remove_contaminants(
    file.path(paths, "kraken_report.txt"),
    quantile = 0.99, exclusive = FALSE
)
microbe_for_plot <- attr(truly_microbe, "truly")[
    order(attr(truly_microbe, "pvalue")[attr(truly_microbe, "truly")])
]
microbe_for_plot <- microbe_for_plot[
    !microbe_for_plot %in% attr(truly_microbe, "exclusive")
]

truly_microbe_density <- ggplot(
    truly_microbe$filter(pl$col("taxid")$is_in(microbe_for_plot))$
        to_data_frame(),
    aes(rpmm),
) +
    geom_density(aes(fill = study), alpha = 0.5) +
    scale_x_log10() +
    facet_wrap(facets = vars(taxa), scales = "free") +
    theme(
        strip.clip = "off",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification.inside = c(1, 0)
    )
ggsave("truly_microbe.pdf",
    plot = truly_microbe_density,
    width = 8, height = 7L
)
```

### Sample level signal denoising

``` r
slsd <- rsahmi::slsd(lapply(sahmi_datasets, `[[`, "kreport"))
real_taxids_slsd <- slsd$filter(
    pl$col("r1")$gt(0),
    pl$col("r2")$gt(0),
    pl$col("r3")$gt(0),
    pl$col("p1")$lt(0.05),
    pl$col("p2")$lt(0.05),
    pl$col("p3")$lt(0.05)
)$get_column("taxid")
```

### Filter contaminants

``` r
umi_list <- lapply(sahmi_datasets, function(dataset) {
    # Barcode level signal denoising (barcode k-mer correlation test)
    blsd <- rsahmi::blsd(dataset$kmer)
    real_taxids <- blsd$filter(pl$col("padj")$lt(0.05))$get_column("taxid")
    # only keep taxids pass Sample level signal denoising
    real_taxids <- real_taxids$filter(real_taxids$is_in(real_taxids_slsd))
    # remove contaminants
    real_taxids <- real_taxids$filter(
        real_taxids$is_in(attr(truly_microbe, "truly"))
    )
    # filter UMI data
    dataset$umi$filter(pl$col("taxid")$is_in(real_taxids))
})
```

### Taxa umi counting

``` r
counts <- rsahmi::taxa_counts(umi_list, basename(names(umi_list)))
counts$write_csv(file.path(odir, "taxa_counts.txt"), separator = "\t")
```

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
#> other attached packages:
#> [1] ggplot2_3.5.2.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5        cli_3.6.3          knitr_1.49         rlang_1.1.4       
#>  [5] xfun_0.49          generics_0.1.3     glue_1.8.0         htmltools_0.5.8.1 
#>  [9] scales_1.4.0       fansi_1.0.6        rmarkdown_2.29     grid_4.4.2        
#> [13] evaluate_1.0.1     tibble_3.2.1       fastmap_1.2.0      yaml_2.3.10       
#> [17] lifecycle_1.0.4    compiler_4.4.2     dplyr_1.1.4        RColorBrewer_1.1-3
#> [21] pkgconfig_2.0.3    farver_2.1.2       digest_0.6.37      R6_2.5.1          
#> [25] tidyselect_1.2.1   utf8_1.2.4         pillar_1.9.0       magrittr_2.0.3    
#> [29] withr_3.0.2        tools_4.4.2        gtable_0.3.6
```
