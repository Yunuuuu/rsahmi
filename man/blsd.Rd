% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blsd.R
\name{blsd}
\alias{blsd}
\title{Barcode level signal denoising}
\usage{
blsd(
  kmer,
  method = "spearman",
  ...,
  min_kmer_len = 3L,
  min_cells = 3L,
  p.adjust = "BH"
)
}
\arguments{
\item{kmer}{kmer file returned by \code{\link[=sckmer]{sckmer()}}.}

\item{method}{A character string indicating which correlation coefficient is
to be used for the test. One of "pearson", "kendall", or "spearman", can be
abbreviated.}

\item{...}{Other arguments passed to \link[stats:cor.test]{cor.test}.}

\item{min_kmer_len}{An integer, the minimal number of kmer to filter taxa.
SAHMI use \code{2}.}

\item{min_cells}{An integer, the minimal number of cell per taxid. SAHMI use
\code{4}.}

\item{p.adjust}{Pvalue correction method, a character string. Can be
abbreviated. Details see \link[stats:p.adjust]{p.adjust}.}
}
\value{
A polars \link[polars:DataFrame_class]{DataFrame}
}
\description{
True taxa are detected on multiple barcodes and with a proprotional number of
total and unique k-mer sequences across barcodes, measured as a significant
Spearman correlation between the number of total and unique k-mers across
barcodes. (\code{padj < 0.05})
}
\seealso{
\url{https://github.com/sjdlabgroup/SAHMI}
}
