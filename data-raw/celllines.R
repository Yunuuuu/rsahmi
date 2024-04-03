## code to prepare `celllines` dataset goes here

library(polars)
# Download file from dropbox
# https://www.dropbox.com/s/r6xvw1589lqyqts/cell.lines.txt?dl=0
cell_lines <- pl$DataFrame(
    read.delim("cell.lines.txt", header = TRUE, sep = " ")
)
cell_lines <- cell_lines$select(
    # study = pl$lit("cell lines"), pl$col("sample"),
    pl$col("taxid"), 
    # pl$col("name")$alias("taxa"), 
    # pl$col("rank"),
    # pl$col("reads")$alias("total_reads"),
    # pl$col("min")$alias("minimizer_len"),
    # pl$col("uniq")$alias("minimizer_n_unique"),
    # pl$col("rpm"),
    pl$col("rpmm")
)
cell_lines$write_parquet(
    "inst/extdata/cell_lines.parquet",
    compression_level = 22L
)

file.remove("cell.lines.txt")
