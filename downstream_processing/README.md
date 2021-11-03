# Downstream Processing

## Differential abundance testing with limma

The R `limma` package is nicely suited to handling differential abundance (or in this case differential intensity) analysis.
Unlike similar DA packages like `DESeq2`, `limma` does not require that the data is provided as count data (i.e. integer values).