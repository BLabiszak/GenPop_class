# Lab 2 — Genome-wide scans for selection

Identify candidate SNPs using pcadapt and OutFLANK, control the false discovery rate, and compare outlier sets across methods and population subsets.

Open [LAB_2.html](LAB_2.html) locally in a web browser and use this folder as your R working directory.

## Files

| File | Contents |
| --- | --- |
| `LAB_2.html` | Polish manual with explanations, code to complete, and report questions |
| `dane.vcf.gz` | Losslessly compressed RADseq dataset containing 3,699 variants and 717 individuals; extract to `dane.vcf` before analysis |
| `metadane.csv` | Sample metadata: `ID`, `SITE`, `LAT`, `LONG`, and `Group`; 24 sites and two groups |
| `LD_Plots.docx` | Supplementary LD plots |

## Starting the exercise

First extract `dane.vcf.gz` with a gzip-compatible archive utility, or run the following once in R from this folder:

```r
install.packages("R.utils") # only if not already installed
R.utils::gunzip("dane.vcf.gz", destname = "dane.vcf", remove = FALSE)
```

This restores the original VCF and retains the compressed copy. Then use `dane.vcf` and `metadane.csv` wherever the manual asks for the input filenames:

```r
vcf.path <- "dane.vcf"
meta.path <- "metadane.csv"
meta <- read.csv(meta.path)
```

The sample order in the supplied metadata matches the sample order in the VCF. Preserve that correspondence when subsetting or reordering either dataset.

The manual uses `vcfR`, `pcadapt`, `qvalue`, `OutFLANK`, `ggVennDiagram`, `ggplot2`, and `pegas`. See the [main README](../README.md#software-and-packages) for installation commands.

Complete the PCA and outlier analyses, compare methods with Venn diagrams, and repeat the OutFLANK analysis for the two selected range-edge populations. Include the code, figures, and answers requested in the manual in your report.
