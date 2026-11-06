# Lab 3 — Genotype–environment associations

Explore Scots pine population structure along the Finland–Baltic–Poland transect, obtain bioclimatic predictors, and use redundancy analysis (RDA) to identify candidate SNP–environment associations. Compare RDA candidates with the pcadapt and OutFLANK results.

Open [LAB_3.html](LAB_3.html) locally in a web browser and use this folder as your R working directory.

## Files

| File | Contents / loading method |
| --- | --- |
| `LAB_3.html` | Polish manual with exercises and report questions |
| `GENIND_data` | Saved `adegenet` genind object; `readRDS("GENIND_data")` |
| `allelic_dosage_pines` | Saved genotype dosage object; `readRDS("allelic_dosage_pines")` |
| `meta_ind.csv` | Individual metadata: `ID`, `pop`, `region` |
| `meta_pops.csv` | Population metadata: `lon`, `lat`, `pop`, `region` |
| `RDA_outliers_finder.R` | Functions to select extreme SNP loadings and annotate correlations with environmental predictors |

The metadata cover 440 individuals and 39 populations. Region codes are `FIN` (Finland), `BAL` (Baltic countries), and `POL` (Poland).

## Practical notes

- Load `adegenet` before working with the genind object. Use `readRDS()` for the extensionless R objects and `read.csv()` for the metadata.
- The manual calls `st_as_sf()` and `corrplot()`; load `sf` and `corrplot` explicitly. For PCA, use `ade4::dudi.pca()` if `dudi.pca()` is not found.
- Where plotting examples use `map_df`, use `meta_pops`, which contains `lon` and `lat`, or define `map_df <- meta_pops` first.
- WorldClim layers are downloaded during the exercise. Save your derived allele-frequency and environmental tables for later use, as directed in the manual.
- The dosage-orientation comments in the manual are inconsistent. Inspect `dim(allelic_dosage)` and sample identifiers before transposing; align the individual dimension with the 440-row metadata table and check the format required by each analysis.
- A scree plot helps explore the importance of RDA axes; it is not a test of statistical significance. The manual also provides the permutation-test call `anova.cca(..., by = "axis")`.

See the [main README](../README.md#software-and-packages) for package installation. Work through the manual section by section and retain the figures and code needed for your report.
