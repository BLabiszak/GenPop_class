# Lab 4 — Genomic vulnerability to future climate

Examine candidate allele-frequency patterns, obtain future climate projections, and use AlleleShift to model allele-frequency changes. Summarize genomic vulnerability across populations and regions and discuss the limitations of the approach.

Open [Lab_4.html](Lab_4.html) locally in a web browser and use this folder as your R working directory. The HTML page retains an earlier `LAB_3` title; its content is the Lab 4 exercise.

## Files

| File | Contents / loading method |
| --- | --- |
| `Lab_4.html` | Polish manual with exercises and report questions |
| `GENIND_data` | Saved Scots pine genind object; load with `readRDS()` |
| `allelic_dosage_pines` | Saved genotype dosage object shared with Lab 3; load with `readRDS()` |
| `allele_freqs` | Saved population allele-frequency table; load with `readRDS()` |
| `env_df` | Saved environmental predictor table; load with `readRDS()` |
| `PAV_df` | Saved candidate-SNP table used in the exercise; load with `readRDS()` |
| `meta_ind.csv` | Individual metadata; load with `read.csv()` |
| `meta_pops.csv` | Population coordinates and region assignments; load with `read.csv()` |
| `calc_gv.R` | Defines `calculate_genomic_vulnerability()` |
| `plot_gv.R` | Defines `plot_genomic_vulnerability()` for spatial visualization |

## Practical notes

Load each helper from its own file. The mapping section in the original manual repeats `source("calc_gv.R")`; it needs `source("plot_gv.R")` for the plotting function:

```r
source("calc_gv.R")
source("plot_gv.R")

meta_pops <- read.csv("meta_pops.csv")
coords_sf <- sf::st_as_sf(meta_pops, coords = c("lon", "lat"), crs = 4326)
```

In addition to the packages listed in the manual, the mapping helper uses `raster` and `sp`. See the [main README](../README.md#software-and-packages) for installation commands.

The helper functions assume consistent population order. Before calculating or plotting vulnerability, check that prediction rows, population names, metadata, and coordinates refer to the same populations in the same order.

The manual uses older RCP labels in several places, including `RCP 4.6`. For the CMIP6 downloads in this exercise, use the following SSP codes:

| Scenario | `ssp` argument in `geodata::cmip6_world()` |
| --- | --- |
| SSP1-2.6 | `"126"` |
| SSP2-4.5 | `"245"` |
| SSP5-8.5 | `"585"` |

The exercise uses model `CNRM-CM6-1` and period `2041-2060`. See the [geodata documentation](https://search.r-project.org/CRAN/refmans/geodata/html/cmip6.html) for the download arguments.

Treat the projected allele-frequency changes as model outputs. A good fit to the current data does not independently validate future predictions. Use the final discussion questions to consider model assumptions, climate uncertainty, marker coverage, and possible conservation applications.
