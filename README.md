# Population Genomics for Bioinformaticians

Practical materials for a five-class course on population genomics, detecting natural selection, and local adaptation. The exercises use R to move from candidate-gene analyses to genome-wide scans, genotype–environment associations, and projections of genomic vulnerability under future climates.

**Instructor:** Bartosz Łabiszak, Adam Mickiewicz University in Poznań.

Each lab folder contains an HTML manual and the data or helper functions needed for the exercises. Materials for **Labs 1–4** are currently available. Lab 1 has English and Polish manuals; Labs 2–4 are in Polish, with R code and package documentation largely in English.

## Course materials

| Lab | Topic | Main analyses | Materials |
| --- | --- | --- | --- |
| 1 | Detecting selection in candidate genes | Genetic diversity, Tajima’s D, linkage disequilibrium, Fay and Wu’s H, and population differentiation (FST) | [Folder](LAB_1/) · [English manual](LAB_1/LAB_1.html) · [Polish manual](LAB_1/LAB_1_PL.html) |
| 2 | Genome-wide scans for selection | PCA, pcadapt, OutFLANK, false discovery rate control, and comparison of candidate SNP sets | [Folder and instructions](LAB_2/) · [Manual](LAB_2/LAB_2.html) |
| 3 | Genotype–environment associations | Population structure, WorldClim data, redundancy analysis (RDA), and comparison with genome-scan outliers | [Folder and instructions](LAB_3/) · [Manual](LAB_3/LAB_3.html) |
| 4 | Genomic vulnerability to climate change | Candidate allele frequencies, CMIP6 climate projections, AlleleShift, and spatial visualization | [Folder and instructions](LAB_4/) · [Manual](LAB_4/Lab_4.html) |

## Getting started

1. Download the repository using **Code → Download ZIP** and extract it, or clone it:

   ```bash
   git clone https://github.com/BLabiszak/GenPop_class.git
   ```

2. Open the selected lab’s `.html` manual in a web browser from your computer. GitHub’s file view may display the HTML source instead of the formatted manual.
3. Start R or RStudio and set the working directory to that lab’s folder. In RStudio, use **Session → Set Working Directory → Choose Directory**.
4. Install the packages required for the lab, then work through the manual in order. Read the folder’s README for file descriptions and practical notes.

The manuals are guided exercises: placeholders such as `XXX`, `...`, and incomplete function arguments are for students to fill in. Run and adapt the code section by section, keeping your own analysis script and the figures requested in the manual.

## Software and packages

The exercises require R; RStudio is a convenient interface. Basic familiarity with R objects, data frames, indexing, and plotting is expected. Internet access is needed to install packages and download climate data in Labs 3–4.

The packages used across the course can be installed as follows. For a single lab, install only the packages it needs.

```r
install.packages(c(
  "ape", "pegas", "adegenet", "hierfstat", "tidyverse",
  "vcfR", "pcadapt", "ggVennDiagram", "ggplot2", "ggpubr", "poppr",
  "rnaturalearth", "rnaturalearthdata", "geodata", "terra", "sf",
  "psych", "adespatial", "vegan", "corrplot", "ade4",
  "AlleleShift", "gstat", "raster", "sp", "RColorBrewer",
  "BiocManager", "remotes"
))

BiocManager::install("qvalue")
remotes::install_github("whitlock/OutFLANK")
remotes::install_github("gavinsimpson/ggvegan")
```

Package versions are not pinned. These commands describe the dependencies used in the materials, rather than a tested environment specification. Keep the output of `sessionInfo()` with your analysis if you need to document your setup.

## Datasets and file formats

| Lab | Dataset | Formats |
| --- | --- | --- |
| 1 | Candidate-gene sequence alignments, including outgroup files | FASTA and R code |
| 2 | RADseq teaching dataset: 717 individuals, 24 sampling localities, two geographic groups, and 3,699 SNPs | Gzip-compressed VCF, CSV metadata, and a supplementary Word document |
| 3–4 | Scots pine (*Pinus sylvestris*) along a Finland–Baltic–Poland transect: metadata for 440 individuals and 39 populations; the manual describes 10,597 SNPs after quality control | RDS objects, CSV metadata, and R helper functions |

For Lab 2, extract `dane.vcf.gz` to `dane.vcf` before analysis; the [Lab 2 instructions](LAB_2/README.md) include an R command for this step.

Files such as `GENIND_data`, `allelic_dosage_pines`, `allele_freqs`, `env_df`, and `PAV_df` are saved R objects **without a `.rds` extension**. Load them with `readRDS()`, for example:

```r
library(adegenet)
genind_pines <- readRDS("GENIND_data")
meta_ind <- read.csv("meta_ind.csv")
meta_pops <- read.csv("meta_pops.csv")
```

Keep the supplied filenames and check sample/population identifiers and their order when combining genetic data with metadata. Labs 3 and 4 include copies of shared inputs so that each folder can be used as its own working directory.

## Interpretation and data sources

These materials are intended for teaching. Outlier SNPs and environmental associations are candidates for further investigation; they do not by themselves establish selection or a causal role in adaptation. Genomic vulnerability exercises illustrate model-based projections and their assumptions.

The Scots pine manuals identify the source study as *Adaptive potential and genomic vulnerability of keystone forest tree species to climate change: a case study in Scots pine*. Follow the references in the manuals for the methods and data sources. Current and projected climate data are obtained from [WorldClim](https://www.worldclim.org/).

**Practical notes:** [Lab 2](LAB_2/README.md), [Lab 3](LAB_3/README.md), and [Lab 4](LAB_4/README.md) describe the supplied files and clarify points needed to follow the original manuals.
