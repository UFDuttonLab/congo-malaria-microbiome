# Gut Microbiome in Infancy Predicts Malaria Susceptibility — Analysis Pipeline

**Dutton, Follis, Munaweera, Maisha, Mulligan, and Moore**
Currently under review at *Frontiers in Microbiology* (revision 1)

This repository contains the full analytical chain — from raw PacBio reads to manuscript figures — for the submitted Rev1 manuscript. Given the two data files archived on Zenodo (see below), re-running this pipeline reproduces every figure, table, and statistical result reported in the manuscript.

---

## Repository structure

- **01_sequence_processing/** — Raw FASTQ → phyloseq object (DADA2 workflow)
  - `Congo_InfantGut_2ndRound_QC.Rmd`
  - `Congo_InfantGut_2ndRound_QC.html`
- **02_analysis/** — phyloseq object → manuscript figures and stats
  - `Congo_Malaria_Rev1.Rmd`
  - `Congo_Malaria_Rev1.html`
  - `Congo_Malaria_Rev1_files/` (supporting images for the HTML)
- **renv/** — Project-local package cache
- `renv.lock` — Exact package versions used for the analysis
- `.Rprofile` — Activates renv on R session start
- `.zenodo.json` — Zenodo metadata for automatic DOI minting
- `CITATION.cff` — Citation metadata for GitHub
- `LICENSE` — MIT license
- `README.md` — This file

---

## Pipeline overview

The pipeline runs in two sequential stages.

### Stage 1: Sequence processing (01_sequence_processing/)

`Congo_InfantGut_2ndRound_QC.Rmd` takes raw PacBio circular consensus sequencing (CCS) FASTQ files from two SMRT cells and runs the DADA2 workflow end-to-end:

1. Primer removal with F27 (AGRGTTYGATYMTGGCTCAG) and R1492 (RGYTACCTTGTTACGACTT)
2. Length filtering (1,000–1,600 bp, maxEE=2, minQ=3)
3. Error learning with the PacBio error model (`PacBioErrfun`)
4. Denoising with DADA2
5. Merging of sequence tables from the two cells
6. Chimera removal (`removeBimeraDenovo`, consensus method)
7. Taxonomy assignment against SILVA v138.1 (species-level training set)
8. Construction of the final phyloseq object (`psCombinedCongo_V5.rds`)

This stage was run on the HiPerGator cluster at the University of Florida (100 GB RAM, 12 cores). The file paths in the Rmd are hardcoded to the original author's environment; reproducers would need to update them to point to their local FASTQ directories.

The compiled HTML provides a full record of the original run, including read counts at each filtering step, error profiles, and chimera statistics.

### Stage 2: Analysis (02_analysis/)

`Congo_Malaria_Rev1.Rmd` takes the phyloseq object from Stage 1 (plus the malaria risk score CSV) and reproduces every figure, table, and statistical result in the manuscript:

| Section | Produces |
|---|---|
| Figure 2 (A–E) | Longitudinal SplinectomeR of top species, PCoA by timepoint, PCoA by malaria status at each timepoint |
| Figure 3 (A–H) | LEfSe differential abundance, SplinectomeR for key taxa |
| Figure 5 (A–D) | ANCOM-BC differential abundance at 6W, 3M, 6M, 1Y |
| Figure 5A sensitivity | Bednet-adjusted ANCOM-BC at 6W (11 taxa remain significant at q < 0.05) |
| Figure 6 (A–D) | Beta and alpha diversity before/during/post malaria |
| Figure S1 | Alpha diversity by timepoint |
| Figure S2 (A–D) | ANCOM-BC heatmaps at each timepoint |
| Classifier (6W) | 72-model grid search, best model (species-level k-NN + Boruta, 90.0% balanced accuracy), feature importance |
| Classifier (1Y) | Parallel analysis at one year |
| Reviewer validation | Repeated stratified CV (10 × 5-fold), LOOCV, permutation testing (1,000 permutations), sensitivity covariate models, post-hoc power analysis, covariate inclusion summary |

This stage is fully reproducible via renv.

---

## Data availability

- **Raw PacBio CCS reads** — deposited at the NCBI Sequence Read Archive (released upon publication of the associated manuscript)
- **Processed data files** (`psCombinedCongo_V5.rds` and `working_malaria_risk_score.csv`) — archived at Zenodo: https://doi.org/10.5281/zenodo.19634458
- **Analysis code** (this repository, versioned snapshot) — archived at Zenodo: https://doi.org/10.5281/zenodo.19634438

---

## Reproducing the analysis

### Prerequisites

- R ≥ 4.3
- Pandoc ≥ 2.0 (bundled with recent RStudio installations)
- System libraries for Bioconductor packages: on Debian/Ubuntu, install `libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev`. Usually already present on macOS and Windows.

### Quick start

1. Clone this repository or download a release tarball from the Releases page on GitHub.

2. Download the data files from Zenodo at https://doi.org/10.5281/zenodo.19634458. Place both files in the `02_analysis/` folder:
   - `psCombinedCongo_V5.rds`
   - `working_malaria_risk_score.csv`

3. Open R in the `02_analysis/` folder (in RStudio, open the Rmd directly — RStudio sets the working directory automatically).

4. Restore the package environment by running `install.packages("renv")` if needed, then `renv::restore()`. This installs every package at the exact version used to produce the submitted results. First-time restoration typically takes 20–60 minutes.

5. Knit the Rmd. In RStudio, click Knit. From the command line, run `Rscript -e 'rmarkdown::render("Congo_Malaria_Rev1.Rmd")'`.

A successful knit produces `Congo_Malaria_Rev1.html` with all figures rendered inline, and populates `figures/`, `SixWeek/`, `SixWeekResults/`, and `FirstYear/` subdirectories with the PNG and PDF files used in the manuscript.

Expected knit time: 15–30 minutes on a modern laptop, dominated by the classifier chunks (72 model configurations trained across six taxonomic levels, plus repeated cross-validation and 1,000 permutation iterations).

---

## Package inventory

The analysis uses 34 R packages spanning CRAN, Bioconductor, and GitHub sources. Exact versions are pinned in `renv.lock`.

**CRAN packages**: plyr, magrittr, dplyr, tidyr, tibble, ggplot2, cowplot, gridExtra, viridis, RColorBrewer, reshape2, circlize, pheatmap, vegan, caret, Boruta, randomForest, pwr, knitr, kableExtra, future

**Bioconductor packages**: phyloseq, microbiome, ANCOMBC, ComplexHeatmap, ALDEx2, DESeq2, metagenomeSeq, Maaslin2, lefser

**GitHub / specialty packages**: microeco, file2meco, splinectomeR, dar

---

## Notes on stochastic operations

All stochastic operations (rarefaction, train/test splits, Boruta feature selection, permutation tests, cross-validation folds) use `set.seed()` explicitly. Primary seeds:

- Global data loading: `set.seed(11)`
- Classifier pipeline: seeds embedded in `trans_classifier$cal_split()` calls
- Repeated cross-validation: `set.seed(123)`
- LOOCV: `set.seed(456)`
- Permutation test: `set.seed(789)`
- Sensitivity covariate analyses: `set.seed(202)`

Given identical package versions (via `renv::restore()`) and identical seeds, the classifier metrics and validation statistics reported in the manuscript are exactly reproducible.

---

## Citation

If you use this pipeline, please cite both the code archive and the associated manuscript.

**Code archive**:
> Dutton, C.L., Follis, M., Munaweera, J., Maisha, F.M., Mulligan, C.J., and Moore, J.M. (2026). Analysis pipeline for "Gut Microbiome in Infancy Predicts Malaria Susceptibility." Zenodo. https://doi.org/10.5281/zenodo.19634438

**Data archive**:
> Dutton, C.L., Follis, M., Munaweera, J., Maisha, F.M., Mulligan, C.J., and Moore, J.M. (2026). Data for: Gut Microbiome in Infancy Predicts Malaria Susceptibility. Zenodo. https://doi.org/10.5281/zenodo.19634458

**Associated manuscript** (citation to be added upon acceptance).

---

## Contact

**Julie M. Moore** — juliemmoore@ufl.edu
Department of Infectious Diseases and Immunology, University of Florida
