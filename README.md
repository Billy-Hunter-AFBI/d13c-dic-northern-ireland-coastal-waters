# Relationships between the isotopic signatures of dissolved inorganic carbon and environmental variability in Northern Ireland coastal waters

## Overview

This repository contains the data and R code used to analyse spatial variability in dissolved inorganic carbon stable isotope signatures (δ¹³C-DIC) across coastal monitoring stations in Northern Ireland.

The study explores relationships between:

- δ¹³C-DIC
- Dissolved inorganic carbon (DIC) concentration
- Chlorophyll-a
- Nutrient concentrations
- Temperature
- Salinity

using simple linear regression models.

Given the limited sample size (n = 12 stations), analysis was restricted to univariate linear models. Multiple regression and interaction models were not considered appropriate.

---

## Repository Structure

```
.
├── data/
│   └── 2020_WFD_del13C_trial.csv
│
├── scripts/
│   ├── WFD_D13C_Analysis.R
│   └── WFD_D13C_Figures.R
│
├── Figures/
│   └── (generated figures)
│
├── Regression_Plots/
│   └── (generated model plots)
│
├── Diagnostics/
│   └── (generated model diagnostics)
│
└── README.md
```

---

## Data

The dataset contains observations from 12 Northern Ireland coastal monitoring stations sampled during July 2020.

Key variables include:

| Variable | Description |
|-----------|-------------|
| STATION | Sampling station identifier |
| Lat | Latitude (decimal degrees) |
| Long | Longitude (decimal degrees) |
| D13c | δ¹³C-DIC signature (‰) |
| DIC_molar | Dissolved inorganic carbon concentration (µmol kg⁻¹) |
| CHL | Chlorophyll-a concentration (mg L⁻¹) |
| INORG.N | Dissolved inorganic nitrogen (µmol L⁻¹) |
| PHOSPHATE | Phosphate concentration (µmol L⁻¹) |
| SILICA | Silicate concentration (µmol L⁻¹) |
| AMMONIA | Ammonia concentration (µmol L⁻¹) |
| Sal | Salinity |
| Temp | Temperature (°C) |

---

## Analysis Workflow

### 1. Statistical Analysis

Script:

```r
WFD_D13C_Analysis.R
```

This script:

- Imports and cleans the dataset
- Calculates correlation matrices
- Fits simple linear regression models
- Produces model summaries
- Generates regression figures
- Generates diagnostic plots
- Repeats all analyses after removal of station WF24 as a sensitivity test

Outputs include:

- `Linear_Model_Summary.csv`
- Regression plots (`*.png`)
- Diagnostic plots (`*.png`)

---

### 2. Figure Generation

Script:

```r
WFD_D13C_Figures.R
```

This script creates publication-ready figures including:

#### Figure 1

Spatial distribution of stations and environmental controls on δ¹³C-DIC.

Panels include:

- Station map
- Temperature
- Chlorophyll-a
- DIN
- Phosphate
- Silicate
- Ammonia
- DIC concentration
- Salinity

#### Supplementary Figure S2

Relationships between DIC concentration and environmental variables.

#### Supplementary Figure S3

Sensitivity analysis excluding station WF24.

---

## Reproducing the Analysis

1. Clone or download this repository.

2. Open the R scripts in RStudio.

3. Install required packages:

```r
install.packages(
  c(
    "tidyverse",
    "broom",
    "sf",
    "patchwork",
    "ggrepel"
  )
)
```

4. Update file paths if required.

5. Run:

```r
source("scripts/WFD_D13C_Analysis.R")
```

followed by:

```r
source("scripts/WFD_D13C_Figures.R")
```

---

## R Version

Analysis was developed using:

```r
sessionInfo()
```

Users are encouraged to record package versions when reproducing analyses.

---

## Data Availability

All data and analytical code used in this study are archived through Zenodo and are publicly available at:

**DOI: [Insert DOI after publication]**

---

## Citation

If using these data or code, please cite:

Hunter, W.R. and Smyth, C.. *Relationships between dissolved inorganic carbon isotopic signatures and environmental variability in Northern Ireland coastal waters*. [Journal].

Dataset and code archive:

DOI: [Insert Zenodo DOI]

---

## License

### Code

MIT License

### Data

Creative Commons Attribution 4.0 International (CC-BY 4.0)
