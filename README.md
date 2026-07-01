# *LRRK2* in Focus: A Global Browser Linking Genetic Diversity to Functional Effects

`GP2 ❤️ Open Science 😍`

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** May 2026

## Summary
This is the online repository for the research project titled ***"LRRK2 in Focus: A Global Browser Linking Genetic Diversity to Functional Effects"***. Here we present an interactive web browser with genomic, clinical, and functional data about *LRRK2* and its relationship with Parkinson's Disease (PD).

### Data Statement 
* All GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson's Disease and are available via application on the website. The GP2 PD case and control data are available via the GP2 website (https://gp2.org; release 11: https://doi.org/10.5281/zenodo.17753486). 
    * Genotyping imputation, quality control, ancestry prediction, and processing were performed using GenoTools (v1.0.0), publicly available on GitHub

### Helpful Links 
- [GP2 Website](https://gp2.org/)
    - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
    - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

# Repository Orientation 
```
THIS_REPO
├── entry.R
├── fetch_data.ipynb
├── global.R
├── ** kinase_activity.tsv
├── ** lrrk2_combined_exome.tsv
├── ** lrrk2_combined_imputed.tsv
├── ** lrrk2_combined_raw.tsv
├── ** lrrk2_combined_wgs.tsv
├── modules
│   ├── annotation_summary_table.R
│   ├── bar_chart.R
│   ├── domain_diagram.R
│   ├── gene_overview.R
│   ├── gene_var_table.R
│   ├── other_resources.R
│   ├── variant_details.R
│   └── world_map.R
├── ** pathogenicity_sources.tsv
├── README.md
├── server.R
├── ui.R
├── util.R
└── www
    ├── clinvar.png
    ├── ensembl.png
    ├── genecards.png
    ├── gnomad.png
    ├── gp2_genome_browser.png
    ├── gp2_gwas_browser.png
    ├── gp2.png
    ├── ldlink.png
    ├── mdsgene.png
    ├── ncbi.png
    ├── omim.png
    ├── uniprot.png
    └── world_map.png

```

**Note: kinase_activity.tsv, pathogenicity_sources.tsv, and lrrk2_combined_[exome/imputed/raw/wgs].tsv are required for the browser but are not hosted in this repository to protect sensitive data. For more details about these files...
- ***kinase_activity.tsv*** includes the columns "Conservation_Score" (Integer between 1-10), "Variant" (ie, A123W), "Mean_pRAB10/RAB10" (Kinase activity measurement), "SD" (Standard deviation), and "Interpretation" ("activating" if that variant is kinase-active, otherwise blank).
- ***pathogenicity_sources.tsv*** includes the columns "Variant" (ie, p.A123W) and "Source" (comma-separated list of resources used to determine disease-association).
- ***lrrk2_combined_[exome/imputed/raw/wgs].tsv*** are each produced by running fetch_data.ipynb with the GP2 datasets.
  
---
### Analysis Scripts
* Languages: Python, R, RShiny, bash

| **File**                           |                  Description                                          |
|------------------------------------|-----------------------------------------------------------------------|
| fetch_data.ipynb                   | Notebook used to fetch and process GP2 data                           |
| global.R                           | Load RShinhy libraries                                                |
| modules/annotation_summary_table.R | Code defining widget for annotation counts across all variants        |
| modules/bar_chart.R                | Code defining widget for kinase activity bar chart                    |
| modules/domain_diagram.R           | Code defining widget for protein and cDNA diagrams                    |
| modules/gene_overview.R            | Code defining widget for positional metadata for LRRK2                |
| modules/gene_var_table.R           | Code defining widget for the main variant data table                  |
| modules/other_resources.R          | Code defining widget for relevant external resources                  |
| modules/variant_details.R          | Code defining widget for variant popup window                         |
| modules/world_map.R                | Code defining widget for world map with variant spectrum donut charts |
| server.R                           | Main server code for RShiny                                           |
| ui.R                               | Main UI code for RShiny                                               |
| util.R                             | Utility functions                                                     |
| www/*                              | .png image files                                                      |

---

# Software 
|               Software              |  Version(s) |            Resource URL             |       RRID      |           Notes            |
|:-----------------------------------:|:-----------:|:-----------------------------------:|:---------------:|:--------------------------:|
| Python Programming Language         | 3.9.19      | http://www.python.org/              | RRID:SCR_008394 | Used for fetching data     |
| R Project for Statistical Computing | 4.4.2       | http://www.r-project.org/           | RRID:SCR_001905 | Used to host RShiny        |
| RShiny                              | 1.9.1       | https://shiny.posit.co/             |       --        | User interface development |
| PLINK                               | 2.0         | http://www.nitrc.org/projects/plink | RRID:SCR_001757 | Used for genetic analyses  |
