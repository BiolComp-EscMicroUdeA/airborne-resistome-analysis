# Resistome profile of Particulate Matter (PM2.5) in Medellín, Colombia

This repository contains the bioinformatic and statistical pipeline to analyze the antimicrobial resistance genes (ARGs) profile (resistome) from PM2.5 particulate matter samples collected over an 18-month period in Medellin, Colombia. The project evaluates variations in ARG abundance, diversity, and composition across different climatic periods (PGE vs. NO PGE).

##  Repository Structure

```text
resistome-pm25-18months/
│
├── data/
│   ├── metadata/
│   │   └── Metadata_PM.csv             # 18-month metadata sample descriptors
│   └── raw/
│       └── .gitkeep                    # Target folder for heavy raw FASTQ files
│
├── results/
│   ├── rgi_output/                     # Raw RGI bwt output directories per sample
│   ├── tables/                         # Processed datasets and statistical outputs
│   └── figures/                        
│
├── scripts/
│   ├── 00_setup_card_database.sh      # Bash script to download/configure CARD DB
│   └── 01_rgi_bwt.sh                   # Bash script for read mapping using RGI bwt
│
└── analysis/                           # R analytical pipeline
    ├── 00_data_loading.R               # Parsing RGI outputs and TPM normalization
    ├── 01_abundance_analysis.R         # Global ARG and Drug Class relative abundance
    ├── 02_diversity_analysis.R         # Alpha/Beta diversity (vegan) and Wilcoxon tests
    ├── 03_circos_plot.R                # Chord diagram generation for AMR families
    └── 04_core_genes_analysis.R        # Shared and persistent ARGs (UpSet/Venn analysis)


```

## Software and versions
| Tool | Version |
|------|---------|
| RGI | 6.0.3 |
| CARD database | 4.0.0 |
| WildCARD database | 4.0.2 |
| Aligner | KMA |
| R | ≥ 4.0 |

**Computational Biology Research Group**  
School of Microbiology, Universidad de Antioquia, Colombia

