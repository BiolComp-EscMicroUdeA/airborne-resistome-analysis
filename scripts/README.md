# Scripts

Pipeline for antimicrobial resistance gene (ARG) identification from metagenomic reads
using [RGI (Resistance Gene Identifier)](https://github.com/arpcard/rgi) and the 
[CARD database](https://card.mcmaster.ca/).

## Requirements
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/)

## Installation

### 1. Create and activate the RGI conda environment

```bash
mamba create --name rgi \
    --channel conda-forge \
    --channel bioconda \
    --channel defaults \
    rgi=6.0.3

conda activate rgi
```

## Usage

### Step 1 — Database setup (only needed once)
Downloads and loads CARD (v4.0.0) and WildCARD (v4.0.2) databases into RGI.

```bash
bash 00_setup_card_database.sh
```

### Step 2 — Run RGI bwt
Runs RGI bwt with KMA aligner over all samples listed in the script.
Before running, open `01_rgi_bwt.sh` and adjust:
- `WORK_DIR`: path to the working directory
- `DATA_DIR`: path to the directory containing the raw reads
- `samples`: list of sample names to process

```bash
bash 01_rgi_bwt.sh
```

## Software versions
| Tool | Version |
|------|---------|
| RGI | 6.0.3 |
| CARD | 4.0.0 |
| WildCARD | 4.0.2 |
| Aligner | KMA |
