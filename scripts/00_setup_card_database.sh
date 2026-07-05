#!/bin/bash

# =============================================================
# Script: CARD and WildCARD database setup for RGI
# RGI version: 6.0.3
# CARD version: 4.0.0
# WildCARD version: 4.0.2
# =============================================================

# --- Paths (adjust according to the environment) ---
WORK_DIR="/path/to/working/directory"

cd "$WORK_DIR"

# --- Download CARD database ---
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json

# --- Download WildCARD database ---
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir -p wildcard
tar -xjf wildcard_data.tar.bz2 -C wildcard
gunzip wildcard/*.gz

# --- Create annotation files ---
rgi card_annotation -i "$WORK_DIR/card.json" > card_annotation.log 2>&1
rgi wildcard_annotation -i wildcard --card_json "$WORK_DIR/card.json" -v 4.0.2 > wildcard_annotation.log 2>&1

# --- Load all data into RGI ---
rgi load \
    --card_json "$WORK_DIR/card.json" \
    --debug --local \
    --card_annotation "$WORK_DIR/localDB/card_database_v4.0.0.fasta" \
    --card_annotation_all_models "$WORK_DIR/localDB/card_database_v4.0.0_all.fasta" \
    --wildcard_annotation wildcard_database_v4.0.2.fasta \
    --wildcard_annotation_all_models wildcard_database_v4.0.2_all.fasta \
    --wildcard_index "$WORK_DIR/wildcard/index-for-model-sequences.txt" \
    --wildcard_version 4.0.2 \
    --amr_kmers "$WORK_DIR/wildcard/all_amr_61mers.txt" \
    --kmer_database "$WORK_DIR/wildcard/61_kmer_db.json" \
    --kmer_size 61

# --- Verify database version ---
echo "Loaded database version:"
rgi database --version --local
