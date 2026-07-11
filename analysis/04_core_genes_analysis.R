# --- Core genes analysis ---

# Step 1: Calculate presence frequency per ARG and period
presencia_frecuencia <- all_data_metadata %>%
  mutate(Presence = ifelse(TPM > 0, 1, 0)) %>%
  group_by(Periodo, `ARO Term`) %>%
  summarise(
    freq = mean(Presence),
    .groups = "drop"
  )

# Step 2: Keep ARGs present in at least 70% of samples per period
genes_estables <- presencia_frecuencia %>%
  filter(freq >= 0.7) %>%
  group_by(Periodo) %>%
  summarise(GENES = list(`ARO Term`))

genes_PGE_estables   <- genes_estables$GENES[[ which(genes_estables$Periodo == "PGE") ]]
genes_NOPGE_estables <- genes_estables$GENES[[ which(genes_estables$Periodo == "NO PGE") ]]

# Step 3: Core genes (present in both periods)
core_genes_estables <- intersect(genes_PGE_estables, genes_NOPGE_estables)

# Step 4: Exclusive genes per period
genes_exclusivos_PGE    <- setdiff(genes_PGE_estables, genes_NOPGE_estables)
genes_exclusivos_NO_PGE <- setdiff(genes_NOPGE_estables, genes_PGE_estables)

# Step 5: Export as a single table
core_genes_df <- data.frame(
  ARG    = c(core_genes_estables,
             genes_exclusivos_PGE,
             genes_exclusivos_NO_PGE),
  Period = c(rep("Both",   length(core_genes_estables)),
             rep("PGE",    length(genes_exclusivos_PGE)),
             rep("NO PGE", length(genes_exclusivos_NO_PGE)))
)

write.csv(core_genes_df, "ruta/core_genes.csv", row.names = FALSE)