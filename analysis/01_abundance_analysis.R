# =============================================================
# Script: Abundance analysis
# Description: Analyzes the relative abundance of AMR genes and
#              drug classes across samples and time periods.
# Input:  all_data, all_data_metadata (from 00_data_loading.R)
# Output: - Abundancia_top10_ARGs_global.png
#         - heatmap_args.png
#         - Abundancia_top10_AMRclass.png
#         - Boxplot_abundance_period.png
# =============================================================


# --- Top 10 globally most abundant ARGs ---

# Identify the 10 ARGs with the highest cumulative TPM across all samples
top_gen_global <- all_data_metadata %>%
  group_by(`ARO Term`) %>%
  summarise(Total_TPM = sum(TPM), .groups = "drop") %>%
  arrange(desc(Total_TPM)) %>%
  slice_head(n = 10) %>%
  pull(`ARO Term`)

# Compute relative abundance (%) per sample for the top 10 ARGs
top_gen_global_data <- all_data_metadata %>%
  filter(`ARO Term` %in% top_gen_global) %>%
  group_by(Sample, `ARO Term`, Periodo) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percent_Abundance = (TPM / sum(TPM)) * 100) %>%
  ungroup() %>%
  mutate(Sample_Periodo = interaction(Periodo, Sample, sep = " - "))


# Stacked bar plot: top 10 ARGs per sample (global ranking)
ggplot(top_gen_global_data,
       aes(x = Sample, y = Percent_Abundance,
           fill = factor(`ARO Term`, levels = top_gen_global))) +
  geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.1) +
  scale_fill_viridis_d(option = "H", name = "ARGs") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Sample", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("results/figures/Abundancia_top10_ARGs_global.png", width = 15, height = 6, dpi = 300)


# --- Presence/absence table for top 10 ARGs ---
# Shows how consistently each ARG appears across samples

top_gen_global_vector <- top_gen_global_data$`ARO Term`

presencia_top10_AMRglobal <- all_data_metadata %>%
  filter(`ARO Term` %in% top_gen_global_vector) %>%
  group_by(Sample, `ARO Term`) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  mutate(Presence = ifelse(TPM > 0, 1, 0)) %>%
  group_by(`ARO Term`) %>%
  summarise(
    Samples_present = sum(Presence),
    Total_samples   = n_distinct(all_data_metadata$Sample)
  ) %>%
  mutate(Proportion = Samples_present / Total_samples)

print(presencia_top10_AMRglobal)


# --- Heatmap: top 30 most abundant ARGs ---

top_gen_heatmap <- all_data_metadata %>%
  group_by(`ARO Term`) %>%
  summarise(Total_TPM = sum(TPM), .groups = "drop") %>%
  arrange(desc(Total_TPM)) %>%
  slice_head(n = 30)

filter_data_heatmap <- all_data_metadata %>%
  filter(`ARO Term` %in% top_gen_heatmap$`ARO Term`)

# Build wide matrix: rows = ARGs, columns = samples
heatmap_data <- filter_data_heatmap %>%
  group_by(`ARO Term`, Sample) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = TPM, values_fill = 0)

matriz_heatmap <- as.data.frame(heatmap_data)
rownames(matriz_heatmap) <- matriz_heatmap$`ARO Term`
matriz_heatmap     <- matriz_heatmap %>% select(-`ARO Term`)
matriz_log_heatmap <- log10(matriz_heatmap + 1)

# Sample annotation by period
anotacion <- all_data_metadata %>%
  distinct(Sample, Periodo) %>%
  column_to_rownames("Sample")

annotation_colors <- list(
  Periodo = c("PGE" = "orange", "NO PGE" = "green")
)

# Heatmap with log10(TPM) values and period annotation
pheatmap(
  matriz_log_heatmap,
  annotation_col          = anotacion,
  annotation_colors       = annotation_colors,
  color                   = colorRampPalette(c("blue", "white", "red"))(100),
  scale                   = "row",
  clustering_distance_rows = "euclidean",
  clustering_method       = "ward.D2",
  cluster_cols            = TRUE,
  cluster_rows            = TRUE,
  show_rownames           = TRUE,
  show_colnames           = TRUE,
  fontsize                = 8,
  angle_col               = 45,
  cellwidth               = 13,
  cellheight              = 13,
  gaps_col                = cumsum(rle(as.character(anotacion$Periodo))$lengths),
  filename                = "results/figures/heatmap_args.png",
  width = 14, height = 9, dpi = 300
)


# --- Top 10 drug classes by relative abundance ---

top_clases <- all_data %>%
  group_by(`Drug Class`) %>%
  summarise(Total_TPM = sum(TPM), .groups = "drop") %>%
  arrange(desc(Total_TPM)) %>%
  slice_head(n = 10) %>%
  pull(`Drug Class`)

top_clases_data <- all_data %>%
  filter(`Drug Class` %in% top_clases) %>%
  group_by(Sample, `Drug Class`) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percent_Abundance = (TPM / sum(TPM)) * 100) %>%
  ungroup()

# Stacked bar plot: top 10 drug classes per sample
ggplot(top_clases_data,
       aes(x = Sample, y = Percent_Abundance,
           fill = factor(`Drug Class`, levels = top_clases))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d(option = "H", name = "Drug Class") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Sample", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 7)
  )

ggsave("results/figures/Abundancia_top10_AMRclass.png", width = 14, height = 6, dpi = 300)


# --- Boxplot: total ARG abundance by period ---

ggplot(all_data_metadata, aes(x = Periodo, y = TPM, fill = Periodo)) +
  geom_boxplot(width = 0.6, alpha = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.09, size = 1.5, alpha = 0.4, aes(color = Periodo)) +
  scale_y_log10(
    labels = scales::scientific_format(),
    breaks = 10^(0:10)
  ) +
  stat_compare_means(
    method  = "wilcox.test",
    mapping = aes(label = paste("p =", after_stat(p.format))),
    label.x = 1.5,
    size    = 5,
    vjust   = 1.5
  ) +
  scale_fill_manual(values  = c("PGE" = "#FF6B00", "NO PGE" = "#4daf4a")) +
  scale_color_manual(values = c("PGE" = "#FF6B00", "NO PGE" = "#4daf4a")) +
  labs(x = "Period", y = "Abundance (log10(TPM))") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )

ggsave("results/figures/Boxplot_abundance_period.png", width = 12, height = 6, dpi = 300)


# --- Normality tests for TPM distribution ---

# Kolmogorov-Smirnov test on raw TPM
ks.test(all_data$TPM, "pnorm",
        mean = mean(all_data$TPM),
        sd   = sd(all_data$TPM))

# Kolmogorov-Smirnov test on log10-transformed TPM
ks.test(all_data$TPM, "pnorm",
        mean = mean(log10(all_data$TPM + 1)),
        sd   = sd(log10(all_data$TPM + 1)))

# Shapiro-Wilk test
# Note: Shapiro-Wilk is recommended for small samples (n < 5000).
# With large datasets it tends to reject normality even for minor deviations.
shapiro.test(all_data$TPM)
shapiro.test(log10(all_data$TPM + 1))

# Log10-transformed TPM column (used in subsequent visualizations)
all_data_metadata$logTPM <- log10(all_data_metadata$TPM + 1)

# Q-Q plots by period
ggqqplot(all_data_metadata, x = "logTPM", facet.by = "Periodo",
         title = "Q-Q plot by period (log10 TPM)")

# Density plots by period
ggdensity(all_data_metadata, x = "logTPM", fill = "Periodo",
          add = "mean", rug = TRUE, facet.by = "Periodo")

# Shapiro-Wilk test by period
shapiro.test(all_data_metadata$logTPM[all_data_metadata$Periodo == "PGE"])
shapiro.test(all_data_metadata$logTPM[all_data_metadata$Periodo == "NO PGE"])