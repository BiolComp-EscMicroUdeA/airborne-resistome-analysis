# =============================================================
# Script: Diversity analysis (alpha and beta)
# Description: Computes alpha diversity indices (Shannon, Simpson,
#              Richness) and beta diversity (Bray-Curtis + NMDS)
#              for AMR gene communities across sampling periods.
#              Includes statistical tests (Wilcoxon, PERMANOVA,
#              betadisper).
# Input:  all_data_metadata (from 00_data_loading.R)
# Output: - results/figures/Alpha_diversity.png
#         - results/figures/nmds_diversity.png
#         - PERMANOVA and betadisper results (printed to console)
# =============================================================


# --- Build abundance matrix ---
# Rows = samples, columns = ARO terms, values = sum of TPM

mtz_diversidad <- all_data_metadata %>%
  select(Sample, `ARO Term`, TPM) %>%
  group_by(Sample, `ARO Term`) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = `ARO Term`, values_from = TPM, values_fill = 0) %>%
  column_to_rownames("Sample")

# Environmental matrix: one row per sample with period information
matriz_ambiental <- all_data_metadata %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample, Periodo) %>%
  column_to_rownames("Sample")



# ALPHA DIVERSITY


# Compute Shannon, Simpson, and species richness indices
diversidad_alfa <- data.frame(
  Sample  = rownames(mtz_diversidad),
  Shannon = diversity(mtz_diversidad, index = "shannon"),
  Simpson = diversity(mtz_diversidad, index = "simpson"),
  Richness = specnumber(mtz_diversidad)
)

# Add period information
diversidad_alfa <- diversidad_alfa %>%
  left_join(
    matriz_ambiental %>% rownames_to_column("Sample"),
    by = "Sample"
  )

# Reshape to long format for plotting
div_largo <- diversidad_alfa %>%
  pivot_longer(
    cols      = c(Shannon, Simpson, Richness),
    names_to  = "Metric",
    values_to = "Value"
  )

# Wilcoxon test per diversity metric, with BH p-value adjustment
stat.test <- div_largo %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ Periodo) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Periodo", dodge = 0.8) %>%
  mutate(label = paste0("p.adj = ", signif(p.adj, 3)))

# Boxplot + jitter with p-value annotations
ggplot(div_largo, aes(x = Periodo, y = Value, fill = Periodo)) +
  geom_boxplot(alpha = 0.1, outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, aes(color = Periodo)) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  geom_text(
    data = data.frame(
      Metric = c("Richness", "Shannon", "Simpson"),
      label  = c(
        paste("p =", signif(stat.test$p[stat.test$Metric == "Richness"], 3)),
        paste("p =", signif(stat.test$p[stat.test$Metric == "Shannon"],  3)),
        paste("p =", signif(stat.test$p[stat.test$Metric == "Simpson"],  3))
      ),
      y_pos = c(
        max(div_largo$Value[div_largo$Metric == "Richness"]) * 1.05,
        max(div_largo$Value[div_largo$Metric == "Shannon"])  * 1.05,
        max(div_largo$Value[div_largo$Metric == "Simpson"])  * 1.05
      )
    ),
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE, size = 4
  ) +
  scale_fill_manual(values  = c("PGE" = "#FF6B00", "NO PGE" = "#4daf4a")) +
  scale_color_manual(values = c("PGE" = "#FF6B00", "NO PGE" = "#4daf4a")) +
  labs(x = "", y = "Observed values") +
  theme_bw(base_size = 14) +
  theme(
    legend.position    = "none",
    strip.background   = element_rect(fill = "gray90"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )

ggsave("results/figures/Alpha_diversity.png", width = 12, height = 6, dpi = 300)



# BETA DIVERSITY


# Bray-Curtis dissimilarity matrix
dist_bray <- vegdist(mtz_diversidad, method = "bray")

# NMDS ordination (k = 2 dimensions)
nmds <- metaMDS(dist_bray, k = 2, trymax = 100)
stressplot(nmds)

# NMDS plot colored by period
nmds_data         <- as.data.frame(nmds$points)
nmds_data$Periodo <- matriz_ambiental$Periodo

ggplot(nmds_data, aes(MDS1, MDS2, color = Periodo)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Periodo), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = c("PGE" = "#2A9D8F", "NO PGE" = "#E76F51")) +
  scale_fill_manual(values  = c("PGE" = "#2A9D8F", "NO PGE" = "#E76F51")) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/figures/nmds_diversity.png", width = 12, height = 6, dpi = 300)


# --- PERMANOVA ---
# Tests whether community composition differs significantly between periods.
# H0: no difference in AMR community composition between PGE and NO PGE.

permanova <- adonis2(dist_bray ~ Periodo, data = matriz_ambiental, permutations = 999)
print(permanova)
# Result: Pr(>F) = 0.361 — no significant difference between periods


# --- Betadisper (homogeneity of dispersion) ---
# Tests whether the within-group variability is similar between periods.
# If p > 0.05: dispersion is homogeneous (PERMANOVA assumption is met).
# If p < 0.05: differences may reflect dispersion rather than composition.

disp <- betadisper(dist_bray, matriz_ambiental$Periodo)
anova(disp)
# Result: p = 0.9679 — dispersion is homogeneous between periods