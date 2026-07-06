# =============================================================
# Script:Circos plot — Top 10 AMR gene families vs. periods
# Description: Visualizes the relationship between the top 10
#              most abundant AMR gene families and sampling
#              periods (PGE / NO PGE) using a Circos plot.
# Input:  all_data_metadata (from 00_data_loading.R)
# Output: results/figures/Chord_family.png
# =============================================================


# --- Build matrix: top 10 AMR gene families vs. periods ---

# Identify the 10 families with the highest cumulative TPM
top_10_family_circos <- all_data_metadata %>%
  group_by(`AMR Gene Family`) %>%
  summarise(total = sum(TPM), .groups = "drop") %>%
  slice_max(total, n = 10) %>%
  pull(`AMR Gene Family`)

# Compute total TPM per family-period combination
chord_data_family <- all_data_metadata %>%
  filter(`AMR Gene Family` %in% top_10_family_circos) %>%
  group_by(`AMR Gene Family`, Periodo) %>%
  summarise(TPM = sum(TPM), .groups = "drop")

# Convert to wide matrix: rows = AMR families, columns = periods
chord_matrix_family <- pivot_wider(
  chord_data_family,
  names_from  = Periodo,
  values_from = TPM,
  values_fill = 0
)

matr           <- as.matrix(chord_matrix_family[, -1])
rownames(matr) <- chord_matrix_family$`AMR Gene Family`


# --- Chord diagram setup ---

families      <- rownames(matr)
letter_labels <- setNames(LETTERS[1:length(families)], families)

# Color palette for AMR families
family_colors <- setNames(
  RColorBrewer::brewer.pal(n = min(10, length(families)), name = "Set3"),
  families
)

# Colors for period sectors
period_colors <- c("PGE" = "#FF6B00", "NO PGE" = "#00A650")

# Combined color vector for all sectors
grid.col <- c(family_colors, period_colors)

# Color matrix for chord links (colored by AMR family)
link_colors <- matrix(
  family_colors[rownames(matr)],
  nrow = nrow(matr), ncol = ncol(matr)
)


# --- Generate chord diagram ---

png("results/figures/Chord_family.png",
    width  = 17,
    height = 6,
    units  = "in",
    res    = 300)

chordDiagram(
  matr,
  grid.col          = grid.col,
  col               = link_colors,
  transparency      = 0.3,
  annotationTrack   = c("grid"),
  preAllocateTracks = list(track.height = 0.1)
)

# Add labels: full names for periods, letters for AMR families
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index <- CELL_META$sector.index
    
    # Period sectors: show full name
    if (sector.index %in% c("PGE", "NO PGE")) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1],
        sector.index,
        facing     = "bending.inside",
        niceFacing = TRUE,
        adj = c(0.5, 0), cex = 1.2, col = "black"
      )
    }
    
    # AMR family sectors: show letter code
    if (sector.index %in% names(letter_labels)) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1],
        letter_labels[sector.index],
        facing     = "inside",
        niceFacing = TRUE,
        adj = c(0.5, 0), cex = 1.2, col = "black"
      )
    }
  },
  bg.border = NA
)

# Legend mapping letters to AMR family names
legend(
  "right",
  legend = paste(letter_labels, "-", names(letter_labels)),
  fill   = family_colors,
  border = NA,
  bty    = "n",
  xpd    = TRUE,
  cex    = 0.8,
  title  = "AMR Gene Family"
)

dev.off()
circos.clear()