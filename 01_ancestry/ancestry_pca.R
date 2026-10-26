# ============================================================
# Ancestry PCA Visualisation
#
# Input:  Ancestry_PCA.eigenval
#         Ancestry_PCA.eigenvec
#         cohorts.txt
# Output: Ancestry_PCA_plot.svg
# ============================================================

library(ggplot2)

# File paths
eigenval_file <- "Ancestry_PCA.eigenval"
eigenvec_file <- "Ancestry_PCA.eigenvec"
cohorts_file  <- "cohorts.txt"

# Colours and shapes
colour_map <- c(
  "NAMA"  = "#9d66a2",
  "YRI"   = "#65b46a",
  "GBR"   = "#6589ab",
  "CHB"   = "#d74c4b",
  "MALAY" = "#fdaa48",
  "UV-HE" = "black",
  "UV-MD" = "black",
  "UV-SD" = "black",
  "V-MD"  = "black"
)

shape_map <- c(
  "NAMA"  = 16,
  "YRI"   = 16,
  "GBR"   = 16,
  "CHB"   = 16,
  "MALAY" = 16,
  "UV-HE" = 15,
  "UV-MD" = 17,
  "UV-SD" = 3,
  "V-MD"  = 8
)

# Eigenvalues: calculate % variance explained
eigenval  <- scan(eigenval_file, quiet = TRUE)
pve       <- round(100 * eigenval / sum(eigenval), 2)
pc1_label <- paste0("PC1 (", pve[1], "% variance explained)")
pc2_label <- paste0("PC2 (", pve[2], "% variance explained)")

# Load eigenvectors
eigenvec <- read.table(
  eigenvec_file,
  header           = TRUE,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

colnames(eigenvec)[1:2] <- c("FID", "IID")

pc_cols <- 3:ncol(eigenvec)
if (!all(grepl("^PC[0-9]+$", colnames(eigenvec)[pc_cols]))) {
  colnames(eigenvec)[pc_cols] <- paste0("PC", seq_along(pc_cols))
}

# Load cohort assignments
cohorts <- read.table(
  cohorts_file,
  header           = FALSE,
  stringsAsFactors = FALSE,
  col.names        = c("IID", "Cohort")
)

# Merge and restore eigenvec row order
pca_data <- merge(eigenvec, cohorts, by = "IID", sort = FALSE)
pca_data <- pca_data[match(eigenvec$IID, pca_data$IID), ]

present_groups   <- unique(pca_data$Cohort)
group_order      <- names(colour_map)[names(colour_map) %in% present_groups]
pca_data$Cohort  <- factor(pca_data$Cohort, levels = group_order)

# PC1 vs PC2 plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Cohort, shape = Cohort)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = colour_map, breaks = group_order) +
  scale_shape_manual(values = shape_map, breaks = group_order) +
  labs(
    x      = pc1_label,
    y      = pc2_label,
    colour = "Population / Cohort",
    shape  = "Population / Cohort"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid   = element_blank()
  )

# Save output
ggsave("Ancestry_PCA_plot.svg", plot = p, width = 8, height = 6)