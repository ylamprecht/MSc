# ==============================================================
# PCA Batch Effect Visualisation
# 
# Input:  merged_batches_autosomes_filtered_pruned_PCA.eigenvec
#         merged_batches_autosomes_filtered_pruned_PCA.eigenval
# Output: batch_pca_PC1_PC2.svg
# ==============================================================

library(ggplot2)
library(dplyr)

# File paths
eigenvec_file <- "merged_batches_autosomes_filtered_pruned_PCA.eigenvec"
eigenval_file <- "merged_batches_autosomes_filtered_pruned_PCA.eigenval"

# Batches
batch1_ids <- c(
  "CVD1008","CVD1133","CVD573","CVD609","CVD614","CVD626","CVD666",
  "CVD684","CVD699","CVD719","CVD755","CVD790","CVD837","CVD838",
  "CVD922","CVD945","CVD953","CVD987","CVD147","CVD190","CVD315",
  "CVD392","CVD415","CVD452","CVD498","CVD555","CVD569","CVD750",
  "CVD789","CVD843","CVD848","CVD930","CVD944","CVD12","CVD1",
  "CVD20","CVD21","CVD22","CVD23","CVD26","CVD27","CVD29","CVD2",
  "CVD30","CVD38","CVD39","CVD3","CVD40","CVD42","CVD43","CVD45",
  "CVD46","CVD4","CVD50","CVD54","CVD59","CVD1004","CVD1128",
  "CVD1158","CVD252","CVD594","CVD595","CVD600","CVD601",
  "CVD605","CVD608","CVD615","CVD622","CVD668","CVD702",
  "CVD721","CVD740","CVD741","CVD744","CVD746","CVD753","CVD781",
  "CVD782","CVD818","CVD863","CVD890","CVD962","CVD967"
)

batch2_ids <- c(
  "1003","1005","1007","1008","1010","1014","1018","1020","1024",
  "1026","1029","1032","1033","1035","1036","1037","1038","1001",
  "1009","1013","1015","1017","1019","1021","1022","1023","1027",
  "1041","1042","1054","1064","1004","1011","1025","1060","1747",
  "1756","549","73","792","899","906"
)

# Load eigenvec 
pca <- read.table(
  eigenvec_file,
  header        = FALSE,
  comment.char  = "",
  stringsAsFactors = FALSE
)

colnames(pca) <- pca[1, ]
colnames(pca)[1] <- "FID"
pca <- pca[-1, ]

pc_cols <- grep("^PC", colnames(pca))
pca[, pc_cols] <- lapply(pca[, pc_cols], as.numeric)

# Assign batch labels
pca <- pca %>%
  mutate(Batch = case_when(
    IID %in% batch1_ids ~ "Batch 1 (30X, DNBSeq-G400)",
    IID %in% batch2_ids ~ "Batch 2 (24X, DNBSeq-T7)"
  ))

# Load eigenvalues → % variance explained
eigenval <- scan(eigenval_file, quiet = TRUE)
pct_var  <- round(100 * eigenval / sum(eigenval), 1)

# Colours
batch_colours <- c(
  "Batch 1 (30X, DNBSeq-G400)" = "#0000FF",
  "Batch 2 (24X, DNBSeq-T7)"   = "#ff7f27"
)

# PC1 vs PC2 plot
p <- ggplot(
  pca %>% dplyr::arrange(Batch == "Batch 2 (24X, DNBSeq-T7)"),
  aes(x = PC1, y = PC2, colour = Batch)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = batch_colours) +
  labs(
    x      = paste0("PC1 (", pct_var[1], "% variance explained)"),
    y      = paste0("PC2 (", pct_var[2], "% variance explained)"),
    colour = "Sequencing batch"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10),
    axis.title      = element_text(size = 11),
    panel.border    = element_rect(colour = "black", fill = NA,
                                   linewidth = 0.6)
  )

# Save outputs
ggsave("batch_pca_PC1_PC2.svg", plot = p, width = 7, height = 5.5)
