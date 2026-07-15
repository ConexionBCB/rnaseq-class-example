# edgeR differential expression + visualisations (volcano plot, heatmap)
# Reads the featureCounts matrix produced by run_pipeline.sh.
# Contrast: mip6-delta vs heat-shocked control (GEO GSE135568).
suppressMessages({
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
})

# ---- 1. Load the count matrix -------------------------------------------
fc <- read.delim("counts/counts.txt", comment.char = "#")
counts <- as.matrix(fc[, 7:ncol(fc)])          # cols 1-6 are annotation
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("_Aligned.*", "", basename(colnames(counts)))

# ---- 2. Experimental groups ---------------------------------------------
# featureCounts orders columns alphabetically: ctrl_1..3 then mip6_1..3.
group <- factor(c("control", "control", "control", "mip6", "mip6", "mip6"),
                levels = c("control", "mip6"))
stopifnot(length(group) == ncol(counts))

# ---- 3. edgeR differential expression -----------------------------------
y <- DGEList(counts = counts, group = group)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)                              # mip6 vs control
res <- topTags(et, n = Inf)$table
res$gene <- rownames(res)
write.csv(res, "counts/de_results.csv", row.names = FALSE)
cat("Wrote counts/de_results.csv (", nrow(res), " genes )\n", sep = "")

dir.create("figures", showWarnings = FALSE)

# ---- 4. Volcano plot ----------------------------------------------------
res$sig <- with(res, ifelse(FDR < 0.05 & abs(logFC) > 1,
                            ifelse(logFC > 0, "up", "down"), "ns"))
ggplot(res, aes(logFC, -log10(FDR), colour = sig)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_colour_manual(values = c(up = "#b2182b", down = "#2166ac", ns = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "mip6 KO vs heat-shock control",
       x = "log2 fold change", y = "-log10 FDR", colour = NULL) +
  theme_bw()
ggsave("figures/volcano.png", width = 7, height = 5, dpi = 150)

# ---- 5. Heatmap of the top differentially expressed genes ---------------
top <- head(rownames(res[order(res$FDR), ]), 40)
logcpm <- cpm(y, log = TRUE)[top, ]
ann <- data.frame(group = group)
rownames(ann) <- colnames(logcpm)
pheatmap(logcpm, scale = "row", annotation_col = ann,
         show_rownames = FALSE, fontsize = 8,
         filename = "figures/heatmap.png", width = 7, height = 8)

cat("Wrote figures/volcano.png and figures/heatmap.png\n")
