## =============================================================================
## 02_diversity.R
## Alpha and beta diversity of the aging mouse model
##
## Part 1. Alpha diversity (Observed genera, Shannon, Pielou)
##           1.3 Abx+ mice over time        1.5 Abx+ groups at Week 8
##           1.7 Abx- mice over time        1.9 Abx- groups at Week 4
## Part 2. Beta diversity (Bray-Curtis, PERMANOVA, PCoA)
##           2.1 Abx+ over time             2.2 Abx+ groups at Week 8
##           2.3 Abx- over time             2.4 Abx- groups at Week 4
##
## Input (written by 01_dada2_phyloseq.R):
##   data/oturare_in_data2_to_phyloseq.rds   genus-level table rarefied to 3,500
##   data/Aging_Miseq_Sample_Information.csv sample metadata
##
## Statistical conventions
##   - Post hoc Dunn tests are reported only for indices whose omnibus
##     Kruskal-Wallis test is significant (p < 0.05). The full Dunn table is
##     exported for reference but is not used for the figures.
##   - Every PERMANOVA uses 999 permutations under a fixed seed (SEED_PERM).
## =============================================================================

library(here)
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(openxlsx)


## ---- Parameters -------------------------------------------------------------
DATA_DIR   <- here("data")
OUT_DIR    <- here("results")
SEED_PERM  <- 123     # seed for every adonis2 call
N_PERM     <- 999     # permutations for every adonis2 call
ALPHA_SIG  <- 0.05    # BH-adjusted omnibus threshold gating the Dunn post hoc

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

diversity_indices <- c("Observed", "Shannon", "Pielou")

group_colors <- c(
  "PBS"   = "#8f8f8f",
  "D-Gal" = "#151515",
  "Live"  = "#539feb",
  "HK"    = "#544ad5",
  "Sup"   = "#f1901b",
  "RCM"   = "#e23b4a"
)
time_colors_abx_plus <- c("Abx 0"  = "#e23b4a",
                          "Abx 7"  = "#FFB366",
                          "Week 8" = "#4CAF50")
time_colors_abx_minus <- c("Week 0" = "#9B59B6",
                           "Week 4" = "#E67E22")
group_levels <- c("PBS", "D-Gal", "Live", "HK", "Sup", "RCM")


## ---- Input ------------------------------------------------------------------
oturare <- readRDS(file.path(DATA_DIR, "oturare_in_data2_to_phyloseq.rds"))
a <- read.csv(file.path(DATA_DIR, "Aging_Miseq_Sample_Information.csv"),
              row.names = 1, check.names = FALSE)


## =============================================================================
## PART 1. ALPHA DIVERSITY
## =============================================================================

## ---- 1.1 Metadata: rename columns and groups --------------------------------
colnames(a)[colnames(a) == "Abx+/-"] <- "Abx_yes_or_no"

a$Group <- dplyr::recode(
  a$Group,
  "NC"      = "PBS",
  "PC"      = "D-Gal",
  "RI Live" = "Live",
  "RI Dead" = "HK",
  "RI Sup"  = "Sup",
  "RCM"     = "RCM"
)

meta_phyloseq <- sample_data(a)
rownames(meta_phyloseq) <- meta_phyloseq$Miseq_Sample_ID


## ---- 1.2 Richness and Pielou's evenness -------------------------------------
alpha.diversity <- estimate_richness(
  oturare,
  measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
)

## oturare has taxa as rows; vegan expects samples as rows.
otu_samples_rows <- t(oturare)
H      <- diversity(otu_samples_rows)                    # Shannon
Pielou <- as.data.frame(H / log(specnumber(otu_samples_rows)))
colnames(Pielou) <- "Pielou"

alpha.diversity_df <- cbind(alpha.diversity, Pielou)
write.csv(alpha.diversity_df, file.path(OUT_DIR, "richness_mouse.csv"))

alpha_diversity_extended <- merge(alpha.diversity_df, meta_phyloseq, by = "row.names")
rownames(alpha_diversity_extended) <- alpha_diversity_extended$Row.names
alpha_diversity_extended$Row.names <- NULL


## ---- Helpers ----------------------------------------------------------------

## Long-format table of the three indices for one grouping variable.
make_alpha_long <- function(df, group_var) {
  df %>%
    dplyr::select(dplyr::all_of(c(group_var, diversity_indices))) %>%
    tidyr::pivot_longer(cols = -dplyr::all_of(group_var),
                        names_to = "Index", values_to = "Value") %>%
    dplyr::filter(!is.na(Value)) %>%
    dplyr::mutate(Index = factor(Index, levels = diversity_indices))
}

## Kruskal-Wallis across the three indices, with effect size.
run_kruskal <- function(alpha_long, group_var) {
  f <- stats::as.formula(paste("Value ~", group_var))
  kw <- alpha_long %>%
    dplyr::group_by(Index) %>%
    rstatix::kruskal_test(f) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p.adj.BH = p.adjust(p, method = "BH"),
                  Significant = p.adj.BH < ALPHA_SIG)
  eff <- alpha_long %>%
    dplyr::group_by(Index) %>%
    rstatix::kruskal_effsize(f) %>%
    dplyr::ungroup() %>%
    dplyr::select(Index, effsize, magnitude)
  dplyr::left_join(kw, eff, by = "Index")
}

## Dunn post hoc with BH correction across the pairwise comparisons of one
## index. Returns both the full table and the gated table; only the gated table
## is used for figures and for the reported statistics.
run_dunn <- function(alpha_long, group_var, kw_res) {
  f <- stats::as.formula(paste("Value ~", group_var))
  dunn_all <- alpha_long %>%
    dplyr::group_by(Index) %>%
    rstatix::dunn_test(f, p.adjust.method = "BH") %>%
    dplyr::ungroup()

  ## The gate uses the BH-adjusted omnibus p, matching the significance that is
  ## reported for the Kruskal-Wallis test. In this dataset the raw p gives the
  ## same gating for every comparison.
  sig_idx <- kw_res %>%
    dplyr::filter(p.adj.BH < ALPHA_SIG) %>%
    dplyr::pull(Index) %>%
    as.character()

  list(all    = dunn_all,
       gated  = dplyr::filter(dunn_all, as.character(Index) %in% sig_idx),
       sig_idx = sig_idx)
}

## Significant pairs for one index, with the y position of the bracket.
## dunn_tbl is passed explicitly so that a figure can never pick up the post hoc
## table of a different comparison.
get_dunn_labels <- function(data, y_var, dunn_tbl) {
  if (is.null(dunn_tbl) || !"p.adj" %in% colnames(dunn_tbl)) return(NULL)

  d <- dunn_tbl %>%
    dplyr::filter(as.character(Index) == y_var, p.adj < 0.05) %>%
    dplyr::mutate(label = "*")
  if (nrow(d) == 0) return(NULL)

  y_max   <- max(data[[y_var]], na.rm = TRUE)
  y_range <- diff(range(data[[y_var]], na.rm = TRUE))
  d$y.position <- y_max + y_range * (0.08 + 0.12 * (seq_len(nrow(d)) - 1))
  d
}

## Boxplot of one index against a grouping variable.
make_box_plot <- function(data, y_var, x_var, palette, dunn_tbl = NULL,
                          title = NULL, x_lab = NULL, box_width = 0.4) {
  dunn_lab <- get_dunn_labels(data, y_var, dunn_tbl)

  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        color = .data[[x_var]])) +
    geom_boxplot(width = box_width, alpha = 0, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    scale_color_manual(values = palette) +
    labs(title = title, x = x_lab, y = y_var) +
    theme_classic() +
    theme(
      plot.title  = element_text(family = "Helvetica", hjust = 0.5, size = 18, color = "black"),
      legend.position = "none",
      axis.text.x = element_text(family = "Helvetica", size = 18, color = "black"),
      axis.text.y = element_text(family = "Helvetica", size = 18, color = "black"),
      axis.title  = element_text(family = "Helvetica", size = 18, color = "black"),
      axis.line   = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
    )

  if (!is.null(dunn_lab)) {
    p <- p + stat_pvalue_manual(dunn_lab, label = "label",
                                y.position = "y.position",
                                tip.length = 0.01, bracket.size = 0.6, size = 7)
  }
  p + scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))
}

add_sheet <- function(wb, name, df) {
  hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(df), headerStyle = hs)
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
  freezePane(wb, name, firstRow = TRUE)
}

## One-way PERMANOVA on a Bray-Curtis distance, always with the same seed and
## permutation count.
run_adonis <- function(dist_obj, metadata, term) {
  f <- stats::as.formula(paste("dist_obj ~", term))
  set.seed(SEED_PERM)
  vegan::adonis2(f, data = metadata, permutations = N_PERM)
}


## ---- 1.3 Alpha diversity over time in Abx+ mice: statistics -----------------
alpha_abx_plus <- alpha_diversity_extended %>%
  dplyr::filter(Abx_yes_or_no == "Abx+") %>%
  mutate(Time = factor(Time, levels = c("Abx 0", "Abx 7", "Week 8")))

alpha_long_abx_plus_time <- make_alpha_long(alpha_abx_plus, "Time")
kw_abx_plus_time   <- run_kruskal(alpha_long_abx_plus_time, "Time")
dunn_abx_plus_time <- run_dunn(alpha_long_abx_plus_time, "Time", kw_abx_plus_time)

print(kw_abx_plus_time)
print(dunn_abx_plus_time$gated)

wb <- createWorkbook()
add_sheet(wb, "Kruskal_Wallis", kw_abx_plus_time)
add_sheet(wb, "Dunn_posthoc",   dunn_abx_plus_time$gated)
add_sheet(wb, "Dunn_all",       dunn_abx_plus_time$all)
saveWorkbook(wb, file.path(OUT_DIR, "alpha_diversity_Abx_plus_time_stats.xlsx"),
             overwrite = TRUE)


## ---- 1.4 Alpha diversity over time in Abx+ mice: figure ---------------------
plots_abx_plus_time <- lapply(
  setNames(diversity_indices, diversity_indices),
  function(idx) make_box_plot(alpha_abx_plus, idx, "Time",
                              time_colors_abx_plus,
                              dunn_tbl = dunn_abx_plus_time$gated,
                              title = "Abx+", x_lab = "Time Point")
)

combined_abx_plus_time <- arrangeGrob(
  plots_abx_plus_time$Observed, plots_abx_plus_time$Shannon, plots_abx_plus_time$Pielou,
  ncol = 3
)
ggsave(file.path(OUT_DIR, "abx_plus_alpha_diversity_Time.svg"),
       combined_abx_plus_time, width = 12, height = 4, dpi = 600)


## ---- 1.5 Alpha diversity across groups at Week 8 (Abx+): statistics ---------
## Week 8 samples belong to the Abx+ arm only.
alpha_wk8 <- alpha_diversity_extended %>%
  dplyr::filter(Time == "Week 8") %>%
  mutate(Group = factor(Group, levels = group_levels))

alpha_long_wk8 <- make_alpha_long(alpha_wk8, "Group")
kw_wk8   <- run_kruskal(alpha_long_wk8, "Group")
dunn_wk8 <- run_dunn(alpha_long_wk8, "Group", kw_wk8)

print(kw_wk8)
print(dunn_wk8$gated)

wb <- createWorkbook()
add_sheet(wb, "Kruskal_Wallis", kw_wk8)
add_sheet(wb, "Dunn_posthoc",   dunn_wk8$gated)
add_sheet(wb, "Dunn_all",       dunn_wk8$all)
saveWorkbook(wb, file.path(OUT_DIR, "alpha_diversity_Abx_plus_groups_w8_stats.xlsx"),
             overwrite = TRUE)


## ---- 1.6 Alpha diversity across groups at Week 8: figure --------------------
plots_wk8 <- lapply(
  setNames(diversity_indices, diversity_indices),
  function(idx) make_box_plot(alpha_wk8, idx, "Group", group_colors,
                              dunn_tbl = dunn_wk8$gated,
                              x_lab = NULL, box_width = 0.8)
)

combined_wk8 <- arrangeGrob(
  plots_wk8$Observed, plots_wk8$Shannon, plots_wk8$Pielou,
  ncol = 3,
  top = textGrob("Abx+ (Week 8)", gp = gpar(fontsize = 18, fontface = "plain"))
)
ggsave(file.path(OUT_DIR, "abx_plus_diversity_groups_week8.svg"),
       combined_wk8, width = 16, height = 6, dpi = 300)


## ---- 1.7 Alpha diversity over time in Abx- mice: statistics -----------------
## Two timepoints, so a Wilcoxon rank-sum test rather than Kruskal-Wallis.
alpha_abx_minus <- alpha_diversity_extended %>%
  dplyr::filter(Abx_yes_or_no == "Abx-") %>%
  mutate(Time = factor(Time, levels = c("Week 0", "Week 4")))

wilcoxon_two_timepoints <- function(data, y_var,
                                    group1 = "Week 0", group2 = "Week 4") {
  x <- data[data$Time == group1, y_var]
  y <- data[data$Time == group2, y_var]
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) >= 2 && length(y) >= 2) {
    res <- wilcox.test(x, y)
    list(comparison = list(c(group1, group2)),
         statistic  = as.numeric(res$statistic),
         p_value    = res$p.value,
         signif     = ifelse(res$p.value < 0.05, "*", "ns"))
  } else {
    list(comparison = list(), statistic = NA, p_value = NA, signif = "ns")
  }
}

wilcoxon_abx_minus <- do.call(rbind, lapply(diversity_indices, function(index) {
  res <- wilcoxon_two_timepoints(alpha_abx_minus, index)
  cat(index, ": W =", ifelse(is.na(res$statistic), "NA", round(res$statistic, 6)),
      ", p =", ifelse(is.na(res$p_value), "NA", round(res$p_value, 6)),
      "(", res$signif, ")\n")
  data.frame(Index = index,
             Statistic = round(res$statistic, 6),
             P_value   = round(res$p_value, 6),
             Significance = res$signif,
             stringsAsFactors = FALSE)
}))

write.xlsx(wilcoxon_abx_minus,
           file.path(OUT_DIR, "Wilcoxon_alpha_Abx_minus_time.xlsx"))


## ---- 1.8 Alpha diversity over time in Abx- mice: figure ---------------------
make_box_plot_wilcoxon <- function(data, y_var) {
  test_info <- wilcoxon_two_timepoints(data, y_var)

  p <- ggplot(data, aes(x = .data[["Time"]], y = .data[[y_var]],
                        color = .data[["Time"]])) +
    geom_boxplot(width = 0.25, alpha = 0, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    scale_color_manual(values = time_colors_abx_minus) +
    labs(title = "Abx-", x = "Time Point", y = y_var) +
    theme_classic() +
    theme(
      plot.title  = element_text(family = "Helvetica", hjust = 0.5, size = 18,
                                 face = "plain", color = "black"),
      legend.position = "none",
      axis.text.x = element_text(family = "Helvetica", angle = 0, hjust = 0.5,
                                 size = 18, color = "black"),
      axis.text.y = element_text(family = "Helvetica", size = 18, color = "black"),
      axis.title  = element_text(family = "Helvetica", size = 18, color = "black"),
      axis.line   = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
    )

  if (!is.na(test_info$p_value) && test_info$signif == "*") {
    p <- p + stat_compare_means(
      comparisons = test_info$comparison, method = "wilcox.test",
      label = "p.signif", hide.ns = TRUE, bracket.size = 0.6,
      tip.length = 0.01, step.increase = 0.3, vjust = -0.5,
      symnum.args = list(cutpoints = c(0, 0.05, 1), symbols = c("*", " "))
    )
  }
  p + scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))
}

plots_abx_minus_time <- lapply(setNames(diversity_indices, diversity_indices),
                               function(idx) make_box_plot_wilcoxon(alpha_abx_minus, idx))

combined_abx_minus_time <- arrangeGrob(
  plots_abx_minus_time$Observed, plots_abx_minus_time$Shannon, plots_abx_minus_time$Pielou,
  ncol = 3,
  top = textGrob("Abx-", gp = gpar(fontsize = 18, fontface = "plain"))
)
ggsave(file.path(OUT_DIR, "abx_minus_alpha_diversity_Time.svg"),
       combined_abx_minus_time, width = 12, height = 6, dpi = 600)


## ---- 1.9 Alpha diversity across groups at Week 4 (Abx-): statistics ---------
## Week 4 samples belong to the Abx- arm only.
alpha_wk4 <- alpha_diversity_extended %>%
  dplyr::filter(Time == "Week 4") %>%
  mutate(Group = factor(Group, levels = group_levels))

alpha_long_wk4 <- make_alpha_long(alpha_wk4, "Group")
kw_wk4   <- run_kruskal(alpha_long_wk4, "Group")
dunn_wk4 <- run_dunn(alpha_long_wk4, "Group", kw_wk4)

print(kw_wk4)
print(dunn_wk4$gated)

wb <- createWorkbook()
add_sheet(wb, "Kruskal_Wallis", kw_wk4)
add_sheet(wb, "Dunn_posthoc",   dunn_wk4$gated)
add_sheet(wb, "Dunn_all",       dunn_wk4$all)
saveWorkbook(wb, file.path(OUT_DIR, "alpha_diversity_Abx_minus_groups_w4_stats.xlsx"),
             overwrite = TRUE)


## ---- 1.10 Alpha diversity across groups at Week 4: figure -------------------
plots_wk4 <- lapply(
  setNames(diversity_indices, diversity_indices),
  function(idx) make_box_plot(alpha_wk4, idx, "Group", group_colors,
                              dunn_tbl = dunn_wk4$gated,
                              x_lab = NULL, box_width = 0.8)
)

combined_wk4 <- arrangeGrob(
  plots_wk4$Observed, plots_wk4$Shannon, plots_wk4$Pielou,
  ncol = 3,
  top = textGrob("Abx- (Week 4)", gp = gpar(fontsize = 18, fontface = "plain"))
)
ggsave(file.path(OUT_DIR, "abx_minus_alpha_diversity_group_week4.svg"),
       combined_wk4, width = 16, height = 6, dpi = 600)


## =============================================================================
## PART 2. BETA DIVERSITY
## =============================================================================

## oturare_t: rows = taxa, columns = samples (the orientation used for subsetting
## samples); each analysis transposes it back before calling vegdist.
oturare_t <- t(otu_samples_rows)

## PCoA scatter plot with 95% confidence ellipses and a PERMANOVA annotation.
make_pcoa_plot <- function(pcoa, df, colour_var, palette, subtitle,
                           annotation_text, annot_x, annot_y,
                           annot_w, annot_h, reverse_x = FALSE) {

  p <- ggplot(df, aes(x = PCoA1, y = PCoA2, color = .data[[colour_var]])) +
    geom_point(aes(fill = .data[[colour_var]]),
               size = 4, alpha = 0.8, shape = 21, color = "white", stroke = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.8, alpha = 0.6, linetype = "solid") +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(
      title = " ",
      subtitle = subtitle,
      x = paste0("PCoA1 (", round(pcoa$eig[1] / sum(pcoa$eig) * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(pcoa$eig[2] / sum(pcoa$eig) * 100, 1), "%)"),
      color = colour_var, fill = colour_var
    ) +
    theme_minimal(base_size = 18, base_family = "Helvetica") +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.8),
      panel.grid.minor = element_line(color = "white", linewidth = 0.4),
      plot.title    = element_text(face = "plain", size = 18, hjust = 0.5,
                                   margin = margin(b = 5), color = "#000000"),
      plot.subtitle = element_text(size = 18, hjust = 0.5,
                                   margin = margin(b = 15), color = "#000000",
                                   face = "plain"),
      axis.title = element_text(size = 18, face = "plain", color = "#000000"),
      axis.text  = element_text(size = 18, color = "#000000"),
      axis.line  = element_line(color = "#000000", linewidth = 0.5),
      legend.title      = element_text(size = 16, face = "plain", color = "#000000"),
      legend.text       = element_text(size = 16, color = "#000000"),
      legend.background = element_rect(fill = "white", color = "#E8E8E8", linewidth = 0.5),
      legend.key        = element_rect(fill = "transparent"),
      legend.position   = "bottom",
      legend.margin     = margin(t = 15),
      panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
      plot.margin  = margin(20, 20, 20, 20)
    ) +
    scale_y_continuous(breaks = pretty_breaks(n = 6), labels = number_format(accuracy = 0.01))

  p <- if (reverse_x) p + scale_x_reverse() else
    p + scale_x_continuous(breaks = pretty_breaks(n = 6),
                           labels = number_format(accuracy = 0.01))

  annotation_bg <- rectGrob(
    x = annot_x, y = annot_y, width = annot_w, height = annot_h,
    hjust = 0, vjust = 1,
    gp = gpar(fill = "white", col = "#E8E8E8", alpha = 0.9, lwd = 1)
  )
  annotation_grob <- textGrob(
    annotation_text, x = annot_x, y = annot_y, hjust = 0, vjust = 1,
    gp = gpar(fontsize = 16, fontface = "plain", col = "#000000",
              fontfamily = "Helvetica")
  )

  p +
    annotation_custom(annotation_bg,   xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    annotation_custom(annotation_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}


## ---- 2.1 Across time points in Abx+ mice ------------------------------------
metadata_abx_plus <- meta_phyloseq[meta_phyloseq$Abx_yes_or_no == "Abx+", ]
common_samples    <- dplyr::intersect(colnames(oturare_t), rownames(metadata_abx_plus))

countable_abx_plus <- oturare_t[, common_samples]
metadata_abx_plus  <- meta_phyloseq[common_samples, ]
stopifnot(all(colnames(countable_abx_plus) == rownames(metadata_abx_plus)))

beta_div_abx_plus <- vegdist(t(countable_abx_plus), method = "bray")

metadata_abx_plus <- as.data.frame(as.matrix(metadata_abx_plus), stringsAsFactors = FALSE)
metadata_abx_plus$Time <- factor(metadata_abx_plus$Time,
                                 levels = c("Abx 0", "Abx 7", "Week 8"))

adonis_abx_plus <- run_adonis(beta_div_abx_plus, metadata_abx_plus, "Time")
print(adonis_abx_plus)
write.xlsx(adonis_abx_plus, file.path(OUT_DIR, "adonis_abx_plus_Time.xlsx"))

p_val  <- adonis_abx_plus$`Pr(>F)`[1]
R2_val <- adonis_abx_plus$R2[1]

## Pairwise comparisons between timepoints, FDR-corrected across the three pairs
group_pairs <- list(c("Abx 0", "Abx 7"),
                    c("Abx 0", "Week 8"),
                    c("Abx 7", "Week 8"))

pairwise_results <- lapply(group_pairs, function(pair) {
  metadata_pair <- metadata_abx_plus[metadata_abx_plus$Time %in% pair, ]
  beta_pair     <- vegdist(t(oturare_t[, rownames(metadata_pair)]), method = "bray")
  res           <- run_adonis(beta_pair, metadata_pair, "Time")
  data.frame(group1 = pair[1], group2 = pair[2],
             F_model = res$F[1], R2 = res$R2[1], p_value = res$`Pr(>F)`[1])
})

pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adjust <- p.adjust(pairwise_df$p_value, method = "fdr")
print(pairwise_df)
write.xlsx(pairwise_df, file.path(OUT_DIR, "pairwise_Time_abx_plus.xlsx"))

pcoa_abx_plus <- cmdscale(beta_div_abx_plus, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  Sample = rownames(as.matrix(beta_div_abx_plus)),
  PCoA1  = pcoa_abx_plus$points[, 1],
  PCoA2  = pcoa_abx_plus$points[, 2],
  Time   = metadata_abx_plus$Time
)

pair_texts <- apply(pairwise_df, 1, function(row) {
  paste0(row["group1"], "/", row["group2"], ": p = ",
         signif(as.numeric(row["p_adjust"]), 3))
})
annotation_text <- paste0("PERMANOVA (Time)\n",
                          "R\u00b2 = ", round(R2_val * 100, 1), "%, ",
                          "p = ", signif(p_val, 3), "\n",
                          paste(pair_texts, collapse = "\n"))

final_plot <- make_pcoa_plot(pcoa_abx_plus, pcoa_df, "Time", time_colors_abx_plus,
                             "Abx+ ", annotation_text,
                             annot_x = 0.64, annot_y = 0.94,
                             annot_w = 0.40, annot_h = 0.40)
ggsave(file.path(OUT_DIR, "abx_plus_beta_diversity_Time.svg"),
       final_plot, width = 12, height = 6, dpi = 600)


## ---- 2.2 Across groups at Week 8 (Abx+ mice) --------------------------------
metadata_abx_plus_w8 <- meta_phyloseq[meta_phyloseq$Time == "Week 8", ]
common_samples       <- dplyr::intersect(colnames(oturare_t), rownames(metadata_abx_plus_w8))

countable_abx_plus_w8 <- oturare_t[, common_samples]
metadata_abx_plus_w8  <- meta_phyloseq[common_samples, ]
stopifnot(all(colnames(countable_abx_plus_w8) == rownames(metadata_abx_plus_w8)))

beta_div_abx_plus_w8 <- vegdist(t(countable_abx_plus_w8), method = "bray")

metadata_abx_plus_w8 <- data.frame(metadata_abx_plus_w8)
adonis_abx_plus_w8 <- run_adonis(beta_div_abx_plus_w8, metadata_abx_plus_w8, "Group")
print(adonis_abx_plus_w8)
write.xlsx(adonis_abx_plus_w8, file.path(OUT_DIR, "adonis_abx_plus_Group_w8.xlsx"))

p_val  <- adonis_abx_plus_w8$`Pr(>F)`[1]
R2_val <- adonis_abx_plus_w8$R2[1]

pcoa_abx_plus_w8 <- cmdscale(beta_div_abx_plus_w8, eig = TRUE, k = 2)
pcoa_df_w8 <- data.frame(
  Sample = rownames(as.matrix(beta_div_abx_plus_w8)),
  PCoA1  = pcoa_abx_plus_w8$points[, 1],
  PCoA2  = pcoa_abx_plus_w8$points[, 2],
  Group  = metadata_abx_plus_w8$Group
)

annotation_text <- paste0("PERMANOVA (Group)\n",
                          "R\u00b2 = ", round(R2_val * 100, 1), "%, ",
                          "p = ", signif(p_val, 3))

final_plot_w8 <- make_pcoa_plot(pcoa_abx_plus_w8, pcoa_df_w8, "Group", group_colors,
                                "Abx+ (Week 8)", annotation_text,
                                annot_x = 0.90, annot_y = 0.15,
                                annot_w = 0.075, annot_h = 0.08)
ggsave(file.path(OUT_DIR, "abx_plus_beta_diversity_group_week8.svg"),
       final_plot_w8, width = 16, height = 6, dpi = 600)


## ---- 2.3 Across time points in Abx- mice ------------------------------------
metadata_abx_minus <- meta_phyloseq[meta_phyloseq$Abx_yes_or_no == "Abx-", ]
common_samples     <- dplyr::intersect(colnames(oturare_t), rownames(metadata_abx_minus))

countable_abx_minus <- oturare_t[, common_samples]
metadata_abx_minus  <- meta_phyloseq[common_samples, ]
stopifnot(all(colnames(countable_abx_minus) == rownames(metadata_abx_minus)))

beta_div_abx_minus <- vegdist(t(countable_abx_minus), method = "bray")

metadata_abx_minus <- as.data.frame(as.matrix(metadata_abx_minus), stringsAsFactors = FALSE)
metadata_abx_minus$Time <- factor(metadata_abx_minus$Time,
                                  levels = c("Week 0", "Week 4"))

adonis_abx_minus <- run_adonis(beta_div_abx_minus, metadata_abx_minus, "Time")
print(adonis_abx_minus)
write.xlsx(adonis_abx_minus, file.path(OUT_DIR, "adonis_abx_minus_Time.xlsx"))

p_val  <- adonis_abx_minus$`Pr(>F)`[1]
R2_val <- adonis_abx_minus$R2[1]

pcoa_abx_minus <- cmdscale(beta_div_abx_minus, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  Sample = rownames(as.matrix(beta_div_abx_minus)),
  PCoA1  = pcoa_abx_minus$points[, 1],
  PCoA2  = pcoa_abx_minus$points[, 2],
  Time   = metadata_abx_minus$Time
)

annotation_text <- paste0("PERMANOVA (Time)\n",
                          "R\u00b2 = ", round(R2_val * 100, 1), "%, ",
                          "p = ", signif(p_val, 3))

final_plot_abx_minus_time <- make_pcoa_plot(
  pcoa_abx_minus, pcoa_df, "Time", time_colors_abx_minus,
  "Abx-", annotation_text,
  annot_x = 0.90, annot_y = 0.15, annot_w = 0.075, annot_h = 0.08
)
ggsave(file.path(OUT_DIR, "abx_minus_beta_diversity_Time.svg"),
       final_plot_abx_minus_time, width = 12, height = 6, dpi = 600)


## ---- 2.4 Across groups at Week 4 (Abx- mice) --------------------------------
metadata_abx_minus_w4 <- meta_phyloseq[meta_phyloseq$Time == "Week 4", ]
common_samples        <- dplyr::intersect(colnames(oturare_t), rownames(metadata_abx_minus_w4))

countable_abx_minus_w4 <- oturare_t[, common_samples]
metadata_abx_minus_w4  <- meta_phyloseq[common_samples, ]
stopifnot(all(colnames(countable_abx_minus_w4) == rownames(metadata_abx_minus_w4)))

beta_div_abx_minus_w4 <- vegdist(t(countable_abx_minus_w4), method = "bray")

metadata_abx_minus_w4 <- data.frame(metadata_abx_minus_w4)
adonis_abx_minus_w4 <- run_adonis(beta_div_abx_minus_w4, metadata_abx_minus_w4, "Group")
print(adonis_abx_minus_w4)
write.xlsx(adonis_abx_minus_w4, file.path(OUT_DIR, "adonis_abx_minus_Group_w4.xlsx"))

p_val  <- adonis_abx_minus_w4$`Pr(>F)`[1]
R2_val <- adonis_abx_minus_w4$R2[1]

pcoa_abx_minus_w4 <- cmdscale(beta_div_abx_minus_w4, eig = TRUE, k = 2)
pcoa_df_w4 <- data.frame(
  Sample = rownames(as.matrix(beta_div_abx_minus_w4)),
  PCoA1  = pcoa_abx_minus_w4$points[, 1],
  PCoA2  = pcoa_abx_minus_w4$points[, 2],
  Group  = metadata_abx_minus_w4$Group
)

annotation_text <- paste0("PERMANOVA (Group)\n",
                          "R\u00b2 = ", round(R2_val * 100, 1), "%, ",
                          "p = ", signif(p_val, 3))

## The x axis is reversed so that the panel orientation matches Fig. 2C.
final_plot_w4 <- make_pcoa_plot(
  pcoa_abx_minus_w4, pcoa_df_w4, "Group", group_colors,
  "Abx- (Week 4)", annotation_text,
  annot_x = 0.90, annot_y = 0.15, annot_w = 0.075, annot_h = 0.08,
  reverse_x = TRUE
)
ggsave(file.path(OUT_DIR, "abx_minus_beta_diversity_Group_week4.svg"),
       final_plot_w4, width = 16, height = 6, dpi = 600)


## ---- Save objects used by the downstream scripts ----------------------------
saveRDS(meta_phyloseq, file.path(DATA_DIR, "meta_phyloseq_metadata.rds"))

sink(file.path(OUT_DIR, "sessionInfo_02.txt")); print(sessionInfo()); sink()
