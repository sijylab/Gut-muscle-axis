## =============================================================================
## 04_differential_abundance.R
## Differential abundance of genera across treatment groups (ALDEx2)
##   Abx+ mice at Week 8  and  Abx- mice at Week 4
##
## Workflow (Week 8)
##   Step 0. Data: genus count table, Week 8 samples
##   Step 1. ALDEx2 Kruskal-Wallis test across the six groups (aldex.kw)
##   Step 2. Select genera with BH-adjusted KW p < 0.05
##   Step 3. Pairwise comparisons (aldex.ttest + aldex.effect) for those genera
##   Step 4. Summary of significant pairwise comparisons (wi.eBH <= 0.05)
##   Step 5. Supplementary table of significant pairwise comparisons
##   Step 6. Boxplots of CLR values
##
## Week 4 repeats Steps 0, 1 and 6. No genus reached BH-adjusted significance in
## the Abx- arm, so the pairwise steps are not run; the figure therefore shows
## the genera with unadjusted KW p < 0.05, as stated in its legend.
##
## Reproduces: 37 significant pairwise comparisons (as reported in the manuscript)
##
## Parameters: ALDEx2 1.38.0, mc.samples = 128, denom = "all", seed = 123 (KW).
##   - Pairwise CLR is computed on the genera selected after the KW test
##     (PAIRWISE_FEATURES = "selected").
##   - SEED_PAIRWISE is NULL because no seed was set for the pairwise step in the
##     original run. Set it to 123 to make that step exactly reproducible; the
##     Monte Carlo p-values then differ marginally from the published ones.
##   - g_Staphylococcus (probable contaminant) is excluded from the reported
##     table (Step 5) and from the figure (Step 6) only.
##
## Two fixes were required to run with ALDEx2 1.38.0 ([FIX 1], [FIX 2] below):
##   [FIX 1] the group vector in the pairwise step is a character vector
##   [FIX 2] rab.win.<group> columns are renamed so tables can be row-bound
##
## Input (written by earlier scripts):
##   data/otutable_in_data2_to_phyloseq.rds   (01_dada2_phyloseq.R)
##   data/aging_model_w8.rds, aging_model_w4.rds  (03_correlation.R)
## =============================================================================

library(here)
library(ALDEx2)
library(dplyr)
library(tibble)
library(openxlsx)


## ---- Parameters --------------------------------------------------------------
DATA_DIR   <- here("data")
OUT_DIR    <- here("results")
MC_SAMPLES <- 128
SEED_KW    <- 123
PAIRWISE_FEATURES <- "selected"   # "selected": pairwise CLR only on the genera
                                  #   selected after the KW test (as submitted)
                                  # "all":      pairwise CLR on all genera,
                                  #   filtered afterwards
SEED_PAIRWISE <- NULL             # see header
CONTAMINANT   <- "g_Staphylococcus"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

group_levels <- c("PBS", "D-Gal", "Live", "HK", "Sup", "RCM")

# Colours for groups (matching GraphPad figures)
group_colors <- c(
  "PBS"   = "#8f8f8f",
  "D-Gal" = "#151515",
  "Live"  = "#539feb",
  "HK"    = "#544ad5",
  "Sup"   = "#f1901b",
  "RCM"   = "#e23b4a"
)

otutable       <- readRDS(file.path(DATA_DIR, "otutable_in_data2_to_phyloseq.rds"))
aging_model_w8 <- readRDS(file.path(DATA_DIR, "aging_model_w8.rds"))
aging_model_w4 <- readRDS(file.path(DATA_DIR, "aging_model_w4.rds"))


## ---- Helpers ------------------------------------------------------------------

## Count table and group vector for one timepoint.
subset_timepoint <- function(counts, aging_model) {
  meta_t <- t(aging_model)
  common <- dplyr::intersect(colnames(counts), colnames(meta_t))
  counts_sub   <- counts[, common, drop = FALSE]
  group_vector <- as.character(meta_t["Group", common])
  names(group_vector) <- common
  stopifnot(identical(colnames(counts_sub), names(group_vector)))
  list(counts = counts_sub, groups = group_vector)
}

## Mean CLR value of one genus for each sample.
get_clr_abundance <- function(genus_name, clr_list, feature_names, sample_names) {
  genus_idx <- which(feature_names == genus_name)
  if (length(genus_idx) == 0) return(NULL)
  clr_values <- sapply(clr_list, function(mc_matrix) mc_matrix[genus_idx, ])
  clr_mean <- rowMeans(clr_values)
  names(clr_mean) <- sample_names
  clr_mean
}

## Strip the rank prefix from a genus label.
clean_genus_name <- function(genus) {
  genus <- gsub("^g_", "", genus)
  genus <- gsub("^uc_f_", "f.", genus)
  genus <- gsub("^uc_o_", "o.", genus)
  genus <- gsub("^uc_c_", "c.", genus)
  genus <- gsub("^uc_k_", "k.", genus)
  gsub("^uc_p_", "p.", genus)
}

## Base-graphics boxplot of CLR values for one genus.
draw_boxplot_clr <- function(clr_abundance, group_factor, title) {
  group_names <- levels(group_factor)
  y_range  <- range(clr_abundance, na.rm = TRUE)
  y_margin <- diff(y_range) * 0.55
  ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)

  plot(1, type = "n",
       xlim = c(0.5, length(group_names) + 0.5), ylim = ylim,
       xlab = "", ylab = "CLR Transformed", main = "",
       cex.lab = 2.1, family = "Helvetica", xaxt = "n", yaxt = "n")

  title(main = title, line = 2, cex.main = 2.1, font.main = 1, family = "Helvetica")
  grid(nx = 0, ny = NULL, col = "lightgray", lty = 3, lwd = 0.5)
  abline(h = 0, col = "red", lty = 2, lwd = 1)   # CLR = 0 is the geometric mean

  for (i in seq_along(group_names)) {
    group_name <- group_names[i]
    group_data <- clr_abundance[group_factor == group_name]
    group_data <- group_data[!is.na(group_data)]
    if (length(group_data) <= 1) next

    q   <- quantile(group_data, probs = c(0.25, 0.5, 0.75))
    iqr <- q[3] - q[1]
    lw  <- max(min(group_data), q[1] - 1.5 * iqr)
    uw  <- min(max(group_data), q[3] + 1.5 * iqr)
    bw  <- 0.4
    col <- group_colors[group_name]

    rect(i - bw, q[1], i + bw, q[3], col = NA, border = col, lwd = 1.2)
    segments(i - bw, q[2], i + bw, q[2], col = col, lwd = 2)
    segments(i, q[1], i, lw, col = col, lwd = 1.2)
    segments(i, q[3], i, uw, col = col, lwd = 1.2)
    segments(i - bw/2, lw, i + bw/2, lw, col = col, lwd = 1.2)
    segments(i - bw/2, uw, i + bw/2, uw, col = col, lwd = 1.2)

    outliers <- group_data[group_data < lw | group_data > uw]
    if (length(outliers) > 0) {
      points(rep(i, length(outliers)), outliers,
             pch = 21, bg = col, col = col, cex = 0.8, lwd = 1.5)
    }
    points(runif(length(group_data), i - 0.025, i + 0.025), group_data,
           pch = 21, bg = col, col = col, cex = 0.5, lwd = 1)
  }

  axis(1, at = seq_along(group_names), labels = FALSE)
  text(seq_along(group_names), par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
       labels = group_names, srt = 45, adj = 1, xpd = TRUE,
       cex = 2.1, family = "Helvetica")
  axis(2, las = 1, cex.axis = 2.1, family = "Helvetica")
}

## Panel of CLR boxplots written to one SVG.
save_clr_panel <- function(file, genera, x_clr, group_factor) {
  clr_list      <- getMonteCarloInstances(x_clr)
  feature_names <- getFeatureNames(x_clr)
  sample_names  <- getSampleIDs(x_clr)

  svg(file, width = 12, height = 8, pointsize = 12, family = "Helvetica")
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3, 7,
                  4, 5, 6, 8), nrow = 2, byrow = TRUE),
         widths = c(1, 1, 1, 0.6), heights = c(1, 1))
  par(mar = c(6.5, 4.5, 3.5, 1.5), family = "Helvetica")

  for (genus in genera) {
    clr_abundance <- get_clr_abundance(genus, clr_list, feature_names, sample_names)
    if (!is.null(clr_abundance)) {
      draw_boxplot_clr(clr_abundance, group_factor, clean_genus_name(genus))
    }
  }
}


## =============================================================================
## WEEK 8 (Abx+ mice)
## =============================================================================

## ---- Step 0. Data -------------------------------------------------------------
w8 <- subset_timepoint(otutable, aging_model_w8)
aldex2_genus_w8 <- w8$counts
group_vector_w8 <- w8$groups
table(group_vector_w8)


## ---- Step 1. Kruskal-Wallis test (ALDEx2) -------------------------------------
set.seed(SEED_KW)
x_clr_w8 <- aldex.clr(aldex2_genus_w8, group_vector_w8,
                      mc.samples = MC_SAMPLES, denom = "all")
kruskal_clr_w8 <- aldex.kw(x_clr_w8)
kruskal_clr_w8 <- kruskal_clr_w8[order(kruskal_clr_w8$kw.ep), ]
kruskal_clr_w8$Genus <- rownames(kruskal_clr_w8)

significant_kw_clr_w8 <- kruskal_clr_w8[kruskal_clr_w8$kw.ep < 0.05, ]
cat("Genera differing across groups (Week 8), unadjusted:",
    nrow(significant_kw_clr_w8), "\n")

significant_kw_clr8_FDR <- kruskal_clr_w8[kruskal_clr_w8$kw.eBH < 0.05, ]
significant_kw_clr8_FDR <- significant_kw_clr8_FDR[order(significant_kw_clr8_FDR$kw.ep), ]
cat("Genera differing across groups (Week 8), BH-adjusted:",
    nrow(significant_kw_clr8_FDR), "\n")
print(significant_kw_clr8_FDR)


## ---- Step 2. Select genera ----------------------------------------------------
significant_genera8 <- rownames(significant_kw_clr8_FDR)
significant_genera8 <- significant_genera8[significant_genera8 %in% rownames(aldex2_genus_w8)]
significant_genera8 <- sort(significant_genera8)
print(significant_genera8)


## ---- Step 3. Pairwise comparisons ---------------------------------------------
aldex_results <- list()

if (length(significant_genera8) > 0) {
  unique_groups <- unique(group_vector_w8)

  for (i in seq_len(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      g1 <- unique_groups[i]
      g2 <- unique_groups[j]
      cat("Comparing", g1, "vs", g2, "...\n")

      selected_samples <- group_vector_w8 %in% c(g1, g2)
      counts_subset <- aldex2_genus_w8[, selected_samples, drop = FALSE]
      if (PAIRWISE_FEATURES == "selected") {
        counts_subset <- counts_subset[rownames(counts_subset) %in% significant_genera8,
                                       , drop = FALSE]
      }

      # [FIX 1] character instead of factor(levels = c(g1, g2))
      group_subset <- as.character(group_vector_w8[selected_samples])

      if (all(table(group_subset) >= 2) && nrow(counts_subset) > 0) {
        tryCatch({
          if (!is.null(SEED_PAIRWISE)) set.seed(SEED_PAIRWISE)
          clr <- aldex.clr(counts_subset, group_subset,
                           mc.samples = MC_SAMPLES, denom = "all", verbose = FALSE)
          clr_tt     <- aldex.ttest(clr, paired.test = FALSE, verbose = FALSE)
          clr_effect <- aldex.effect(clr, verbose = FALSE)

          clr_out <- data.frame(Genus = rownames(clr_tt), clr_tt, clr_effect,
                                check.names = FALSE)

          # [FIX 2] standardise the rab.win.<group> column names
          rab_idx   <- grep("^rab\\.win\\.", colnames(clr_out))
          grp_names <- sub("^rab\\.win\\.", "", colnames(clr_out)[rab_idx])
          colnames(clr_out)[rab_idx] <- c("rab.win.grp1", "rab.win.grp2")
          clr_out$Group1 <- grp_names[1]
          clr_out$Group2 <- grp_names[2]

          clr_out <- clr_out[clr_out$Genus %in% significant_genera8, , drop = FALSE]
          clr_out$Comparison <- paste(g1, "vs", g2, sep = "_")
          aldex_results[[paste(g1, "vs", g2, sep = "_")]] <- clr_out
        }, error = function(e) {
          cat("Error in comparison", g1, "vs", g2, ":", conditionMessage(e), "\n")
        })
      } else {
        cat("Skipping", g1, "vs", g2, "- insufficient samples per group\n")
      }
    }
  }
  cat("Completed", length(aldex_results), "pairwise comparisons (Week 8)\n")
}


## ---- Step 4. Significant pairwise comparisons (wi.eBH <= 0.05) ----------------
all_significant_clr <- NULL

if (length(aldex_results) > 0) {
  significant_hits_clr <- list()

  for (name in names(aldex_results)) {
    sig_df <- aldex_results[[name]] %>%
      dplyr::filter(wi.eBH <= 0.05) %>%
      dplyr::mutate(
        Comparison   = name,
        Overall_KW_p = ifelse(Genus %in% significant_genera8,
                              significant_kw_clr8_FDR[Genus, "kw.eBH"], NA)
      ) %>%
      dplyr::arrange(wi.eBH)
    if (nrow(sig_df) > 0) significant_hits_clr[[name]] <- sig_df
  }

  if (length(significant_hits_clr) > 0) {
    all_significant_clr <- dplyr::bind_rows(significant_hits_clr)
    cat("Total significant pairwise comparisons (Week 8):",
        nrow(all_significant_clr), "\n")
    cat("Unique genera:", length(unique(all_significant_clr$Genus)), "\n")

    summary_table <- all_significant_clr %>%
      dplyr::group_by(Genus) %>%
      dplyr::summarise(
        n_significant_comparisons = dplyr::n(),
        min_pairwise_p  = min(wi.eBH),
        max_effect_size = max(abs(effect)),
        overall_kw_p    = dplyr::first(Overall_KW_p),
        .groups = "drop"
      ) %>%
      dplyr::arrange(overall_kw_p, min_pairwise_p)
    print(summary_table)
  } else {
    cat("No significant pairwise differences found (Week 8).\n")
  }
}


## ---- Step 5. Supplementary table ----------------------------------------------
if (!is.null(all_significant_clr)) {
  table_pairwise_w8 <- all_significant_clr %>%
    transmute(Variables = Comparison, Genus = Genus,
              we.ep = we.ep, we.eBH = we.eBH,
              wi.ep = wi.ep, wi.eBH = wi.eBH, overlap = overlap) %>%
    arrange(match(Variables, names(aldex_results)), wi.eBH) %>%
    filter(Genus != CONTAMINANT)

  write.xlsx(table_pairwise_w8,
             file.path(OUT_DIR, "Aldex2_pairwise_by_groups_Abx_plus_week8.xlsx"))
}
write.xlsx(kruskal_clr_w8,
           file.path(OUT_DIR, "Aldex2_kruskal_clr_Abx_plus_week8.xlsx"))


## ---- Step 6. CLR boxplots -----------------------------------------------------
genera_plot_w8  <- dplyr::setdiff(significant_genera8, CONTAMINANT)
group_factor_w8 <- factor(group_vector_w8, levels = group_levels)

save_clr_panel(file.path(OUT_DIR, "aldex2_clr_boxplots_w8.svg"),
               genera_plot_w8, x_clr_w8, group_factor_w8)


## =============================================================================
## WEEK 4 (Abx- mice)
## =============================================================================

## ---- Step 0. Data -------------------------------------------------------------
w4 <- subset_timepoint(otutable, aging_model_w4)
aldex2_genus_w4 <- w4$counts
group_vector_w4 <- w4$groups
table(group_vector_w4)


## ---- Step 1. Kruskal-Wallis test (ALDEx2) -------------------------------------
set.seed(SEED_KW)
x_clr_w4 <- aldex.clr(aldex2_genus_w4, group_vector_w4,
                      mc.samples = MC_SAMPLES, denom = "all")
kruskal_clr_w4 <- aldex.kw(x_clr_w4)
kruskal_clr_w4 <- kruskal_clr_w4[order(kruskal_clr_w4$kw.ep), ]
kruskal_clr_w4$Genus <- rownames(kruskal_clr_w4)
write.xlsx(kruskal_clr_w4,
           file.path(OUT_DIR, "Aldex2_kruskal_clr_Abx_minus_week4.xlsx"))

significant_kw_clr_w4 <- kruskal_clr_w4[kruskal_clr_w4$kw.ep < 0.05, ]
cat("Genera differing across groups (Week 4), unadjusted:",
    nrow(significant_kw_clr_w4), "\n")

significant_kw_clr4_FDR <- kruskal_clr_w4[kruskal_clr_w4$kw.eBH < 0.05, ]
significant_kw_clr4_FDR <- significant_kw_clr4_FDR[order(significant_kw_clr4_FDR$kw.ep), ]
cat("Genera differing across groups (Week 4), BH-adjusted:",
    nrow(significant_kw_clr4_FDR), "\n")
print(significant_kw_clr4_FDR)


## ---- Step 6. CLR boxplots -----------------------------------------------------
## No genus is significant after BH correction in the Abx- arm, so the figure
## shows the genera with unadjusted KW p < 0.05 (stated in the figure legend).
## The pairwise steps are therefore not run for Week 4.
significant_genera4 <- rownames(significant_kw_clr_w4)
significant_genera4 <- significant_genera4[significant_genera4 %in% rownames(aldex2_genus_w4)]
significant_genera4 <- sort(significant_genera4)

has_prefix_4 <- grepl("\\.", significant_genera4)
genera_plot_w4 <- c(sort(significant_genera4[has_prefix_4]),
                    sort(significant_genera4[!has_prefix_4]))
print(genera_plot_w4)

group_factor_w4 <- factor(group_vector_w4, levels = group_levels)

save_clr_panel(file.path(OUT_DIR, "aldex2_clr_boxplots_w4.svg"),
               genera_plot_w4, x_clr_w4, group_factor_w4)


sink(file.path(OUT_DIR, "sessionInfo_04.txt")); print(sessionInfo()); sink()
