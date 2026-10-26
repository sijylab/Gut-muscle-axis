## =============================================================================
## 03_correlation.R
## GMPR normalisation and Spearman correlation between genera and host
## phenotypes (Abx+ Week 8, Abx- Week 4)
##
## Part 1. GMPR size factors and normalised genus table
## Part 2. Genus-phenotype correlation, all mice at the endpoint
## Part 3. Genus-phenotype correlation within each treatment group
##
## Input (written by earlier scripts):
##   data/meta_phyloseq_metadata.rds                      (02_diversity.R)
##   data/mouse_silva_table_cleaned_genus_glom.csv        (01_dada2_phyloseq.R)
##   data/<prefix>_gmpr_wk8_dominant.csv / _wk4_dominant.csv
##     GMPR-normalised tables restricted to the dominant genera, merged with the
##     phenotype columns. See README.md for how these are produced.
##
## Significance on the heatmaps is marked at FDR < 0.1, corrected within each
## phenotype (Part 2) or within each Group x phenotype pair (Part 3). This
## threshold is stated in every figure title.
## =============================================================================

library(here)
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(openxlsx)


## ---- Parameters -------------------------------------------------------------
DATA_DIR  <- here("data")
OUT_DIR   <- here("results")
FDR_MARK  <- 0.1      # threshold used to mark cells on the heatmaps
MIN_N_ALL <- 3        # minimum complete pairs for the pooled correlation
MIN_N_GRP <- 5        # minimum complete pairs for the within-group correlation
FONT      <- "Helvetica"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

group_vector_levels <- c("PBS", "D-Gal", "Live", "HK", "Sup", "RCM")

## Genus label cleanup: drop the "g_" prefix and shorten unclassified ranks.
replace_patterns <- c(
  "^g_"    = "",
  "^uc_f_" = "f.",
  "^uc_o_" = "o.",
  "^uc_c_" = "c.",
  "^uc_k_" = "k.",
  "^uc_p_" = "p."
)

## Non-genus columns that may be present in the GMPR tables.
META_COLS <- c("Miseq_Sample_ID", "Index.Primer", "Animal.set",
               "Abx_yes_or_no", "Time", "Group")


## =============================================================================
## PART 1. GMPR SIZE FACTORS
## =============================================================================
## GMPR() comes from the reference implementation of Chen et al. (2018),
## vendored in R/lib/GMPR.R. Run this part once to regenerate
## mouse_silva_table_cleaned_genus_glom_gmpr.csv.

run_gmpr <- function() {
  source(here("R", "lib", "GMPR.R"))

  otu.tab <- read.csv(file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom.csv"),
                      row.names = 1, check.names = FALSE)
  otu.tab <- as.matrix(otu.tab)              # rows = genera, columns = samples

  gmpr.size.factor <- GMPR(otu.tab, intersect.no = 4)
  otu.tab.norm <- t(t(otu.tab) / gmpr.size.factor)

  write.csv(otu.tab.norm,
            file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom_gmpr.csv"))
  invisible(otu.tab.norm)
}

# run_gmpr()   # uncomment to recompute the GMPR-normalised table


## ---- Helpers ----------------------------------------------------------------

meta_phyloseq <- readRDS(file.path(DATA_DIR, "meta_phyloseq_metadata.rds"))

## Genus columns of a GMPR table: everything that is neither metadata nor a
## phenotype column, and that is numeric.
genus_columns <- function(tbl, clinical_vars) {
  candidates <- dplyr::setdiff(colnames(tbl), c(META_COLS, clinical_vars))
  candidates[sapply(tbl[, candidates, drop = FALSE], is.numeric)]
}

## Spearman correlation of every phenotype against every genus, pooled over all
## mice at one timepoint. FDR is corrected within each phenotype.
correlate_pooled <- function(tbl, clinical_vars, genus_vars) {
  dims <- list(clinical_vars, genus_vars)
  cor_results  <- matrix(NA, length(clinical_vars), length(genus_vars), dimnames = dims)
  pval_results <- cor_results

  for (cln in clinical_vars) {
    for (gen in genus_vars) {
      ok <- !is.na(tbl[[cln]]) & !is.na(tbl[[gen]])
      if (sum(ok) >= MIN_N_ALL && sd(tbl[[cln]][ok]) > 0 && sd(tbl[[gen]][ok]) > 0) {
        tmp <- suppressWarnings(
          cor.test(tbl[[cln]][ok], tbl[[gen]][ok], method = "spearman", exact = FALSE)
        )
        cor_results[cln, gen]  <- tmp$estimate
        pval_results[cln, gen] <- tmp$p.value
      }
    }
  }

  fdr_results <- pval_results
  for (cln in clinical_vars) {
    pv <- pval_results[cln, ]
    ok <- !is.na(pv)
    if (any(ok)) {
      fdr_results[cln, ok] <- p.adjust(pv[ok], method = "BH")
      cat("Phenotype:", cln,
          "| tests:", sum(ok),
          "| p range:", paste(signif(range(pv[ok]), 3), collapse = "-"),
          "| FDR range:", paste(signif(range(fdr_results[cln, ok]), 3), collapse = "-"),
          "\n")
    }
  }

  genus_clean <- str_replace_all(genus_vars, replace_patterns)
  colnames(cor_results) <- colnames(pval_results) <- colnames(fdr_results) <- genus_clean

  long <- expand.grid(Clinical_Var = clinical_vars, Genus = genus_clean,
                      stringsAsFactors = FALSE) %>%
    mutate(Correlation     = cor_results[cbind(Clinical_Var, Genus)],
           P_value         = pval_results[cbind(Clinical_Var, Genus)],
           FDR_By_Clinical = fdr_results[cbind(Clinical_Var, Genus)]) %>%
    filter(!is.na(P_value))

  list(cor = cor_results, p = pval_results, fdr = fdr_results, long = long)
}

## Heatmap of the pooled correlations: x = phenotype, y = genus.
plot_pooled_heatmap <- function(res, title) {
  cor_df <- reshape2::melt(res$cor, varnames = c("Clinical", "Genus"),
                           value.name = "Correlation")
  fdr_df <- reshape2::melt(res$fdr, varnames = c("Clinical", "Genus"),
                           value.name = "FDR")
  heatmap_data <- merge(cor_df, fdr_df, by = c("Clinical", "Genus"))
  heatmap_data$significant <- ifelse(!is.na(heatmap_data$FDR) & heatmap_data$FDR < FDR_MARK,
                                     "*", "")

  ggplot(heatmap_data, aes(x = Clinical, y = Genus)) +
    geom_tile(aes(fill = Correlation), color = "black", linewidth = 0.15) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#313695",
                         midpoint = 0, na.value = "black",
                         name = "Spearman\nCorrelation", limits = c(-1, 1)) +
    geom_text(aes(label = significant), size = 5, color = "black",
              fontface = "bold", family = FONT) +
    theme_minimal(base_family = FONT, base_size = 18) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1, size = 18,
                                  face = "plain", color = "black"),
      axis.text.y  = element_text(size = 17, face = "plain", color = "black"),
      axis.title   = element_text(size = 18, face = "plain", color = "black"),
      plot.title   = element_text(size = 18, hjust = 0.5, face = "plain", color = "black"),
      legend.title = element_text(size = 18, face = "plain", color = "black"),
      legend.text  = element_text(size = 18, color = "black"),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.85)
    ) +
    labs(title = paste0(title, "\n(*: FDR < ", FDR_MARK, ")"),
         x = "Host phenotype", y = "Genus")
}

## Spearman correlation within each treatment group. FDR is corrected within
## each Group x phenotype pair.
correlate_by_group <- function(tbl, meta_df, clinical_vars, genus_vars) {
  common_ids <- dplyr::intersect(rownames(meta_df), rownames(tbl))

  clin_df <- meta_df[common_ids, clinical_vars, drop = FALSE]
  clin_df[] <- lapply(clin_df, function(v) suppressWarnings(as.numeric(v)))

  dat <- cbind(
    data.frame(Group = factor(meta_df[common_ids, "Group"]), row.names = common_ids),
    clin_df[common_ids, , drop = FALSE],
    tbl[common_ids, genus_vars, drop = FALSE]
  )

  genus_map <- setNames(str_replace_all(genus_vars, replace_patterns), genus_vars)

  res_list <- list()
  k <- 1
  for (grp in levels(dat$Group)) {
    sub <- dat[dat$Group == grp, ]
    for (cln in clinical_vars) {
      for (gen in genus_vars) {
        x <- sub[[cln]]; y <- sub[[gen]]
        ok <- !is.na(x) & !is.na(y)
        if (sum(ok) >= MIN_N_GRP && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
          tmp <- suppressWarnings(
            cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
          )
          res_list[[k]] <- data.frame(
            Group = grp, Clinical = cln, Genus = unname(genus_map[gen]),
            n = sum(ok), rho = unname(tmp$estimate), p = tmp$p.value
          )
          k <- k + 1
        }
      }
    }
  }

  res_group <- bind_rows(res_list) %>%
    group_by(Group, Clinical) %>%
    mutate(FDR = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(sig_FDR = FDR < FDR_MARK)

  list(res = res_group, genus_map = genus_map)
}

## Heatmap of the within-group correlations: x = group, y = genus, one panel per
## phenotype.
plot_group_heatmap <- function(res_group, genus_map, clinical_vars, title) {
  unique_genus <- unique(unname(genus_map))
  has_prefix   <- grepl("\\.", unique_genus)
  all_genus    <- c(sort(unique_genus[has_prefix],  decreasing = TRUE),
                    sort(unique_genus[!has_prefix], decreasing = TRUE))

  plot_df <- res_group %>%
    mutate(Group    = factor(Group, levels = group_vector_levels),
           Clinical = factor(Clinical, levels = clinical_vars)) %>%
    complete(Group    = factor(group_vector_levels, levels = group_vector_levels),
             Clinical = factor(clinical_vars, levels = clinical_vars),
             Genus    = all_genus) %>%
    mutate(Group    = factor(Group, levels = group_vector_levels),
           Clinical = factor(Clinical, levels = clinical_vars),
           Genus    = factor(Genus, levels = rev(all_genus)))

  p <- ggplot(plot_df, aes(x = Group, y = Genus)) +
    geom_tile(aes(fill = rho), color = "black", linewidth = 0.15) +
    scale_fill_gradient2(low = "#313695", mid = "white", high = "#d73027",
                         midpoint = 0, na.value = "grey85",
                         name = "Spearman\nCorrelation", limits = c(-1, 1)) +
    geom_point(data = function(d) d[d$sig_FDR %in% TRUE, ],
               shape = 8, size = 2.5, stroke = 1, color = "black") +
    facet_wrap(~ Clinical, nrow = 1) +
    coord_fixed(ratio = 1) +
    theme_minimal(base_family = FONT, base_size = 18) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1, size = 18,
                                  face = "plain", color = "black"),
      axis.text.y  = element_text(size = 17, face = "plain", color = "black"),
      axis.title   = element_text(size = 18, face = "plain", color = "black"),
      plot.title   = element_text(size = 18, hjust = 0.5, face = "plain", color = "black"),
      strip.text   = element_text(size = 18, face = "plain", color = "black"),
      legend.title = element_text(size = 18, face = "plain", color = "black"),
      legend.text  = element_text(size = 18, color = "black"),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.85)
    ) +
    labs(title = paste0(title, "\n(*: FDR < ", FDR_MARK, ")"),
         x = "Group", y = "Genus")

  list(plot = p, data = plot_df)
}


## =============================================================================
## PART 2 + 3. ABX+ MICE AT WEEK 8
## =============================================================================
gmpr_w8 <- read.csv(
  file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom_gmpr_wk8_dominant.csv"),
  row.names = 1, check.names = FALSE
)
aging_model_w8 <- meta_phyloseq[meta_phyloseq$Time == "Week 8", ]

## The phenotype columns are renamed by position. Print the original names and
## confirm they match before running on a regenerated metadata file.
print(colnames(aging_model_w8)[c(11, 15, 16, 17)])
colnames(aging_model_w8)[11] <- "grip_wk8"
colnames(aging_model_w8)[15] <- "speed_wk8"
colnames(aging_model_w8)[16] <- "time_wk8"
colnames(aging_model_w8)[17] <- "distance_wk8"

clinical_vars_w8 <- c("grip_wk8", "speed_wk8", "time_wk8", "distance_wk8")
genus_vars_w8    <- genus_columns(gmpr_w8, clinical_vars_w8)

res_w8 <- correlate_pooled(gmpr_w8, clinical_vars_w8, genus_vars_w8)
cat("\n=== SUMMARY WEEK 8 ===\n")
cat("Number of tests:", nrow(res_w8$long), "\n")
cat("Significant (FDR <", FDR_MARK, "):",
    sum(res_w8$long$FDR_By_Clinical < FDR_MARK, na.rm = TRUE), "\n")
write.xlsx(res_w8$long, file.path(OUT_DIR, "correlation_results_w8.xlsx"))

p_w8 <- plot_pooled_heatmap(res_w8, "Abx+ (Week 8)")
ggsave(file.path(OUT_DIR, "abx_plus_week8_correlation.svg"),
       p_w8, width = 12, height = 10, dpi = 600)

aging_df_w8 <- as.data.frame(as.matrix(aging_model_w8), stringsAsFactors = FALSE)
grp_w8 <- correlate_by_group(gmpr_w8, aging_df_w8, clinical_vars_w8, genus_vars_w8)
heat_w8 <- plot_group_heatmap(grp_w8$res, grp_w8$genus_map, clinical_vars_w8,
                              "Genus-phenotype correlation by group (Week 8)")
ggsave(file.path(OUT_DIR, "abx_plus_week8_correlation_group.svg"),
       heat_w8$plot, width = 12, height = 10, dpi = 600)
write.xlsx(heat_w8$data, file.path(OUT_DIR, "correlation_results_w8_group.xlsx"))


## =============================================================================
## PART 2 + 3. ABX- MICE AT WEEK 4
## =============================================================================
gmpr_w4 <- read.csv(
  file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom_gmpr_wk4_dominant.csv"),
  row.names = 1, check.names = FALSE
)
aging_model_w4 <- meta_phyloseq[meta_phyloseq$Time == "Week 4", ]

print(colnames(aging_model_w4)[c(9, 12, 13, 14)])
colnames(aging_model_w4)[9]  <- "grip_wk4"
colnames(aging_model_w4)[12] <- "speed_wk4"
colnames(aging_model_w4)[13] <- "time_wk4"
colnames(aging_model_w4)[14] <- "distance_wk4"

clinical_vars_w4 <- c("grip_wk4", "speed_wk4", "time_wk4", "distance_wk4")
genus_vars_w4    <- genus_columns(gmpr_w4, clinical_vars_w4)

res_w4 <- correlate_pooled(gmpr_w4, clinical_vars_w4, genus_vars_w4)
cat("\n=== SUMMARY WEEK 4 ===\n")
cat("Number of tests:", nrow(res_w4$long), "\n")
cat("Significant (FDR <", FDR_MARK, "):",
    sum(res_w4$long$FDR_By_Clinical < FDR_MARK, na.rm = TRUE), "\n")
write.xlsx(res_w4$long, file.path(OUT_DIR, "correlation_results_w4.xlsx"))

p_w4 <- plot_pooled_heatmap(res_w4, "Abx- (Week 4)")
ggsave(file.path(OUT_DIR, "abx_minus_week4_correlation.svg"),
       p_w4, width = 12, height = 10, dpi = 600)

aging_df_w4 <- as.data.frame(as.matrix(aging_model_w4), stringsAsFactors = FALSE)
grp_w4 <- correlate_by_group(gmpr_w4, aging_df_w4, clinical_vars_w4, genus_vars_w4)
heat_w4 <- plot_group_heatmap(grp_w4$res, grp_w4$genus_map, clinical_vars_w4,
                              "Genus-phenotype correlation by group (Week 4)")
ggsave(file.path(OUT_DIR, "abx_minus_week4_correlation_group.svg"),
       heat_w4$plot, width = 12, height = 10, dpi = 600)
write.xlsx(heat_w4$data, file.path(OUT_DIR, "correlation_results_w4_group.xlsx"))


## ---- Save objects used by 04_differential_abundance.R -----------------------
saveRDS(aging_model_w8, file.path(DATA_DIR, "aging_model_w8.rds"))
saveRDS(aging_model_w4, file.path(DATA_DIR, "aging_model_w4.rds"))

sink(file.path(OUT_DIR, "sessionInfo_03.txt")); print(sessionInfo()); sink()
