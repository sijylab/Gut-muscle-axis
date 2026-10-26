## =============================================================================
## 01_dada2_phyloseq.R
## 16S rRNA amplicon processing: DADA2 -> phyloseq -> genus table -> rarefaction
##
## Part 1. DADA2 (demultiplexed paired-end reads, R1 = forward, R2 = reverse)
## Part 2. DADA2 output -> phyloseq object, merge with sample metadata
## Part 3. Genus-level agglomeration and rarefaction
##
## All paths are relative to the project root (the directory containing this
## script's parent folder). Raw fastq files and the SILVA reference database are
## not distributed with this repository; set RAW_READS_DIR and SILVA_* below to
## their location on your machine. See README.md for the expected layout.
##
## Outputs written to data/: track.tsv, sequence_table_SV.tsv, taxa_SV.tsv,
## mouse_silva_table_cleaned_genus_glom.csv and the rarefied tables used by
## 02_diversity.R.
## =============================================================================

library(here)
library(dada2)
library(phyloseq)
library(data.table)

packageVersion("dada2")
packageVersion("phyloseq")


## ---- Parameters -------------------------------------------------------------
## Everything that affects the result is set here, so the whole pipeline can be
## reproduced from this block alone.

RAW_READS_DIR <- here("data", "raw_reads")          # demultiplexed *.fastq.gz
SILVA_TRAIN   <- here("data", "reference", "silva_nr_v138_train_set.fa.gz")
SILVA_SPECIES <- here("data", "reference", "silva_species_assignment_v138.fa.gz")
DATA_DIR      <- here("data")

# Filtering and trimming (dada2::filterAndTrim)
TRIM_LEFT   <- c(20, 20)   # primer length removed from the 5' end of R1 / R2
MIN_LEN     <- 100
MAX_N       <- 0
MAX_EE      <- c(2, 2)
TRUNC_Q     <- 15
RM_PHIX     <- TRUE

# Error model
ERR_NREADS  <- 1e6
SEED_ERR    <- 12345

# Sample filtering
MIN_SAMPLE_READS <- 1000   # samples with fewer total reads are dropped

# Rarefaction
RAREFY_DEPTH <- 3500       # depth used for the diversity analyses in the paper
RAREFY_SEED  <- 711

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)


## =============================================================================
## PART 1. DADA2
## =============================================================================

## ---- 1.1 Input files --------------------------------------------------------
fns     <- list.files(RAW_READS_DIR)
fastqs  <- sort(fns[grepl("\\.fastq\\.gz$", fns)])   # sorting keeps R1/R2 aligned

fnFs <- fastqs[grepl("_1.fastq.gz", fastqs)]   # forward reads
fnRs <- fastqs[grepl("_2.fastq.gz", fastqs)]   # reverse reads

# Sample names = the part of the forward read file name before "_1.fastq.gz"
sample.names <- sapply(strsplit(fnFs, "_1.fastq.gz"), `[`, 1)

fnFs <- file.path(RAW_READS_DIR, fnFs)
fnRs <- file.path(RAW_READS_DIR, fnRs)

stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)


## ---- 1.2 Quality profiles (raw reads) ---------------------------------------
pdf(file.path(DATA_DIR, "plotQualityProfile.pdf"), onefile = TRUE)
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])
dev.off()


## ---- 1.3 Filtering and trimming ---------------------------------------------
## Primers are removed with trimLeft; low-quality regions are removed by the
## expected-error (maxEE) and quality (truncQ) filters.
filt_path <- file.path(RAW_READS_DIR, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      trimLeft = TRIM_LEFT, minLen = MIN_LEN, maxN = MAX_N,
                      maxEE = MAX_EE, truncQ = TRUNC_Q, rm.phix = RM_PHIX,
                      compress = TRUE, multithread = TRUE)
head(out2)

pdf(file.path(DATA_DIR, "plotQualityProfile.filt.pdf"), onefile = TRUE)
plotQualityProfile(filtFs[1:15])
plotQualityProfile(filtRs[1:15])
dev.off()


## ---- 1.4 Error models -------------------------------------------------------
set.seed(SEED_ERR)

errF <- learnErrors(filtFs, nread = ERR_NREADS, multithread = TRUE)
errR <- learnErrors(filtRs, nread = ERR_NREADS, multithread = TRUE)

pdf(file.path(DATA_DIR, "plotErrors_F.pdf"), onefile = TRUE)
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(DATA_DIR, "plotErrors_R.pdf"), onefile = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()


## ---- 1.5 Dereplication, sample inference and merging ------------------------
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)


## ---- 1.6 Sequence table and chimera removal ---------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))   # distribution of sequence lengths

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)     # fraction of reads retained


## ---- 1.7 Read tracking (Supplementary Table S15) ----------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out2,
               sapply(dadaFs, getN),
               sapply(mergers, getN),
               rowSums(seqtab),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)


## ---- 1.8 Taxonomy assignment (SILVA v138) -----------------------------------
taxHS <- assignTaxonomy(seqtab.nochim, SILVA_TRAIN, multithread = TRUE)
taxHS <- addSpecies(taxHS, SILVA_SPECIES)
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(taxHS))


## ---- 1.9 Write to disk ------------------------------------------------------
write.table(track,         file = file.path(DATA_DIR, "track.tsv"),             quote = FALSE)
write.table(seqtab.nochim, file = file.path(DATA_DIR, "sequence_table_SV.tsv"), quote = FALSE)
write.table(taxHS,         file = file.path(DATA_DIR, "taxa_SV.tsv"),           quote = FALSE)


## =============================================================================
## PART 2. DADA2 -> PHYLOSEQ
## =============================================================================

## ---- 2.1 Function: DADA2 output -> phyloseq ---------------------------------
dada2_to_phyloseq <- function(abundance.file, tax.file, sep.otus = " ", sep.tax = " ") {

  data <- fread(abundance.file, header = TRUE, fill = TRUE, dec = ".", sep = sep.otus)
  data.colnames <- colnames(data)[1:length(colnames(data)) - 1]
  data <- data.frame(data, row.names = 1)
  colnames(data) <- data.colnames

  data <- t(data)   # rows = taxa, columns = samples

  tax <- read.table(tax.file, header = TRUE, row.names = 1, dec = ".", sep = sep.tax)
  tax <- as.matrix(tax)

  percNA     <- table(is.na(tax[, 6]))
  percent_NA <- (percNA["TRUE"] / sum(percNA)) * 100
  cat(abundance.file, "\n")
  cat("%NA genus level: ", percent_NA, "\n")

  # Prefix each rank, and label unresolved ranks as "uc_<last known rank>"
  for (i in rownames(tax)) {
    temp.line <- tax[i, ]
    temp.line[1] <- paste("k", temp.line[1], sep = "_")
    temp.line[2] <- paste("p", temp.line[2], sep = "_")
    temp.line[3] <- paste("c", temp.line[3], sep = "_")
    temp.line[4] <- paste("o", temp.line[4], sep = "_")
    temp.line[5] <- paste("f", temp.line[5], sep = "_")
    temp.line[6] <- paste("g", temp.line[6], sep = "_")

    if (length(grep("_NA", temp.line)) > 0) {
      na.ind   <- grep("_NA", temp.line)     # positions of the missing labels
      na.lev   <- min(na.ind)                # first missing rank
      last.tax <- temp.line[na.lev - 1]      # last known rank
      temp.line[na.ind] <- rep(paste("uc", last.tax, sep = "_"), length(temp.line[na.ind]))
    }
    tax[i, ] <- temp.line
  }

  OTU <- otu_table(data, taxa_are_rows = TRUE)
  TAX <- tax_table(tax[rownames(OTU), ])
  phyloseq(OTU, TAX)
}


## ---- 2.2 Build phyloseq object and merge sample data ------------------------
dada2_phyloseq <- dada2_to_phyloseq(
  abundance.file = file.path(DATA_DIR, "sequence_table_SV.tsv"),
  tax.file       = file.path(DATA_DIR, "taxa_SV.tsv"),
  sep.otus = " ", sep.tax = " "
)

a <- read.csv(file.path(DATA_DIR, "Aging_Miseq_Sample_Information.csv"),
              row.names = 1, check.names = FALSE)
meta_phyloseq   <- sample_data(a)
dada2_phyloseq2 <- merge_phyloseq(dada2_phyloseq, meta_phyloseq)

save(dada2_phyloseq2, file = file.path(DATA_DIR, "dada2_phyloseq2.RData"))


## ---- 2.3 Filtering ----------------------------------------------------------
## Remove non-target and contaminant taxa
in_physeq <- subset_taxa(
  dada2_phyloseq2,
  Kingdom != "k_Eukaryota" &
  Kingdom != "uc_k_Bacteria" &
  Kingdom != "uc_" &
  Phylum  != "NA" &
  Order   != "o_Chloroplast" &
  Species != "s_lentus"
)

## Remove taxa with zero total abundance
taxSums    <- taxa_sums(in_physeq)
in_physeq2 <- prune_taxa(names(taxSums[taxSums > 0]), in_physeq)

## Remove samples with fewer than MIN_SAMPLE_READS total reads (1,000)
in_physeq3 <- prune_samples(sample_sums(in_physeq2) >= MIN_SAMPLE_READS, in_physeq2)

save(in_physeq3, file = file.path(DATA_DIR, "phyloseq_filtered.RData"))

## Agglomerate to genus level. All diversity analyses in this study are
## performed on the genus-level table produced here.
in_physeq3_genus <- tax_glom(in_physeq3, taxrank = "Genus")
save(in_physeq3_genus, file = file.path(DATA_DIR, "phyloseq_filtered_genus_glom.RData"))

## Genus-level count table (rows = genera, columns = samples)
otu_genus <- as(otu_table(in_physeq3_genus), "matrix")
rownames(otu_genus) <- as.character(tax_table(in_physeq3_genus)[, "Genus"])
write.csv(otu_genus, file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom.csv"))

write.csv(as(sample_data(in_physeq3_genus), "data.frame"),
          file.path(DATA_DIR, "sample_metadata_filtered.csv"))


## =============================================================================
## PART 3. RAREFACTION
## =============================================================================
## The diversity analyses in 02_diversity.R use the table rarefied to
## RAREFY_DEPTH. The additional depths are written for the rarefaction-depth
## sensitivity analysis reported in Supplementary Tables S16-S18.

otutable <- read.csv(file.path(DATA_DIR, "mouse_silva_table_cleaned_genus_glom.csv"),
                     row.names = 1, check.names = FALSE)
min(colSums(otutable))

otutable_phyloseq <- phyloseq(otu_table(as.matrix(otutable), taxa_are_rows = TRUE))

rarefy_to <- function(physeq, depth, seed = RAREFY_SEED) {
  rarefy_even_depth(physeq, sample.size = depth, rngseed = seed,
                    replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
}

## Primary table used throughout the paper
oturare <- rarefy_to(otutable_phyloseq, RAREFY_DEPTH)
write.csv(oturare,
          file.path(DATA_DIR, sprintf("mouse_silva_table_cleaned_genus_glom_%d.csv", RAREFY_DEPTH)))
saveRDS(oturare,  file.path(DATA_DIR, "oturare_in_data2_to_phyloseq.rds"))
saveRDS(otutable, file.path(DATA_DIR, "otutable_in_data2_to_phyloseq.rds"))

## Additional depths for the sensitivity analysis. These are written to their
## own files and must NOT overwrite oturare_in_data2_to_phyloseq.rds, which is
## the object every downstream script reads.
for (depth in c(2000, 6000)) {
  write.csv(rarefy_to(otutable_phyloseq, depth),
            file.path(DATA_DIR, sprintf("mouse_silva_table_cleaned_genus_glom_%d.csv", depth)))
}

sink(file.path(DATA_DIR, "sessionInfo_01.txt")); print(sessionInfo()); sink()
