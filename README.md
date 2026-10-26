# Gut-muscle-axis
Context-dependent modulation of gut–muscle axis phenotypes by Roseburia intestinalis–derived interventions

# 16S rRNA Gene Amplicon Data

16S rRNA gene amplicon data from mice in the aging model were analysed in R.

## Experimental Groups

- PBS
- D-Gal
- Live
- HK
- Sup
- RCM

## Sampling Time Points

- **Abx−:** Week 0, Week 4
- **Abx+:** Abx 0, Abx 7, Week 8

## Phenotypic Variables

Measured at Week 4 (Abx−) and Week 8 (Abx+):

- Grip strength
- Speed
- Time
- Distance

| Script | Content |
|---|---|
| `01_dada2_phyloseq.R` | DADA2 (SILVA v138) -> phyloseq object -> genus table -> rarefaction (3500 reads/sample) |
| `02_diversity.R` | Alpha diversity (Observed, Shannon, Pielou; Kruskal-Wallis, Dunn), beta diversity (Bray-Curtis, PERMANOVA, PCoA) |
| `03_correlation.R` | GMPR normalisation; Spearman correlation between genera and clinical variables (all mice, and within each group) |
| `04_differential_abundance.R` | ALDEx2 Kruskal-Wallis across groups at Week 4, Week 8, pairwise comparisons, CLR boxplots |


## Data Availability

- **Raw sequencing reads:** European Nucleotide Archive (ENA), project [PRJEB95113](https://www.ebi.ac.uk/ena/browser/view/PRJEB95113).
- **GC–MS metabolomics raw data:** MassIVE, accession [MSV000103652](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000103652).
- Processed GC–MS metabolomics data:
1) https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=02ba2a3af7f44207b84a3f253f129ca7
2) https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=c7be457a77c54dc694a5b456392b84cd
3) https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=08c980729dda4e7baf8968085496b3ca
4) https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=fc4b924d89a7496fb883caf366904235
