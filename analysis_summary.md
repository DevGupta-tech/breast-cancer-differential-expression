# Public Dataset Capstone: Breast Cancer Differential Expression

## Project question
Which Affymetrix GPL570 probe measurements differ between breast tumour samples and normal breast tissue samples in GEO series GSE7904?

## Dataset
- GEO accession: **GSE7904**
- Organism: *Homo sapiens*
- Platform: Affymetrix Human Genome U133 Plus 2.0 (GPL570)
- Source: NCBI Gene Expression Omnibus series matrix included in `data/`
- Samples in series: 62
- Tumour group: 43 samples (Basal, BRCA1 and non-BLC labels)
- Normal breast group: 7 samples (NB label)
- Excluded group: 12 samples labelled NO; these are not pooled with normal breast tissue.

The GEO sample records describe the measurements as dChip signal values and state that the data were processed using dChip with default normalization settings.

## Analysis and statistical methods
1. Used the processed GEO series-matrix signal values.
2. Kept the 43 tumour and 7 normal-breast samples; excluded the 12 NO samples.
3. Applied `log2(max(signal, 1))` to make the signal scale suitable for fold-change interpretation. The floor at 1 is explicit because some dChip values are zero or negative; this handling is a pragmatic limitation and should not be mistaken for reprocessing raw CEL files.
4. For each of 54,675 probes, calculated the difference between mean log2 signal in tumour and normal samples.
5. Tested the group difference using a two-sided Welch two-sample t-test (does not assume equal group variances).
6. Controlled multiple testing across probes using the Benjamini-Hochberg false-discovery-rate procedure.
7. Defined statistical significance as BH-adjusted p-value < 0.05. No additional fold-change cutoff was imposed.
8. Generated a volcano plot and PCA plot of the 500 most variable probes for exploratory visualization.

## Results
- Probes tested: **54,675**
- Probes significant at BH FDR < 0.05: **18,250**
- Significant probes with higher mean signal in tumour: **6,883**
- Significant probes with lower mean signal in tumour: **11,367**
- PCA variance explained: PC1 **19.20%**, PC2 **12.04%**

The complete ranked table is `results/all_probe_results.csv`; the FDR-significant subset is `results/significant_probes_FDR_0.05.csv`. Probe IDs are retained as probe identifiers; they have not been silently treated as unique gene symbols.

## Interpretation
The analysis identifies a set of probe-level measurements whose mean log2 signal differs between the selected tumour and normal-breast samples after BH correction. Positive log2 fold change means higher signal in the tumour group; negative log2 fold change means lower signal in the tumour group. The volcano plot summarizes effect direction and adjusted statistical evidence. PCA is an exploratory view of sample-level variation and is not a hypothesis test.

These results are associations in this dataset. They do not show that altered expression causes cancer, and they do not establish diagnostic or therapeutic utility. Probe-level results are not necessarily one-to-one with genes.

## Limitations and possible confounders
- **Unequal sample sizes:** 43 tumour samples are compared with only 7 normal-breast samples, limiting precision in the normal group.
- **Tumour heterogeneity:** Basal, BRCA1-associated and non-BLC samples are pooled into one tumour group. Subtype-specific patterns may be hidden.
- **Tissue composition:** Bulk tissue expression can reflect differences in epithelial, stromal, immune and other cell proportions, not only tumour-cell regulation.
- **Clinical and technical confounding:** Patient characteristics, sample handling, batch effects and other unrecorded variables may contribute to observed differences.
- **Processed data:** This analysis uses GEO's processed dChip signal values, not raw CEL files. A full reprocessing with a documented Affymetrix workflow (for example, RMA from raw CEL files) would be a useful sensitivity analysis.
- **Negative signal values:** Some processed dChip signals are non-positive. Flooring values below 1 before log2 transformation can affect low-intensity probes and is a material preprocessing choice.
- **Probe annotation:** Results are reported by GPL570 probe ID. Probe-to-gene mapping can be ambiguous or many-to-one and should be performed and documented before gene-level biological claims.
- **Observational design:** Statistical association does not establish causality.

## Reproducibility
Run from the project root in R/RStudio:

```r
source("scripts/run_analysis.R")
```

The script uses base R only and recreates the sample classification, probe-level Welch tests, BH correction, result tables, volcano plot, PCA plot and summary report. It reads the supplied GEO series matrix from `data/GSE7904_series_matrix.txt.gz`.

## Data source
NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904
