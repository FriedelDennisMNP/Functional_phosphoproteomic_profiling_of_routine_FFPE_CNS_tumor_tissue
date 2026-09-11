
### Functional proteomic and phosphoproteomic profiling of routine FFPE CNS tumor tissue complements genomic and epigenomic characterization

Genomic and DNA methylation-based analyses have transformed the classification of central nervous system tumors by enabling 
robust molecular diagnosis and subclassification. 
However, these approaches primarily define tumor identity and do not directly capture the functional state of intracellular 
signaling networks. Proteomic and phosphoproteomic profiling provide complementary molecular layers by measuring protein
abundance and site-specific phosphorylation. Recent methodological advances have made comprehensive analysis of formalin-fixed,
paraffin-embedded (FFPE) tissue increasingly feasible, raising the question of how functional proteomic information can 
complement established molecular diagnostics. 
Here, we applied an FFPE-compatible workflow for integrated proteomic and phosphoproteomic profiling of ten comprehensively characterized glioblastomas, 
including five EGFR-amplified and five non-amplified tumors. 
Mass spectrometry generated robust proteomic and phosphoproteomic coverage across all cases. Global proteomic profiling 
revealed group-associated protein abundance and pathway differences,including increased EGFR abundance in EGFR-amplified tumors,
while showing substantial intertumoral overlap. Global phosphoproteomic profiles similarly showed partial group-level 
separation with substantial intertumoral overlap. Differential phosphoproteomic analysis identified 443 phosphosites based on predefined exploratory statistical criteria, with only one phosphosite remaining significant after correction for multiple testing. Multiple EGFR-associated phosphosites showed increased abundance in EGFR-amplified tumors, while kinase-substrate enrichment further revealed coordinated EGFR- and SRC-associated signaling alongside considerable heterogeneity within both molecular groups.
These findings demonstrate that routine FFPE tissue retains biologically 
coherent functional information at both the protein abundance and phosphosignaling levels that can be interpreted alongside 
genomic and epigenomic data. Phosphoproteomics therefore represents an orthogonal functional layer that may complement established molecular characterization of CNS tumors.



## How to reproduce the results from the paper using the information from this repository 

# 1) Create an environment with Micromamba using the provided .yaml file:

```bash
micromamba create -f Functional_phosphoproteomic_profiling_of_routine_FFPE_CNS_tumor_tissue.yaml
micromamba activate Functional_phosphoproteomic_profiling_of_routine_FFPE_CNS_tumor_tissue

```

# 2) Execute the code provided in R/

## Overview of the R scripts

The scripts in `R/` are intended to be executed in numerical order. Scripts 1 and 2 prepare the whole-proteome and phosphoproteome data, respectively. Scripts 3 to 5 use these processed results to generate the figures and perform kinase-activity inference. Script 6 generates single-sample reports from the phosphoproteomic data.

- `1_WP_PRC_DEA_20260817_RR1296_PCF.R`: Loads DIA-NN whole-proteome results, annotates the samples, performs quality control and preprocessing, and runs differential abundance analysis between EGFR-amplified and non-amplified tumors.
- `2_PTM_PRC_DEA_20260813_RRS_1296_PCF_Phospho.R`: Loads the MaxQuant phosphoproteomic results, filters and annotates phosphosites, performs quality control and preprocessing, and runs differential phosphoproteomic analysis between EGFR-amplified and non-amplified tumors.
- `3_PTM_PRC_ADJ_DEA_RRS_1296_PCF_Phospho.R`: Loads the processed intensities and corrects phosphoproteome by proteome and performs differential phosphoproteomic analysis between EGFR-amplified and non-amplified tumors.
- `4_Figure1.R`: Loads the processed whole-proteome and phosphoproteome results and creates the protein/phosphosite identification and data-quality panels for Figure 1 and Supplemental Figure 1.
- `5_Figure_Prot2_Phos3.R`: Creates PCA plots, differential-abundance volcano plots, and pathway-enrichment visualizations for the whole-proteome and phosphoproteome results, including Figures 2 and 3.
- `6_Figure3c_Kinase_inference.R`: Performs kinase-substrate enrichment and kinase-activity inference using OmniPath enzyme-substrate resources, then generates kinase-activity plots and Supplemental Table 3.
- `7_SingleSampleReport_v2_MBR.R`: Generates a single-sample report for each phosphoproteomic sample, including sample metadata, quality-control metrics, kinase activity, and pathway-activity summaries.

## Questions

The provided environment uses R version 4.5.2. The workflow uses the following R packages: `arrow`, `assertthat`, `circlize`, `clusterProfiler`, `ComplexHeatmap`, `cowplot`, `data.table`, `decoupler`, `dplyr`, `EnhancedVolcano`, `enrichplot`, `factoextra`, `fgsea`, `furrr`, `ggally`, `ggplot2`, `ggpubr`, `ggsci`, `here`, `imputeLCMD`, `janitor`, `limma`, `lubridate`, `magrittr`, `matrixStats`, `msigdbr`, `OmnipathR`, `openxlsx`, `patchwork`, `pathview`, `plyr`, `purrr`, `remotes`, `reshape2`, `rtsne`, `r.utils`, `S4Vectors`, `shiny`, `stringr`, `SummarizedExperiment`, `sva`, `tibble`, `tidyr`, `umap`, and `vsn`.

For questions or issues, use the issue section in this repository. I will try to
respond as soon as possible.



