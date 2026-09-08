
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
separation with substantial intertumoral overlap. Differential phosphoproteomic analysis identified 483 
phosphosites based on predefined exploratory statistical criteria, with multiple EGFR phosphosites, including Y1110, Y1197, and T693,
among the strongest changes in EGFR-amplified tumors. 
Kinase-substrate enrichment further revealed coordinated EGFR- and SRC-associated signaling alongside considerable 
heterogeneity within both molecular groups. These findings demonstrate that routine FFPE tissue retains biologically 
coherent functional information at both the protein abundance and phosphosignaling levels that can be interpreted alongside 
genomic and epigenomic data.
Phosphoproteomics therefore represents an orthogonal functional layer that may complement established molecular characterization of CNS tumors.

## How to reproduce the results from the paper using the information from this repository 

# 1) Create an environment with Micromamba using the provided .yaml file:

```bash
micromamba create -f Functional_phosphoproteomic_profiling_Friedel_et_al.yaml
micromamba activate Functional_phosphoproteomic_profiling_Friedel_et_al

```

# 2) Download data form PRoteomics Identification Database (PRIDE) 
https://www.ebi.ac.uk/pride/ using the accession IDS 

# 3) Exectue the code provided in R/analysis/figures

## Questions ? 

For questions or issues use the issue section in this repository, I will try to
respond as soon as possible.



