# Redefining Housekeeping Genes 
Source code, datasets, and instructions for reproducing the results of our paper 'Redefining Housekeeping Genes in the Context of Vertebrate Development' 

## Reproducing the Results

This repository is organized into two main folders:

- **`result_codes/`** – Python scripts used to reproduce each analysis presented in the paper.
- **`files/`** – Required datasets and annotation files (**provided as compressed files**).  

---

## Overview of Scripts and Outputs

| Script(s) | Purpose | Output (Figures/Tables) | Required Files (from `files/`) |
|---|---|---|---|
| `dev_classification.py` | Identifies developmental gene classes | Figures **1**, **S1** | Bulk embryonic RNA-seq dataset (White et al., 2017) |
| `bulk_tissue.py` | Classifies genes by tissue specificity and compares to developmental classes | Figures **3**, **S2** | Tissue RNA-seq dataset (Hu et al., 2015; GEO GSE62221) |
| `TF_cofactors.py` | Transcription factors and cofactors analysis | Figures **2B, 2C** | `Danio_rerio_TF`, `Danio_rerio_cof` (AnimalTFDB4) |
| `figure_GO_terms.py` | Visualization of GO terms enrichment | Figure **2A** | Output of GO enrichment scripts |
| `zebrafish_ontology.py`, `anatomical_pleio_subclasses.py`, `ontology_genes_phen_per_stage_and_lethality.py`, `new_pleio_forms.py` | Anatomical pleiotropy and lethality analysis | Figures **4**, **S3**, Table **S4** | ZFIN ontology and phenotype files |
| `gene_origin.py` | Evolutionary gene age analysis | Figures **5**, **S4**, **S5** | `gene_origin.csv` (GenOrigin) |
| `homology.py` | Paralog and ortholog homology analysis | Figures **6**, **S6**, Table **S6** | `paralog.txt`, `human_orthos.txt` |
| `gene_ontology.py` | GO enrichment for gene classes | Tables **S2**, **S3** | `go.obo`, `zfin.gaf` *(run `great_library_phenotype.py` first)* |
| `enrichment_U_hourglass.py`, `young_genes_enrichment.py` | GO enrichment of evolutionary subsets (U-old, S/hS-young classes) | Table **S5** | `go.obo`, `zfin.gaf` *(run `great_library_phenotype.py` first)* |
| `tau_cv_tissue_dev.py` | Analysis with continuous metrics | Figure **S7**, Table **S7** | Bulk embryonic and tissue RNA-seq datasets |
| `different_thresholds.py` | Evaluates the robustness of the gene expression thresholds | Figure **S8** | Bulk embryonic RNA-seq dataset (White et al., 2017) |
| `Wang_conversion.R`, `create_data_sc_adult_h5ad.py`, `comparison_single_cell_data.py` | Analysis of ubiquitous genes in single-cell transcriptomic data | Figure **S9** | Developmental scRNA-seq dataset (Lange et al., 2024; Zebrahub portal). Adult scRNA-seq dataset (Wang et al., 2023; Zebrafish Cell Landscape v2.0) |

---

## References

White, R. J., et al. (2017). A high-resolution mRNA expression time course of embryonic development in zebrafish. eLife, 6, e30860.
https://doi.org/10.7554/eLife.30860

Hu, P., et al. (2015). Global identification of the genetic networks and cis-regulatory elements of the cold response in zebrafish. Nucleic Acids Research, 43(19), 9198–9213.
https://doi.org/10.1093/nar/gkv780

Lange, M., et al. (2024). A multimodal zebrafish developmental atlas reveals the state-transition dynamics of late-vertebrate pluripotent axial progenitors. Cell, 187(23), 6742-6759.e17.
https://doi.org/10.1016/j.cell.2024.09.047

Wang, R., et al. (2023). Construction of a cross-species cell landscape at single-cell level. Nucleic Acids Research, 51(2), 501–516.
https://doi.org/10.1093/nar/gkac633

