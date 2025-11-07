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
| `dev_classification.py` | Identifies developmental gene classes and generates visualizations | Figures **1**, **S1** | Bulk embryonic RNA-seq dataset (White et al., 2017) |
| `bulk_tissue.py` | Classifies genes by tissue specificity and compares to developmental classes | Tissue-specific classifications | Tissue RNA-seq dataset (Hu et al., 2015; GEO GSE62221) |
| `TF_cofactors.py` | Transcription factor and cofactor class composition analysis | Figures **2B–2C** | `Danio_rerio_TF`, `Danio_rerio_cof` (AnimalTFDB4) |
| `figure_GO_terms.py` | Summary visualization of GO terms enrichment | Figure **2A** | Output of GO enrichment scripts |
| `zebrafish_ontology.py`, `anatomical_pleio_subclasses.py`, `ontology_genes_phen_per_stage_and_lethality.py` | Anatomical pleiotropy and developmental lethality analysis | Figures **4**, **S3**, Table **S4** | ZFIN ontology and phenotype files |
| `gene_origin.py` | Evolutionary gene age determination | Figures **5**, **S4** | `gene_origin.csv` (GenOrigin) |
| `homology.py` | Paralog and ortholog homology analysis | Figures **6**, **S5** | `paralog.txt`, `human_orthos.txt` |
| `gene_ontology.py` | GO enrichment for developmental gene classes | Tables **S2**, **S3** | `go.obo`, `zfin.gaf` *(run `great_library_phenotype.py` first)* |
| `enrichment_U_hourglass.py`, `young_genes_enrichment.py` | GO enrichment of evolutionary subsets (U-old, S/hS-young classes) | Table **S5** | Same as `gene_ontology.py` + *(run `great_library_phenotype.py` first)* |

---

## References

White, R. J., et al. (2017). A high-resolution mRNA expression time course of embryonic development in zebrafish. eLife, 6, e30860.
https://doi.org/10.7554/eLife.30860

Hu, P., et al. (2015). Global identification of the genetic networks and cis-regulatory elements of the cold response in zebrafish. Nucleic Acids Research, 43(19), 9198–9213.
https://doi.org/10.1093/nar/gkv780

