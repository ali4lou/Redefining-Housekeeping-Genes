"""
Created on Mon Oct 27 12:16:42 2025

@author: Alicia
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
# ------------------------------------------------
#GO terms chosen from Table S2
# ------------------------------------------------
data = {
    "Set": ["U"] * 5 + ["S"] * 5 + ["hS"] * 5,
    "GO_term": [
        "Metabolic process", "Biosynthetic process", "Gene expression", "RNA metabolic process", "Protein-containing complex assembly",
        "Regulation of transcription by RNA Pol II", "Cell differentiation", "Anatomical structure morphogenesis", "Neurogenesis", "Cellular developmental process",
        "Cell–cell signaling", "Monoatomic ion transport", "Synaptic signaling", "Chemical synaptic transmission", "Nervous system process"
    ],
    "n_genes_subset": [
        4506, 3043, 2411, 1767, 944,   # U genes
        825, 1159, 1051, 576, 1160,     # S genes
        246, 428, 201, 198, 248         # hS genes
    ],
    "pval_corr": [
        1e-64, 1e-48, 1e-42, 1e-29, 1e-26,   # U genes
        1e-16, 1e-22, 1e-19, 1e-15, 1e-22,   # S genes
        1e-21, 1e-32, 1e-19, 1e-19, 1e-16    # hS genes
    ]
}

df = pd.DataFrame(data)
df["-log10(pval)"] = -np.log10(df["pval_corr"])

df["Set"] = pd.Categorical(df["Set"], categories=["U", "S", "hS"], ordered=True)
df = df.sort_values(["Set", "pval_corr"], ascending=[True, True]).reset_index(drop=True)

colors = {"U": "royalblue", "S": "mediumturquoise", "hS": "tomato"}

plt.figure(figsize=(11, 9))

y_pos = np.arange(len(df))

plt.barh(
    y=y_pos,
    width=df["-log10(pval)"],
    color=[colors[s] for s in df["Set"]],
    alpha=0.85,
    edgecolor=None
)

for i, row in df.iterrows():
    plt.text(
        row["-log10(pval)"] + 0.3,   # valor x del texto
        i,                           # posición y
        str(int(row["n_genes_subset"])),
        va="center",
        fontsize=14,
        color="black"
    )

plt.yticks(y_pos, df["GO_term"], fontsize=18)
plt.xlabel("-log10(corrected p-value)", fontsize=20, labelpad=12, fontweight='bold')
plt.ylabel("GO term", fontsize=25, labelpad=12, fontweight='bold')
plt.xticks(fontsize=20)

plt.grid(False)

handles = [plt.Rectangle((0,0),1,1, color=colors[s], label=f"{s} genes") for s in colors]
plt.legend(handles=handles, fontsize=18, loc="lower right", frameon=False)

plt.gca().invert_yaxis()

plt.tight_layout()
plt.savefig(path_save_data+'GO_relevant_terms_gene_class.png', dpi=600, bbox_inches='tight')

plt.show()
