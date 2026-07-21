import pandas as pd
import numpy as np
import re
import glob
import matplotlib.pyplot as plt 
import os

path_save_data = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'

# =============================================================================
# 1. GTF to map TCONS to Ensembl (ENSDARG)
# =============================================================================
print("Parseando el archivo GTF...")
tcons_to_ens = {}
with open(path_save_data+'merge.annotated.gtf', 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        if '\ttranscript\t' in line:
            tcons_match = re.search(r'transcript_id "([^"]+)"', line)
            ens_match = re.search(r'ref_gene_id "([^"]+)"', line)
            xloc_match = re.search(r'gene_id "([^"]+)"', line)
            
            if tcons_match:
                tcons_id = tcons_match.group(1)
                if ens_match:
                    tcons_to_ens[tcons_id] = ens_match.group(1)
                elif xloc_match:
                    tcons_to_ens[tcons_id] = xloc_match.group(1)

print(f"Éxito: Se mapearon {len(tcons_to_ens)} transcriptos (TCONS).")

# =============================================================================
# 2. Read TPM count matrix
# =============================================================================
search_path = os.path.join(path, '*isoforms_expression.xls')
files = glob.glob(search_path)

stages_to_ignore = ['unfertilized', 'oblong', 'day6']

print(f"n files: {len(files)}") 

stage_data = {}

for file in files:
    basename = os.path.basename(file)
    stage = re.sub(r'_[0-9]+\.isoforms_expression\.xls', '', basename)
    
    if stage in stages_to_ignore:
        continue
    
    df = pd.read_csv(file, sep=r'\s+', engine='python', skiprows=1, 
                     names=['Gene', 'Isoform', 'Length', 'Count', 'TPM', 'FPKM'])
    
    df['Ensembl_ID'] = df['Isoform'].map(tcons_to_ens)
    df = df.dropna(subset=['Ensembl_ID'])
    
    gene_tpm = df.groupby('Ensembl_ID')['TPM'].sum().reset_index()
    gene_tpm = gene_tpm.rename(columns={'TPM': basename})
    
    if stage not in stage_data:
        stage_data[stage] = []
    stage_data[stage].append(gene_tpm)

print(f"Extracted data: {list(stage_data.keys())}")

# =============================================================================
# 3. Compute average TPM per stage and re-TPM to sum up 1 million
# =============================================================================
print("Compute avrage per stage and normalize...")
mean_tpm_per_stage = pd.DataFrame()

for stage, df_list in stage_data.items():
    merged_stage = df_list[0]
    for i in range(1, len(df_list)):
        merged_stage = pd.merge(merged_stage, df_list[i], on='Ensembl_ID', how='outer').fillna(0)
    
    cols_to_mean = [col for col in merged_stage.columns if col != 'Ensembl_ID']
    merged_stage[stage] = merged_stage[cols_to_mean].mean(axis=1)
    
    suma_etapa = merged_stage[stage].sum()
    if suma_etapa > 0:
        merged_stage[stage] = (merged_stage[stage] / suma_etapa) * 1000000
    
    if mean_tpm_per_stage.empty:
        mean_tpm_per_stage = merged_stage[['Ensembl_ID', stage]]
    else:
        mean_tpm_per_stage = pd.merge(mean_tpm_per_stage, merged_stage[['Ensembl_ID', stage]], on='Ensembl_ID', how='outer').fillna(0)

mean_tpm_per_stage = mean_tpm_per_stage.set_index('Ensembl_ID')

print("\nSum TPM per stage")
print(mean_tpm_per_stage.sum().round(2))

mean_bulk_matrix = mean_tpm_per_stage.values
genes = mean_tpm_per_stage.index.values
genes_name = genes

stages_names = mean_tpm_per_stage.columns.values
total_stages = len(stages_names) # Deberían ser 21

print(f"Final matrix: {len(genes)} genes x {total_stages} etapas.")

# =============================================================================
# 4. Classify genes
# =============================================================================
#TPM>1
mean_bulk_matrix_clean = np.where(mean_bulk_matrix < 1, 0, mean_bulk_matrix)

np.savetxt(path_save_data+'mean_bulk_matrix.txt', mean_bulk_matrix, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'mean_bulk_matrix_clean.txt', mean_bulk_matrix_clean, fmt='%s', delimiter=',')

plt.figure(figsize=(4, 3), dpi=600)
nonzero_data = mean_bulk_matrix[mean_bulk_matrix > 0]
plt.hist(nonzero_data, bins=100, color='lightseagreen', log=True)
plt.xlabel('TPM per gene and stage', fontsize=14, fontweight='bold')
plt.ylabel('freq', fontsize=14, fontweight='bold')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(False)
plt.savefig(path_save_data+'total_mean_gene_counts_per_gene_stage.png', dpi=600, bbox_inches='tight')
plt.show()

act_stage_per_gene = np.zeros(len(genes))
for i in range(len(genes)):
    non_null = len(np.where(mean_bulk_matrix_clean[i, :] >= 1)[0]) # Cambiado a >= 1 por el threshold
    act_stage_per_gene[i] = non_null

n_stages_arr, n_times = np.unique(act_stage_per_gene, return_counts=True)
n_stages_arr = np.array(n_stages_arr, dtype=int)

plt.figure(figsize=(8, 5), dpi=600)
plt.bar(n_stages_arr, n_times, color='darkorange', edgecolor='black')
plt.xticks(ticks=n_stages_arr, labels=n_stages_arr, rotation=90, fontsize=15)
plt.yticks(fontsize=20)
plt.xlabel('# stages expressed', fontsize=25, fontweight='bold')
plt.ylabel('# genes (TPM>1)', fontsize=25, fontweight='bold')
plt.savefig(path_save_data+'act_stages_per_gene.png', dpi=600, bbox_inches='tight')
plt.show()

# 4.) We divided the genes in specific, high-specific and ubiqutouly expressed
specific = []
high_specific = []
ubiq_genes = []
null_genes = []
null_genes_id = []
hS_id = []
U_id = []
S_id = []
act_stages_S = []
act_stages_hS = []
act_stages_U = []

ubiq_cut = total_stages
# hs_cut = total_stages // 2 # <= 10 stages
hs_cut=8

for i in range(len(genes)):
    if act_stage_per_gene[i] == 0:
        null_genes.append(genes_name[i])
        null_genes_id.append(genes[i])
    elif 0 < act_stage_per_gene[i] <= hs_cut:
        high_specific.append(genes_name[i])
        hS_id.append(genes[i])
        act_stages_hS.append(act_stage_per_gene[i])
    elif hs_cut < act_stage_per_gene[i] < ubiq_cut:
        specific.append(genes_name[i])
        S_id.append(genes[i])
        act_stages_S.append(act_stage_per_gene[i])
    elif act_stage_per_gene[i] == ubiq_cut:
        ubiq_genes.append(genes_name[i])
        U_id.append(genes[i])
        act_stages_U.append(act_stage_per_gene[i])

print(f"Total genes: {len(U_id)+len(S_id)+len(hS_id)+len(null_genes_id)}")

np.savetxt(path_save_data+'hS_id_bulk.txt', hS_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'S_id_bulk.txt', S_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'U_id_bulk.txt', U_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'null_genes.txt', null_genes_id, fmt='%s', delimiter=',')

# Boxplot
data = [act_stages_U, act_stages_S, act_stages_hS]
fig, ax = plt.subplots(figsize=(4, 3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"], widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("# expressed stages", fontsize=15, fontweight='bold')
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=16)
plt.savefig(path_save_data+'act_stages_per_class_dev.png', dpi=600, bbox_inches='tight')
plt.show()

df = pd.DataFrame({
    'gene': genes_name,
    'stages_active': act_stage_per_gene
})
df.to_csv(path_save_data+'genes_stage_count.csv', index=False)


# =============================================================================
# 5. Comparison with White et al
# =============================================================================
print("\n" + "="*60)
print("COMPARISON")
print("="*60)

old_path = '/media/alicia/TOSHIBA EXT/zebrafish/zebrafish_code/bulk/'

# 1. Files
def load_old_file(filename):
    filepath = os.path.join(old_path, filename)
    try:
        return set(np.loadtxt(filepath, dtype=str))
    except Exception as e:
        print(f"⚠️ No se pudo cargar {filename}: {e}")
        return set()

old_hS = load_old_file('hS_id_bulk.txt')
old_S  = load_old_file('S_id_bulk.txt')
old_U  = load_old_file('U_id_bulk.txt')
old_null = load_old_file('null_genes.txt')

all_old_genes = old_hS.union(old_S).union(old_U).union(old_null)

new_hS = set(hS_id)
new_S  = set(S_id)
new_U  = set(U_id)
new_null = set(null_genes_id)

def trace_new_class(class_name, new_set, old_hS, old_S, old_U, old_null, all_old):
    total_new = len(new_set)
    if total_new == 0:
        return
    
    # 1. Conserved genes
    conservados = new_set.intersection(all_old)
    total_conservados = len(conservados)
    
    # 2. New genes
    nuevos_reales = new_set.difference(all_old)
    
    print(f"\n--- ANÁLISIS DE LA NUEVA CLASE: {class_name} (Total actual: {total_new} genes) ---")
    print(f"🔹 Genes conservados (estaban en el dataset antiguo): {total_conservados} genes ({(total_conservados/total_new)*100:.2f}%)")
    print(f"🔸 Genes totalmente nuevos (NO estaban en el antiguo): {len(nuevos_reales)} genes ({(len(nuevos_reales)/total_new)*100:.2f}%)")
    
    if total_conservados > 0:
        print(f"\n   De los {total_conservados} genes conservados de la nueva clase {class_name}, en el dataset antiguo eran:")
        
        # 3. From the conserved ones, in which stage are they
        came_from_hS = conservados.intersection(old_hS)
        came_from_S  = conservados.intersection(old_S)
        came_from_U  = conservados.intersection(old_U)
        came_from_null = conservados.intersection(old_null)
        
        if len(came_from_hS) > 0:
            print(f"      -> Eran hS: {len(came_from_hS)} genes ({(len(came_from_hS)/total_conservados)*100:.2f}%)")
        if len(came_from_S) > 0:
            print(f"      -> Eran S:  {len(came_from_S)} genes ({(len(came_from_S)/total_conservados)*100:.2f}%)")
        if len(came_from_U) > 0:
            print(f"      -> Eran U:  {len(came_from_U)} genes ({(len(came_from_U)/total_conservados)*100:.2f}%)")
        if len(came_from_null) > 0:
            print(f"      -> Eran Nulos: {len(came_from_null)} genes ({(len(came_from_null)/total_conservados)*100:.2f}%)")

# 4. Trace new classes
trace_new_class('U (Ubiquitous)', new_U, old_hS, old_S, old_U, old_null, all_old_genes)
trace_new_class('S (Specific)', new_S, old_hS, old_S, old_U, old_null, all_old_genes)
trace_new_class('hS (High Specific)', new_hS, old_hS, old_S, old_U, old_null, all_old_genes)
trace_new_class('Nulos (TPM<1)', new_null, old_hS, old_S, old_U, old_null, all_old_genes)


# =============================================================================
# 6. DOTPLOT: figure
# =============================================================================

old_classes = [old_U, old_S, old_hS, old_null]
new_classes = [new_U, new_S, new_hS, new_null]
class_names = ['U', 'S', 'hS', 'NE'] 

matrix_counts = np.zeros((4, 4))
matrix_fracs = np.zeros((4, 4))
for i, new_set in enumerate(new_classes):
    conservados = new_set.intersection(all_old_genes)
    total_conservados = len(conservados)
    for j, old_set in enumerate(old_classes):
        intersection_count = len(conservados.intersection(old_set))
        matrix_counts[i, j] = intersection_count
        if total_conservados > 0:
            matrix_fracs[i, j] = intersection_count / total_conservados
        else:
            matrix_fracs[i, j] = 0

fig, ax = plt.subplots(figsize=(4.1, 4.1), dpi=600)
colors = ['royalblue', 'mediumturquoise', 'tomato', 'gray']
for i in range(4): # Eje Y (Nuevo Dataset - Bo et al.)
    for j in range(4): # Eje X (Antiguo Dataset - White et al.)
        count = matrix_counts[i, j]
        frac = matrix_fracs[i, j]
        plt.scatter(j + 0.5, i + 0.5, s=frac * 2500, color=colors[i], alpha=0.8)
        if count > 0:
            plt.text(j + 0.5, i + 0.5, int(count), color='black',
                     fontsize=10, ha='center', va='center', fontweight='bold')
for k in range(3):
    plt.axhline(y=k+1, color='gray', linestyle='--', alpha=0.5, lw=1)
    plt.axvline(x=k+1, color='gray', linestyle='--', alpha=0.5, lw=1)
ax.set_xticks(np.arange(4) + 0.5)
ax.set_xticklabels(class_names, fontsize=16, fontweight='bold')
ax.set_yticks(np.arange(4) + 0.5)
ax.set_yticklabels(class_names, fontsize=16, fontweight='bold')
plt.xlim(0, 4)
plt.ylim(4, 0) # Al invertir (4 a 0 en lugar de 0 a 4), la coordenada 0 queda arriba.
plt.ylabel("Bo et al. (2025)", fontsize=18, fontweight='bold', color='black')
plt.xlabel("White et al. (2017)", fontsize=18, fontweight='bold', color='black')
ax.xaxis.set_label_position('top') 
ax.xaxis.tick_top()
row_totals = np.sum(matrix_counts, axis=1)
for i, total in enumerate(row_totals):
    plt.text(4.2, i + 0.5, str(int(total)), va='center', ha='left', fontsize=14)
col_totals = np.sum(matrix_counts, axis=0)
for j, total in enumerate(col_totals):
    plt.text(j + 0.5, 4.3, str(int(total)), va='top', ha='center', fontsize=14)
plt.savefig(path_save_data + 'dataset_comparison_dotplot_U_first.png', dpi=600, bbox_inches='tight')
plt.show()
