# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 09:36:54 2025

GO ontology:
    format-version: 1.2
    data-version: releases/2025-02-06

@author: Alicia
"""

import obonet
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import great_library_phenotypes as glp
from scipy.stats import fisher_exact, false_discovery_control
from scipy import stats



def flatten_and_unique(nested_list):
    result = []
    for sublist in nested_list:
        unique_items = set()
        for item in sublist:
            unique_items.update(item)  
        result.append(list(unique_items)) 
    return result


def create_association_matrix_gene_go_term(gene_subset, GO_specific_terms, pathlist_specific_terms, GO_specific_terms_descrip):
    
    common_genes=np.intersect1d(gene_subset, genes_id_type)
    print(len(common_genes))
    GO_specific_terms=np.array(GO_specific_terms)
    
    matrix=np.zeros((len(common_genes), len(GO_specific_terms)))
    
    for i in range(len(common_genes)):
        ind_gene=np.where(genes_id_type==common_genes[i])[0]
        for j in range(len(go_type[int(ind_gene)])):
            
            go_ind=np.where(GO_specific_terms==go_type[int(ind_gene)][j])[0]
            if len(go_ind)>0:
                ind_matrix = np.isin(GO_specific_terms, pathlist_specific_terms[int(go_ind)])
                matrix[i, ind_matrix] = 1
    return matrix
    



def enrichement_go(big_matrix_all_genes, submatrix_gene_go, go_list, go_list_term, label):

    #numberof times that a go_term is associated with a gene
    go_n_times_all_genes=np.zeros(len(go_list))
    for i in range(len(go_list)):
        go_n_times_all_genes[i]=np.sum(big_matrix_all_genes[:, i])

    odd_ratio_enrich=np.zeros(len(go_list))
    p_value_enrich=np.zeros(len(go_list))
    n_genes_subset=len(submatrix_gene_go[:, 0])
    go_enrich_fisher_genes_subset=[]
    go_term_enrich=[]
    p_value_enriched_go_term=[]
    n_genes=[]
    n_genes_subset_associated_go=[]
    subset_analyzed=[]
    #For each phenotype we compute a score that indicates if the phenotypes is enriched
    for fen in range(len(go_list)):
        go_n_times_subset=np.sum(submatrix_gene_go[:, fen])
        tabla=[[go_n_times_subset, n_genes_subset-go_n_times_subset],[go_n_times_all_genes[fen], len(big_matrix_all_genes[:, 0])-go_n_times_all_genes[fen]]]
        odd_ratio_enrich[fen], p_value_enrich[fen] = fisher_exact(tabla, alternative="greater") 
        if p_value_enrich[fen]<0.001:
            go_enrich_fisher_genes_subset.append(go_list[fen])
            go_term_enrich.append(go_list_term[fen])
            p_value_enriched_go_term.append(p_value_enrich[fen])
            n_genes.append(go_n_times_all_genes[fen])
            n_genes_subset_associated_go.append(go_n_times_subset)
            subset_analyzed.append(label)
            
    return np.array(subset_analyzed), np.array(go_enrich_fisher_genes_subset), np.array(go_term_enrich), np.array(p_value_enriched_go_term), np.array(n_genes), np.array(n_genes_subset_associated_go)




path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'



#1) ANALYSIS OF WORM PHENOTYPE ONTOLOGY (WPO)

#We read the data sheet with the Worm Phenotype Ontology
# url='https://downloads.wormbase.org/releases/current-production-release/ONTOLOGY/phenotype_ontology.WS290.obo'
url=path_save_data+'go.obo'
graph = obonet.read_obo(url)

#create mappings
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
gene_description = {data['def']: id_ for id_, data in graph.nodes(data=True) if 'def' in data}

del url


#We keep nodes -> id phenotype and name phenotype
GO_term=glp.convert_dictionary_to_array(id_to_name)
go_description=glp.convert_dictionary_to_array(gene_description)

#We find all the possible paths for each phenotype
pathlist_bio_process=[]
GO_bio_process=[]
GO_bio_process_descrip=[]
for i in range(len(GO_term)):
    start=id_to_name[GO_term[i][0]]
    
    paths = nx.all_simple_paths(
        graph,
        source=name_to_id[start],
        target=name_to_id['biological_process']
    )
    innerlist = []
    for path in paths:
        innerlist.append(path)
    
    if len(innerlist)>0:
        pathlist_bio_process.append(innerlist)
        GO_bio_process.append(GO_term[i][0])
        GO_bio_process_descrip.append(GO_term[i][1])
        
pathlist_cell_comp=[]
GO_cell_comp=[]
GO_cell_comp_descrip=[]
for i in range(len(GO_term)):
    start=id_to_name[GO_term[i][0]]
    
    paths = nx.all_simple_paths(
        graph,
        source=name_to_id[start],
        target=name_to_id['cellular_component']
    )
    innerlist = []
    for path in paths:
        innerlist.append(path)
        
    if len(innerlist)>0:
        pathlist_cell_comp.append(innerlist)
        GO_cell_comp.append(GO_term[i][0])
        GO_cell_comp_descrip.append(GO_term[i][1])
    
pathlist_molecular=[]
GO_molecular=[]
GO_molecular_descrip=[]
for i in range(len(GO_term)):
    start=id_to_name[GO_term[i][0]]
    
    paths = nx.all_simple_paths(
        graph,
        source=name_to_id[start],
        target=name_to_id['molecular_function']
    )
    innerlist = []
    for path in paths:
        innerlist.append(path)
    
    if len(innerlist)>0:
        pathlist_molecular.append(innerlist)
        GO_molecular.append(GO_term[i][0])
        GO_molecular_descrip.append(GO_term[i][1])




#2.) We read the gene-GOterm association dataset 
file_path = path_save_data + "zfin.gaf"

df = pd.read_csv(
    file_path, 
    sep="\t", 
    comment="!",  # Ignore comment lines starting with '!'
    header=None,  # No predefined header
    on_bad_lines="skip",  # Skip problematic lines (pandas >= 1.3.0)
    engine="python"  # More robust parsing
)

print(df.head())  # Show the first few rows


gene_id=np.array(df.iloc[:, 1] , dtype=str)
gene_name=df.iloc[:, 2]
GOterm_association=np.array(df.iloc[:, 4], dtype=str)
type_association=df.iloc[:, 3]



genes_id_type=np.unique(gene_id)

# genes_id_type2=list(set(gene_id))
# print(genes_id_type[5420])


#3.) In this loop we find the associatied go_terms with each unique gene
go_type=[]
#go_terms corresponing genes_id_tipo 
for i in range(len(genes_id_type)):
    ind_genes=np.where(gene_id==genes_id_type[i])
    ind_genes=ind_genes[0]
    # inner_list_bio_process=[]
    # inner_list_cell_comp=[]
    # inner_list_molecular=[]
    inner_list=[]
    for j in range(len(ind_genes)):
        index=ind_genes[j]
        inner_list.append(GOterm_association[index])
    go_type.append(inner_list)
del ind_genes, inner_list, index



#4.) Zfin gene related to ENS gene
f=open(path_save_data+'zfin_to_ENS.txt', 'r')
txt = f.read()
zfin_to_ENS = txt.split('\n')
del txt, f
zfin_to_ENS=np.delete(zfin_to_ENS, len(zfin_to_ENS)-1)
zfin_to_ENS=np.array(zfin_to_ENS)

zfin_all = []
ENS_all=[]

for elemento in zfin_to_ENS:
    partes = elemento.split("\t")

    for i in range(len(partes)):
    
        if partes[i].startswith("ZDB"):
            zfin_all.append(partes[i])
        
        if partes[i].startswith("ENS"):
            ENS_all.append(partes[i])

ENS_all=np.array(ENS_all)
zfin_all=np.array(zfin_all)

del zfin_to_ENS, partes, elemento


#4.1.) We find the final genes and their associations
common_genes_zfin, ind_zfin_go, ind_ens=np.intersect1d(genes_id_type, zfin_all, return_indices=True)


genes=ENS_all[ind_ens]
genes_unique=np.unique(genes)


go_type_list_all_genes=[]
for i in range(len(ind_zfin_go)):
    go_type_list_all_genes.append(go_type[int(ind_zfin_go[i])])
    
    
#4.2.) We find the final GO list with the unique GO terms
go_type_list=[]
count=0
for i in range(len(genes_unique)):
    ind_gene=np.where(genes==genes_unique[i])[0]
    if len(ind_gene)>1:
        go_type_list.append(np.unique(np.concatenate((go_type_list_all_genes[int(ind_gene[0])], go_type_list_all_genes[int(ind_gene[1])]))))
    else:
        go_type_list.append(go_type_list_all_genes[int(ind_gene)])

del genes, go_type_list_all_genes

del go_type, genes_id_type, zfin_all

go_type=go_type_list
genes_id_type=genes_unique

del genes_unique, go_type_list





#ENRICHMENT
#5.) Gene-GO term associaciation matrix  and enrichment (for each subset of genes to analyze)
pathlist_cell_comp_unique = flatten_and_unique(pathlist_cell_comp)
pathlist_bio_process_unique=flatten_and_unique(pathlist_bio_process)
pathlist_molecular_unique=flatten_and_unique(pathlist_molecular)


# GO_bio_process_descrip=np.array(GO_bio_process_descrip)
# np.where(GO_bio_process_descrip=='biological_process')[0]

#5.0.) We read the U, hS and S genes
f=open(path_save_data+'hS_id_bulk.txt', 'r')
txt = f.read()
hS_id_bulk = txt.split('\n')
del txt, f
hS_id_bulk=np.delete(hS_id_bulk, len(hS_id_bulk)-1)
hS_id_bulk=np.array(hS_id_bulk)

f=open(path_save_data+'S_id_bulk.txt', 'r')
txt = f.read()
S_id_bulk = txt.split('\n')
del txt, f
S_id_bulk=np.delete(S_id_bulk, len(S_id_bulk)-1)
S_id_bulk=np.array(S_id_bulk)

f=open(path_save_data+'U_id_bulk.txt', 'r')
txt = f.read()
U_id_bulk = txt.split('\n')
del txt, f
U_id_bulk=np.delete(U_id_bulk, len(U_id_bulk)-1)
U_id_bulk=np.array(U_id_bulk)

genes_analyze_list=[U_id_bulk, S_id_bulk, hS_id_bulk]
genes_label_list=['U_id_bulk', 'S_id_bulk', 'hS_id_bulk']
genes=genes_id_type

for i in range(3):

    genes_analyze=genes_analyze_list[i]
    genes_label=genes_label_list[i]
    label_analyzed=genes_label_list[i]
    
    # #5.1.) BIOLOGICAL PROCESS
    
    big_matrix_bio_process=create_association_matrix_gene_go_term(genes, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
    
    
    matrix_subset_bio_process=create_association_matrix_gene_go_term(genes_analyze, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
    subset, enrich_subset_bio_process, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_bio_process, matrix_subset_bio_process, GO_bio_process_descrip, GO_bio_process,  label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_bio_process=pd.DataFrame()
    df_bio_process['subset']=subset
    df_bio_process['GO term']=enrich_go_term[index_sorted]
    df_bio_process['GO description']=enrich_subset_bio_process[index_sorted]
    df_bio_process['n genes']=n_genes[index_sorted]
    df_bio_process['n genes subset']=n_genes_subset[index_sorted]
    df_bio_process['p-value']=p_value[index_sorted]
    df_bio_process['corected p-value']=corrected_p_value[index_sorted]
    
    df_bio_process.to_csv(path_save_data+'enrich_go_bio_process_%s.csv' %genes_label, sep='\t')
    np.savetxt(path_save_data+"matrix_bio_process_%s.txt" %genes_label, matrix_subset_bio_process, fmt="%f")
    
    del big_matrix_bio_process, matrix_subset_bio_process, enrich_subset_bio_process, df_bio_process
    
    # #5.2.) MOLECULAR FUNCTION
    
    big_matrix_molecular_function=create_association_matrix_gene_go_term(genes, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
    
    
    matrix_subset_molecular_func=create_association_matrix_gene_go_term(genes_analyze, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
    subset, enrich_subset_molecular_func, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_molecular_function, matrix_subset_molecular_func, GO_molecular_descrip, GO_molecular, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_molecular=pd.DataFrame()
    df_molecular['subset']=subset
    df_molecular['GO term']=enrich_go_term[index_sorted]
    df_molecular['GO description']=enrich_subset_molecular_func[index_sorted]
    df_molecular['n genes']=n_genes[index_sorted]
    df_molecular['n genes subset']=n_genes_subset[index_sorted]
    df_molecular['p-value']=p_value[index_sorted]
    df_molecular['corected p-value']=corrected_p_value[index_sorted]
    
    df_molecular.to_csv(path_save_data+'enrich_go_molecular_func_%s.csv' %genes_label, sep='\t')
    np.savetxt(path_save_data+"matrix_molecular_func_%s.txt" %genes_label, matrix_subset_molecular_func, fmt="%f")

    
    del big_matrix_molecular_function, matrix_subset_molecular_func, enrich_subset_molecular_func, df_molecular
    
    
    # #5.3.) CELL COMPONENT
    
    big_matrix_cell_comp=create_association_matrix_gene_go_term(genes, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
    
    
    matrix_subset_cell_comp=create_association_matrix_gene_go_term(genes_analyze, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
    subset, enrich_subset_cell_comp, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_cell_comp, matrix_subset_cell_comp, GO_cell_comp_descrip, GO_cell_comp, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_cell_comp=pd.DataFrame()
    df_cell_comp['subset']=subset
    df_cell_comp['GO term']=enrich_go_term[index_sorted]
    df_cell_comp['GO description']=enrich_subset_cell_comp[index_sorted]
    df_cell_comp['n genes']=n_genes[index_sorted]
    df_cell_comp['n genes subset']=n_genes_subset[index_sorted]
    df_cell_comp['p-value']=p_value[index_sorted]
    df_cell_comp['corected p-value']=corrected_p_value[index_sorted]
    
    df_cell_comp.to_csv(path_save_data+'enrich_go_cell_comp_%s.csv' %genes_label, sep='\t')
    np.savetxt(path_save_data+"matrix_cell_comp_%s.txt" %genes_label, matrix_subset_cell_comp, fmt="%f")
    
    del big_matrix_cell_comp, matrix_subset_cell_comp, enrich_subset_cell_comp, df_cell_comp
    
    print(i)


#6.) Pleiotropy

#6.1.) Biological process pleiotropy

matrix_U_bio_process = np.loadtxt(path_save_data+"matrix_bio_process_U_id_bulk.txt", dtype=float)
matrix_S_bio_process = np.loadtxt(path_save_data+"matrix_bio_process_S_id_bulk.txt", dtype=float)
matrix_hS_bio_process = np.loadtxt(path_save_data+"matrix_bio_process_hS_id_bulk.txt", dtype=float)

pleio_U_bp=np.sum(matrix_U_bio_process, axis=1)
pleio_S_bp=np.sum(matrix_S_bio_process, axis=1)
pleio_hS_bp=np.sum(matrix_hS_bio_process, axis=1)

#figure boxplot
data=[pleio_U_bp, pleio_S_bp, pleio_hS_bp]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

plt.ylabel("Pleiotropy\n(Biological process)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 2.2*y_max, y_max + 3.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.1 * y_max, label, ha='center', fontsize=10, c='red')
    
plt.savefig(path_save_data+'pleiotropy_bio_process_boxplot_log.png', dpi=600, bbox_inches='tight')

plt.show()

#KS test
#6.1.1.) ks test between pleio distrib
ks_U_S, pvalue_U_S = stats.ks_2samp(pleio_U_bp, pleio_S_bp)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(pleio_U_bp, pleio_hS_bp)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(pleio_S_bp, pleio_hS_bp)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between hS and S, U pleiotropy distrib')
    
print(np.median(pleio_U_bp), np.median(pleio_S_bp), np.median(pleio_hS_bp))

#6.2.) Molecular function pleiotropy

matrix_U_mf = np.loadtxt(path_save_data+"matrix_molecular_func_U_id_bulk.txt", dtype=float)
matrix_S_mf = np.loadtxt(path_save_data+"matrix_molecular_func_S_id_bulk.txt", dtype=float)
matrix_hS_mf = np.loadtxt(path_save_data+"matrix_molecular_func_hS_id_bulk.txt", dtype=float)

pleio_U_mf=np.sum(matrix_U_mf, axis=1)
pleio_S_mf=np.sum(matrix_S_mf, axis=1)
pleio_hS_mf=np.sum(matrix_hS_mf, axis=1)

#figure boxplot
data=[pleio_U_mf, pleio_S_mf, pleio_hS_mf]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

plt.ylabel("Pleiotropy\n(Molecular function)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 2.2*y_max, y_max + 3.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.1 * y_max, label, ha='center', fontsize=10, c='red')

plt.savefig(path_save_data+'pleiotropy_molecular_function_boxplot_log.png', dpi=600, bbox_inches='tight')

plt.show()

#KS test
#6.2.1.) ks test between pleio distrib
ks_U_S, pvalue_U_S = stats.ks_2samp(pleio_U_mf, pleio_S_mf)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(pleio_U_mf, pleio_hS_mf)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(pleio_S_mf, pleio_hS_mf)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between hS and S, U pleiotropy distrib')
    
print(np.median(pleio_U_mf), np.median(pleio_S_mf), np.median(pleio_hS_mf))


#6.2.) Celular component pleiotropy

matrix_U_cc = np.loadtxt(path_save_data+"matrix_cell_comp_U_id_bulk.txt", dtype=float)
matrix_S_cc = np.loadtxt(path_save_data+"matrix_cell_comp_S_id_bulk.txt", dtype=float)
matrix_hS_cc = np.loadtxt(path_save_data+"matrix_cell_comp_hS_id_bulk.txt", dtype=float)

pleio_U_cc=np.sum(matrix_U_cc, axis=1)
pleio_S_cc=np.sum(matrix_S_cc, axis=1)
pleio_hS_cc=np.sum(matrix_hS_cc, axis=1)

#figure boxplot
data=[pleio_U_cc, pleio_S_cc, pleio_hS_cc]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

plt.ylabel("Pleiotropy\n(Cellular component)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 2.2*y_max, y_max + 3.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.1 * y_max, label, ha='center', fontsize=10, c='red')

plt.savefig(path_save_data+'pleiotropy_cell_comp_boxplot_log.png', dpi=600, bbox_inches='tight')

plt.show()

#KS test
#6.2.1.) ks test between pleio distrib
ks_U_S, pvalue_U_S = stats.ks_2samp(pleio_U_cc, pleio_S_cc)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(pleio_U_cc, pleio_hS_cc)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(pleio_S_cc, pleio_hS_cc)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between hS and S, U pleiotropy distrib')
    
print(np.median(pleio_U_cc), np.median(pleio_S_cc), np.median(pleio_hS_cc))






#7.) GO enrichment for each developmental category subclassified in tissues

genes=genes_id_type

pathlist_cell_comp_unique = flatten_and_unique(pathlist_cell_comp)
pathlist_bio_process_unique=flatten_and_unique(pathlist_bio_process)
pathlist_molecular_unique=flatten_and_unique(pathlist_molecular)


big_matrix_bio_process=create_association_matrix_gene_go_term(genes, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
big_matrix_molecular_function=create_association_matrix_gene_go_term(genes, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
big_matrix_cell_comp=create_association_matrix_gene_go_term(genes, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)


dfs_bio_process={}
dfs_molecular_func={}
dfs_cell_comp={}
#9.1.) ENRICHMENT U genes
label=['U_dev_U_tis', 'U_dev_S_tis', 'U_dev_hS_tis', 'U_dev_null_tis']
for k in range(len(label)):

    genes_analyze=np.loadtxt(path_save_data+"%s.txt" %label[k], dtype=str)
    genes_label=label[k]
    label_analyzed=label[k]

    # #9.1.1.) BIOLOGICAL PROCESS
    
    matrix_subset_bio_process=create_association_matrix_gene_go_term(genes_analyze, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
    subset, enrich_subset_bio_process, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_bio_process, matrix_subset_bio_process, GO_bio_process_descrip, GO_bio_process,  label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_bio_process=pd.DataFrame()
    df_bio_process['subset']=subset
    df_bio_process['GO term']=enrich_go_term[index_sorted]
    df_bio_process['GO description']=enrich_subset_bio_process[index_sorted]
    df_bio_process['n genes']=n_genes[index_sorted]
    df_bio_process['n genes subset']=n_genes_subset[index_sorted]
    df_bio_process['p-value']=p_value[index_sorted]
    df_bio_process['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_bio_process[f'U_genes_cluster_{k}'] = df_bio_process
    
    del matrix_subset_bio_process, enrich_subset_bio_process, df_bio_process
    
    
    # #9.1.2.) MOLECULAR FUNCTION
    
    matrix_subset_molecular_func=create_association_matrix_gene_go_term(genes_analyze, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
    subset, enrich_subset_molecular_func, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_molecular_function, matrix_subset_molecular_func, GO_molecular_descrip, GO_molecular, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_molecular=pd.DataFrame()
    df_molecular['subset']=subset
    df_molecular['GO term']=enrich_go_term[index_sorted]
    df_molecular['GO description']=enrich_subset_molecular_func[index_sorted]
    df_molecular['n genes']=n_genes[index_sorted]
    df_molecular['n genes subset']=n_genes_subset[index_sorted]
    df_molecular['p-value']=p_value[index_sorted]
    df_molecular['corected p-value']=corrected_p_value[index_sorted]

    dfs_molecular_func[f'U_genes_cluster_{k}'] = df_molecular        

    del matrix_subset_molecular_func, enrich_subset_molecular_func, df_molecular
    
    
    # #9.1.3.) CELL COMPONENT
    
    matrix_subset_cell_comp=create_association_matrix_gene_go_term(genes_analyze, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
    subset, enrich_subset_cell_comp, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_cell_comp, matrix_subset_cell_comp, GO_cell_comp_descrip, GO_cell_comp, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_cell_comp=pd.DataFrame()
    df_cell_comp['subset']=subset
    df_cell_comp['GO term']=enrich_go_term[index_sorted]
    df_cell_comp['GO description']=enrich_subset_cell_comp[index_sorted]
    df_cell_comp['n genes']=n_genes[index_sorted]
    df_cell_comp['n genes subset']=n_genes_subset[index_sorted]
    df_cell_comp['p-value']=p_value[index_sorted]
    df_cell_comp['corected p-value']=corrected_p_value[index_sorted]

    dfs_cell_comp[f'U_genes_cluster_{k}'] = df_cell_comp        
    
    del matrix_subset_cell_comp, enrich_subset_cell_comp, df_cell_comp
    
    print(k)
        
        
#Save bio process
output_excel = path_save_data + 'enriched_U_genes_cluster_bio_process.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_bio_process.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
#save molecular func
output_excel = path_save_data + 'enriched_U_genes_cluster_molecular_func.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_molecular_func.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

#save cell comp
output_excel = path_save_data + 'enriched_U_genes_cluster_cell_comp.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_cell_comp.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)



dfs_bio_process={}
dfs_molecular_func={}
dfs_cell_comp={}
#9.2.) ENRICHMENT S genes
label=['S_dev_U_tis', 'S_dev_S_tis', 'S_dev_hS_tis', 'S_dev_null_tis']
for k in range(len(label)):
    
    genes_analyze=np.loadtxt(path_save_data+"%s.txt" %label[k], dtype=str)
    genes_label=label[k]
    label_analyzed=label[k]


    # #9.2.1.) BIOLOGICAL PROCESS    
    
    matrix_subset_bio_process=create_association_matrix_gene_go_term(genes_analyze, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
    subset, enrich_subset_bio_process, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_bio_process, matrix_subset_bio_process, GO_bio_process_descrip, GO_bio_process,  label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_bio_process=pd.DataFrame()
    df_bio_process['subset']=subset
    df_bio_process['GO term']=enrich_go_term[index_sorted]
    df_bio_process['GO description']=enrich_subset_bio_process[index_sorted]
    df_bio_process['n genes']=n_genes[index_sorted]
    df_bio_process['n genes subset']=n_genes_subset[index_sorted]
    df_bio_process['p-value']=p_value[index_sorted]
    df_bio_process['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_bio_process[f'S_genes_cluster_{k}'] = df_bio_process
    
    del matrix_subset_bio_process, enrich_subset_bio_process, df_bio_process
    
    
    # #9.2.2.) MOLECULAR FUNCTION    
    
    matrix_subset_molecular_func=create_association_matrix_gene_go_term(genes_analyze, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
    subset, enrich_subset_molecular_func, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_molecular_function, matrix_subset_molecular_func, GO_molecular_descrip, GO_molecular, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_molecular=pd.DataFrame()
    df_molecular['subset']=subset
    df_molecular['GO term']=enrich_go_term[index_sorted]
    df_molecular['GO description']=enrich_subset_molecular_func[index_sorted]
    df_molecular['n genes']=n_genes[index_sorted]
    df_molecular['n genes subset']=n_genes_subset[index_sorted]
    df_molecular['p-value']=p_value[index_sorted]
    df_molecular['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_molecular_func[f'S_genes_cluster_{k}'] = df_molecular        

    del matrix_subset_molecular_func, enrich_subset_molecular_func, df_molecular
    
    
    # #9.2.3.) CELL COMPONENT    
    
    matrix_subset_cell_comp=create_association_matrix_gene_go_term(genes_analyze, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
    subset, enrich_subset_cell_comp, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_cell_comp, matrix_subset_cell_comp, GO_cell_comp_descrip, GO_cell_comp, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_cell_comp=pd.DataFrame()
    df_cell_comp['subset']=subset
    df_cell_comp['GO term']=enrich_go_term[index_sorted]
    df_cell_comp['GO description']=enrich_subset_cell_comp[index_sorted]
    df_cell_comp['n genes']=n_genes[index_sorted]
    df_cell_comp['n genes subset']=n_genes_subset[index_sorted]
    df_cell_comp['p-value']=p_value[index_sorted]
    df_cell_comp['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_cell_comp[f'S_genes_cluster_{k}'] = df_cell_comp        
    
    del matrix_subset_cell_comp, enrich_subset_cell_comp, df_cell_comp
    
    print(k)
    
        
#Save bio process
output_excel = path_save_data + 'enriched_S_genes_cluster_bio_process.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_bio_process.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
#save molecular func
output_excel = path_save_data + 'enriched_S_genes_cluster_molecular_func.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_molecular_func.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

#save cell comp
output_excel = path_save_data + 'enriched_S_genes_cluster_cell_comp.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_cell_comp.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)



dfs_bio_process={}
dfs_molecular_func={}
dfs_cell_comp={}
#9.3.) ENRICHMENT hS genes
label=['hS_dev_U_tis', 'hS_dev_S_tis', 'hS_dev_hS_tis', 'hS_dev_null_tis']
for k in range(len(label)):
    
    genes_analyze=np.loadtxt(path_save_data+"%s.txt" %label[k], dtype=str)
    genes_label=label[k]
    label_analyzed=label[k]


    # #9.3.1.) BIOLOGICAL PROCESS    
    
    matrix_subset_bio_process=create_association_matrix_gene_go_term(genes_analyze, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
    subset, enrich_subset_bio_process, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_bio_process, matrix_subset_bio_process, GO_bio_process_descrip, GO_bio_process,  label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_bio_process=pd.DataFrame()
    df_bio_process['subset']=subset
    df_bio_process['GO term']=enrich_go_term[index_sorted]
    df_bio_process['GO description']=enrich_subset_bio_process[index_sorted]
    df_bio_process['n genes']=n_genes[index_sorted]
    df_bio_process['n genes subset']=n_genes_subset[index_sorted]
    df_bio_process['p-value']=p_value[index_sorted]
    df_bio_process['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_bio_process[f'hS_genes_cluster_{k}'] = df_bio_process
    
    del matrix_subset_bio_process, enrich_subset_bio_process, df_bio_process
    
    
    # #9.3.2.) MOLECULAR FUNCTION    
    
    matrix_subset_molecular_func=create_association_matrix_gene_go_term(genes_analyze, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
    subset, enrich_subset_molecular_func, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_molecular_function, matrix_subset_molecular_func, GO_molecular_descrip, GO_molecular, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_molecular=pd.DataFrame()
    df_molecular['subset']=subset
    df_molecular['GO term']=enrich_go_term[index_sorted]
    df_molecular['GO description']=enrich_subset_molecular_func[index_sorted]
    df_molecular['n genes']=n_genes[index_sorted]
    df_molecular['n genes subset']=n_genes_subset[index_sorted]
    df_molecular['p-value']=p_value[index_sorted]
    df_molecular['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_molecular_func[f'hS_genes_cluster_{k}'] = df_molecular        

    del matrix_subset_molecular_func, enrich_subset_molecular_func, df_molecular
    
    
    # #9.3.3.) CELL COMPONENT    
    
    matrix_subset_cell_comp=create_association_matrix_gene_go_term(genes_analyze, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
    subset, enrich_subset_cell_comp, enrich_go_term, p_value, n_genes, n_genes_subset=enrichement_go(big_matrix_cell_comp, matrix_subset_cell_comp, GO_cell_comp_descrip, GO_cell_comp, label_analyzed)
    
    index_sorted=np.argsort(p_value)
    
    corrected_p_value=false_discovery_control(p_value, method='bh')
    n_genes=np.array(n_genes, dtype=int)
    n_genes_subset=np.array(n_genes_subset, dtype=int)
    
    
    df_cell_comp=pd.DataFrame()
    df_cell_comp['subset']=subset
    df_cell_comp['GO term']=enrich_go_term[index_sorted]
    df_cell_comp['GO description']=enrich_subset_cell_comp[index_sorted]
    df_cell_comp['n genes']=n_genes[index_sorted]
    df_cell_comp['n genes subset']=n_genes_subset[index_sorted]
    df_cell_comp['p-value']=p_value[index_sorted]
    df_cell_comp['corected p-value']=corrected_p_value[index_sorted]
    
    dfs_cell_comp[f'hS_genes_cluster_{k}'] = df_cell_comp        
    
    del matrix_subset_cell_comp, enrich_subset_cell_comp, df_cell_comp
    
    print(k)
        
        
#Save bio process
output_excel = path_save_data + 'enriched_hS_genes_cluster_bio_process.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_bio_process.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
#save molecular func
output_excel = path_save_data + 'enriched_hS_genes_cluster_molecular_func.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_molecular_func.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

#save cell comp
output_excel = path_save_data + 'enriched_hS_genes_cluster_cell_comp.xlsx'
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    for sheet_name, df in dfs_cell_comp.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)















