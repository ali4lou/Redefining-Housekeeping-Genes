# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:31:55 2025

Comparison of marker genes, differential genes 

@author: Alicia
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.stats import hypergeom

def test_hipergeometrico(total_poblacion, exitos_totales, tamaño_muestra, exitos_observados, verbose=True):

    N = total_poblacion
    K = exitos_totales
    n = tamaño_muestra
    k = exitos_observados

    # P-valor (cola superior)
    p_valor = hypergeom.sf(k-1, N, K, n)
    
    print(p_valor)
    return p_valor


path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'

#1.) We read the bulk data
f=open(path_save_data+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)

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


#======================================================================================
#======================================================================================
#TRANSCRIPTION FACTORS
#2.) We read the transcription factors
f=open(path_save_data+'Danio_rerio_TF', 'r')
txt = f.read()
df_TF = txt.split('\n')
del txt, f
df_TF=np.delete(df_TF, len(df_TF)-1)
df_TF=np.array(df_TF)

genes_TF=[]
for elemento in df_TF:
    partes = elemento.split("\t")

    if partes[2].startswith('ENS'):
        genes_TF.append(partes[2])

genes_TF=np.array(genes_TF)

print('there are genes TF', len(genes_TF), 'in the datasheet')

#3.) We find the common genes
common_gene_names, ind_total, ind_df=np.intersect1d(genes_bulk, genes_TF, return_indices=True)
print('We just find', len(common_gene_names), 'in our data')

#4.) We have to classify them
U_genes_marker=np.intersect1d(U_id_bulk, common_gene_names)    
S_genes_marker=np.intersect1d(S_id_bulk, common_gene_names)    
hS_genes_marker=np.intersect1d(hS_id_bulk, common_gene_names)    

#figure
frac_genes=np.array([(len(common_gene_names)-(len(U_genes_marker)+len(S_genes_marker)+len(hS_genes_marker)))/len(common_gene_names), 
                     len(U_genes_marker)/len(common_gene_names), 
                     len(S_genes_marker)/len(common_gene_names), 
                     len(hS_genes_marker)/len(common_gene_names)])


frac_genes_TF=np.array([len(U_genes_marker)/len(common_gene_names), 
                     len(S_genes_marker)/len(common_gene_names), 
                     len(hS_genes_marker)/len(common_gene_names), 
                     (len(common_gene_names)-(len(U_genes_marker)+len(S_genes_marker)+len(hS_genes_marker)))/len(common_gene_names)])

groups=['NE', 'U', 'S', 'hS']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=15, fontweight='bold')
plt.yticks(fontsize=14)
# plt.title('Single Cell Base', fontsize=16, fontweight='bold')
plt.ylabel('% TF', fontsize=16, fontweight='bold')
for i in range(1, len(groups)):
    plt.text(i-0.25, frac_genes[i]-0.07, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
plt.text(-0.25, 0.45, "# genes = %d" %len(common_gene_names), fontsize=9, color="black",
        bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"))
# --- NUEVO: Asterisco encima de la barra 'S' ---
if "S" in groups:
    idx_S = groups.index("S")  # posición de la barra 'S'
    plt.text(idx_S, frac_genes[idx_S]-0.03, "*",
             ha='center', va='bottom', fontsize=18, color='red')
plt.savefig(path_save_data+'gene_frac_TF.png', dpi=600, bbox_inches='tight')
plt.show()



print('U genes', len(U_genes_marker), len(U_id_bulk))
p_valor_U_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(U_id_bulk),
    exitos_observados=len(U_genes_marker)
)

print('S genes', len(S_genes_marker), len(S_id_bulk))
p_valor_S_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(S_id_bulk),
    exitos_observados=len(S_genes_marker)
)

print('hS genes', len(hS_genes_marker), len(hS_id_bulk))
p_valor_hS_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(hS_id_bulk),
    exitos_observados=len(hS_genes_marker)
)


#COFACTORS
#5.) We read the cofactors of TF
f=open(path_save_data+'Danio_rerio_Cof', 'r')
txt = f.read()
df_TF = txt.split('\n')
del txt, f
df_TF=np.delete(df_TF, len(df_TF)-1)
df_TF=np.array(df_TF)

genes_TF=[]
for elemento in df_TF:
    partes = elemento.split("\t")

    if partes[2].startswith('ENS'):
        genes_TF.append(partes[2])

genes_TF=np.array(genes_TF)
        

print('there are cofactors TF', len(genes_TF), 'in the datasheet')

#3.) We find the common genes
common_gene_names, ind_total, ind_df=np.intersect1d(genes_bulk, genes_TF, return_indices=True)
print('We just find', len(common_gene_names), 'in our data')

#4.) We have to classify them
U_genes_marker=np.intersect1d(U_id_bulk, common_gene_names)    
S_genes_marker=np.intersect1d(S_id_bulk, common_gene_names)    
hS_genes_marker=np.intersect1d(hS_id_bulk, common_gene_names)    

#figure
frac_genes_cof=np.array([len(U_genes_marker)/len(common_gene_names), 
                     len(S_genes_marker)/len(common_gene_names), 
                     len(hS_genes_marker)/len(common_gene_names), 
                     (len(common_gene_names)-(len(U_genes_marker)+len(S_genes_marker)+len(hS_genes_marker)))/len(common_gene_names)])

frac_genes=np.array([(len(common_gene_names)-(len(U_genes_marker)+len(S_genes_marker)+len(hS_genes_marker)))/len(common_gene_names), 
                     len(U_genes_marker)/len(common_gene_names), 
                     len(S_genes_marker)/len(common_gene_names), 
                     len(hS_genes_marker)/len(common_gene_names)])
groups=['NE', 'U', 'S', 'hS']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=15, fontweight='bold')
plt.yticks(fontsize=14)
# plt.title('Single Cell Base', fontsize=16, fontweight='bold')
plt.ylabel('% cofactors', fontsize=16, fontweight='bold')
for i in range(1, len(groups)):
    if i==3:
        plt.text(i-0.25, frac_genes[i], np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='black')
    else:
        plt.text(i-0.25, frac_genes[i]-0.07, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
plt.text(1.85, 0.65, "# genes = %d" %len(common_gene_names), fontsize=9, color="black",
        bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"))
if "U" in groups:
    idx_U = groups.index("U")  # posición de la barra 'S'
    plt.text(idx_U, frac_genes[idx_U]-0.04, "*",
             ha='center', va='bottom', fontsize=18, color='red')
plt.savefig(path_save_data+'gene_frac_cofactors_TF.png', dpi=600, bbox_inches='tight')
plt.show()



print('U genes', len(U_genes_marker), len(U_id_bulk))
p_valor_U_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(U_id_bulk),
    exitos_observados=len(U_genes_marker)
)

print('S genes', len(S_genes_marker), len(S_id_bulk))
p_valor_S_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(S_id_bulk),
    exitos_observados=len(S_genes_marker)
)

print('hS genes', len(hS_genes_marker), len(hS_id_bulk))
p_valor_hS_TF = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(common_gene_names),
    tamaño_muestra=len(hS_id_bulk),
    exitos_observados=len(hS_genes_marker))


