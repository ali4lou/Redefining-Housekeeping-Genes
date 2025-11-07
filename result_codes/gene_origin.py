# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:10:27 2025

Gene origin

@author: Alicia
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 


def age_vs_max_exp(genes_analyze):
    genes_analyze_age, i, ind_age=np.intersect1d(genes_analyze, genes_df, return_indices=True)
    age_genes_analyze=age_df[ind_age]
    age_genes_analyze = ['4290' if x == '>4290' else x for x in age_genes_analyze]
    age_genes_analyze=np.array(age_genes_analyze, dtype=int)
    
    g, ind_bulk, i=np.intersect1d(genes_bulk, genes_analyze_age, return_indices=True)
    bulk_genes_analyze=mean_bulk_matrix[ind_bulk, :]
    
    stage_of_max_exp=np.zeros(len(genes_analyze_age))
    for i in range(len(genes_analyze_age)):
        ind_max=np.argmax(bulk_genes_analyze[i, :])
        stage_of_max_exp[i]=ind_max
    
    #sliding window
    sw=200
    # sw=int(len(age_genes_analyze)/8)
    stage_of_max_exp_sorted=np.sort(stage_of_max_exp)
    index_stage_sort=np.argsort(stage_of_max_exp)
    age_per_gen_sorted=np.zeros(len(stage_of_max_exp))
    for i in range(len(index_stage_sort)):
        age_per_gen_sorted[i]=age_genes_analyze[int(index_stage_sort[i])]
    

    serie1 = pd.Series(stage_of_max_exp_sorted)
    serie2 = pd.Series(age_per_gen_sorted)
    
    sw_stage = serie1.rolling(window=sw, center=False).mean()
    sw_age = serie2.rolling(window=sw, center=False).mean()
    
    return sw_stage, sw_age, age_genes_analyze




path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'




#1.) We read the bulk data
mean_bulk_matrix= np.genfromtxt(path+'mean_bulk_matrix.txt', delimiter=',', dtype=None, encoding=None)

f=open(path+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)


f=open(path+'time_bulk.txt', 'r')
txt = f.read()
time_bulk = txt.split('\n')
del txt, f
time_bulk=np.delete(time_bulk, len(time_bulk)-1)
time_bulk=np.array(time_bulk)



f=open(path+'hS_id_bulk.txt', 'r')
txt = f.read()
hS_id_bulk = txt.split('\n')
del txt, f
hS_id_bulk=np.delete(hS_id_bulk, len(hS_id_bulk)-1)
hS_id_bulk=np.array(hS_id_bulk)

f=open(path+'S_id_bulk.txt', 'r')
txt = f.read()
S_id_bulk = txt.split('\n')
del txt, f
S_id_bulk=np.delete(S_id_bulk, len(S_id_bulk)-1)
S_id_bulk=np.array(S_id_bulk)

f=open(path+'U_id_bulk.txt', 'r')
txt = f.read()
U_id_bulk = txt.split('\n')
del txt, f
U_id_bulk=np.delete(U_id_bulk, len(U_id_bulk)-1)
U_id_bulk=np.array(U_id_bulk)


#2.) We read the gene origin datasheet
df_gene_origin=pd.read_csv(path+"gene_origin.csv")
genes_df=list(df_gene_origin['ensembl_gene_id'])
age_df=list(df_gene_origin['gene_age'])

genes_df=np.array(genes_df)
age_df=np.array(age_df)

age_df_unique=np.unique(age_df)

#3.) We classify them
U_genes_age, i, ind_age_U=np.intersect1d(U_id_bulk, genes_df, return_indices=True)    
age_U=age_df[ind_age_U]
age_U = ['4290' if x == '>4290' else x for x in age_U]
age_U=np.array(age_U, dtype=int)

S_genes_age, i, ind_age_S=np.intersect1d(S_id_bulk, genes_df, return_indices=True)   
age_S=age_df[ind_age_S]
age_S = ['4290' if x == '>4290' else x for x in age_S]
age_S=np.array(age_S, dtype=int)
 

hS_genes_age, i, ind_age_hS=np.intersect1d(hS_id_bulk, genes_df, return_indices=True)
age_hS=age_df[ind_age_hS]
age_hS = ['4290' if x == '>4290' else x for x in age_hS]
age_hS=np.array(age_hS, dtype=int)

print('median:', np.median(age_U), np.median(age_S), np.median(age_hS))


#3.) Boxplot
data=[age_U, age_S, age_hS]
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

# plt.yscale(log)
plt.ylabel("Age (myr)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.02 * y_max  
line_heights = [y_max + 0.1*y_max, y_max + 0.2*y_max, y_max + 0.3*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
plt.savefig(path+'age_boxplot.png', dpi=600, bbox_inches='tight')

plt.show()


#3.) Boxplot WITHOUT OUTLIERS
data=[age_U, age_S, age_hS]
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

# plt.yscale(log)
plt.ylabel("Age (myr)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.ylim(0, 2800)

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
# y_max = max([max(d) for d in data])
y_max=2000
h = 0.02 * y_max  
line_heights = [y_max + 0.1*y_max, y_max + 0.2*y_max, y_max + 0.3*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
plt.savefig(path+'age_boxplot_with_lim_y.png', dpi=600, bbox_inches='tight')

plt.show()



#4.) KS test
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(age_U, age_S)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(age_U, age_hS)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(age_S, age_hS)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')
    

#5.) Age related to developmental expression
g, ind_bulk_hS, i=np.intersect1d(genes_bulk, hS_genes_age, return_indices=True)
g, ind_bulk_S, i=np.intersect1d(genes_bulk, S_genes_age, return_indices=True)
g, ind_bulk_U, i=np.intersect1d(genes_bulk, U_genes_age, return_indices=True)

bulk_U=mean_bulk_matrix[ind_bulk_U, :]
bulk_S=mean_bulk_matrix[ind_bulk_S, :]
bulk_hS=mean_bulk_matrix[ind_bulk_hS, :]

#5.1.) For each gene we find the stage of maximal exp 
stage_of_max_exp_U=np.zeros(len(U_genes_age))
U_old_genes_hourglass=[]
for i in range(len(U_genes_age)):
    ind_max=np.argmax(bulk_U[i, :])
    stage_of_max_exp_U[i]=ind_max
    
    #we search the U genes with max expression in stages 10.3, 16, 19 and 24 hpf
    if (ind_max==8) & (age_U[i]==4290): #time 10.3hpf
        # print(ind_max, age_U[i])
        U_old_genes_hourglass.append(U_genes_age[i])
    if (ind_max==9) & (age_U[i]==4290): #time 16 hpf
        U_old_genes_hourglass.append(U_genes_age[i])
    if (ind_max==10) & (age_U[i]==4290): #time 19 hpf
        U_old_genes_hourglass.append(U_genes_age[i])
    if (ind_max==11) & (age_U[i]==4290): #time 24 hpf
        U_old_genes_hourglass.append(U_genes_age[i])
        
U_old_genes_hourglass=np.array(U_old_genes_hourglass)
np.savetxt(path+'U_old_genes_hourglass.txt', U_old_genes_hourglass, fmt='%s') 


#5.2.) For each gene we find the stage of maximal exp 
stage_of_max_exp_S=np.zeros(len(S_genes_age))
S_genes_gastrulation=[]
for i in range(len(S_genes_age)):
    ind_max=np.argmax(bulk_S[i, :])
    stage_of_max_exp_S[i]=ind_max
    
    #we search the U genes with max expression in stages 10.3, 16, 19 and 24 hpf
    if (ind_max==5) & (age_S[i]<600): #time 10.3hpf
        # print(ind_max, age_U[i])
        S_genes_gastrulation.append(S_genes_age[i])
    if (ind_max==6) & (age_S[i]<600): #time 16 hpf
        S_genes_gastrulation.append(S_genes_age[i])
    
        
S_genes_gastrulation=np.array(S_genes_gastrulation)
np.savetxt(path+'S_genes_gastrulation.txt', S_genes_gastrulation, fmt='%s') 


#5.3.) For each gene we find the stage of maximal exp 
stage_of_max_exp_hS=np.zeros(len(hS_genes_age))
hS_genes_gastrulation=[]
for i in range(len(hS_genes_age)):
    ind_max=np.argmax(bulk_hS[i, :])
    stage_of_max_exp_hS[i]=ind_max
    
    #we search the U genes with max expression in stages 10.3, 16, 19 and 24 hpf
    if (ind_max==5) & (age_hS[i]<600): #time 10.3hpf
        # print(ind_max, age_U[i])
        hS_genes_gastrulation.append(hS_genes_age[i])
    if (ind_max==6) & (age_hS[i]<600): #time 16 hpf
        hS_genes_gastrulation.append(hS_genes_age[i])
    
        
hS_genes_gastrulation=np.array(hS_genes_gastrulation)
np.savetxt(path+'hS_genes_gastrulation.txt', hS_genes_gastrulation, fmt='%s') 


stage_of_max_exp_S=np.zeros(len(S_genes_age))
for i in range(len(S_genes_age)):
    ind_max=np.argmax(bulk_S[i, :])
    stage_of_max_exp_S[i]=ind_max
    
stage_of_max_exp_hS=np.zeros(len(hS_genes_age))
for i in range(len(hS_genes_age)):
    ind_max=np.argmax(bulk_hS[i, :])
    stage_of_max_exp_hS[i]=ind_max



#6.) Sliding window
sw=400
#U genes
stage_of_max_exp_sorted_U=np.sort(stage_of_max_exp_U)
index_stage_sort_U=np.argsort(stage_of_max_exp_U)
age_per_gen_sorted_U=np.zeros(len(stage_of_max_exp_U))
for i in range(len(index_stage_sort_U)):
    age_per_gen_sorted_U[i]=age_U[int(index_stage_sort_U[i])]
    

serie1U = pd.Series(stage_of_max_exp_sorted_U)
serie2U = pd.Series(age_per_gen_sorted_U)

sw_stage_U = serie1U.rolling(window=sw, center=False).mean()
sw_pleio_U = serie2U.rolling(window=sw, center=False).mean()

#S genes
stage_of_max_exp_sorted_S=np.sort(stage_of_max_exp_S)
index_stage_sort_S=np.argsort(stage_of_max_exp_S)
age_per_gen_sorted_S=np.zeros(len(stage_of_max_exp_S))
for i in range(len(index_stage_sort_S)):
    age_per_gen_sorted_S[i]=age_S[int(index_stage_sort_S[i])]

serie1S = pd.Series(stage_of_max_exp_sorted_S)
serie2S = pd.Series(age_per_gen_sorted_S)


sw_stage_S = serie1S.rolling(window=sw, center=False).mean()
sw_pleio_S = serie2S.rolling(window=sw, center=False).mean()

#hS genes
stage_of_max_exp_sorted_hS=np.sort(stage_of_max_exp_hS)
index_stage_sort_hS=np.argsort(stage_of_max_exp_hS)
age_per_gen_sorted_hS=np.zeros(len(stage_of_max_exp_hS))
for i in range(len(index_stage_sort_hS)):
    age_per_gen_sorted_hS[i]=age_hS[int(index_stage_sort_hS[i])]

serie1hS = pd.Series(stage_of_max_exp_sorted_hS)
serie2hS = pd.Series(age_per_gen_sorted_hS)

sw_stage_hS = serie1hS.rolling(window=sw, center=False).mean()
sw_pleio_hS = serie2hS.rolling(window=sw, center=False).mean()


# # figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(sw_stage_U, sw_pleio_U, s=1, color='royalblue', label='U')
plt.scatter(sw_stage_S, sw_pleio_S, s=1, color='mediumturquoise', label='S')
plt.scatter(sw_stage_hS, sw_pleio_hS, s=1, color='tomato', label='hS')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.savefig(path+'age_vs_max_stage_exp_sw_%d.png' %sw, dpi=600, bbox_inches='tight')
plt.show()


#6.1.) Scatter plot with density points
from scipy.stats import gaussian_kde

#U genes
x=stage_of_max_exp_sorted_U
y=age_per_gen_sorted_U

# Density computation of each point
xy = np.vstack([x, y])
density = gaussian_kde(xy)(xy)

#Figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(x, y, c=density, cmap="gnuplot_r", s=10)
plt.colorbar(label="Density")
plt.title('U genes', fontsize=15, fontweight='bold')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.savefig(path+'age_vs_max_stage_exp_U_genes.png', dpi=600, bbox_inches='tight')
plt.show()

#S genes
x=stage_of_max_exp_sorted_S
y=age_per_gen_sorted_S

# Density computation of each point
xy = np.vstack([x, y])
density = gaussian_kde(xy)(xy)

#Figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(x, y, c=density, cmap="gnuplot_r", s=10)
plt.colorbar(label="Density")
plt.title('S genes', fontsize=15, fontweight='bold')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.savefig(path+'age_vs_max_stage_exp_S_genes.png', dpi=600, bbox_inches='tight')
plt.show()

#hS genes
x=stage_of_max_exp_sorted_hS
y=age_per_gen_sorted_hS

# Density computation of each point
xy = np.vstack([x, y])
density = gaussian_kde(xy)(xy)

#Figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(x, y, c=density, cmap="gnuplot_r", s=10)
plt.colorbar(label="Density")
plt.title('hS genes', fontsize=15, fontweight='bold')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.savefig(path+'age_vs_max_stage_exp_hS_genes.png', dpi=600, bbox_inches='tight')
plt.show()



#===================================================================================
#7.) Sliding window for each subclass

#7.1.) We read the subclasses
label=['U_dev_U_tis', 'U_dev_S_tis', 'U_dev_hS_tis', 
       'S_dev_U_tis', 'S_dev_S_tis', 'S_dev_hS_tis',
       'hS_dev_U_tis', 'hS_dev_S_tis', 'hS_dev_hS_tis']
for k in (label):
    genes_analyze=np.loadtxt(path+"%s.txt" %k, dtype=str)
    genes_label=k
    globals()[f"{k}_sw_stage"], globals()[f"{k}_sw_age"], globals()[f"{k}_age"]=age_vs_max_exp(genes_analyze)
    


plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(U_dev_U_tis_sw_stage, U_dev_U_tis_sw_age, s=1, color='#0D47A1', label='U-Ut')
plt.scatter(U_dev_S_tis_sw_stage, U_dev_S_tis_sw_age, s=1, color='#1976D2', label='U-St')
plt.scatter(U_dev_hS_tis_sw_stage, U_dev_hS_tis_sw_age, s=1, color='#5C6BC0', label='U-hSt')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.title('U genes', fontsize=13, fontweight='bold')
plt.savefig(path+'age_vs_max_stage_exp_sw_U_genes_dev.png', dpi=600, bbox_inches='tight')
plt.show()


plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(S_dev_U_tis_sw_stage, S_dev_U_tis_sw_age, s=1, color='#7E57C2', label='S-Ut')
plt.scatter(S_dev_S_tis_sw_stage, S_dev_S_tis_sw_age, s=1, color='#AB47BC', label='S-St')
plt.scatter(S_dev_hS_tis_sw_stage, S_dev_hS_tis_sw_age, s=1, color='#EC407A', label='S-hSt')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.title('S genes', fontsize=13, fontweight='bold')
plt.savefig(path+'age_vs_max_stage_exp_sw_S_genes_dev.png', dpi=600, bbox_inches='tight')
plt.show()


plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(hS_dev_U_tis_sw_stage, hS_dev_U_tis_sw_age, s=1, color='mediumvioletred', label='hS-Ut')
plt.scatter(hS_dev_S_tis_sw_stage, hS_dev_S_tis_sw_age, s=1, color='hotpink', label='hS-St')
plt.scatter(hS_dev_hS_tis_sw_stage, hS_dev_hS_tis_sw_age, s=1, color='pink', label='hS-hSt')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.title('hS genes', fontsize=13, fontweight='bold')
plt.savefig(path+'age_vs_max_stage_exp_sw_hS_genes_dev_sw.png', dpi=600, bbox_inches='tight')
plt.show()


#figure of gene age distrib
data=[U_dev_U_tis_age, U_dev_S_tis_age, U_dev_hS_tis_age, 
      S_dev_U_tis_age, S_dev_S_tis_age, S_dev_hS_tis_age, 
      hS_dev_U_tis_age, hS_dev_S_tis_age, hS_dev_hS_tis_age]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U-Ut", "U-St", "U-hSt", 
                               "S-Ut", "S-St", "S-hSt",
                               "hS-Ut", "hS-St", "hS-hSt"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = [
    "#0D47A1",  # Azul profundo
    "#1976D2",  # Azul medio
    "#5C6BC0",  # Azul violáceo
    "#7E57C2",  # Violeta suave
    "#AB47BC",  # Púrpura
    "#EC407A",  # Fucsia fuerte
    "mediumvioletred",  # Rosa medio
    "hotpink",  # Rosa claro
    "pink" ]  # Rosa pálido

for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

# plt.yscale(log)
plt.ylabel("Age (myr)", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold', rotation=90)
plt.yticks(fontsize=13)

plt.ylim(0, 2800)

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '**']   
# y_max = max([max(d) for d in data])
y_max=2000
h = 0.02 * y_max  
line_heights = [y_max + 0.1*y_max, y_max + 0.2*y_max, y_max + 0.3*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
comparisons = [(3,4), (3,5), (4,5)]  
labels = ['***', '****', '**']   
h = 0.02 * y_max  
line_heights = [y_max + 0.1*y_max, y_max + 0.2*y_max, y_max + 0.3*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
comparisons = [(6,7), (6,8), (7,8)]  
labels = ['', '**', '**']   
h = 0.02 * y_max  
line_heights = [y_max + 0.1*y_max, y_max + 0.2*y_max, y_max + 0.3*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    

plt.savefig(path+'age_boxplot_subclasses_with_lim_y.png', dpi=600, bbox_inches='tight')

plt.show()




#4.) KS test
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(U_dev_U_tis_age, U_dev_S_tis_age)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(U_dev_U_tis_age, U_dev_hS_tis_age)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(U_dev_S_tis_age, U_dev_hS_tis_age)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')
    
ks_U_S, pvalue_U_S = stats.ks_2samp(S_dev_U_tis_age, S_dev_S_tis_age)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(S_dev_U_tis_age, S_dev_hS_tis_age)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(S_dev_S_tis_age, S_dev_hS_tis_age)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')

ks_U_S, pvalue_U_S = stats.ks_2samp(hS_dev_U_tis_age, hS_dev_S_tis_age)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(hS_dev_U_tis_age, hS_dev_hS_tis_age)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(hS_dev_S_tis_age, hS_dev_hS_tis_age)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')








