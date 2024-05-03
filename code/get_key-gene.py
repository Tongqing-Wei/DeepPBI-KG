import os
import re
import pandas as pd
import numpy as np
from collections import Counter
import seaborn as sns
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	The path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name


if __name__ == '__main__':
    # file output path
    walkfile = walkFile('/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein')
    name = []
    gene = []
    for ls in walkfile:
        if ls[:3] == 'cds':
            name.append(ls[15:-6])
            res_dir = "/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/" + ls
            fo = open(res_dir, 'r+')
            flist = fo.readlines()
            for i, line in enumerate(flist):
                if line[0] == '>':
                    l = re.findall(r"\[(.+?)\]",flist[i])
                    if l[0][:4] == 'gene':
                        gene.append(l[0])
                    else:
                        gene.append(l[1]) 
                        
    gene_total = list(set(gene))
    KP_df = pd.DataFrame(index = name,columns = gene_total)
    KP_df = KP_df.replace(np.nan,0)
    ### Fill KP_df with the gene value of KP
    for ls in walkfile:
        gene_single = []
        if ls[:3] == 'cds':
            res_dir = "/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/" + ls
            fo = open(res_dir, 'r+')
            flist = fo.readlines()
            for i, line in enumerate(flist):
                if line[0] == '>':
                    l = re.findall(r"\[(.+?)\]",flist[i])
                    if l[0][:4] == 'gene':
                        gene_single.append(l[0])
                    else:
                        gene_single.append(l[1])  
            gene_single = Counter(gene_single) 
            gene_single = pd.Series(gene_single)
            for index in gene_single.index: 
                KP_df.loc[ls[15:-6], index] = gene_single[index]  
    KP_df.to_csv("/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/phage_KP_key-gene.csv")
    KP_df = np.array(KP_df).astype(float)
    fig_name = 'phage_KP_key-gene_heatmap.png'
    fig_path = '/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/' + fig_name
    fig = sns.heatmap(KP_df)
    heatmap = fig.get_figure()
    heatmap.savefig(fig_path, dpi = 400)

# Look at the unique genes of each sample
KP_dff = pd.DataFrame(KP_df, index = name, columns = gene_total )
KP_dff.loc["Column summation"] = KP_dff.apply(lambda x:x.sum())
condidate_gene = []
for i,j in enumerate(KP_dff.loc["Column summation"]):
    if j == 1.0 or j == 1:
        condidate_gene.append(KP_dff.loc["Column summation"].index[i])
KP_dff[condidate_gene].to_csv("/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/phage_KP_condidate-gene.csv")
KP_df_condidate = np.array(KP_dff[condidate_gene]).astype(float)
fig_name = 'phage_KP_condidate-gene_heatmap.png'
fig_path = '/bios-store1/home/TongqingWei/KP_phage/phage_KP_marker_protein/' + fig_name
fig = sns.heatmap(KP_df_condidate)
heatmap = fig.get_figure()
heatmap.savefig(fig_path, dpi = 400)