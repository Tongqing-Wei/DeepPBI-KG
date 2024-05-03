import os
import re
import pandas as pd
import numpy as np
from collections import Counter
import argparse

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

def get_gene(path):
    walkfile = walkFile(path)
    name = []
    gene = []
    for ls in walkfile:
        if ls[:3] == 'cds':
            name.append(ls[4:-6])
            res_dir = path + os.sep + ls
            fo = open(res_dir, 'r+')
            flist = fo.readlines()
            for i, line in enumerate(flist):
                if line[0] == '>' and 'hypothetical protein' not in line:
                    l = re.findall(r"\[(.+?)\]",flist[i])
                    if l[0][:4] == 'gene':
                        gene.append(l[0])
                    else:
                        gene.append(l[1])     
    return name, gene

if __name__ == '__main__':
    # The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help = 'input phage(host) cds dna fasta file path')
    parser.add_argument('--type', type=str, help = 'process phage or host')
    parser.add_argument('--output', type=str, help = 'output phage(host) key gene vector hash table')
    opt = parser.parse_args()    
        
    # Gets the key_gene list for phage(host)
    name_list, gene_list = get_gene(opt.input)
    gene_total = list(set(gene_list))
    KP_df = pd.DataFrame(index = name_list,columns = gene_total)
    KP_df = KP_df.replace(np.nan,0)
    ### Fill KP_df with the gene value of KP
    walkfile = walkFile(opt.input)
    for ls in walkfile:
        gene_single = []
        if ls[:3] == 'cds':
            res_dir = opt.input + os.sep + ls
            fo = open(res_dir, 'r+')
            flist = fo.readlines()
            for i, line in enumerate(flist):
                if line[0] == '>' and 'hypothetical protein' not in line:
                    l = re.findall(r"\[(.+?)\]",flist[i])
                    if l[0][:4] == 'gene':
                        gene_single.append(l[0])
                    else:
                        gene_single.append(l[1])  
            gene_single = Counter(gene_single) 
            gene_single = pd.Series(gene_single)
            for index in gene_single.index: 
                KP_df.loc[ls[4:-6], index] = gene_single[index]  
    KP_df.to_csv(opt.output + os.sep + opt.type + "_key_gene_vector.csv", index = True)
