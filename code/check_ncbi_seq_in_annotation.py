# Check whether the RBP found from NCBI is in phage_protein and extract the ID of the corresponding config
import pandas as pd
import os
import re
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt

def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	the path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files in this folder
        # traversal file
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name



tail = pd.read_csv('.../Pipeline/raw_Data/receptor.txt', sep = '\t', header = None)
tail = tail.drop_duplicates()
xl = pd.read_excel('.../Pipeline/raw_Data/RF_all_sample.xlsx', sheet_name = 'ID', header = 0)
dic = {}
for i in range(len(tail)):
    p = xl[xl['host_name'] == tail.iloc[i,0]].index.tolist()[0]
    tail.iloc[i,0] = xl.iloc[p,1]
for i in range(len(tail)):
    if tail.iloc[i,0] not in dic.keys():
        dic[tail.iloc[i,0]] = [tail.iloc[i,1]]
    else:
        dic[tail.iloc[i,0]].append(tail.iloc[i,1])

        
# 文件输出路径
walkfile = walkFile('.../Pipeline/raw_Data/host_protein')
f = open('.../Pipeline/raw_Data/rbp_host.txt', 'w')
name = []
for k in dic.keys():
    for ls in walkfile:
        if ls[13:-6][ls[13:-6].find('_') + 1:] == k:
            res_dir = ".../Pipeline/raw_Data/host_protein/" + ls
            records = (r for r in SeqIO.parse(res_dir, "fasta"))
            test = []
            for i in records:
                test.append(i)
            for v in dic[k]:
                for j in test:
                    if j.seq == v:
                        f.writelines('%s\t%s\n' % (ls[13:-6][ls[13:-6].find('_') + 1:], j.seq))
                        name.append(j.name)
name = pd.DataFrame(name)
name.to_csv('.../Pipeline/raw_Data/rbp_host.csv', index = False)
f.close()