import re
import pandas as pd
import numpy as np

def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root    	the path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	the list of files in this folder
        # traversal file
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name


walkfile = walkFile('.../Pipeline/raw_Data/host_blast')
for i in range(0, len(walkfile), 2):
    ls = walkfile[i]
    lt = walkfile[i+1]
    data = pd.read_csv('.../raw_Data/host_blast' + os.sep + lt, sep = '\t', header = None)
    ref = pd.read_csv('.../Pipeline/raw_Data/host_blast' + os.sep + ls, sep = '\t', header = None)

    ind = []
    for i in range(len(ref)):
        if ref.iloc[i,0][:3] == 'CDD':
            ind.append(ref.iloc[i,0])

    dic = {}
    for i in range(len(ind)):
        dic[ind[i][:10]] = ' '.join(re.findall(r'(.+?)\  ', ind[i])[:2])

    for i in range(len(data)):
        if data.iloc[i,1] in dic.keys():
            data.iloc[i,10] = dic[data.iloc[i,1]]
    data.iloc[:,[0,1,6,7,8,9,10]].to_csv('.../Pipeline/raw_Data/host_blast' + os.sep + lt, sep = '\t', index = False)