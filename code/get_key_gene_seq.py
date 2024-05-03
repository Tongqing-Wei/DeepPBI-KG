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


import re
import os
from Bio import SeqIO
import argparse

#The parameters are defined and encapsulated
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help = 'input phage(host) cds protein fasta file path')
parser.add_argument('--type', type=str, help = 'phage or host')
parser.add_argument('--output', type=str, help = 'output phage(host) all key gene seq fasta file')
opt = parser.parse_args()     
walkfile = walkFile(opt.input) 
lt = []
for ls in walkfile:
    if ls[-5:] == 'fasta':
        gb = list(SeqIO.parse(opt.input + os.sep + ls, 'fasta'))
        for i in gb:
            if re.findall(r'\[(.+?)\]', i.description)[0][:4] == 'gene':
                lt.append(re.findall(r'\[(.+?)\]', i.description)[0])
            elif re.findall(r'\[(.+?)\]', i.description)[1] != 'protein=hypothetical protein':
                lt.append(re.findall(r'\[(.+?)\]', i.description)[1])

ll = list(set(lt))

import re
import os
from Bio import SeqIO
walkfile = walkFile(opt.input) 
#f = open('/bios-store1/home/TongqingWei/database/host.txt', 'w')
dic = {}
for k in ll:
    for ls in walkfile:
        if ls[-5:] == 'fasta':
            gb = list(SeqIO.parse(opt.input + os.sep + ls, 'fasta'))
            for i in gb:
                if re.findall(r'\[(.+?)\]', i.description)[0][:4] == 'gene' and re.findall(r'\[(.+?)\]', i.description)[0] == k:
                    dic[k] = str(i.seq)
                    #f.writelines('>%s\n' % re.findall(r'\[(.+?)\]', i.description)[0])
                    #f.writelines('%s\n' % i.seq)
                elif re.findall(r'\[(.+?)\]', i.description)[1] == k:
                    dic[k] = str(i.seq)
                    #f.writelines('>%s\n' % re.findall(r'\[(.+?)\]', i.description)[1])
                    #.writelines('%s\n' % i.seq)
                
#f.close()

f = open(opt.output + os.sep + opt.type + '_key_gene.txt', 'w')
for k, v in dic.items():
    f.writelines('>%s\n' % k)
    f.writelines('%s\n' % v)   
f.close()