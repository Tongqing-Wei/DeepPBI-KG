import os
import pandas as pd
import numpy as np
import re
from Bio import SeqIO

def format_fasta(ana, seq, num):
    """
    Format the text in fasta format
    :param ana: annotation information
    :param seq: sequence
    :param num: The number of characters when the sequence is wrapped
    :return: fasta format text
    """
    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq += "\n"
    return '>' + ana + '\n' + format_seq + "\n"


def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	the path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name



if __name__ == '__main__':
    rbp = pd.read_csv('.../Pipeline/raw_Data/rbp_phage.csv', header = 0)
    dic = {}
    for i in range(len(rbp)):
        pre = rbp.iloc[i,0].split('.')[0]
        if pre in dic:
            dic[pre].append(rbp.iloc[i,0])
        else:
            dic[pre] = [rbp.iloc[i,0]]


    walkfile = walkFile('.../Pipeline/raw_Data/phage_protein')
    gene = []
    for ls in walkfile:
        if ls[14:-6][ls[14:-6].find('_') + 1:] in dic:
            target = dic[ls[14:-6][ls[14:-6].find('_') + 1:]]
        else:
            continue
        f = open('.../Pipeline/tmp_file/phage_rbp_protein/' + os.sep + ls, 'w')
        cds_fasta = ''
        res_dir = ".../Pipeline/raw_Data/phage_protein/" + ls
        records = (r for r in SeqIO.parse(res_dir, "fasta"))  
        for i in records:
            if i.description.split(' ')[0] in target:
                cds_fasta += format_fasta(i.description, i.seq, 70)
                if re.findall(r'\[(.+?)\]', i.description)[0][:4] == 'gene':
                    gene.append(re.findall(r'\[(.+?)\]', i.description)[0])
                elif re.findall(r'\[(.+?)\]', i.description)[1] != 'protein=hypothetical protein':
                    gene.append(re.findall(r'\[(.+?)\]', i.description)[1])
        f.write(cds_fasta)
        f.close()

    gene = pd.DataFrame(list(set(gene))).to_csv('.../Pipeline/tmp_file/phage_gene.csv', index = False)