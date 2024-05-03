from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import re

def NCBI_seq(key_words):
    
    Entrez.email = "20210700125@fudan.edu.cn"     # Always tell NCBI who you are
    # find corresponding ID    
    handle = Entrez.esearch(db="protein", term="%s receptor" % key_words)
    record = Entrez.read(handle)
    gi_list = record["IdList"]

    # get Genbank information
    test = []
    seqindex = []
    seq = []
    label = []
    dic = []
    gi_str = ",".join(gi_list)
    handle1 = Entrez.efetch(db="protein", id=gi_str, rettype="gb", retmode="text")
    records = SeqIO.parse(handle1, "gb")
    for i in records:
        test.append(i)
    for i in test:
        if 'receptor' in re.findall(r'(.+?)\ ', i.description):
            label.append(i.description)
    label = list(set(label))
    for i in label:
        dic.append({j:len(ls.seq) for j,ls in enumerate(test) if ls.description == i})
    for i in dic:
        seqindex.append(max(i, key = i.get))
    for i in seqindex:
        seq.append(test[i].seq)
        
    return seq
        
    
if __name__ == '__main__':
    xl = pd.read_excel('.../Pipeline/raw_Data/PHI_final.xlsx', sheet_name = 'ID', header = 0)
    file = open('.../Pipeline/raw_Data/receptor.txt', 'w')
    phage = []
    for i in range(len(xl)):
        phage.append(xl.iloc[i][3])
    phage = list(set(phage))
    
    for i in range(len(phage)):
        if ',' in phage[i]:
            result = NCBI_seq(phage[i].split(',')[0])
        elif 'genome' in phage[i]:
            result = NCBI_seq(phage[i].split('genome')[0])
        elif 'complete genome' in phage[i]:
            result = NCBI_seq(phage[i].split('complete genome')[0])
        elif 'chromosome' in phage[i]:
            result = NCBI_seq(phage[i].split('chromosome')[0])
        else:
            result = NCBI_seq(phage[i])
        for j in result:
            file.writelines('%s\t%s\n' % (phage[i], j))
    file.close()