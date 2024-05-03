import os
from Bio import SeqIO
import argparse


def format_fasta(ana, seq, num):
    """
    Format the text in fasta format
    :param ana: Annotation information
    :param seq: sequence
    :param num: The number of characters when the sequence is wrapped
    :return: fasta format text
    """
    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq += "\n"
    return ana + format_seq + "\n"


def get_complete(gb_file):
    """
    Extract the cds sequence and its complete sequence from the genbank file
    :param gb_file: genbank file path
    :param f_cds: Whether to obtain only one CDS sequence
    :return: CDS sequence in fasta format, complete sequence in fasta format
    """
    # Extract the complete sequence and format it in fasta
    gb_seq = SeqIO.parse(open(gb_file), "genbank")
    complete_total = []
    complete_seq_total = []
    for seq_record in gb_seq:
        complete_seq = str(seq_record.seq)
        complete_seq_total.append(complete_seq)
        complete_ana = ">" + seq_record.id + ":" + " " + seq_record.description + "\n"
        complete_fasta = format_fasta(complete_ana, complete_seq, 70)
        complete_total.append(complete_fasta)
        
    return complete_total, complete_seq_total    


def get_cds_dna(gb_file, complete_seq_total):

    gb_seq = SeqIO.parse(open(gb_file), "genbank")
    
    # The CDS dna sequence was extracted and formatted as fasta
    cds_num = 1
    cds_fasta = ""
    for i, seq_record in enumerate(gb_seq):
        for ele in seq_record.features:
            if ele.type == "CDS":
                cds_seq = ""
                if list(ele.qualifiers)[0] == 'gene':
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                              " [gene=" + ele.qualifiers['gene'][0] + "]" + "[inference=" + ele.qualifiers['inference'][1] + "]" + \
                              " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                              " [gbkey=CDS]\n"
                else:
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                              " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                              " [gbkey=CDS]\n"
                cds_num += 1
                for ele1 in ele.location.parts:
                    cds_seq += complete_seq_total[i][ele1.start:ele1.end]
                cds_fasta += format_fasta(cds_ana, cds_seq, 70)
    return cds_fasta


def get_cds_protein(gb_file, f_cds):
    """
    Extract the cds sequence and its complete sequence from the genbank file
    :param gb_file: genbank file path
    :param f_cds: Whether to obtain only one CDS sequence
    :return: CDS sequence in fasta format, complete sequence in fasta format
    """
    gb_seq = SeqIO.parse(open(gb_file), "genbank")

    # The CDS protein sequence was extracted and the format was fasta
    cds_num = 1
    cds_fasta = ""
    for i, seq_record in enumerate(gb_seq):
        for ele in seq_record.features:
            if ele.type == "CDS":
                cds_seq = ""
                if list(ele.qualifiers)[0] == 'gene':
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                            " [gene=" + ele.qualifiers['gene'][0] + "]" + "[inference=" + ele.qualifiers['inference'][1] + "]" + \
                            " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                            " [gbkey=CDS]\n"
                else:
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                            " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                            " [gbkey=CDS]\n"
                cds_num += 1
                cds_fasta += format_fasta(cds_ana, ele.qualifiers['translation'][0], 70)
                if (f_cds):
                    break
    return cds_fasta


def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	the path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in dirs:
            phage_name.append(os.path.join(f))
    return phage_name


if __name__ == '__main__':
    # The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help = 'input file path')
    parser.add_argument('--dna_output', type=str, help = 'dna cds output file path')
    parser.add_argument('--protein_output', type=str, help='protein cds output file path')
    opt = parser.parse_args()

    # file output path
    walkfile = walkFile(opt.input)
    #Handle bp Spaces in gbk annotation files
    for ls in walkfile:
        res_dir = opt.input + os.sep + ls
        for file in os.listdir(res_dir):
            if file[-3:] == 'gbk':
                fo = open(res_dir + os.sep + file, 'r+')
                flist = fo.readlines()
                for k, line in enumerate(flist):
                    if line[:5] == 'LOCUS' and 'length' in line :
                        line = line[:line.find('length')] + line[line.find('length'):line.find('cov')].replace(' ', '') + line[line.find('cov'):]
                        res = line.split('_')
                        num = res[res.index('length') + 1]

                        index_list = []
                        index = line.find(num)
                        while index != -1:
                            index_list.append(index)
                            index = line.find(num, index + 1)
                        
                        line = line[:index_list[0]] + line[index_list[0]:index_list[-1]].replace(" ", "") + ' ' + line[index_list[-1]:]
                        flist[k] = line

                fo = open(res_dir + os.sep + file, 'w+')
                fo.writelines(flist)
                fo.close()
    
    
    
    for ls in walkfile:
        #cds_dna_file = opt.dna_output + "cds_" + ls[7:] + ".fasta"
        #cds_protein_file = opt.protein_output + "cds_pro_" + ls[7:] + ".fasta"
        cds_dna_file = opt.dna_output + os.sep + "cds_" + ls[7:] + ".fasta"
        cds_protein_file = opt.protein_output + os.sep + "cds_pro_" + ls[7:] + ".fasta"
        # genbank file path
        res_dir = opt.input + os.sep + ls
        cds_dna_file_obj = open(cds_dna_file, "w")
        cds_protein_file_obj = open(cds_protein_file, 'w')
        
        for file in os.listdir(res_dir):
            if file[-3:] == 'gbk':
                complete_fasta, complete_seq_fasta = get_complete(res_dir + os.sep + file)
                cds_dna_fasta = get_cds_dna(res_dir + os.sep + file, complete_seq_fasta)
                cds_dna_file_obj.write(cds_dna_fasta)   
                
                cds_protein_fasta = get_cds_protein(res_dir + os.sep + file, False)
                cds_protein_file_obj.write(cds_protein_fasta)                
                
        cds_dna_file_obj.close()
        cds_protein_file_obj.close()