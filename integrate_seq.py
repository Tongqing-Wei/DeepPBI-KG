import os 

def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name

path = os.getcwd()
phage = walkFile(path + os.sep + 'phage_raw_data')
f = open(path + os.sep + 'all_phage_seq.fasta', 'w')
pp = []
for i in range(len(phage)):
    ff = open(path + os.sep + 'phage_raw_data' + os.sep + phage[i], 'r')
    line = ff.readlines()
    ff.close()
    pp.extend(line)
f.writelines(pp)
f.close()

host = walkFile(path + os.sep + 'host_raw_data')
f = open(path + os.sep + 'all_host_seq.fasta', 'w')
hh = []
for i in range(len(host)):
    ff = open(path + os.sep + 'host_raw_data' + os.sep + host[i], 'r')
    line = ff.readlines()
    ff.close()
    hh.extend(line)
f.writelines(hh)
f.close()