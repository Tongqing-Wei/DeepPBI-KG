import sys
import gzip

USAGE = "\nusage: python %s idmapping.tb.gz blastout outputfile\n" % sys.argv[0]

if len(sys.argv) != 4:
    print(USAGE)
    sys.exit()

def parseIDmapping(filename):
    UniProt_GO = {}
    with gzip.open(filename, 'r') as f:
        for line in f:
            lsplit = line.decode().rstrip().split("\t")
            if lsplit[7]:
                UniProt_GO[lsplit[0]] = lsplit[6]
    return UniProt_GO

def parseBlastOut(filename):
    tab_res = {}
    with open(filename, 'r') as f:
        for line in f:
            lsplit = line.split()
            tab_res.setdefault(lsplit[0], set()).add(lsplit[1])
    return tab_res


UniProtKB_GO = parseIDmapping(sys.argv[1])
BlastOut = parseBlastOut(sys.argv[2])
OUT = open(sys.argv[3], 'w')


for i in BlastOut:
    temp = []
    for j in BlastOut[i]:
        if j in UniProtKB_GO:
            go = UniProtKB_GO[j].split("; ")
            temp = temp + go
        else:
            continue
    if temp:
        OUT.write(i + "\t" +  ",".join(set(temp)) + "\n")

OUT.close()
