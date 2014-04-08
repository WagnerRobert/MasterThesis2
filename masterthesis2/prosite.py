import re


def readProsite():
    prosite = {}
    prositePath = "./euka_small/ProSite.txt"
    f = open(prositePath, 'r')

    for line in f:
        if line.startswith(">"):
            protein = line.rstrip().split('>')[1]
        if protein not in prosite:
            prosite[protein] = []
        match = re.search( r"(\d+)\s\-\s(\d+)", line)
        if match:
            prosite[protein].append( (int(match.group(1)), int(match.group(2))) )
            #print (int(match.group(1)), int(match.group(2)))
    return prosite

prosite = readProsite()