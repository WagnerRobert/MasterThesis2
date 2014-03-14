import os

__author__ = 'wagnerr'

def getLoc2Prot(paths):
    path = paths["kmer_dir"]

    Loc2Prot = {}

    f = open( os.path.join(path, "result.txt"), 'r')
    for line in f:
        if not line.startswith("#"):
            tmp = line.rstrip().split('\t')
            if tmp[1] not in Loc2Prot:
                Loc2Prot[tmp[1]] = []
            Loc2Prot[tmp[1]].append(tmp[0])
    f.close()
    return Loc2Prot