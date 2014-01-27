import os

__author__ = 'delur'


def read_resultfile(paths):
    path = paths["kmer_dir"]

    result_info = {}

    f = open( os.path.join(path, "result.txt"), 'r')
    for line in f:
        if not line.startswith("#"):
            tmp = line.rstrip().split('\t')
            result_info[tmp[0]] = (tmp[1], tmp[2])
    f.close()
    return result_info