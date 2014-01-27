import os

__author__ = 'delur'

def read_treefile(paths):
    path = paths["working_dir"]
    infile = open(os.path.join(path, "tree.txt"), 'r')

    tree = {}
    for line in infile:
        tmp = line.rstrip().split('\t')
        tree[tmp[0]] = tmp[1].split(':')
    infile.close()

    return tree

