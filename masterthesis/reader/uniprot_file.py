import os

__author__ = 'delur'

def read_uniprotfile(name, paths):
    entry = ""
    filelocation = os.path.join(paths["uniprot"],  name + ".txt")

    f = open(filelocation, 'r')
    for line in f:
        entry += line
    f.close()

    entry = entry.split('\n')

    return entry