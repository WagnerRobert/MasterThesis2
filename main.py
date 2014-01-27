import os

__author__ = 'delur'

from masterthesis import *

paths = {}
paths["working_dir"] = "/home/delur/Desktop/master/bacteria"
paths["kmer_dir"] = os.path.join(paths["working_dir"], "kmers")

quant = 0.1

# sets up the directory in which all calculations will be done
setUp(paths["working_dir"])
result = reader.read_resultfile(paths)
tree = reader.read_treefile(paths)
#reader.read_kmerfiles(paths, quant) #reads and prepares kmerweights
#files, saves complete result (dict[svm][protein]) kmers.pkl in pickles dir

kmerlist = reader.read_picklefile("kmers", paths)
#checkOrder(kmerlist, result) # to checkOrder you need to outcomment the doQuant call in read_kmerfiles

for svm in sorted(kmerlist):
    print svm
    for protein in sorted(kmerlist[svm]):
        print "\t" + protein
        processProtein(protein, kmerlist[svm][protein], result[protein], tree[svm])
