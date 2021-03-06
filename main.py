import os
import sys

__author__ = 'delur'

from masterthesis import *

constants = {}
constants["working_dir"] = "/home/delur/Desktop/master/archaea"
constants["kmer_dir"] = os.path.join(constants["working_dir"], "kmers")
constants["uniprot"] = os.path.join(constants["working_dir"], "uniprot")
constants["fasta"] = os.path.join(constants["working_dir"], "fasta")
constants["blast_dir"] = os.path.join(constants["working_dir"], "blast")
constants["blast_tool"] = "/usr/bin/blastpgp"
constants["big_80"] = "/home/delur/Desktop/master/big/big_80"
constants["mfasta"] = os.path.join(constants["working_dir"], "mfasta")
constants["uniqueprot"] = "uniqueprot"
constants["msa"] = os.path.join(constants["working_dir"], "msa")
constants["clustalo"] = "/home/delur/Desktop/master/test/clustalo"
constants["num_cores"] = "1"
constants["pdf"] = os.path.join(constants["working_dir"], "pdf")
constants["needle"] = "needle"
constants["needle_dir"] = os.path.join(constants["working_dir"], "needle")
quant = 1.0

# sets up the directory in which all calculations will be done
setUp(constants)
result = reader.read_resultfile(constants)
tree = reader.read_treefile(constants)
reader.read_kmerfiles(constants, quant) #reads and prepares kmerweights
#files, saves complete result (dict[svm][protein]) kmers.pkl in pickles dir

kmerlist = reader.read_picklefile("kmers", constants)
#checkOrder(kmerlist, result) # to checkOrder you need to outcomment the doQuant call in read_kmerfiles

for svm in sorted(kmerlist):
     print svm
     for protein in sorted(kmerlist[svm]):
         print "\t" + protein
         processProtein(protein, kmerlist[svm][protein], result[protein], tree[svm], constants)
