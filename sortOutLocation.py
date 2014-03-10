#!/usr/bin/env python
import math
import os
import sys
from masterthesis import reader
import masterthesis
from masterthesis import *

__author__ = 'delur'

constants = {}
constants["working_dir"] = "/mnt/project/locbloc-ha/studs/robert/euka/5"
constants["kmer_dir"] = os.path.join(constants["working_dir"], "kmers")
constants["uniprot"] = os.path.join(constants["working_dir"], "uniprot")
constants["fasta"] = os.path.join(constants["working_dir"], "fasta")
constants["blast_dir"] = os.path.join(constants["working_dir"], "blast")
constants["blast_tool"] = "/usr/bin/blastpgp"
constants["big_80"] = "/var/tmp/rost_db/data/big/big_80"
constants["mfasta"] = os.path.join(constants["working_dir"], "mfasta")
constants["uniqueprot"] = "uniqueprot"
constants["msa"] = os.path.join(constants["working_dir"], "msa")
constants["clustalo"] = "/home/delur/Desktop/master/test/clustalo"
constants["num_cores"] = "1"
constants["pdf"] = os.path.join(constants["working_dir"], "pdf")
constants["needle"] = "needle"
constants["needle_dir"] = os.path.join(constants["working_dir"], "needle")
constants["qsub"] = ['qsub', '-o', '/dev/null', '-e', '/dev/null', '-b', 'y']
quant = 0.1


result = reader.read_resultfile(constants)


def kmer_file(kmer_file_path):
    kmer_list = []
    f = open(kmer_file_path,'r')
    for line in f:
        tmp = line.rstrip().split(':')
        if len(tmp) > 1:
            kmer_list.append( (tmp[0], float(tmp[1])) )
        else:
            factor = math.sqrt(float(tmp[0]))
    f.close()

    for i in range(len(kmer_list)):
        kmer_list[i] = (kmer_list[i][0] , kmer_list[i][1] / factor)
    kmer_list = sorted(kmer_list, key=lambda  x:x[1])

    return  kmer_list


def kmer_dir(kmer_svm_path, paths):
    results = masterthesis.reader.read_resultfile(paths)
    svm_kmer_dict = {}
    for root, dirs, files in os.walk(kmer_svm_path):
        for kmer in sorted(files):
            proteinname = kmer.split('.')[0]
            if kmer.endswith(".kmerweights.txt"):
                svm_kmer_dict[proteinname] = kmer_file(os.path.join(kmer_svm_path, kmer))

    return svm_kmer_dict

svm_dict = kmer_dir("/mnt/project/locbloc-ha/studs/robert/euka/5/kmers/kmerweights/SVM_14", constants)

def sortOutLocation(location, svm_dict):
    for protein in svm_dict:
        if result[protein][0] == location:
            print protein +".kmerweights.txt"

sortOutLocation("cytopla", svm_dict)