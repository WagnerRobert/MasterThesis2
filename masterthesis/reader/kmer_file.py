import os
import math
from masterthesis import doQuant
import masterthesis.reader
import masterthesis.writer
__author__ = 'delur'

paths = {}

def kmer_file(kmer_file_path, quant):
    kmer_list = []
    f = open(kmer_file_path,'r')
    for line in f:
        tmp = line.rstrip().split(' ')
        if len(tmp) > 1:
            kmer_list.append( (tmp[0], float(tmp[1])) )
        else:
            factor = math.sqrt(float(tmp[0]))
    f.close()

    for i in range(len(kmer_list)):
        kmer_list[i] = (kmer_list[i][0] , kmer_list[i][1] / factor)
    kmer_list = sorted(kmer_list, key=lambda  x:x[1])
    kmer_list = doQuant(kmer_list,quant)
    return  kmer_list


def kmer_dir(kmer_svm_path, tree, paths, quant):
    results = masterthesis.reader.read_resultfile(paths)
    svm_kmer_dict = {}
    for root, dirs, files in os.walk(kmer_svm_path):
        for kmer in sorted(files):
            proteinname = kmer.split('.')[0]
            if kmer.endswith(".kmerweights.txt") and (results[proteinname][0] in tree[0] or results[proteinname][0] in tree[1]):
                print "\t"+proteinname
                svm_kmer_dict[proteinname] = kmer_file(os.path.join(kmer_svm_path, kmer), quant)

    return svm_kmer_dict


def read_kmerfiles(paths, quant):
    tree = masterthesis.reader.read_treefile(paths)
    kmerweights_dir_path = os.path.join(paths["kmer_dir"] , "kmerweights")
    kmer_dict = {}
    for root, dirs, files in os.walk(kmerweights_dir_path):
        for svm in sorted(dirs):
            print svm
            if svm.startswith("SVM"):
                kmer_dict[svm] = kmer_dir(os.path.join(kmerweights_dir_path, svm), tree[svm], paths, quant)

    masterthesis.writer.write_picklefile(kmer_dict,"kmers",  paths)
