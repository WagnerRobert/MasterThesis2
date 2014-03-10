#!/usr/bin/env python
import copy
import math
import os
import sys
import operator
from masterthesis import reader
import masterthesis
from masterthesis import *

__author__ = 'delur'

constants = {}
constants["working_dir"] = "/mnt/project/locbloc-ha/studs/robert/euka_small/"
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


#result = reader.read_resultfile(constants)

def doQuant(kmerlist, quant):
    positive_list = []
    negative_list = []

    for kmer, value in kmerlist:
        if value < 0 :
            negative_list.append( (kmer, value) )
        else:
            positive_list.append( (kmer, value) )
    positive_list.reverse()

    total = 0.0
    for kmer, value in positive_list :
        total += value
    posquant = quant * total

    while total > posquant and len(positive_list) > 0:
        total -= positive_list.pop()[1]

    total = 0.0
    for kmer, value in negative_list:
        total += value
    negquant = quant * total

    while total < negquant and len(negative_list) > 0:
        total -= negative_list.pop()[1]

    return (positive_list, negative_list)


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
    kmer_list = sorted(kmer_list, key=lambda  x:x[1], reverse=True)

    kmer_list = doQuant(kmer_list,quant)
    return  kmer_list


def kmer_dir(kmer_svm_path, paths):
    svm_kmer_dict = {}
    for root, dirs, files in os.walk(kmer_svm_path):
        for kmer in sorted(files):
            proteinname = kmer.split('.')[0]
            if kmer.endswith(".kmerweights.txt"):
                svm_kmer_dict[proteinname] = kmer_file(os.path.join(kmer_svm_path, kmer))

    return svm_kmer_dict

cytopla_dict = kmer_dir("/mnt/project/locbloc-ha/studs/robert/euka_small/cytopla", constants)
write_picklefile(cytopla_dict,"cytopla_dict", constants)
nucleus_dict = kmer_dir("/mnt/project/locbloc-ha/studs/robert/euka_small/nucleus", constants)
write_picklefile(nucleus_dict,"nucleus_dict", constants)

cytopla_dict = read_picklefile("cytopla_dict", constants)
nucleus_dict = read_picklefile("nucleus_dict", constants)

import numpy as np
from  scipy import stats

def zscore(svm_dict, posOrNeg):
    zscores_location = np.array([])
    kmers_location = []
    for protein in svm_dict:
        kmers = []
        values = []
        for kmer,value in svm_dict[protein][posOrNeg]:
            kmers.append(kmer)
            values.append(value)
        zscores_protein = stats.zscore(values)
        kmers_location += kmers
        zscores_location = np.append(zscores_location, zscores_protein)


    zscoreLocList = zscores_location.tolist()
    #print zscoreLocList
    dict_with_zscores = {}
    while len(dict_with_zscores) < 30:
        max_index, max_value = max(enumerate(zscoreLocList), key=operator.itemgetter(1))
        if kmers_location[max_index] not in dict_with_zscores:
            dict_with_zscores[kmers_location[max_index]] = (max_value, 1)
        else:
            dict_with_zscores[kmers_location[max_index]] = (dict_with_zscores[kmers_location[max_index]][0], dict_with_zscores[kmers_location[max_index]][1] +1)
        kmers_location.pop(max_index)
        zscoreLocList.pop(max_index)
    for key, value_tuple in sorted(dict_with_zscores.iteritems(), key=operator.itemgetter(1), reverse=True):
        print key + "\t" + '%.2f' % value_tuple[0] + "\t" + str(value_tuple[1])

    return dict_with_zscores




print "starting zscore stuff now"
pos = 0
neg = 1
cytopla_dict = zscore(cytopla_dict, neg)
nucleus_dict = zscore(nucleus_dict, pos)

removelist = []
for index1, element1 in enumerate(sorted(cytopla_dict.iteritems(), key=operator.itemgetter(1), reverse=True)):
    for index2, element2 in enumerate(sorted(nucleus_dict.iteritems(), key=operator.itemgetter(1), reverse=True)):
        if element1[0] == element2[0]:
            print element1[0] + " found in both"
            removelist.append(element1[0])

#print removelist
clean_cytopla_dict =  copy.deepcopy(cytopla_dict)
for index, element in enumerate(sorted(cytopla_dict.iteritems(), key=operator.itemgetter(1), reverse=True)):
    #print element
    if element[0] in removelist:
        #print element[0]
        clean_cytopla_dict.pop(element[0])

#print removelist
clean_nucleus_dict =  copy.deepcopy(nucleus_dict)
for index, element in enumerate(sorted(nucleus_dict.iteritems(), key=operator.itemgetter(1), reverse=True)):
    #print element
    if element[0] in removelist:
        #print element[0]
        clean_nucleus_dict.pop(element[0])

print "!!Cleared List!!"

print "cytoplas"
for key, value_tuple in sorted(clean_cytopla_dict.iteritems(), key=operator.itemgetter(1), reverse=True):
        print "\t" + key + "\t" + '%.2f' % value_tuple[0] + "\t" + str(value_tuple[1])
print "nucleus"
for key, value_tuple in sorted(clean_nucleus_dict.iteritems(), key=operator.itemgetter(1), reverse=True):
        print "\t" + key + "\t" + '%.2f' % value_tuple[0] + "\t" + str(value_tuple[1])
