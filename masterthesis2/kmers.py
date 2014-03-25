import os
import math
from masterthesis import reader

__author__ = 'wagnerr'
def readKmers(SVM, quant,loc2Prot, constants):
    tree = reader.read_treefile(constants)

    locKmerList = {}
    for location in loc2Prot:
        #determin if kmers supporting this Localization have positive or negative values
        locationPosNeg = "?"
        if location in tree[SVM][0] or location in tree[SVM][1]:
            if location not in locKmerList:
                locKmerList[location] = {}
        if location in tree[SVM][0]:
            locationPosNeg = "neg"
        elif location in tree[SVM][1]:
            locationPosNeg = "pos"
        else:
            continue

        print location

        #read appropriate kmers
        for protein in sorted(loc2Prot[location]):
            locKmerList[location][protein] = []
            print "\t" + protein
            svmdir = os.path.join(constants["kmer_dir"],"kmerweights/"+ SVM)
            proteinPath = os.path.join(svmdir , protein + ".kmerweights.txt")

            kmer_list = []
            # this reads the appropriate kmers, and turnes the negative kmervalues positive so further comparison is easier
            f = open(proteinPath, 'r')
            for line in f:
                tmp = line.rstrip().split(':')
                if len(tmp) > 1:
                    if locationPosNeg == "neg" and float(tmp[1]) < 0:
                        kmer_list.append( (tmp[0], -float(tmp[1])) )
                    elif locationPosNeg == "pos" and float(tmp[1]) > 0:
                        kmer_list.append( (tmp[0], float(tmp[1])) )
                else:
                    factor = math.sqrt(float(tmp[0]))
            f.close()

            for i in range(len(kmer_list)):
                kmer_list[i] = (kmer_list[i][0] , kmer_list[i][1] / factor)
            kmer_list = sorted(kmer_list, key=lambda  x:x[1])

            quanted_list = []

            #determin quant
            total = 0.0
            for kmer, value in kmer_list :
                total += value
            posquant = quant * total

            #apply quant
            quanted_sum = 0.0
            while quanted_sum < posquant:
                kmer = kmer_list.pop()
                quanted_sum += kmer[1]
                quanted_list.append(kmer)

            #add quantedList to locKmerList
            locKmerList[location][protein] = quanted_list

    return locKmerList