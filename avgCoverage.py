import os
import math
from masterthesis import reader

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


#todo
#get dict with locations as keys and lists of according proteines as values
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

loc2Prot = getLoc2Prot(constants)

# for each protein in localization for a specifies SVM level
def readKmers(SVM, quant, constants):
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
        for protein in loc2Prot[location]:
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
                print kmer + "\t" + str(value)
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

locKmerList = readKmers("SVM_14", 0.1, constants)


    # match kmers on protein and on profile proteines
    # calculate coverage for each amino acid
    # read prosite file
    # calculate which amino acids are inside prosite regions, and which are outside
    # average the coverage over all inside and outside regions
# plot each protein in localization as a point in a plot, using avg inside for one axis and avg outside for the other
