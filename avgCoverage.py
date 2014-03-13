import os
import math
import re
import sys
from masterthesis import write_picklefile, reader, build_pairwise_alignments, get_fasta, get_uniprot
from masterthesis import match_kmers_pairwise

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

if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList.pkl")):
    locKmerList = reader.read_picklefile("locKmerList", constants)
else:
    locKmerList = readKmers("SVM_14", 0.1, constants)
    write_picklefile(locKmerList, "locKmerList", constants)

# read prosite file
def readProsite():
    prosite = {}
    prositePath = os.path.join(constants["working_dir"], "ProSite.txt")
    f = open(prositePath, 'r')

    for line in f:
        if line.startswith(">"):
            protein = line.rstrip().split('>')[1]
        if protein not in prosite:
            prosite[protein] = []
        match = re.search( r"(\d+)\s\-\s(\d+)", line)
        if match:
            prosite[protein].append( (int(match.group(1)), int(match.group(2))) )
            #print (int(match.group(1)), int(match.group(2)))
    return prosite

prosite = readProsite()

locAvgCoverage = {}
for location in locKmerList:
    locAvgCoverage[location]= []
    print location
    for protein in locKmerList[location]:
        if protein not in prosite:
            continue
        print "\t" + protein + "\t" + str(len(locKmerList[location][protein]))
        # match kmers on protein and on profile proteines
        foundUniprot, entry = get_uniprot(protein, constants, False)
        sequence = get_fasta(protein, entry, constants, False)
        pairwise_alignments = build_pairwise_alignments(protein, constants, False)
        if len(pairwise_alignments) == 0:
            continue
        pro_matches = match_kmers_pairwise(protein, sequence, pairwise_alignments, locKmerList[location][protein])

        # calculate coverage for each amino acid position
        pos_count = None
        pos_count_noGaps = [0] * len(sequence)

        for prot in pro_matches:
            for match_seq , start, end in pro_matches[prot]:
                pos_count = [0] * len(match_seq)
                for j in range(start, end):
                    if j < len(match_seq):
                        pos_count[j] += 1
                x = 0
                for i in range(len(match_seq)):
                    if x == len(sequence):
                        break
                    if match_seq[i] == '-':
                        pass
                    else:
                        pos_count_noGaps[x] += pos_count[i]
                        x += 1

        numProfileProteins = len(pairwise_alignments)
        for i in range(len(sequence)):
            pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins



        # calculate which amino acids are inside prosite regions, and which are outside
        prosite_regions = [0] * len(pos_count_noGaps)
        for start, end in prosite[protein]:
            for i in range(start,end):
                prosite_regions[i] = 1

        no_prosite_coverage_list = []
        prosite_coverage_list = []
        for i in range(len(pos_count_noGaps)):
            if prosite_regions[i] == 0:
                no_prosite_coverage_list.append(pos_count_noGaps[i])
            elif prosite_regions[i] == 1:
                prosite_coverage_list.append(pos_count_noGaps[i])
        # average the coverage over all inside and outside regions
        print "\t\t" +'%.2f' % (math.fsum(prosite_coverage_list)/float(len(prosite_coverage_list))) + "\t" + '%.2f' % (math.fsum(no_prosite_coverage_list)/float(len(no_prosite_coverage_list))) + "\t" + str( len(prosite_coverage_list)) + "\t" + str( len(no_prosite_coverage_list))
        locAvgCoverage[location].append( (protein, math.fsum(prosite_coverage_list)/float(len(prosite_coverage_list)), math.fsum(no_prosite_coverage_list)/float(len(no_prosite_coverage_list))))

# plot each protein in localization as a point in a plot, using avg inside for one axis and avg outside for the other
import numpy as np
import matplotlib.pyplot as plt
for location in locAvgCoverage:
    locInPrositeRegion = []
    locOutPrositeRegion = []
    for entry in locAvgCoverage[location]:
        locInPrositeRegion.append(entry[1])
        locOutPrositeRegion.append(entry[2])
    plt.plot(locInPrositeRegion, locOutPrositeRegion, marker='.')
    x1,x2,n,m,b = min(locInPrositeRegion),max(locInPrositeRegion),1000,1.,0.
    x = np.r_[x1:x2:n*1j] #http://docs.scipy.org/doc/numpy/reference/generated/numpy.r_.html
    plt.plot(x,m*x + b); plt.grid();
    plt.show()
