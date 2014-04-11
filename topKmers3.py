import os
import operator
import itertools
import sys
import timeit
import time
import masterthesis.reader.pickle_file
import masterthesis2.kmers
import masterthesis2.zscore
import masterthesis2.loc2prot
import masterthesis2.locSeq

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


# for each protein in localization
print "\nreading result file to get a location to protein mapping"
loc2prot = masterthesis2.loc2prot.getCorrectPredictedLoc2Prot(constants)

print "\nreading the kmers and filtering to the quantile"
if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList2.pkl")):
    locKmerList = masterthesis.reader.pickle_file.read_picklefile("locKmerList2", constants)
else:
    locKmerList = masterthesis2.kmers.readKmers("SVM_14", 0.1, loc2prot, constants)
    masterthesis.writer.pickle_file(locKmerList, "locKmerList2", constants)

locSeqDict = masterthesis2.locSeq.getlocSeqDict("/mnt/project/locbloc-ha/studs/robert/euka_small/eukaryota.1682.fa")

if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList2Clean.pkl")):
    locKmerList = masterthesis.reader.read_picklefile("locKmerList2Clean", constants)
else:
    locKmerList = masterthesis2.cleanKmers.cleanKmers(locKmerList, locSeqDict)
    masterthesis.writer.write_picklefile(locKmerList, "locKmerList2Clean", constants)

locZscore = masterthesis2.zscore.calc_zscoreDict(locKmerList)
locKmerList = None
print "\ncounting all the kmers for each location"
locKmerDict = {}
for location in locZscore:
    locKmerDict[location] = {}
    for protein in locZscore[location]:
        for kmer in locZscore[location][protein]:
            if kmer not in locKmerDict[location]:
                locKmerDict[location][kmer] = []
            locKmerDict[location][kmer].append(locZscore[location][protein][kmer])


print "\ngetting the highest scoring Kmers"
topKmer_dict = {}
for location in locKmerDict:
    topKmer_dict[location] = {}
    locKmerList = []
    for kmer in locKmerDict[location]:
        for value in locKmerDict[location][kmer]:
            locKmerList.append( (kmer, value) )
    locKmerList = sorted(locKmerList, key=operator.itemgetter(1), reverse=True)
    for kmer, value in locKmerList:
            if kmer not in topKmer_dict[location]:
                if len(topKmer_dict[location]) < 15000:
                    topKmer_dict[location][kmer] = []
                else:
                    break
            topKmer_dict[location][kmer].append(value)

print"\n Top x Kmers for:"
for location in topKmer_dict:
    print location
    for kmer in sorted(topKmer_dict[location], key=lambda x: topKmer_dict[location][x][0] , reverse=True):
        print "\t" + kmer + "\t" + str(len(topKmer_dict[location][kmer])) + "\t" + '%.2f' % (topKmer_dict[location][kmer][0])

    print ""
    print "same output again for sequence logo"
    for kmer in topKmer_dict[location]:
        print ">" + kmer
        print kmer

f = open("NLS_clear_experimental.txt", 'r')
motifs = []
for line in f:
    motifs.append(line.rstrip())
f.close()
f = open("NLS_clear_experimental_partial.txt", 'r')
partial_motifs = []
for line in f:
    partial_motifs.append(line.rstrip())

print "\nChecking top kmers against NLSdb motifs:"
print str(len(topKmer_dict["nucleus"])) + " kmers in nucleus"
count = 0
for kmer in topKmer_dict["nucleus"]:
    #print kmer
    for motif in motifs:
        if len(kmer) <= len(motif):
            if kmer in motif:
                count += 1
                print str(count) + "\t found " + kmer + " Kmer in " + motif
        else:
            if motif in kmer:
                count += 1
                print str(count) + "\tmotif " + motif + " is substring of Kmer " + kmer
    for motif in partial_motifs:
        if len(kmer) <= len(motif):
            if kmer in motif:
                count += 1
                print str(count) + "\t\t found " + kmer + " Kmer in " + motif
        else:
            if motif in kmer:
                count += 1
                print str(count) + "\t\tmotif " + motif + " is substring of Kmer " + kmer

sys.exit()
locKmerDict = None


locTree2Uniprot = {}
locTree2Uniprot["cytopla"] = "cytoplasm"
locTree2Uniprot["nucleus"] = "nucleus"
locTree2Uniprot["cellmemb"] = "cellmembrane"
locTree2Uniprot["memmitoc"] = "memmitochondria"
locTree2Uniprot["peroxis"] = "peroxisome"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["er"] = "er"
locTree2Uniprot["secrete"] = "secreted"
locTree2Uniprot["chloropl"] = "chloroplast"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["mitochon"] = "mitochondria"


print "\nCounting all found kmers for eatch protein and building all kombinations of kmers found"
counts = {}
for location in locZscore:
    counts[location] = {}
    for protein in locZscore[location]:
        foundKmers = []
        for kmer in locZscore[location][protein]:
            if kmer in locSeqDict[locTree2Uniprot[location]][protein]:
                foundKmers.append(kmer)
        #build all combinations from the found kmers
        for L in range(1, len(foundKmers)+1):
            for subset in itertools.combinations(sorted(foundKmers), L):
                if subset in counts[location]:
                    counts[location][subset] += 1
                else:
                    counts[location][subset] = 1

print "\nremoving kombinations that are subsets of larger kombinations if they dont have higher counts than the larger combinations"
for location in counts:
    print location + "\t" + str(len(counts[location]))
    removelist = []
    i = 0
    keepGoing = True
    workedList = []



    for index, subsetA in enumerate(sorted(counts[location], key=lambda x: counts[location][x])):
        if index % 1000 == 0:
            print str(index) + "/" + str(len(counts[location])) #+ "\t" + str(subsetA)
        sys.stdout.flush()
        for subsetB in counts[location]:
            if subsetA != subsetB:
                if counts[location][subsetA] <= counts[location][subsetB]:
                    totalySubSet = True
                    for kmer in subsetA:
                        if kmer not in subsetB:
                            totalySubSet = False
                            break
                    if totalySubSet:
                        #print str(subsetA) + "\t" + str(counts[location][subsetA]) + " is subset of " + str(subsetB)+ "\t" + str(counts[location][subsetB])
                        del counts[location][subsetA]
                        break

    # while keepGoing:
    #     keepGoing = False
    #     oldLength = len(counts[location])
    #     print oldLength
    #     print "Length workedlist: " + str(len(workedList))
    #     for subsetA in sorted(counts[location], key= lambda x : len(x), reverse=True):
    #         if subsetA in workedList:
    #             continue
    #

    #             print str(subsetA) + "\t" + str(subsetB)
    #             sys.stdout.flush()
    #             if counts[location][subsetB] < counts[location][subsetA]:
    #                 if set(subsetB).issubset(set(subsetA)) :
    #                     print "\t" + str(subsetB) + "\t" + str(subsetA)
    #
    #                     del counts[location][subsetB]
    #         newLength = len(counts[location])
    #         print newLength
    #         if oldLength != newLength:
    #             workedList.append(subsetA)
    #             keepGoing = True
    #             break
    #     newLength = len(counts[location])
    #     print newLength
    #     if oldLength != newLength:
    #         keepGoing = True



f = open( "top4", 'w')
for location in counts:
    print location
    f.write(location + "\n")
    for subset in sorted(sorted(counts[location], key= lambda x : counts[location][x], reverse=True), key= lambda  x : len(x)):
        f.write(str(subset) + "\t" + str(counts[location][subset]) + "\n")
        print str(subset) + "\t" + str(counts[location][subset])