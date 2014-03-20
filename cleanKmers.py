import os
import operator
from sets import Set
import sys
import masterthesis2.loc2prot
import masterthesis2.kmers
import masterthesis2.zscore
import masterthesis.writer.pickle_file

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
# for each protein in localization
print "reading result file to get a location to protein mapping"
loc2prot = masterthesis2.loc2prot.getLoc2Prot(constants)

    # read appropriate kmers
print "reading the kmers and filtering to the quantile"
if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList2.pkl")):
    locKmerList = masterthesis.reader.read_picklefile("locKmerList2", constants)
else:
    locKmerList = masterthesis2.kmers.readKmers("SVM_14", 0.1, loc2prot, constants)
    masterthesis.writer.pickle_file(locKmerList, "locKmerList2", constants)

    # calculate the zscore
print "calculating the zscores on the kmer vaules"
locZscore = masterthesis2.zscore.calc_zscore(locKmerList)
locKmerList = None

    # get the kmers of each location, count how often which one occurs and the highest(all?) zscore(s)
print "counting all the kmers for each location"
locKmerDict = {}
for location in locZscore:
    locKmerDict[location] = {}
    for protein in locZscore[location]:
        for kmer, value in locZscore[location][protein]:
            if kmer not in locKmerDict[location]:
                locKmerDict[location][kmer] = []
            locKmerDict[location][kmer].append(value)
locZscore = None


print "reading all the annotated sequences from swissprot"
f = open("/mnt/project/locbloc-ha/sp042011/SP13_11/eukaryotes.SP13_11.expSL.50.Before05_11.new.fa", 'r')
name = ""
localisation = ""
locSeqDict= {}
i = 0
for line in f:
    if line.startswith(">"):
        tmp = line.split(' ')
        name = tmp[0][1:]
        localisation = tmp[1].rstrip()
        if localisation not in locSeqDict:
            locSeqDict[localisation] = {}
        if name not in locSeqDict[localisation]:
            locSeqDict[localisation][name] = ""
        else:
            print name + " is already in " + localisation
    else:
        locSeqDict[localisation][name] = line.rstrip()
f.close()

locTree2Uniprot = {}
locTree2Uniprot["cytopla"] = "cytoplasm"
locTree2Uniprot["nucleus"] = "nucleus"

removelist = []
for location in locKmerDict:
    print location
    i = 0
    lenLoc = len(locKmerDict[location].keys())
    for loc in locSeqDict:
        if loc != locTree2Uniprot[location]:
            print "\t" + loc
            for kmer in locKmerDict[location].keys():
                for protein in locSeqDict[loc]:
                    # print locSeqDict[loc][protein]
                    if kmer in locSeqDict[loc][protein]:
                        i += 1
                        #print str(i) + "/" + str(lenLoc)  + ": \tfound " + location + " kmer " + kmer + " \tin " + loc + " protein " + protein + " \t- removing it from all location lists"

                        for tmplocation in locKmerDict:
                            if kmer in locKmerDict[tmplocation]:
                                del locKmerDict[tmplocation][kmer]
                        break
            print "\t\t" + str(i) + "/" + str(lenLoc)



for loc in locSeqDict:
    print loc + "\t" + str(len(locSeqDict[loc]))