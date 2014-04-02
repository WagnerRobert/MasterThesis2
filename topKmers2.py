import os
import operator
from sets import Set
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
loc2prot = masterthesis2.loc2prot.getCorrectPredictedLoc2Prot(constants)

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

    # compare kmers across localizations list everything that appears on multiple localizations extra
print "building the remove list, this can take some time"
removeList = []
i = 0
for location1 in locKmerDict:
    i += 1
    print "location1: " + str(i) + "\t" + str(len(locKmerDict))
    for location2 in locKmerDict:
        if location1 != location2:
            loc1Keys = locKmerDict[location1].keys()
            loc2Keys = locKmerDict[location2].keys()
            loc2KeySet = set(loc2Keys)
            localRemovelist = [val for val in loc1Keys if val in loc2KeySet]
            #for kmer in loc1Keys:
            #    if kmer in loc2Keys:
            #        removeList.append(kmer)
            for kmer in localRemovelist:
                removeList.append(kmer)
for location in locKmerDict:
    locKeyList = locKmerDict[location].keys()
    for kmer in Set(removeList):
        if kmer in locKeyList:
            #print "removing " + kmer + " from " + location + " :"
            #print locKmerDict[location][kmer]
            del locKmerDict[location][kmer]

    # print the reduced list
print "outputting the cleaned up top30list"

topKmer_dict = {}
for location in locKmerDict:
    topKmer_dict[location] = {}
    print location
    locKmerList = []
    for kmer in locKmerDict[location]:
        for value in locKmerDict[location][kmer]:
            locKmerList.append( (kmer, value) )
    locKmerList = sorted(locKmerList, key=operator.itemgetter(1), reverse=True)
    for kmer, value in locKmerList:
            if kmer not in topKmer_dict[location]:
                if len(topKmer_dict[location]) < 45:
                    topKmer_dict[location][kmer] = []
                else:
                    break
            topKmer_dict[location][kmer].append(value)

for location in topKmer_dict:
    print location
    for kmer in sorted(topKmer_dict[location], key=lambda x: topKmer_dict[location][x][0] , reverse=True):
        print "\t" + kmer + "\t" + str(len(topKmer_dict[location][kmer])) + "\t" + '%.2f' % (topKmer_dict[location][kmer][0])

    print ""
    print "same output again for sequence logo"
    for kmer in topKmer_dict[location]:
        print ">" + kmer
        print kmer


locKmerDict = None


