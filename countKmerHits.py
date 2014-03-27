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

__author__ = 'delur'

locKmerDict = masterthesis.reader.read_picklefile("cleanlocKmerDict", constants)

locKmerDictList = {}
for location in locKmerDict:
    if location not in locKmerDictList:
        locKmerDictList[location] = []
    for kmer in locKmerDict[location]:
        for value in locKmerDict[location][kmer]:
            locKmerDictList[location].append( (kmer, value) )
    locKmerDictList[location] = sorted(locKmerDictList[location], key=operator.itemgetter(1), reverse=True)
    locKmerDictList[location] = [i[0] for i in locKmerDictList[location]]
print "reading all the annotated sequences from swissprot"
f = open("/mnt/project/locbloc-ha/studs/robert/euka_small/eukaryota.1682.fa", 'r')
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

Overwrite = True
if os.path.exists(os.path.join(constants["working_dir"], "pickles/locCountDict.pkl")) and Overwrite == False:
    locCountDict = masterthesis.reader.read_picklefile("locCountDict", constants)
else:
    locCountDict = {}
    #print locKmerDict
    for location in locKmerDict:
        locCountDict[location] = {}

        print location
        i = 0
        lenLoc = len(locKmerDict[location].keys())
        for loc in locSeqDict:
            if loc == locTree2Uniprot[location]:
                print "\t" + loc
                for kmer in locKmerDict[location].keys():
                    if kmer not in locCountDict[location]:
                        locCountDict[location][kmer] = []
                    for protein in locSeqDict[loc]:
                        # print locSeqDict[loc][protein]
                        if kmer in locSeqDict[loc][protein]:
                            locCountDict[location][kmer].append(protein)
                            #print "found kmer " + kmer + " \tin " + loc + " protein " + protein


    masterthesis.writer.write_picklefile(locCountDict, "locCountDict", constants)


for location in locCountDict:
    i = 0
    print location + "\t" + str(len(locSeqDict[locTree2Uniprot[location]])) + "\t" + str(len(locKmerDictList[location]))
    for k in sorted(locCountDict[location], key=lambda k: len(locCountDict[location][k]), reverse=True):
        i += 1

        print k + "\t" + str(len(locCountDict[location][k])) + "\t" + str(len(Set(locCountDict[location][k]))) + "\t" + '%.2f' % (locKmerDictList[location].index(k)*100/float(len(locKmerDictList[location])))
        if i >= 45:
            break



