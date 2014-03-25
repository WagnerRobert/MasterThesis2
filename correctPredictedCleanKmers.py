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
locTree2Uniprot["cellmemb"] = "plasma_membrane"
locTree2Uniprot["memmitoc"] = "mitochondria_membrane"
locTree2Uniprot["peroxis"] = "peroxisome"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["er"] = "er"
locTree2Uniprot["secrete"] = "secreted"
locTree2Uniprot["chloropl"] = "chloroplast"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["mitochon"] = "mitochondria"

#cross check loc2Prot against locSeqDict, to find Proteines that are wrongfully predicted to be in a location
i = 0
for localisation in loc2prot:
    for protein in loc2prot[localisation]:
        if protein not in locSeqDict[locTree2Uniprot[localisation]]:
            i +=1
            print i
            del loc2prot[localisation][protein]
            foundLoc = None
            for loc in locSeqDict:
                if protein in locSeqDict[loc]:
                    foundLoc = loc
                    print protein + "\tpredicted localisation: " + localisation + "\tfound localization: " + loc
                    break
            if foundLoc == None:
                print protein + " predicted localisation: " + localisation + " not found!"


#     # read appropriate kmers
# print "reading the kmers and filtering to the quantile"
# if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList2.pkl")):
#     locKmerList = masterthesis.reader.read_picklefile("locKmerList2", constants)
# else:
#     locKmerList = masterthesis2.kmers.readKmers("SVM_14", 0.1, loc2prot, constants)
#     masterthesis.writer.write_picklefile(locKmerList, "locKmerList2", constants)


