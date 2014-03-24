
import operator
import os
from sets import Set
import masterthesis.writer.pickle_file

__author__ = 'wagnerr'

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

locKmerDict = masterthesis.reader.read_picklefile("cleanlocKmerDict", constants)

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


protCountKmerDict = {}
for location in locKmerDict:
    loc = locTree2Uniprot[location]
    if loc not in protCountKmerDict:
        protCountKmerDict[loc] = {}
    locKmerList = []
    for kmer in locKmerDict[location]:
        for value in locKmerDict[location][kmer]:
            locKmerList.append( (kmer, value) )
    locKmerList = sorted(locKmerList, key=operator.itemgetter(1), reverse=True)
    top30_dict = {}
    for kmer, value in locKmerList:
            if kmer not in top30_dict:
                if len(top30_dict) < 45:
                    top30_dict[kmer] = []
                else:
                    break
            top30_dict[kmer].append(value)
    print "Top 45 kmers in " + location
    for kmer in top30_dict:
        print "\t" + kmer + "\t" + str(len(top30_dict[kmer]))

    for kmer in top30_dict:
        protCountKmerDict[loc][kmer] = []
        for sequence in locSeqDict[loc]:
            if kmer in locSeqDict[loc][sequence]:
                protCountKmerDict[loc][kmer].append(sequence)


print "hits in sequences"
for location in protCountKmerDict:
    print location
    for k in sorted(protCountKmerDict[location], key=lambda k: len(protCountKmerDict[location][k]), reverse=True):
        print k + "\t" + str(len(protCountKmerDict[location][k]))
