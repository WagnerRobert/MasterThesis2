import os
import re
import operator
import sys
import masterthesis2.loc2prot
import masterthesis2.kmers
import masterthesis2.locSeq
import masterthesis2.cleanKmers
import masterthesis.reader.pickle_file
import masterthesis.writer.pickle_file
import masterthesis.writer.getUniprot
import masterthesis.writer.getFasta
import masterthesis.writer.build_pairwise_alignments
import masterthesis2.matchKmers2
import masterthesis.reader.result_file
import masterthesis.reader.tree_file
import masterthesis2.plot
import masterthesis2.prosite
import masterthesis2.evalNoPlot
import numpy as np

__author__ = 'delur'

def getConstants():
    constants = {}
    #constants["working_dir"] = "/home/delur/Desktop/euka_small"
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
    return constants

def getCorrectPredictedKmers(constants):
    locKmerList = masterthesis.reader.read_picklefile("locKmerList2Clean", constants)
    return  locKmerList

def getSequences(constants):
    locSeqDict = masterthesis2.locSeq.getlocSeqDict(os.path.join(constants["working_dir"], "eukaryota.1682.fa"))
    return locSeqDict

def getTopKmers(locKmerList, location):
    locZscore = masterthesis2.zscore.calc_zscoreDict(locKmerList)
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
        locKmerListA = []
        for kmer in locKmerDict[location]:
            for value in locKmerDict[location][kmer]:
                locKmerListA.append( (kmer, value) )
        locKmerListA = sorted(locKmerListA, key=operator.itemgetter(1), reverse=True)
        for kmer, value in locKmerListA:
                if kmer not in topKmer_dict[location]:
                    if len(topKmer_dict[location]) < 50 or False:
                        topKmer_dict[location][kmer] = []
                    else:
                        break
                topKmer_dict[location][kmer].append(value)
    top_kmerlist = sorted(topKmer_dict[location], key=lambda x: topKmer_dict[location], reverse=True)


    return top_kmerlist

def matchKmers(protein, constants, location, locKmerList):
    overwrite = False
    clean_name = protein.split('#')[0]
    if os.path.exists(os.path.join(constants["needle_dir"], clean_name + ".needle")):
        pass
    else:
        print "no needle file for " + protein
        return None
    foundUniprot, entry = masterthesis.writer.getUniprot.get_uniprot(clean_name, constants, overwrite)
    sequence = masterthesis.writer.getFasta.get_fasta(clean_name, entry, constants, overwrite)

    pairwise_alignments = masterthesis.writer.build_pairwise_alignments.build_pairwise_alignments(clean_name, constants, overwrite)
    kmerlist = locKmerList[location][protein].keys()
    print "Protein Kmerlist has: " + str(len(kmerlist)) + " elements"

    matches = masterthesis2.matchKmers2.match_kmers_pairwise(clean_name, sequence, pairwise_alignments, kmerlist)

    return matches

def getNLSdbPatternMatches(constants, locSeqDict, location, protein):
    patternMatches = []

    f = open(os.path.join(constants["working_dir"], "NLSdb.txt"), 'r')
    motifs = []
    for line in f:
        motifs.append(line.replace("x", ".").rstrip())
    f.close()

    f = open(os.path.join(constants["working_dir"], "NLSdb_potential.txt"), 'r')
    potential_motifs = []
    for line in f:
        potential_motifs.append(line.replace("x", ".").rstrip())
    f.close()

    for motif in motifs:
        regex = re.compile(motif)
        i = 0
        while True:
            match = regex.search(locSeqDict[location][protein], i)
            if  match:
                start, end = match.span()
                end = end -1

                patternMatches.append( (start,end) )
                i = match.start()+1
            else:
                break

    for motif in potential_motifs:
        regex = re.compile(motif)
        i = 0
        while True:
            match = regex.search(locSeqDict[location][protein], i)
            if  match:
                start, end = match.span()
                end = end -1

                patternMatches.append( (start,end) )
                i = match.start()+1
            else:
                break

    return patternMatches

def getValidNESMatches(constants, locSeqDict, location, protein):
    patternMatches = []

    f = open(os.path.join(constants["working_dir"], "validNES_filtered.txt"), 'r')
    motifs = []
    for line in f:
        motifs.append(line.replace("x", ".").rstrip())
    f.close()


    for motif in motifs:
        regex = re.compile(motif)
        i = 0
        while True:
            match = regex.search(locSeqDict[location][protein], i)
            if  match:
                start, end = match.span()
                end = end -1

                patternMatches.append( (start,end) )
                i = match.start()+1
            else:
                break

    return patternMatches

def getPrositeMatches(protein):
    prosite = masterthesis2.prosite.readProsite(constants)
    prositeFeatures = []
    if protein in prosite:
        prositeFeatures = prosite[protein]
    else:
        prositeFeatures = []
    return  prositeFeatures

def evaulutePerAminoAcid(patternMatches, zscore_count, cutoff):
    def findHit(patternMatches, i):
        for start, end in patternMatches:
            if i in range(start, end):
                return True
        return False

    TruePositives = 0
    FalsePositives = 0
    TrueNegatives = 0
    FalseNegatives = 0

    for i in range(len(zscore_count)):
        if zscore_count[i] > cutoff:
            if findHit(patternMatches, i):
                TruePositives += 1
            else:
                FalsePositives += 1
        else:
            if findHit(patternMatches, i):
                FalseNegatives += 1
            else:
                TrueNegatives += 1

    print "\t\t|Positives\t|Negatives"
    print "True\t|" + str(TruePositives) + "\t\t\t|" + str(TrueNegatives)
    print "False\t|" + str(FalsePositives) + "\t\t\t|" + str(FalseNegatives)

def evaluatePerSegment(protein, patternMatches, pro_matches, cutoff):
    overwrite = False
    clean_name = protein.split('#')[0]
    foundUniprot, entry = masterthesis.writer.getUniprot.get_uniprot(clean_name, constants, overwrite)
    sequence = masterthesis.writer.getFasta.get_fasta(clean_name, entry, constants, overwrite)
    pairwise_alignments = masterthesis.writer.build_pairwise_alignments.build_pairwise_alignments(clean_name, constants, overwrite)
    answer = masterthesis2.evalNoPlot.eval_segment((clean_name,sequence), pro_matches, len(pairwise_alignments), patternMatches, [], cutoff)

    return answer

def evaluatePerAminoacid(protein, patternMatches, pro_matches, cutoff):
    overwrite = False
    clean_name = protein.split('#')[0]
    foundUniprot, entry = masterthesis.writer.getUniprot.get_uniprot(clean_name, constants, overwrite)
    sequence = masterthesis.writer.getFasta.get_fasta(clean_name, entry, constants, overwrite)
    pairwise_alignments = masterthesis.writer.build_pairwise_alignments.build_pairwise_alignments(clean_name, constants, overwrite)
    answer = masterthesis2.evalNoPlot.eval_aminoacid((clean_name,sequence), pro_matches, len(pairwise_alignments), patternMatches, [], cutoff)

    return answer

def getLocKmerDict(locZscore, location):
    locKmerDict = {}
    locKmerDict[location] = {}
    for protein in locZscore[location]:
        for kmer in locZscore[location][protein]:
            if kmer not in locKmerDict[location]:
                locKmerDict[location][kmer] = []
            locKmerDict[location][kmer].append(locZscore[location][protein][kmer])
    return locKmerDict

location = "nucleus"
#get constants required for work
constants = getConstants()

#get correct predicted kmers
locKmerList = getCorrectPredictedKmers(constants)

locZscore = masterthesis2.zscore.calc_zscoreDict(locKmerList)

locKmerDict = getLocKmerDict(locZscore, location)

#get sequences of correct predicted kmers
locSeqDict = getSequences(constants)

#get top kmers for localization
top_kmerlist = getTopKmers(locKmerList, location)

def theEvaluation(type, cutoff, proMatchesDict, top_kmerlist):
    i = 0
    precisionList = []
    recallList = []
    for protein in sorted(locKmerList[location]):
        #get the kmer matches
        pro_matches = proMatchesDict[location][protein]

        if pro_matches is None:
            continue
        #todo this is bugged, needs fixing, top_kmerlist has the wrong format
        def dotopMatches():
            top_matches = {}
            top_matches[location] = {}
            top_matches[location][protein] = {}
            for kmer in top_kmerlist:
                top_matches[location][protein][kmer] = True
            top_matches = matchKmers(protein, constants, location, top_matches)
            return top_matches
        #get the pattern matches
        #prositeMatches = getPrositeMatches()

        if type == "NLSdbAA":
            patterns = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerAminoacid(protein, patterns, pro_matches, cutoff)
        if type == "NLSdbSeg":
            patterns = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerSegment(protein, patterns, pro_matches, cutoff)

        if type == "PrositeAA":
            patterns = getPrositeMatches(protein)
            answer = evaluatePerAminoacid(protein, patterns, pro_matches, cutoff)
        if type == "PrositeSeg":
            patterns = getPrositeMatches(protein)
            answer = evaluatePerSegment(protein, patterns, pro_matches, cutoff)

        if type == "ValidNESAA":
            patterns = getValidNESMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerAminoacid(protein, patterns, pro_matches, cutoff)
        if type == "VaildNESSeg":
            patterns = getValidNESMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerSegment(protein, patterns, pro_matches, cutoff)

        if type == "NLSdbAATop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerAminoacid(protein, patterns, top_matches, cutoff)
        if type == "NLSdbSegTop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerSegment(protein, patterns, top_matches, cutoff)

        if type == "PrositeAATop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getPrositeMatches(protein)
            answer = evaluatePerAminoacid(protein, patterns, top_matches, cutoff)
        if type == "PrositeSegTop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getPrositeMatches(protein)
            answer = evaluatePerSegment(protein, patterns, top_matches, cutoff)

        if type == "ValidNESAATop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getValidNESMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerAminoacid(protein, patterns, top_matches, cutoff)
        if type == "VaildNESSegTop":
            top_matches = None
            top_matches = dotopMatches()
            patterns = getValidNESMatches(constants, locSeqDict, location, protein)
            answer = evaluatePerSegment(protein, patterns, top_matches, cutoff)

        if type == "mergedAA":
            nldPatt = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            proPatt = getPrositeMatches(protein)
            valPatt = getValidNESMatches(constants, locSeqDict, location, protein)
            if nldPatt is None:
                nldPatt = []
            if proPatt is None:
                proPatt = []
            if valPatt is None:
                valPatt = []
            patterns = nldPatt + proPatt + valPatt
            answer = evaluatePerAminoacid(protein, patterns, pro_matches, cutoff)
        if type == "mergedSeg":
            nldPatt = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            proPatt = getPrositeMatches(protein)
            valPatt = getValidNESMatches(constants, locSeqDict, location, protein)
            if nldPatt is None:
                nldPatt = []
            if proPatt is None:
                proPatt = []
            if valPatt is None:
                valPatt = []
            patterns = nldPatt + proPatt + valPatt
            answer = evaluatePerSegment(protein, patterns, pro_matches, cutoff)

        if type == "mergedAATop":
            top_matches = None
            top_matches = dotopMatches()
            nldPatt = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            proPatt = getPrositeMatches(protein)
            valPatt = getValidNESMatches(constants, locSeqDict, location, protein)
            if nldPatt is None:
                nldPatt = []
            if proPatt is None:
                proPatt = []
            if valPatt is None:
                valPatt = []
            patterns = nldPatt + proPatt + valPatt
            answer = evaluatePerAminoacid(protein, patterns, top_matches, cutoff)
        if type == "mergedSegTop":
            top_matches = None
            top_matches = dotopMatches()
            nldPatt = getNLSdbPatternMatches(constants, locSeqDict, location, protein)
            proPatt = getPrositeMatches(protein)
            valPatt = getValidNESMatches(constants, locSeqDict, location, protein)
            if nldPatt is None:
                nldPatt = []
            if proPatt is None:
                proPatt = []
            if valPatt is None:
                valPatt = []
            patterns = nldPatt + proPatt + valPatt
            answer = evaluatePerSegment(protein, patterns, top_matches, cutoff)

        if answer is None:
            continue
        else:
            precision, recall = answer
        i += 1
        print "\t" + str(i) + "\t" + protein
        precisionList.append(precision)
        recallList.append(recall)
    return precisionList, recallList


def doFMeasure(index, precision, recall):
    return ((1+index)*precision*recall) / (index*precision + recall)


def getMatches(constants, location, locKmerList):
    proMatchesDict = {}
    proMatchesDict[location] = {}
    for protein in sorted(locKmerList[location]):
        pro_matches = matchKmers(protein, constants, location, locKmerList)
        proMatchesDict[location][protein] = pro_matches
    return proMatchesDict


proMatchesDict = getMatches(constants, location, locKmerList)
topMatchesDict = top_kmerlist
f = open(os.path.join(constants["working_dir"], "patternEval.txt"), 'w' )

f.write("Stats for full quant:\n")
print("Stats for full quant:\n")
# f.write("NLSdbAA\tprecision\trecall\tF1\tF0.5\n")
# print("NLSdbAA\tprecision\trecall\tF1\tF0.5\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("NLSdbAA",i, proMatchesDict, topMatchesDict)
#     f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
#     print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
# f.write("NLSdbSeg\tprecision\trecall\tF1\tF0.5\n")
# print("NLSdbSeg\tprecision\trecall\tF1\tF0.5\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("NLSdbSeg",i , proMatchesDict, topMatchesDict)
#     f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
#     print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
#
f.write("ValidNESAA\tprecision\trecall\tF1\tF0.5\n")
print("ValidNESAA\tprecision\trecall\tF1\tF0.5\n")
for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
    precisionList, recallList = theEvaluation("ValidNESAA",i, proMatchesDict, topMatchesDict)
    f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
    print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
f.write("VaildNESSeg\tprecision\trecall\tF1\tF0.5\n")
print("VaildNESSeg\tprecision\trecall\tF1\tF0.5\n")
for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
    precisionList, recallList = theEvaluation("VaildNESSeg",i, proMatchesDict, topMatchesDict)
    f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
    print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")

# f.write("PrositeAA\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("PrositeAA",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
# f.write("PrositeSeg\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("PrositeSeg",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")

# f.write("mergedAA\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("mergedAA",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
# f.write("mergedSeg\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("mergedSeg",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
#
f.write("Stats for top quant:\n")
print("Stats for top quant:\n")
# f.write("NLSdbAA\tprecision\trecall\tF1\tF0.5\n")
# print("NLSdbAA\tprecision\trecall\tF1\tF0.5\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("NLSdbAATop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
#     print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
# f.write("NLSdbSeg\tprecision\trecall\tF1\tF0.5\n")
# print("NLSdbSeg\tprecision\trecall\tF1\tF0.5\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("NLSdbSegTop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
#     print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")

f.write("ValidNESAA\tprecision\trecall\tF1\tF0.5\n")
print("ValidNESAA\tprecision\trecall\tF1\tF0.5\n")
for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
    precisionList, recallList = theEvaluation("ValidNESAATop",i, proMatchesDict, topMatchesDict)
    f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
    print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
f.write("VaildNESSeg\tprecision\trecall\tF1\tF0.5\n")
print("VaildNESSeg\tprecision\trecall\tF1\tF0.5\n")
for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
    precisionList, recallList = theEvaluation("VaildNESSegTop",i, proMatchesDict, topMatchesDict)
    f.write("\t" +"%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")
    print("\t" + "%.2f" % np.average(precisionList) +"\t" + "%.2f" % np.average(recallList) +"\t" + "%.2f" % doFMeasure(1,np.average(precisionList), np.average(recallList) )+"\t" + "%.2f" % doFMeasure(0.5,np.average(precisionList), np.average(recallList) )+"\n")

# f.write("PrositeAA\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("PrositeAATop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
# f.write("PrositeSeg\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("PrositeSegTop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
#
# f.write("mergedAA\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("mergedAATop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")
# f.write("mergedSeg\tprecision\trecall\n")
# for i in [1.0, 0.5, 0.0, -0.5, -1.0]:
#     precisionList, recallList = theEvaluation("mergedSegTop",i, proMatchesDict, topMatchesDict)
#     f.write("\t" + str(np.average(precisionList)) +"\t" + str(np.average(recallList)) +"\t" + str(doFMeasure(1,np.average(precisionList), np.average(recallList) ))+"\t" + str(doFMeasure(0.5,np.average(precisionList), np.average(recallList) ))+"\n")

f.close()