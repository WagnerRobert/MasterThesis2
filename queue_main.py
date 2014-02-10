#!/usr/bin/env python
import os
import sys
from masterthesis import *
from masterthesis.writer import create_plot

__author__ = 'wagnerr'

constants = {}
constants["working_dir"] = "/mnt/project/locbloc-ha/studs/robert/Bact/"
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
quant = 0.1

# sets up the directory in which all calculations will be done
setUp(constants)
result = reader.read_resultfile(constants)
tree = reader.read_treefile(constants)
#reader.read_kmerfiles(constants, quant) #reads and prepares kmerweights
#files, saves complete result (dict[svm][protein]) kmers.pkl in pickles dir

kmerlist = reader.read_picklefile("kmers", constants)
#checkOrder(kmerlist, result) # to checkOrder you need to outcomment the doQuant call in read_kmerfile

overwrite = False
queue = True
def queue_blast():
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):
            print "\t" + protein
            clean_name = protein.split('#')[0]
            foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
            sequence = get_fasta(clean_name, entry, constants, overwrite)
            blastProtein(clean_name, constants, overwrite, queue)
            #processProtein(protein, kmerlist[svm][protein], result[protein], tree[svm], constants)


def get_fasta_files():
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):
            print "\t" + protein
            clean_name = protein.split('#')[0]
            foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
            sequence = get_fasta(clean_name, entry, constants, overwrite)
            profileProteines = blastProtein(clean_name, constants, overwrite, queue)
            print "\t\tcollecting fasta sequences to build mfasta file from"
            wget_fasta(profileProteines, constants, overwrite)

def queue_uniqueprot():
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):
            print "\t" + protein
            clean_name = protein.split('#')[0]
            foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
            sequence = get_fasta(clean_name, entry, constants, overwrite)
            profileProteines = blastProtein(clean_name, constants, overwrite, queue)
            print "\t\tbuild mfasta file and clean it with uniqueprot"
            build_mfasta(clean_name, sequence, profileProteines, constants, overwrite, queue)

def pairwise():
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):
            print "\t" + protein
            clean_name = protein.split('#')[0]
            pairwise_alignments = build_pairwise_alignments(clean_name, constants, overwrite)

def doPlots():
    i = 0
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):

            print "\t" + protein + "\t" + str(i)
            clean_name = protein.split('#')[0]
            if os.path.exists(os.path.join(constants["needle_dir"], clean_name + ".needle")):
                pass
            else:
                print "no needle file for " + protein
                continue
            foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
            sequence = get_fasta(clean_name, entry, constants, overwrite)

            pairwise_alignments = build_pairwise_alignments(clean_name, constants, overwrite)
            kmerlists = kmerlist[svm][protein]
            pro_kmerlist = []
            con_kmerlist = []
            if result[protein][0] in tree[svm][0]:
                pro_kmerlist = kmerlists[1]
                con_kmerlist = kmerlists[0]
            elif result[protein][0] in tree[svm][1]:
                pro_kmerlist = kmerlists[0]
                con_kmerlist = kmerlists[1]
            pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
            con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
            create_plot((clean_name,sequence), pro_matches, con_matches, entry, len(pairwise_alignments), result, constants, svm)
            i += 1

def doQuantCountPlots():
    kmersPerQuant = {}
    kmersPerQuantLocationBased = {}

    if os.path.exists(os.path.join(os.path.join(constants["working_dir"] , "pickles"), "kmercount.pkl")):
        kmersPerQuant = read_picklefile("kmercount", constants)
    else:
        print "fail"
        for svm in sorted(kmerlist):
            kmersPerQuant[svm] = {}
            print svm
            for protein in sorted(kmerlist[svm]):
                kmersPerQuant[svm][protein] = {}
                print "\t" + protein
                clean_name = protein.split('#')[0]
                path = os.path.join(os.path.join(os.path.join(constants["kmer_dir"], "kmerweights"), svm), protein + ".kmerweights.txt")
                for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
                    kmersPerQuant[svm][protein][i] = None
                    kmerlisting = reader.kmer_file(path, i)
                    kmersPerQuant[svm][protein][i] =  len(kmerlisting[0]), len(kmerlisting[1])
                    #print str(i) + "\t" + str(len(kmerlisting[0])) + "\t" + str(len(kmerlisting[1]))
        write_picklefile(kmersPerQuant, "kmercount", constants)

    for svm in kmersPerQuant:
        kmersPerQuantLocationBased[svm] = {}
        for location in tree[svm]:
            tmp = location.split(',')
            for split_location in tmp:
                kmersPerQuantLocationBased[svm][split_location] = {}
        for protein in kmersPerQuant[svm]:
            kmersPerQuantLocationBased[svm][result[protein][0]][protein] = {}
            if result[protein][0] in tree[svm][0]:
                for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
                    kmersPerQuantLocationBased[svm][result[protein][0]][protein][i] = kmersPerQuant[svm][protein][i][1], kmersPerQuant[svm][protein][i][0]
            elif result[protein][0] in tree[svm][1]:
                for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
                    kmersPerQuantLocationBased[svm][result[protein][0]][protein][i] = kmersPerQuant[svm][protein][i][0], kmersPerQuant[svm][protein][i][1]

    svm_location_dict = {}
    for svm in kmersPerQuantLocationBased:
        svm_location_dict[svm] = {}
        for location in kmersPerQuantLocationBased[svm]:
            svm_location_dict[svm][location] = None
            pos_count = []
            neg_count = []
            for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
                    pos_count.append(0)
                    neg_count.append(0)
            for protein in kmersPerQuantLocationBased[svm][location]:
                j = 0
                for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
                    pos_count[j] += kmersPerQuantLocationBased[svm][location][protein][i][0]
                    neg_count[j] += kmersPerQuantLocationBased[svm][location][protein][i][1]
                    j += 1
            for i in range(len(pos_count)):
                pos_count[i] = pos_count[i] / float(len(kmersPerQuantLocationBased[svm][location]) )
                neg_count[i] = neg_count[i] / float(len(kmersPerQuantLocationBased[svm][location]) )

            svm_location_dict[svm][location] = (pos_count, neg_count)

    print svm_location_dict.keys()
    for svm in svm_location_dict:
        print svm
        for location in svm_location_dict[svm]:
            print "\t" + str(location)
            print svm_location_dict[svm][location][0]
            print svm_location_dict[svm][location][1]

            x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            import matplotlib.pyplot as plt
            plt.clf()
            plt.cla()
            fig, (ax0, ax1) = plt.subplots(nrows=2)
            ax0.plot(x, svm_location_dict[svm][location][0], 'o', color=(0.0, 0.0, 1.0))
            ax0.plot(x, svm_location_dict[svm][location][1], 'o', color=(1.0, 0.0, 0.0))
            ax0.set_xlim(0.0, 1.0)
            ax0.yaxis.grid(True)
            ax0.xaxis.grid(True)

            ax1.plot(x[0:6], svm_location_dict[svm][location][0][0:6],'o', color=(0.0, 0.0, 1.0))
            ax1.plot(x[0:6], svm_location_dict[svm][location][1][0:6], 'o', color=(1.0, 0.0, 0.0))
            ax1.set_xlim(0.0, 0.5)
            ax1.yaxis.grid(True)
            ax1.xaxis.grid(True)

            plt.savefig(os.path.join(constants["pdf"], svm+location+".pdf"))

def calcHitWidth():
    for svm in sorted(kmerlist):
        print svm
        for protein in sorted(kmerlist[svm]):
            print "\t" + protein
            clean_name = protein.split('#')[0]
            if os.path.exists(os.path.join(constants["needle_dir"], clean_name + ".needle")):
                pass
            else:
                break
            foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
            sequence = get_fasta(clean_name, entry, constants, overwrite)

            pairwise_alignments = build_pairwise_alignments(clean_name, constants, overwrite)
            kmerlists = kmerlist[svm][protein]
            pro_kmerlist = []
            con_kmerlist = []
            if result[protein][0] in tree[svm][0]:
                pro_kmerlist = kmerlists[1]
                con_kmerlist = kmerlists[0]
            elif result[protein][0] in tree[svm][1]:
                pro_kmerlist = kmerlists[0]
                con_kmerlist = kmerlists[1]
            #pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
            pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
            con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
            #print pro_matches

            pos_count = None
            pos_count_noGaps = [0] * len(sequence)
            neg_count = None
            neg_count_noGaps = [0] * len(sequence)

            print pro_matches

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
                            print x
                            print len(sequence)
                            pos_count_noGaps[x] += pos_count[i]
                            x += 1

            for prot in con_matches:
                for match_seq , start, end in con_matches[prot]:
                    neg_count = [0] * len(match_seq)
                    for j in range(start, end):
                        if j < len(match_seq):
                            neg_count[j] += 1
                    x = 0
                    for i in range(len(match_seq)):
                        if x == len(sequence):
                            break
                        if match_seq[i] == '-':
                            pass
                        else:
                            neg_count_noGaps[x] += neg_count[i]
                            x += 1

            numProfileProteins =  float(len(pairwise_alignments))
            print numProfileProteins
            for i in range(len(sequence)):
                pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins
                neg_count_noGaps[i] = neg_count_noGaps[i] * 100 / numProfileProteins

            print pos_count_noGaps[i]
            sys.exit()

def countNumProfProteines():
    f = open(os.path.join(constants["working_dir"],"numProfProts.txt" ) , 'r')
    numProfProt = {}
    for line in f:
        tmp = line.split(':')
        numProfProt[tmp[0].split('.')[0]] = int(tmp[1])
    f.close()

    numProfProtPerLoc = {}
    for protein in numProfProt:
        #print protein
        if result[protein][0] not in numProfProtPerLoc:
            numProfProtPerLoc[result[protein][0]] = []
        numProfProtPerLoc[result[protein][0]].append( numProfProt[protein])

    for location in numProfProtPerLoc:
        print location
        print sorted(numProfProtPerLoc[location])
        import matplotlib.pyplot as plt
        plt.clf()
        plt.cla()
        # maximum = max(numProfProtPerLoc[location])
        # minimum = min(numProfProtPerLoc[location])
        # L = maximum - minimum
        # plt.hist(numProfProtPerLoc[location], L/5)

        plt.hist(numProfProtPerLoc[location], [0, 10, 100, 250, 500, 750, 1000, 1250, 1500])
        plt.xticks([0, 10, 100, 250, 500, 750, 1000, 1250, 1500])
        plt.title(location + " Total: " + str(len(numProfProtPerLoc[location])))
        plt.ylabel('Number of Query Proteins')
        plt.xlabel('Number of ProfileProteins')
        plt.savefig(os.path.join(constants["pdf"], location+".pdf"))

#queue_blast()
#get_fasta_files()
#queue_uniqueprot()
#pairwise()

#doPlots()
#doQuantCountPlots()
calcHitWidth()
#countNumProfProteines()
