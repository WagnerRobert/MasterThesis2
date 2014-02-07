#!/usr/bin/env python
import os
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
#checkOrder(kmerlist, result) # to checkOrder you need to outcomment the doQuant call in read_kmerfiles

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
            pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
            con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
            create_plot((clean_name,sequence), pro_matches, con_matches, entry, len(pairwise_alignments), result, constants)

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
            plt.plot(x, svm_location_dict[svm][location][0],'o', color=(0.0, 0.0, 1.0), x, svm_location_dict[svm][location][1], 'o', color=(1.0, 0.0, 0.0))
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
            pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
            con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
            print pro_matches



#queue_blast()
#get_fasta_files()
#queue_uniqueprot()
#pairwise()

#doPlots()
doQuantCountPlots()
#calcHitWidth()

