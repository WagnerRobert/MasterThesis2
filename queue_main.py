#!/usr/bin/env python
import os
from masterthesis import *
from masterthesis.writer import create_plot

__author__ = 'wagnerr'

constants = {}
constants["working_dir"] = "/mnt/home/wagnerr/master/Bact"
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
#queue_blast()
#get_fasta_files()
#queue_uniqueprot()
#pairwise()

doPlots()

