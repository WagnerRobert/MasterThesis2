import os
import datetime
from masterthesis.match_kmers import match_kmers, match_kmers_pairwise
from masterthesis.writer.blastProtein import blastProtein
from masterthesis.writer import get_uniprot
from masterthesis.writer import get_fasta
from masterthesis.writer.build_mfasta import build_mfasta
from masterthesis.writer.build_msa import build_msa
from masterthesis.writer.build_pairwise_alignments import build_pairwise_alignments
from masterthesis.writer.createPlot import create_plot

__author__ = 'delur'
import sys




def processProtein(name, kmerlists, result_info, tree, constants):
    pro_kmerlist = []
    con_kmerlist = []
    if result_info[0] in tree[0]:
        pro_kmerlist = kmerlists[1]
        con_kmerlist = kmerlists[0]
    elif result_info[0] in tree[1]:
        pro_kmerlist = kmerlists[0]
        con_kmerlist = kmerlists[1]

    clean_name = name.split('#')[0]
    overwrite = False
    queue = True
    print "\t\tgetting Uniprot entry..."
    foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
    print "\t\tfinished"
    if not foundUniprot:
        print "!" + name + "has no uniprot entry"
        return
    sequence = get_fasta(clean_name, entry, constants, overwrite)
    print "\t\tblasting Protein against Big_80..." + str(datetime.datetime.now().time())
    profileProteines = blastProtein(clean_name, constants, overwrite, queue)
    print "\t\tfinished "  + str(datetime.datetime.now().time())
    print "\t\tloading fasta Files from Uniprot, building mfasta file, then cleaning it with uniqueprot"
    build_mfasta(clean_name, sequence, profileProteines, constants, overwrite, queue)
    print "\t\tfinished " + str(datetime.datetime.now().time())
    print "\t\tcalculating pairwise sequence alignments between query and profile protein sequences with needle"
    #msa = build_msa(clean_name, constants, overwrite)
    pairwise_alignments = build_pairwise_alignments(clean_name, constants, overwrite)
    print "\t\tfinished"

    pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
    con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
#    pro_matches = match_kmers(pro_kmerlist, msa)
#    con_matches = match_kmers(con_kmerlist, msa)
    create_plot((name,sequence), pro_matches, con_matches, entry, len(pairwise_alignments), result_info, constants)