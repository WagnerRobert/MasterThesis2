import os
import masterthesis2.loc2prot
import masterthesis2.kmers
import masterthesis2.locSeq
import masterthesis2.cleanKmers

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

loc2prot = masterthesis2.loc2prot.getCorrectPredictedLoc2Prot(constants)
locKmerList = masterthesis2.kmers.readKmers("SVM_14", 0.1, loc2prot, constants)
locSeqDict = masterthesis2.locSeq.getlocSeqDict("/mnt/project/locbloc-ha/studs/robert/euka_small/eukaryota.1682.fa")

locTree2Uniprot = {}
locTree2Uniprot["cytopla"] = "cytoplasm"
locTree2Uniprot["nucleus"] = "nucleus"
locTree2Uniprot["cellmemb"] = "cellmembrane"
locTree2Uniprot["memmitoc"] = "memmitochondria"
locTree2Uniprot["peroxis"] = "peroxisome"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["er"] = "er"
locTree2Uniprot["secrete"] = "secreted"
locTree2Uniprot["chloropl"] = "chloroplast"
locTree2Uniprot["mitochon"] = "mitochondria"
locTree2Uniprot["mitochon"] = "mitochondria"

locKmerList = masterthesis2.cleanKmers.cleanKmers(locKmerList, locSeqDict)
