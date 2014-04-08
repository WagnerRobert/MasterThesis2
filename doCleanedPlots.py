import os
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
print "reading the kmers and filtering to the quantile"
if os.path.exists(os.path.join(constants["working_dir"], "pickles/locKmerList2.pkl")):
    locKmerList = masterthesis.reader.read_picklefile("locKmerList2", constants)
else:
    locKmerList = masterthesis2.kmers.readKmers("SVM_14", 0.1, loc2prot, constants)
    masterthesis.writer.write_picklefile(locKmerList, "locKmerList2", constants)
locSeqDict = masterthesis2.locSeq.getlocSeqDict("/mnt/project/locbloc-ha/studs/robert/euka_small/eukaryota.1682.fa")
prosite = masterthesis2.prosite.readProsite()

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
overwrite = False
result = masterthesis.reader.result_file.read_resultfile(constants)
tree = masterthesis.reader.tree_file.read_treefile(constants)
i = 0
for location in locKmerList:
    print location
    for protein in sorted(locKmerList[location]):
        if protein not in locSeqDict[locTree2Uniprot[location]]:
            print "protein not sound in uniprotfile " + protein
            continue

        print "\t" + protein + "\t" + str(i)
        clean_name = protein.split('#')[0]
        if os.path.exists(os.path.join(constants["needle_dir"], clean_name + ".needle")):
            pass
        else:
            print "no needle file for " + protein
            continue
        foundUniprot, entry = masterthesis.writer.getUniprot.get_uniprot(clean_name, constants, overwrite)
        sequence = masterthesis.writer.getFasta.get_fasta(clean_name, entry, constants, overwrite)

        pairwise_alignments = masterthesis.writer.build_pairwise_alignments.build_pairwise_alignments(clean_name, constants, overwrite)
        kmerlist = locKmerList[location][protein].keys()
        pro_matches = masterthesis2.matchKmers2.match_kmers_pairwise(clean_name, sequence, pairwise_alignments, kmerlist)

        prositeFeatures = []
        if protein in prosite:
            prositeFeatures = prosite[protein]
        else:
            prositeFeatures = []

        masterthesis2.plot.createPlotWithProsite((clean_name,sequence), pro_matches, entry, len(pairwise_alignments), result, constants, prositeFeatures)
        i += 1
