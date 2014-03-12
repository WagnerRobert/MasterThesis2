import os

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
#get dict with locations as keys and lists of according proteines as values
def getLoc2Prot(paths):
    path = paths["kmer_dir"]

    Loc2Prot = {}

    f = open( os.path.join(path, "result.txt"), 'r')
    for line in f:
        if not line.startswith("#"):
            tmp = line.rstrip().split('\t')
            if tmp[1] not in Loc2Prot:
                Loc2Prot[tmp[1]] = []
            Loc2Prot[tmp[1]].append(tmp[0])
    f.close()
    return Loc2Prot
loc2Prot = getLoc2Prot(constants)
for location in loc2Prot:
    print location
    for protein in loc2Prot[location]:
        print "\t" + protein
# for each protein in localization
    # read appropriate kmers
    # match kmers on protein and on profile proteines
    # calculate coverage for each amino acid
    # read prosite file
    # calculate which amino acids are inside prosite regions, and which are outside
    # average the coverage over all inside and outside regions
# plot each protein in localization as a point in a plot, using avg inside for one axis and avg outside for the other
