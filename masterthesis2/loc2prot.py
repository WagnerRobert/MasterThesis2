import os

__author__ = 'wagnerr'

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


def getCorrectPredictedLoc2Prot(paths):
    loc2prot = getLoc2Prot(paths)

    print "reading all the annotated sequences from swissprot"
    f = open(os.path.join(paths["working_dir"], "eukaryota.1682.fa"), 'r')
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
    locTree2Uniprot["cellmemb"] = "cellmembrane"
    locTree2Uniprot["memmitoc"] = "memmitochondria"
    locTree2Uniprot["peroxis"] = "peroxisome"
    locTree2Uniprot["mitochon"] = "mitochondria"
    locTree2Uniprot["er"] = "er"
    locTree2Uniprot["secrete"] = "secreted"
    locTree2Uniprot["chloropl"] = "chloroplast"
    locTree2Uniprot["mitochon"] = "mitochondria"
    locTree2Uniprot["mitochon"] = "mitochondria"

    #cross check loc2Prot against locSeqDict, to find Proteins that are wrongfully predicted to be in a location
    i = 0
    removeList = []
    for localisation in loc2prot:
        for protein in sorted(loc2prot[localisation]):
            if protein not in locSeqDict[locTree2Uniprot[localisation]]:
                i +=1
                print i
                removeList.append(protein)
                foundLoc = None
                for loc in locSeqDict:
                    if protein in locSeqDict[loc]:
                        foundLoc = loc
                        print protein + "\tpredicted localisation: " + localisation + "\tfound localization: " + loc
                        break
                if foundLoc is None:
                    print protein + " predicted localisation: " + localisation + " not found!"
    for localisation in loc2prot:
        for protein in removeList:
            loc2prot[localisation] = filter(lambda a: a != protein, loc2prot[localisation])
    return  loc2prot