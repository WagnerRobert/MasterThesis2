import operator
import masterthesis2.zscore

__author__ = 'delur'

def cleanKmers(locKmerList, locSeqDict):

    locTree2Uniprot = {}
    locTree2Uniprot["cytopla"] = "cytoplasm"
    locTree2Uniprot["nucleus"] = "nucleus"

    locKmerDict = {}
    for location in locKmerList:
        locKmerDict[location] = {}
        for protein in locKmerList[location]:
            for kmer, value in locKmerList[location][protein]:
                if kmer not in locKmerDict[location]:
                    locKmerDict[location][kmer] = []
                locKmerDict[location][kmer].append(value)

    print "\nSearching and removing kmers that match in multiple locations"
    removelist = []
    for location in locKmerDict:
        print location
        i = 0
        lenLoc = len(locKmerDict[location].keys())
        for loc in locSeqDict:
            if loc != locTree2Uniprot[location]:
                print "\t" + loc
                for kmer in locKmerDict[location].keys():
                    for protein in locSeqDict[loc]:
                        # print locSeqDict[loc][protein]
                        if kmer in locSeqDict[loc][protein]:
                            i += 1
                            #print str(i) + "/" + str(lenLoc)  + ": \tfound " + location + " kmer " + kmer + " \tin " + loc + " protein " + protein + " \t- removing it from all location lists"
                            removelist.append(kmer)
                            for tmplocation in locKmerDict:
                                if kmer in locKmerDict[tmplocation]:
                                    del locKmerDict[tmplocation][kmer]
                            break
                print "\t\t" + str(i) + "/" + str(lenLoc)
    locKmerDict = None




    locProtKmerDict = {}
    for location in locKmerList:
         locProtKmerDict[location] = {}
         for protein in locKmerList[location]:
             locProtKmerDict[location][protein] = {}
             for kmer, value in locKmerList[location][protein]:
                 if kmer not in locProtKmerDict[location]:
                     locProtKmerDict[location][protein][kmer] = []
                 locProtKmerDict[location][protein][kmer].append(value)

    locKmerList = None

    print "\nFound the following kmers in multiple locations:"
    count = 0
    for kmer in removelist:
        count += 1
        print str(count) + "\t" + kmer
        for location in locProtKmerDict:
            for protein in locProtKmerDict[location]:
                 if kmer in locProtKmerDict[location][protein]:
                     del locProtKmerDict[location][protein][kmer]



    print "Searching for kmers that LocTree deems important to multiple regions and building the remove list, this can take some time"
    removeList = []
    i = 0
    for location1 in locProtKmerDict:
        i += 1
        print "location1: " + str(i) + "\t" + str(len(locProtKmerDict))
        for location2 in locProtKmerDict:
            if location1 != location2:
                loc1KmerList = []
                print "LOCATION 1 KMERLIST:"
                for protein in locProtKmerDict[location1]:
                    for kmer in locProtKmerDict[location1][protein].keys():
                        #print kmer
                        loc1KmerList.append(kmer)
                loc2KmerSet = set()
                print "LOCATION 2 KMERLIST:"
                for protein in locProtKmerDict[location2]:
                    for kmer in locProtKmerDict[location2][protein].keys():
                        #print kmer
                        loc2KmerSet.add(kmer)
                localRemovelist = [val for val in loc1KmerList if val in loc2KmerSet]
                #for kmer in loc1Keys:
                #    if kmer in loc2Keys:
                #        removeList.append(kmer)
                print "REMOVE KMERLIST:"
                for kmer in localRemovelist:
                #    print kmer
                    removeList.append(kmer)

    print "\nFound the following kmers in multiple locations:"
    count = 0
    for kmer in removeList:
        count += 1
        print str(count) + "\t" + kmer
        for location in locProtKmerDict:
            for protein in locProtKmerDict[location]:
                 if kmer in locProtKmerDict[location][protein]:
                     del locProtKmerDict[location][protein][kmer]


    return locProtKmerDict