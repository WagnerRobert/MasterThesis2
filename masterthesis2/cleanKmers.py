import operator
import masterthesis2.zscore

__author__ = 'delur'

def cleanKmers(locKmerList, locSeqDict):

         # get the kmers of each location, count how often which one occurs and the highest(all?) zscore(s)
    # print "counting all the kmers for each location"
    locKmerDict = {}
    for location in locKmerList:
         locKmerDict[location] = {}
         for protein in locKmerList[location]:
             for kmer, value in locKmerList[location][protein]:
                 if kmer not in locKmerDict[location]:
                     locKmerDict[location][protein] = []
                 locKmerDict[location][protein].append(kmer)

    locKmerList = None


    locTree2Uniprot = {}
    locTree2Uniprot["cytopla"] = "cytoplasm"
    locTree2Uniprot["nucleus"] = "nucleus"

    # removelist = []
    # for location in locKmerList :
    #     for protein in locKmerList[location]:
    #         for kmer, value in locKmerList[location][protein]:
    #             if kmer not in removelist:
    #                 for loc in locSeqDict:
    #                     if kmer in removelist:
    #                         break
    #                     elif locTree2Uniprot[location] != loc :
    #                         for protein in locSeqDict[loc]:
    #                             if kmer in locSeqDict[loc][protein]:
    #                                 removelist.append(kmer)
    #                                 break


    print "building the remove list, this can take some time"
    removeList = []
    i = 0
    for location1 in locKmerDict:
        i += 1
        print "location1: " + str(i) + "\t" + str(len(locKmerDict))
        for location2 in locKmerDict:
            if location1 != location2:
                loc1KmerList = []
                print "LOCATION 1 KMERLIST:"
                for protein in locKmerDict[location1]:
                    for kmer in locKmerDict[location1][protein].keys():
                        print kmer
                        loc1KmerList.append(kmer)
                loc2KmerSet = set()
                print "LOCATION 2 KMERLIST:"
                for protein in locKmerDict[location2]:
                    for kmer in locKmerDict[location2][protein].keys():
                        print kmer
                        loc2KmerSet.add(kmer)
                localRemovelist = [val for val in locKmerDict if val in loc2KmerSet]
                #for kmer in loc1Keys:
                #    if kmer in loc2Keys:
                #        removeList.append(kmer)
                print "REMOVE KMERLIST:"
                for kmer in localRemovelist:
                    print kmer
                    removeList.append(kmer)

    print "\nFound the following kmers in multiple locations:"
    count = 0
    for kmer in removeList:
        count += 1
        print str(count) + "\t" + kmer
        for location in locKmerDict:
            for protein in locKmerDict[location]:
                 if kmer in locKmerDict[location][protein]:
                     del locKmerDict[location][protein][kmer]


    return locKmerDict