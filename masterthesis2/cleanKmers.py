import operator
import masterthesis2.zscore

__author__ = 'delur'

def cleanKmers(locKmerList, locSeqDict):

    #     # get the kmers of each location, count how often which one occurs and the highest(all?) zscore(s)
    # print "counting all the kmers for each location"
    # locKmerDict = {}
    # for location in locZscore:
    #     locKmerDict[location] = {}
    #     for protein in locZscore[location]:
    #         for kmer, value in locZscore[location][protein]:
    #             if kmer not in locKmerDict[location]:
    #                 locKmerDict[location][kmer] = []
    #             locKmerDict[location][kmer].append(value)
    # locZscore = None

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
    for location1 in locKmerList:
        i += 1
        print "location1: " + str(i) + "\t" + str(len(locKmerList))
        for location2 in locKmerList:
            if location1 != location2:
                loc1KmerList = []
                print "LOCATION 1 KMERLIST:"
                for protein in locKmerList[location1]:
                    for kmer in locKmerList[location1][protein]:
                        print kmer
                        loc1KmerList.append(kmer)
                loc2KmerSet = set()
                print "LOCATION 2 KMERLIST:"
                for protein in locKmerList[location2]:
                    for kmer in locKmerList[location2][protein]:
                        print kmer
                        loc2KmerSet.add(kmer)
                localRemovelist = [val for val in loc1KmerList if val in loc2KmerSet]
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
        for location in locKmerList:
            for protein in locKmerList[location]:
                 if kmer in locKmerList[location][protein]:
                     del locKmerList[location][protein][kmer]


    return locKmerList