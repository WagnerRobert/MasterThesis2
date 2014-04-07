import operator
import masterthesis2.zscore

__author__ = 'delur'

def cleanKmers(locKmerList, locSeqDict):
        # calculate the zscore
    print "calculating the zscores on the kmer vaules"
    locKmerList = None

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

    removelist = []
    for location in locKmerList :
        for protein in locKmerList[location]:
            for kmer, value in locKmerList[location][protein]:
                if kmer not in removelist:
                    for loc in locSeqDict:
                        if kmer in removelist:
                            break
                        elif locTree2Uniprot[location] != loc :
                            for protein in locSeqDict[loc]:
                                if kmer in locSeqDict[loc][protein]:
                                    removelist.append(kmer)
                                    break
    print "\nFound the following kmers in multiple locations:"
    count = 0
    for kmer in removelist:
        count += 1
        print str(count) + "\t" + kmer
        for location in locKmerList:
            for protein in locKmerList[location]:
                 if kmer in locKmerList[location][protein]:
                     del locKmerList[location][protein][kmer]


    return locKmerList