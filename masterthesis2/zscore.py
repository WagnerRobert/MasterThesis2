import numpy as np
from scipy import stats

__author__ = 'wagnerr'

def calc_zscore(locKmerList):

    for location in locKmerList:
        for protein in locKmerList[location]:
            kmers = []
            values = []
            for kmer,value in locKmerList[location][protein]:
                kmers.append(kmer)
                values.append(value)
            values = stats.zscore(values)
            locKmerList[location][protein] = []
            for i in range(len(kmers)):
                locKmerList[location][protein].append( (kmers[i], values[i]))
    return locKmerList