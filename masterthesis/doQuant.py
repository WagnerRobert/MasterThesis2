import sys

__author__ = 'delur'

def doQuant(kmerlist, quant):
    positive_list = []
    negative_list = []

    for kmer, value in kmerlist:
        if value < 0 :
            negative_list.append( (kmer, value) )
        else:
            positive_list.append( (kmer, value) )
    positive_list.reverse()

    total = 0.0
    for kmer, value in positive_list:
        total += value
    posquant = quant * total

    while total > posquant:
        total -= positive_list.pop()[1]

    total = 0.0
    for kmer, value in negative_list:
        total += value
    negquant = quant * total

    while total < negquant:
        total -= negative_list.pop()[1]

    return (positive_list, negative_list)