import re
import sys

__author__ = 'delur'

def match_kmers(kmerlisting, msa):
    if len(kmerlisting) == 0:
        return {}
    pattern = ""
    matches = {}

    for kmer_tuple in kmerlisting:
        kmer = kmer_tuple[0]
        for letter in kmer:
            pattern += letter + "[-]*"
        pattern = pattern[0:len(pattern)-4]
        pattern = pattern + "|"
    pattern = pattern[0:len(pattern)-1]
    regex = re.compile(pattern)

    for name, sequence in msa:
        i = 0
        while True:
            match = regex.search(sequence, i)
            if  match:
                if name in matches:
                    pass
                else:
                    matches[name] = []
                matches[name].append(match.span())
                i = match.start()+1
            else:
                break

    return matches


def do_mapping(sequence, gapped_sequence):
    gaps = 0
    mapping = []
    for i in range(len(gapped_sequence)):
        if gapped_sequence[i] is '-':
            gaps += 1
            mapping.append(-1)
        else:
            mapping.append(i-gaps)
    print "Len Mapping: " + str(len(mapping))
    print "Len sequence: " + str(len(sequence))
    print "Len gapped_sequence: " + str(len(gapped_sequence))

    print sequence
    print gapped_sequence
    for i in range(len(mapping)):
        print i
        print mapping[i]
        if mapping[i] == -1:
            pass
        elif sequence[mapping[i]] is not gapped_sequence[i]:
            # print "Fail!!!"
            # print sequence
            # print gapped_sequence
            # print i
            # print sequence[mapping[i]]
            # print gapped_sequence[i]
            # sys.exit()
            pass
    return mapping



def match_kmers_pairwise(clean_name, sequence, pairwise_alignments, kmerlisting):
    if len(kmerlisting) == 0:
        return {}
    pattern = ""
    matches = {}

    for kmer_tuple in kmerlisting:
        kmer = kmer_tuple[0]
        for letter in kmer:
            pattern += letter + "[-]*"
        pattern = pattern[0:len(pattern)-4]
        pattern = pattern + "|"
    pattern = pattern[0:len(pattern)-1]
    regex = re.compile(pattern)

    for first, second in pairwise_alignments:
        mapping = do_mapping(sequence, pairwise_alignments[first,second][0])
        i = 0
        while True:
            match = regex.search(pairwise_alignments[first,second][1], i)
            if  match:
                if second in matches:
                    pass
                else:
                    matches[second] = []
                span = match.span()
                text = match.group()
                start, end = span
                end = end -1
                for j in range(span[0],span[1]):
                    if mapping[start] == -1:
                        start += 1
                for j in range(start,span[1]):
                    if mapping[end] == -1:
                        end -= 1
                if start > end or mapping[end] == -1:
                    break
                start = mapping[start]
                end = mapping[end]
                matches[second].append( (start,end) )
                i = match.start()+1
            else:
                break

    return matches