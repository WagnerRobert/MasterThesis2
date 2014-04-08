import re

__author__ = 'wagnerr'


def match_kmers_pairwise(clean_name, sequence, pairwise_alignments, kmerlisting):
    if len(kmerlisting) == 0:
        return {}
    pattern = ""
    matches = {}

    for kmer in kmerlisting:
        for letter in kmer:
            pattern += letter + "[-]*"
        pattern = pattern[0:len(pattern)-4]
        pattern = pattern + "|"
    pattern = pattern[0:len(pattern)-1]
    regex = re.compile(pattern)


    for first, second in pairwise_alignments:
        i = 0
        while True:
            match = regex.search(pairwise_alignments[first,second][1], i)
            if  match:
                if second in matches:
                    pass
                else:
                    matches[second] = []

                start, end = match.span()
                end = end -1

                matches[second].append( (pairwise_alignments[first,second][0], start,end) )
                i = match.start()+1
            else:
                break

        while True:
            match = regex.search(sequence, i)
            if  match:
                if second in matches:
                    pass
                else:
                    matches[second] = []

                start, end = match.span()
                end = end -1

                matches[second].append( (pairwise_alignments[first,second][0], start,end) )
                i = match.start()+1
            else:
                break

    return matches