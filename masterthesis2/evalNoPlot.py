import os
import sys
import masterthesis2.features
import masterthesis2.evaluation
from scipy import stats



def eval_without_plot(query_protein_sequence, pro_matches, numProfileProteins, patternMatches, top_matches, cutoff):
    name = query_protein_sequence[0]
    #print name
    sequence = query_protein_sequence[1]
    #print sequence

    pos_count = None
    pos_count_noGaps = [0] * len(sequence)
    top_count_noGaps = [0] * len(sequence)


    #print pro_matches

    for prot in pro_matches:
        for match_seq , start, end in pro_matches[prot]:
            pos_count = [0] * len(match_seq)
            for j in range(start, end):
                if j < len(match_seq):
                    pos_count[j] += 1
            x = 0
            for i in range(len(match_seq)):
                if x == len(sequence):
                    break
                if match_seq[i] == '-':
                    pass
                else:
                    pos_count_noGaps[x] += pos_count[i]
                    x += 1
    for prot in top_matches:
        for match_seq , start, end in top_matches[prot]:
            top_count = [0] * len(match_seq)
            for j in range(start, end):
                if j < len(match_seq):
                    top_count[j] += 1
            x = 0
            for i in range(len(match_seq)):
                if x == len(sequence):
                    break
                if match_seq[i] == '-':
                    pass
                else:
                    top_count_noGaps[x] += top_count[i]
                    x += 1


    zscore_count = stats.zscore(pos_count_noGaps)


    for i in range(len(sequence)):
        pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins
        top_count_noGaps[i] = top_count_noGaps[i] * 100 / numProfileProteins

    patternMatches = masterthesis2.evaluation.prefilterPatternSegments(patternMatches,zscore_count)
    if len(patternMatches) == 0:
        print "protein has no nlsDB motif matches and is therefore not considered for the evaluation"
        return None
    precision, recall = masterthesis2.evaluation.doSegmentEvaluation(patternMatches,zscore_count, cutoff)


    return precision, recall