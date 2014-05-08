import sys

__author__ = 'delur'


def findHit(patternMatches, i):
    for start, end in patternMatches:
        if i in range(start, end):
            return True
    return False


def doEvaluation(patternMatches, zscore_count, cutoff):

    TruePositives = 0
    FalsePositives = 0
    TrueNegatives = 0
    FalseNegatives = 0

    for i in range(len(zscore_count)):
        hits = 0
        if zscore_count[i] > cutoff:
            if findHit(patternMatches, i):
                TruePositives += 1
            else:
                FalsePositives += 1
        else:
            if findHit(patternMatches, i):
                FalseNegatives += 1
            else:
                TrueNegatives += 1

    print "\t|Positives\t|Negatives"
    print "True\t|" + str(TruePositives) + "\t\t|" + str(TrueNegatives)
    print "False\t|" + str(FalsePositives) + "\t\t|" + str(FalseNegatives)

    if TruePositives == 0:
        precision = 0.0
        recall = 0.0
    else:
        precision = TruePositives / (float(TruePositives) + float(FalsePositives))
        recall = TruePositives / (float(TruePositives) + float(FalseNegatives))

    return precision, recall


def getKmerSegments(zscore_count, cutoff):
    kmerSegments = []
    start = 0
    end = 0
    run = False
    for i in range(len(zscore_count)):
        if run == False:
            if zscore_count[i] > cutoff:
                run = True
                start = i
                end = i
        else:
            if zscore_count[i] > cutoff:
                end = i
            else:
                run = False
                kmerSegments.append( (start, end))
    if zscore_count[len(zscore_count)-1] > cutoff:
        kmerSegments.append( (start, end))
    return  kmerSegments


def getPatternSegmentSequence(patternSegments, size):
    patternSegmentSequence = [False] * size
    for start, end in patternSegments:
        for i in range(start, end):
            patternSegmentSequence[i] = True
    return patternSegmentSequence

def prefilterPatternSegments(patternSegments, zscore_count):
    filteredPatternSegments = []
    for start, end in patternSegments:
        if start - end >= 0.5* len(zscore_count):
            continue
        else:
            filteredPatternSegments.append((start, end))
    return filteredPatternSegments

def doSegmentEvaluation(patternSegments, zscore_count, cutoff):
    patternSegmentSequence = getPatternSegmentSequence(patternSegments, len(zscore_count))

    kmerSegments = getKmerSegments(zscore_count, cutoff)

    print patternSegments
    print kmerSegments
    PatternSegmentsHit = 0
    PatternSegmentsMissed = 0
    KmerSegmentsHit = 0
    KmerSegmentsMissed = 0


    #filtering out too long segments this is allready done in evalNoPlot.py
    # patternSegments = prefilterPatternSegments(patternSegments, zscore_count)

    for start, end in patternSegments:
        if start - end >= 0.5* len(zscore_count):
            continue
        hits = 0
        for i in range(start, end+1):
            if zscore_count[i] > cutoff:
                hits += 1

        if hits >= 5 or hits >= len(range(start, end+1)) :
            PatternSegmentsHit += 1
        else:
            PatternSegmentsMissed += 1

    for start, end in kmerSegments:
        hits = 0
        for i in range(start, end):
            if patternSegmentSequence[i]:
                hits += 1

        if hits >= 5 or hits >= len(range(start, end+1)) :
            KmerSegmentsHit += 1
        else:
            KmerSegmentsMissed += 1

    print "PatternSegments Hit: " + str(PatternSegmentsHit)
    print "PatternSegments Missed: " + str(PatternSegmentsMissed)
    print "KmerSegments Hit: " + str(KmerSegmentsHit)
    print "KmerSegments Missed: " + str(KmerSegmentsMissed)

    if PatternSegmentsHit == 0:
        precision = 0.0
        recall = 0.0
    else:
        precision = PatternSegmentsHit / (float(PatternSegmentsHit) + float(KmerSegmentsMissed))
        recall = PatternSegmentsHit / (float(PatternSegmentsHit) + float(PatternSegmentsMissed))

    return precision, recall


def doCombinedEvaluation(patternMatches, prosite, zscore_count, cutoff):

    TruePositives = 0
    FalsePositives = 0
    TrueNegatives = 0
    FalseNegatives = 0

    for i in range(len(zscore_count)):
        hits = 0
        if zscore_count[i] > cutoff:
            if findHit(patternMatches, i):
                TruePositives += 1
            elif findHit(prosite, i):
                TruePositives += 1
            else:
                FalsePositives += 1
        else:
            if findHit(patternMatches, i):
                FalseNegatives += 1
            elif findHit(prosite, i):
                FalseNegatives += 1
            else:
                TrueNegatives += 1

    print "\t\t|Positives\t|Negatives"
    print "True\t|" + str(TruePositives) + "\t\t\t|" + str(TrueNegatives)
    print "False\t|" + str(FalsePositives) + "\t\t\t|" + str(FalseNegatives)

    if TruePositives == 0:
        precision = 0.0
        recall = 0.0
    else:
        precision = TruePositives / (float(TruePositives) + float(FalsePositives))
        recall = TruePositives / (float(TruePositives) + float(FalseNegatives))

    return precision, recall


