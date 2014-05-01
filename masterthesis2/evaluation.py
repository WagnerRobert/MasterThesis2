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
    print "True\t|" + str(TruePositives) + "\t|" + str(TrueNegatives)
    print "False\t|" + str(FalsePositives) + "\t|" + str(FalseNegatives)


