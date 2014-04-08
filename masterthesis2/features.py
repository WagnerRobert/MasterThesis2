def getFeatures(entry):
    features = {}

    for line in entry:
        if line.startswith("FT   "):
            #print line
            tmp = line.rstrip().split()
            if len(tmp) > 3 :
                if tmp[1] in features:
                    pass
                else:
                   features[tmp[1]] = []

                try:
                    start = int(tmp[2])
                    end = int(tmp[3])
                    if "By similarity" in line or "Potential" in line or "Propable" in line :
                        features[tmp[1]].append( (start, end, False) )
                    else:
                        features[tmp[1]].append( (start, end, True) )
                except ValueError:
                    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!Encountered slight Problem"
                    continue

    return features