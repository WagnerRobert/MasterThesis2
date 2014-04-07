__author__ = 'delur'
def getlocSeqDict(path):

    f = open(path, 'r')
    name = ""
    localisation = ""
    locSeqDict= {}
    i = 0
    for line in f:
        if line.startswith(">"):
            tmp = line.split(' ')
            name = tmp[0][1:]
            localisation = tmp[1].rstrip()
            if localisation not in locSeqDict:
                locSeqDict[localisation] = {}
            if name not in locSeqDict[localisation]:
                locSeqDict[localisation][name] = ""
            else:
                print name + " is already in " + localisation
        else:
            locSeqDict[localisation][name] = line.rstrip()
    f.close()
    return  locSeqDict
