import os
import urllib2
from masterthesis.reader import read_uniprotfile
__author__ = 'delur'

def get_uniprot(protein, paths, overwrite):
    if not os.path.exists(paths["uniprot"]):
        os.makedirs(paths["uniprot"])

    def loadUni(protein, paths):
        response = urllib2.urlopen("http://www.uniprot.org/uniprot/" + protein + ".txt")
        entry = response.read()
        if str(entry) == "":
            print "Did not find uniprot entry for: " + protein + " !"
            return False, ""
        else:
            f = open(os.path.join(paths["uniprot"], protein+".txt"), 'w')
            f.write(entry)
            f.close()
            entry = entry.split('\n')
            return True, entry

    if os.path.isfile(os.path.join(paths["uniprot"], protein + ".txt")):
        if overwrite:
            return loadUni(protein,paths)
        else:
            return True, read_uniprotfile(protein, paths)
    else:
        return loadUni(protein, paths)
