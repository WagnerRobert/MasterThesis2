import os
import subprocess

__author__ = 'delur'

def build_msa(protein, constants, overwrite):
    if not os.path.exists(constants["msa"]):
        os.makedirs(constants["msa"])

    def runClustalO(protein, paths):
        subprocess.call([paths["clustalo"], '-i' ,os.path.join(paths["mfasta"], protein + ".clean"),'-o',os.path.join(paths["msa"], protein + ".msa"), '--outfmt=clu', '--force', '--wrap=9999', "--threads=" + constants["num_cores"]])

    if os.path.isfile(os.path.join(constants["msa"], protein + ".msa")):
        if overwrite:
            runClustalO(protein, constants)
    else:
        runClustalO(protein, constants)

    msa = []
    f = open (os.path.join(constants["msa"], protein + ".msa"), 'r')

    #skipp clustalO header
    i = 0
    for line in f:
        line = line.rstrip()
        i+=1
        if i > 3:
            if not line == "" and not line.startswith(" "):
                name = line.split(' ', 1)[0]
                sequence = line.rsplit(' ', 1)[1]
            #print name + " : " + sequence
            msa.append( (name, sequence))
    f.close()
    return msa