import os
import subprocess

__author__ = 'delur'

def wget_fasta(profileProteins, constants, overwrite):
    if not os.path.exists(constants["mfasta"]):
        os.makedirs(constants["mfasta"])
    def grabFasta(prot, paths):
            subprocess.Popen(["wget", "--no-verbose", "http://www.uniprot.org/uniprot/" + prot + ".fasta", "-O", os.path.join(target_dir, prot+".fa")])

    for prot in profileProteins:
        tmp = prot[0]
        target_dir = os.path.join(constants["fasta"], tmp)
        if not os.path.exists(target_dir):
                os.makedirs(target_dir)

        if os.path.isfile(os.path.join(target_dir, prot + ".fa")):
            if overwrite:
                grabFasta(prot, constants)
        else:
            grabFasta(prot, constants)