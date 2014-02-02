import os
import subprocess

__author__ = 'delur'

def wget_fasta(profileProteins, constants, overwrite):
    if not os.path.exists(constants["mfasta"]):
        os.makedirs(constants["mfasta"])
    def grabFasta(prot, paths):
            subprocess.Popen(["wget", "http://www.uniprot.org/uniprot/" + prot + ".fasta", "-O", os.path.join(constants["fasta"], prot+".fa")])
            if str(fasta) == "":
                print "\t\t\tDid not find uniprot entry for ProfileProtein: " + prot + " !"
            else:
                fastafile = open(os.path.join(paths["fasta"], prot+".fa"), 'w')
                fastafile.write(fasta)
                fastafile.close()
                return fasta

    for prot in profileProteins:
        fasta = ""
        if os.path.isfile(os.path.join(constants["fasta"], prot + ".fa")):
            if overwrite:
                grabFasta(prot, constants)
        else:
            fasta = grabFasta(prot, constants)