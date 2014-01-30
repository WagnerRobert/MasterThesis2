import os
import urllib2
import subprocess

__author__ = 'delur'

def build_mfasta(protein, sequence, profileProteins, paths, overwrite):
    if not os.path.exists(paths["mfasta"]):
        os.makedirs(paths["mfasta"])
    def grabAndWriteFasta(prot, paths):
            response = urllib2.urlopen("http://www.uniprot.org/uniprot/" + prot + ".fasta")
            fasta = response.read()
            if str(fasta) == "":
                print "\t\t\tDid not find uniprot entry for ProfileProtein: " + prot + " !"
            else:
                fastafile = open(os.path.join(paths["fasta"], prot+".fa"), 'w')
                fastafile.write(fasta)
                fastafile.close()
                return fasta

    def loadFasta(protein, profileProteins, paths):
        mfastafile = open(os.path.join(paths["mfasta"], protein+".mfasta"), 'w')
        mfastafile.write("")
        mfastafile.close()

        mfastafile = open(os.path.join(paths["mfasta"], protein+".mfasta"), 'a')
        for prot in profileProteins:
            fasta = ""
            if os.path.isfile(os.path.join(paths["fasta"], prot + ".fa")):
                if overwrite:
                    fasta = grabAndWriteFasta(prot, paths)
                else:
                    fastafile = open(os.path.join(paths["fasta"], prot+".fa"), 'r')
                    fasta = fastafile.read()
                    fastafile.close()
            else:
                fasta = grabAndWriteFasta(prot, paths)
            if fasta is not None:
                mfastafile.write(str(fasta))
        mfastafile.close()

    if os.path.isfile(os.path.join(paths["mfasta"], protein + ".mfasta")):
        if overwrite:
            loadFasta(protein, profileProteins, paths)
    else:
        loadFasta(protein, profileProteins, paths)

    def clean_mfasta(protein, paths):
        subprocess.call([paths["uniqueprot"], '-i' ,os.path.join(paths["mfasta"], protein + ".mfasta"),'-o',os.path.join(paths["mfasta"], protein + ".clean"), '-t', '20'])
        # with open(os.path.join(paths["mfasta"], protein + ".clean"), "r+") as f:
        #         old = f.read()
        #         f.seek(0)
        #         f.write(">" +  protein + "\n" + sequence + "\n" + old)

    if os.path.isfile(os.path.join(paths["mfasta"], protein + ".clean")):
        if overwrite:
            clean_mfasta(protein, paths)
    else:
        clean_mfasta(protein, paths)
