import os
import urllib2
import subprocess

__author__ = 'delur'

def build_mfasta(protein, sequence, profileProteins, paths, overwrite, queue):
    uniqueprot_call = [paths["uniqueprot"], '-i' ,os.path.join(paths["mfasta"], protein + ".mfasta"),'-o',os.path.join(paths["mfasta"], protein + ".clean"), '-t', '20']
    if not os.path.exists(paths["mfasta"]):
        os.makedirs(paths["mfasta"])
    def grabAndWriteFasta(prot, target_dir):
            response = urllib2.urlopen("http://www.uniprot.org/uniprot/" + prot + ".fasta")
            fasta = response.read()
            if str(fasta) == "":
                print "\t\t\tDid not find uniprot entry for ProfileProtein: " + prot + " !"
            else:
                fastafile = open(os.path.join(target_dir, prot+".fa"), 'w')
                fastafile.write(fasta)
                fastafile.close()
                return fasta

    def loadFasta(protein, profileProteins, paths):
        mfastafile = open(os.path.join(paths["mfasta"], protein+".mfasta"), 'w')
        mfastafile.write("")
        mfastafile.close()

        mfastafile = open(os.path.join(paths["mfasta"], protein+".mfasta"), 'a')
        for prot in profileProteins:
            print prot
            fasta = ""
            tmp = prot[0]
            target_dir = os.path.join(paths["fasta"], tmp)
            if os.path.isfile(os.path.join(target_dir, prot + ".fa")):
                if overwrite:
                    fasta = grabAndWriteFasta(prot, target_dir)
                else:
                    fastafile = open(os.path.join(target_dir, prot+".fa"), 'r')
                    fasta = fastafile.read()
                    fastafile.close()
            else:
                fasta = grabAndWriteFasta(prot, target_dir)
            if fasta is not None:
                mfastafile.write(str(fasta))
        mfastafile.close()

    if os.path.isfile(os.path.join(paths["mfasta"], protein + ".mfasta")):
        if overwrite:
            loadFasta(protein, profileProteins, paths)
    else:
        loadFasta(protein, profileProteins, paths)

    def clean_mfasta(protein, paths):
        subprocess.call(paths["qsub"] + ["-l", "long=TRUE"] + uniqueprot_call if queue else uniqueprot_call)
        # with open(os.path.join(paths["mfasta"], protein + ".clean"), "r+") as f:
        #         old = f.read()
        #         f.seek(0)
        #         f.write(">" +  protein + "\n" + sequence + "\n" + old)

    if os.path.isfile(os.path.join(paths["mfasta"], protein + ".clean")):
        if overwrite:
            clean_mfasta(protein, paths)
    else:
        clean_mfasta(protein, paths)
