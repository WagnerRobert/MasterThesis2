import os
import re
import subprocess
import sys

__author__ = 'delur'

def blastProtein(protein, constants, overwrite, queue):
    blast_call = [constants["blast_tool"], '-F', 'F', '-a', '1', '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat']
    if not os.path.exists(constants["blast_dir"]):
        os.makedirs(constants["blast_dir"])
    if os.path.isfile(os.path.join(constants["blast_dir"], protein + ".blast")):
        if overwrite:
            subprocess.call(constants["qsub"]+blast_call if queue else blast_call)
            if queue:
                return
            #subprocess.call([constants["blast_tool"], '-F', 'F', '-a', constants["num_cores"], '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-m', '8', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])
    else:
        subprocess.call(constants["qsub"]+blast_call if queue else blast_call)
        if queue:
            return
        #subprocess.call([constants["blast_tool"], '-F', 'F', '-a', constants["num_cores"], '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-m', '8', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])

    ProfileProteines = []

    def text_blast():
        f = open (os.path.join(constants["blast_dir"], protein + ".blast"), 'r')
        run = False


        inSequence = False
        SequenceName = ""
        SequenceEntry = ""
        curRound = 0
        RoundList = [None,{},{},{}]
        importantRound = []
        for line in f:
            if line.startswith("Results from round"):
                curRound += 1

            if line.startswith(">tr") or line.startswith(">sp"):
                if inSequence:
                    RoundList[curRound][SequenceName] = SequenceEntry
                    SequenceEntry = ""
                    inSequence = False
                if not inSequence:
                    RoundList[curRound][line.split('|')[1]] = ""
                    SequenceName = line.split('|')[1]
                    inSequence = True
            elif inSequence:
                SequenceEntry += line.rstrip()

            if "round 3" in line:
                importantRound = RoundList[2]
                break
            if "CONVERGED!" in line:
                if curRound == 3:
                    importantRound = RoundList[2]
                elif curRound == 2:
                    importantRound = RoundList[1]
                elif curRound == 1:
                    importantRound = []
        f.close()
        regex = re.compile(r"Expect = (\S+),")
        for sequence in importantRound:
            #print importantRound[sequence]
            #print sequence
            match = re.search(regex, importantRound[sequence])
            print sequence + "\t" + match.group(1)
        sys.exit


    def tab_blast():
        f = open (os.path.join(constants["blast_dir"], protein + ".blast"), 'r')
        for line in  f:
            tmp = line.split('\t')
            if tmp[1].startswith("tr|") or tmp[1].startswith("sp|"):
                if tmp[1].split('|')[1] not in ProfileProteines:
                    ProfileProteines.append(tmp[1].split('|')[1])
            else:
                print "\t\t\tignoring pdb entry " + tmp[1].rstrip()
        f.close()

    #tab_blast()
    text_blast()

    ProfileProteines_copy = []
    for entry in ProfileProteines:
        ProfileProteines_copy.append(entry)

    for i in range(len(ProfileProteines_copy)):
        if ProfileProteines_copy[i] == protein:
            ProfileProteines.pop(i)

    return ProfileProteines
