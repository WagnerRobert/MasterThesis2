import os
import subprocess

__author__ = 'delur'

def blastProtein(protein, constants, overwrite):
    if not os.path.exists(constants["blast_dir"]):
        os.makedirs(constants["blast_dir"])
    if os.path.isfile(os.path.join(constants["blast_dir"], protein + ".blast")):
        if overwrite:
            subprocess.call([constants["blast_tool"], '-F', 'F', '-a', '1', '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])
            #subprocess.call([constants["blast_tool"], '-F', 'F', '-a', constants["num_cores"], '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-m', '8', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])
    else:
        subprocess.call([constants["blast_tool"], '-F', 'F', '-a', '1', '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])
        #subprocess.call([constants["blast_tool"], '-F', 'F', '-a', constants["num_cores"], '-j', '3', '-b', '3000', '-e', '1', '-h', '1e-3', '-m', '8', '-d', constants["big_80"], '-i', os.path.join(constants["fasta"], protein + ".fa"), '-o', os.path.join(constants["blast_dir"], protein + ".blast"), '-C', 'tmpfile.chk', '-Q', 'tmpfile.blastPsiMat'])

    ProfileProteines = []

    def text_blast():
        f = open (os.path.join(constants["blast_dir"], protein + ".blast"), 'r')
        run3 = False

        for line in f:
            if "round 3" in line:
                run3 = True
            if run3:
                if line.startswith(">tr") or line.startswith(">sp"):
                    ProfileProteines.append(line.split('|')[1])
                elif line.startswith(">"):
                    print "\t\t\tignoring entry " + line
        f.close()

    def tab_blast():
        f = open (os.path.join(constants["blast_dir"], protein + ".blast"), 'r')
        for line in  f:
            tmp = line.split('\t')
            if tmp[1].startswith("tr|") or tmp[1].startswith("sp|"):
                if tmp[1].split('|')[1] not in ProfileProteines:
                    ProfileProteines.append(tmp[1].split('|')[1])
            else:
                print "\t\t\tignoring pdb entry " + tmp[1]
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
