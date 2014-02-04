import os
import subprocess

__author__ = 'delur'

def build_pairwise_alignments(protein, constants, overwrite):
    if not os.path.exists(constants["needle_dir"]):
        os.makedirs(constants["needle_dir"])

    def runNeedle(protein, paths):
        if os.path.isfile(os.path.join(constants["mfasta"], protein + ".clean")):
            f = open(os.path.join(constants["mfasta"], protein + ".clean"), 'r')
            text = f.read()
            if text == "":
                pass
            else:
                subprocess.call([paths["needle"], '-asequence' ,os.path.join(paths["fasta"], protein + ".fa"),'-bsequence',os.path.join(paths["mfasta"], protein + ".clean"), '-outfile', os.path.join(paths["needle_dir"], protein + ".needle"), '-gapopen', '10.0', '-gapextend', '0.5'])

    if os.path.isfile(os.path.join(constants["needle_dir"], protein + ".needle")):
        if overwrite:
            runNeedle(protein, constants)
    else:
        runNeedle(protein, constants)

    f = open(os.path.join(constants["needle_dir"], protein + ".needle"), 'r')
    pairwise = {}
    first = "-"
    second = "-"
    for line in f:
        if line.startswith("# 1: "):
            first = line.split(':')[1].strip()
            #print first
        elif line == "\n" or line.startswith(' '):
            pass
        elif line.startswith("# 2: "):
            second = line.split(':')[1].strip()
            #print second
            if first is not second:
                pairwise[(first,second)] = ["",""]
        elif first is not second:
            if line[0:13].rstrip() in first:
                tmp = line.split()
                pairwise[(first,second)][0] += tmp[2]
            if line[0:13].rstrip() in second:
                tmp = line.split()
                pairwise[(first,second)][1] += tmp[2]
    f.close()

    # for first,second in pairwise:
        # print first
        # print pairwise[(first,second)][0]
        # print pairwise[(first,second)][1]
        # print second
        # print ""

    return pairwise

