import os

__author__ = 'delur'

def get_fasta(protein, entry, paths, overwrite):
    if not os.path.exists(paths["fasta"]):
        os.makedirs(paths["fasta"])

    start = -1
    end = -1
    sequence = ""
    for i in range(len(entry)):
        if entry[i].startswith("SQ"):
            start = i+1
        if entry[i].startswith("//"):
            end = i

    for i in range(start, end):
        sequence += entry[i]

    sequences = sequence.split(' ')
    sequence = ""
    for line in sequences:
        sequence += line

    def write_fasta(protein, sequence, paths):
        f = open(os.path.join(paths["fasta"], protein+".fa"), 'w')
        f.write(">" + protein + "\n" + sequence)
        f.close()

    if os.path.isfile(os.path.join(paths["fasta"], protein + ".fa")):
        if overwrite:
            write_fasta(protein, sequence, paths)
    else:
        write_fasta(protein, sequence, paths)
    return sequence