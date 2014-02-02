import os
import datetime
from masterthesis.match_kmers import match_kmers, match_kmers_pairwise
from masterthesis.writer.blastProtein import blastProtein
from masterthesis.writer import get_uniprot
from masterthesis.writer import get_fasta
from masterthesis.writer.build_mfasta import build_mfasta
from masterthesis.writer.build_msa import build_msa
from masterthesis.writer.build_pairwise_alignments import build_pairwise_alignments

__author__ = 'delur'
import sys


def getFeatures(entry):
    features = {}

    for line in entry:
        if line.startswith("FT   "):
            print line
            tmp = line.rstrip().split()
            if len(tmp) > 3 :
                if tmp[1] in features:
                    pass
                else:
                   features[tmp[1]] = []

                try:
                    start = int(tmp[2])
                    end = int(tmp[3])
                    if "By similarity" in line or "Potential" in line or "Propable" in line :
                        features[tmp[1]].append( (start, end, False) )
                    else:
                        features[tmp[1]].append( (start, end, True) )
                except ValueError:
                    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!Encountered slight Problem"
                    continue

    return features


def create_plot(query_protein_sequence, pos_matches, neg_matches, entry, numProfileProteins, resultfile_info, paths):
    name = query_protein_sequence[0]
    print name
    sequence = query_protein_sequence[1]
    print sequence

    pos_count = [0] * len(sequence)
    for protein in pos_matches:
        for start,end in pos_matches[protein]:
            if start == -1:
                print "start = " + str(start)
                sys.exit()
            if end == -1:
                print "end = " + str(end)
                sys.exit()

            for j in range(start, end):
                pos_count[j] += 1

    neg_count = [0] * len(sequence)
    for protein in neg_matches:
        for start,end in neg_matches[protein]:
            for j in range(start, end):
                neg_count[j] += 1

    pos_count_noGaps = []
    neg_count_noGaps = []
    seq_noGap = ""
    text = ""
    for i in range(len(sequence)):
        if sequence[i] == '-':
            pass
        else:
            text += str(pos_count[i])
            pos_count_noGaps += [pos_count[i]]
            neg_count_noGaps += [neg_count[i]]
            seq_noGap += sequence[i]

    for i in range(len(pos_count_noGaps)):
        pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins
        neg_count_noGaps[i] = neg_count_noGaps[i] * 100 / numProfileProteins

    # positions_above_50 = []
    # for i in range(len(pos_count_noGaps)):
    #     if pos_count_noGaps[i] >= 50:
    #         positions_above_50.append(i)
    #
    # length = 1
    # position_length_list = []
    # for i in range(1,len(positions_above_50)):
    #     if positions_above_50[i-1] +1 == positions_above_50[i]:
    #         length += 1
    #     else:
    #         position_length_list.append(length)
    #         length = 1
    #
    # avg_len = 0.0
    # for i in range(len(position_length_list)):
    #     avg_len += position_length_list[i]
    # avg_len = avg_len / float(len(position_length_list))
    #
    # print "Average signal length is: " + str(avg_len)


    import matplotlib.pyplot as plt

    x = range(1, len(seq_noGap)+1)
    plt.clf()
    plt.cla()
    #print len(x)
    plt.plot(x,pos_count_noGaps, color='#336699')
    plt.plot(x, neg_count_noGaps, color='#CC0000')
    plt.ylabel('Coverage')

    long_name = ""
    location = ""
    confidence = ""
    for protein in resultfile_info:
        if name in protein:
            long_name = protein
            location = resultfile_info[protein][0]
            confidence = resultfile_info[protein][1]
    #plt.title(long_name + "|" + location + "|"+ confidence +"|" + str(numProfileProteins) + " PPs|" + "%.1f" % avg_len +" avg len pos")
    plt.title(long_name + "|" + location + "|"+ confidence +"|" + str(numProfileProteins) + " PPs|")
    ax = plt.gca()


    plt.xticks(x, seq_noGap)




    ax.set_yticks( [0,50,100,150,200])
    ax.set_yticklabels(["0", "50", "100", "150", "200"])
    ax.yaxis.grid(True)
    #ax.xaxis.grid(True)

    ypos = -20
    height = 20
    #rect.set_alpha(0.5)

    #rect = plt.Rectangle((20 - 0.5, -30), 20, 10, facecolor="#0000FF")
    #rect.set_alpha(0.5)
    #plt.gca().add_patch(rect)

    features = getFeatures(entry)

    for feature in features:
        if feature == "TURN":
            color = "#4169e1"
        elif feature == "STRAND":
            color = "#8470ff"
        elif feature == "HELIX":
            color = "#20b2aa"
        elif feature == "CHAIN":
            color = "#000000"
            continue
        elif feature == "METAL":
            color = "#708090"
        elif feature == "NP_BIND":
            color = "#eedd82"
        elif feature == "BINDING":
            color = "#ff8c00"
        elif feature == "ACT_SITE":
            color = "#ff0000"
        elif feature == "DOMAIN":
            color = "#32cd32"
        elif feature == "INIT_MET":
            color = "#a52a2a"
        elif feature == "TRANSMEM":
            color = "#FFFF00"
        elif feature == "TOPO_DOM":
            color = "#FF4500"
        elif feature == "CONFLICT":
            color = "#40e0d0"
        elif feature == "SIGNAL":
            color = "#800080"
        else:
            print feature
            color = "#FF1493"
        ypos = ypos - height

        for start,end, experimental in features[feature]:
            if experimental:
                rect = plt.Rectangle((start - 0.5, ypos), end-start + 1, height, facecolor=color)
            else:
                rect = plt.Rectangle((start - 0.5, ypos), end-start + 1, height, facecolor=color, hatch='//')
            plt.gca().add_patch(rect)


    #for i in range (0,6):
    fig = plt.gcf()
    fig.set_size_inches(2+ (len(x)/10),4)
    plt.ylim( -10 +ypos, 200)
    plt.xlim( plt.xlim()[0], len(seq_noGap)+1)
    #plt.xlim( (i*200,i*200 + 200))
    plt.tight_layout()
    plt.savefig(os.path.join(paths["pdf"], name + ".pdf"))


def processProtein(name, kmerlists, result_info, tree, constants):
    pro_kmerlist = []
    con_kmerlist = []
    if result_info[0] in tree[0]:
        pro_kmerlist = kmerlists[1]
        con_kmerlist = kmerlists[0]
    elif result_info[0] in tree[1]:
        pro_kmerlist = kmerlists[0]
        con_kmerlist = kmerlists[1]

    clean_name = name.split('#')[0]
    overwrite = False
    queue = True
    print "\t\tgetting Uniprot entry..."
    foundUniprot, entry = get_uniprot(clean_name, constants, overwrite)
    print "\t\tfinished"
    if not foundUniprot:
        print "!" + name + "has no uniprot entry"
        return
    sequence = get_fasta(clean_name, entry, constants, overwrite)
    print "\t\tblasting Protein against Big_80..." + str(datetime.datetime.now().time())
    profileProteines = blastProtein(clean_name, constants, overwrite, queue)
    print "\t\tfinished "  + str(datetime.datetime.now().time())
    print "\t\tloading fasta Files from Uniprot, building mfasta file, then cleaning it with uniqueprot"
    build_mfasta(clean_name, sequence, profileProteines, constants, overwrite, queue)
    print "\t\tfinished " + str(datetime.datetime.now().time())
    print "\t\tcalculating pairwise sequence alignments between query and profile protein sequences with needle"
    #msa = build_msa(clean_name, constants, overwrite)
    pairwise_alignments = build_pairwise_alignments(clean_name, constants, overwrite)
    print "\t\tfinished"

    pro_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, pro_kmerlist)
    con_matches = match_kmers_pairwise(clean_name, sequence, pairwise_alignments, con_kmerlist)
#    pro_matches = match_kmers(pro_kmerlist, msa)
#    con_matches = match_kmers(con_kmerlist, msa)
    create_plot((name,sequence), pro_matches, con_matches, entry, len(pairwise_alignments), result_info, constants)