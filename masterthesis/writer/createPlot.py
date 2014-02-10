import os
import sys

__author__ = 'delur'
def getFeatures(entry):
    features = {}

    for line in entry:
        if line.startswith("FT   "):
            #print line
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

def create_plot(query_protein_sequence, pro_matches, con_matches, entry, numProfileProteins, resultfile_info, paths, svm):
    name = query_protein_sequence[0]
    print name
    sequence = query_protein_sequence[1]
    #print sequence

    pos_count = None
    pos_count_noGaps = [0] * len(sequence)
    neg_count = None
    neg_count_noGaps = [0] * len(sequence)


    for prot in pro_matches:
        for match_seq , start, end in pro_matches[prot]:
            pos_count = [0] * len(match_seq)
            for j in range(start, end):
                if j < len(match_seq):
                    pos_count[j] += 1
            x = 0
            for i in range(len(match_seq)):
                if x == len(sequence):
                    break
                if match_seq[i] == '-':
                    pass
                else:
                    pos_count_noGaps[x] += pos_count[i]
                    x += 1

    for prot in con_matches:
        for match_seq , start, end in con_matches[prot]:
            neg_count = [0] * len(match_seq)
            for j in range(start, end):
                if j < len(match_seq):
                    neg_count[j] += 1
            x = 0
            for i in range(len(match_seq)):
                if x == len(sequence):
                    break
                if match_seq[i] == '-':
                    pass
                else:
                    neg_count_noGaps[x] += neg_count[i]
                    x += 1

    for i in range(len(sequence)):
        pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins
        neg_count_noGaps[i] = neg_count_noGaps[i] * 100 / numProfileProteins

    # pos_count = [0] * len(sequence)
    # for protein in pos_matches:
    #     for start,end in pos_matches[protein]:
    #         if start == -1:
    #             print "start = " + str(start)
    #             sys.exit()
    #         if end == -1:
    #             print "end = " + str(end)
    #             sys.exit()
    #         for j in range(start, end):
    #             if j < len(sequence):
    #                 pos_count[j] += 1
    #
    # neg_count = [0] * len(sequence)
    # for protein in neg_matches:
    #     for start,end in neg_matches[protein]:
    #         for j in range(start, end):
    #             if j < len(sequence):
    #                 neg_count[j] += 1
    #
    # pos_count_noGaps = []
    # neg_count_noGaps = []
    # seq_noGap = ""
    # text = ""
    # for i in range(len(sequence)):
    #     if sequence[i] == '-':
    #         pass
    #     else:
    #         text += str(pos_count[i])
    #         pos_count_noGaps += [pos_count[i]]
    #         neg_count_noGaps += [neg_count[i]]
    #         seq_noGap += sequence[i]
    #
    # for i in range(len(pos_count_noGaps)):
    #     pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / float(numProfileProteins)
    #     neg_count_noGaps[i] = neg_count_noGaps[i] * 100 / float(numProfileProteins)

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

    x = range(1, len(sequence)+1)
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


    plt.xticks(x, sequence)




    ax.set_yticks( [0,50,100])
    ax.set_yticklabels(["0", "50", "100"])
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
    plt.ylim( -10 +ypos, 100)
    plt.xlim( plt.xlim()[0], len(sequence)+1)
    #plt.xlim( (i*200,i*200 + 200))
    plt.tight_layout()
    pdfpath = os.path.join(os.path.join(paths["pdf"], svm), location)
    if not os.path.exists(pdfpath):
        os.makedirs(pdfpath)
    #print os.path.join(pdfpath, name + ".pdf")
    plt.savefig(os.path.join(pdfpath, name + ".pdf"))
