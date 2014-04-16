import os
import sys
import masterthesis2.features
from scipy import stats



def create_plot(query_protein_sequence, pro_matches, entry, numProfileProteins, resultfile_info, paths, prosite, patternMatches):
    name = query_protein_sequence[0]
    #print name
    sequence = query_protein_sequence[1]
    #print sequence

    pos_count = None
    pos_count_noGaps = [0] * len(sequence)



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

    zscore_count = stats.zscore(pos_count_noGaps)


    for i in range(len(sequence)):
        pos_count_noGaps[i] = pos_count_noGaps[i] * 100 / numProfileProteins

    import matplotlib.pyplot as plt

    x = range(1, len(sequence)+1)
    plt.clf()
    plt.cla()
    #print len(x)
    plt.plot(x,pos_count_noGaps, color='#336699', linestyle='-')
    plt.ylabel('Coverage')

    long_name = ""
    location = ""
    confidence = ""
    for protein in resultfile_info:
        if name in protein:
            long_name = protein
            location = resultfile_info[protein][0]
            confidence = resultfile_info[protein][1]
    plt.title(long_name + "|" + location + "|"+ confidence +"|" + str(numProfileProteins) + " PPs|")
    ax = plt.gca()


    plt.xticks(x, sequence)
    #[i.set_color("red") for i in plt.gca().get_xticklabels()]
    colorCode = {"H": "blue", "K": "blue", "R": "blue", "D": "red", "E": "red", "S": "green", "T": "green", "N": "green", "Q": "green", "C" : "yellow", "F" : "#AA00FF", "Y" : "#AA00FF", "W" : "#AA00FF" }
    xticklabels = plt.gca().get_xticklabels()
    for i in xticklabels:
        aminoacid = i.get_text()
        if aminoacid in colorCode:
            color = colorCode[aminoacid]
        else:
            color = "black"
        i.set_color(color)



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

    features = masterthesis2.features.getFeatures(entry)

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
            plt.text(start, ypos + (height/3.0), feature)
            plt.gca().add_patch(rect)

    for start, end in prosite:
        rect = plt.Rectangle((start - 0.5, 0), end-start + 1, 100, facecolor="#0000FF")
        rect.set_alpha(0.5)
        plt.gca().add_patch(rect)

    for start, end in patternMatches:
        rect = plt.Rectangle((start - 0.5, 0), end-start + 1, 100, facecolor="#00FF00")
        rect.set_alpha(0.5)
        plt.gca().add_patch(rect)



    for i in range(len(zscore_count)):
        if zscore_count[i] > 0.5:
            plt.plot([i+1, i+1], [-1,1], color='red')


    #for i in range (0,6):
    fig = plt.gcf()
    fig.set_size_inches(2+ (len(x)/10),4)
    plt.ylim( -10 +ypos, 100)
    plt.xlim( plt.xlim()[0], len(sequence)+1)
    #plt.xlim( (i*200,i*200 + 200))
    plt.tight_layout()
    pdfpath = os.path.join(os.path.join(paths["pdf"], "prosite"),location)
    if not os.path.exists(pdfpath):
        os.makedirs(pdfpath)
    #print os.path.join(pdfpath, name + ".pdf")
    plt.savefig(os.path.join(pdfpath, name + ".pdf"))