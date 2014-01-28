__author__ = 'delur'
import sys

def processProtein(name, kmerlists, result_info, tree):
    pro_kmerlist = []
    con_kmerlist = []
    if result_info[0] in tree[0]:
        print tree[0]
        pro_kmerlist = kmerlists[1]
        con_kmerlist = kmerlists[0]
    elif result_info[0] in tree[1]:
        print tree[1]
        pro_kmerlist = kmerlists[0]
        con_kmerlist = kmerlists[1]

    clean_name = name.split('#')[0]
    getUniprot(clean_name)
    sys.exit()

    pass