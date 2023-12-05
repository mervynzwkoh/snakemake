import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("source", help="The source location of the MSA file")
args = parser.parse_args()

# src is the ABSOLUTE path to MSA file (eg. /Users/xyz/covid/dataset.msa)
src = args.source
PATH = os.path.dirname(src)


# input: mafft msa file src
# output: dictionary of samples
def generateDict(src):
    f1 = open(src, "r")
    lines = f1.readlines()
    out = {}
    for i in range(0, len(lines)):
        line = lines[i].strip()
        if line and line[0] == ">":
            key = line[1:20]  # 9 for OU/OV, 11 for ERR
            # print(key)
            if i < len(lines):
                i += 1
                temp = lines[i].strip()
                while (
                    i + 1 < len(lines)
                    and lines[i].strip()
                    and str(lines[i].strip())[0] != ">"
                ):
                    temp += lines[i].strip()
                    i += 1
                temp = temp.strip()
                if True:
                    out[key] = temp
                    continue
    f1.close()
    return out  ## out is dictionary of terms in fasta


with open(PATH + "/msadict.txt", "w") as f:
    dictSamples = generateDict(src)
    f.write(str(dictSamples))
# print(os.path(dictSamples))

import ast

# if msadict already made, input here
msadictpath = PATH + "/msadict.txt"
with open(msadictpath, "r") as f:
    dictSamples = ast.literal_eval(f.read())
PATH = os.path.dirname(msadictpath)


#################################
#            Distances          #
#################################
def estimate_nucleotide_frequencies(seq):
    # seq = seq.replace('-','').upper()
    # seq = seq.replace('N','').upper() #remove all N
    seq = seq.upper()
    A = seq.count("A")
    C = seq.count("C")
    G = seq.count("G")
    T = seq.count("T")
    length = float(len(seq))
    return [x / length for x in [A, C, G, T]]


def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1, seq2):
        pairs.append(x)
    # for (x,y) in zip(seq1,seq2):
    for x, y in pairs:
        if x != y:
            p += 1
    # length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    return float(p) / length


def levendistance(seq1, seq2):
    import Levenshtein as lv

    # seq1 = seq1.replace('-','').upper()
    # seq1 = seq1.replace('N','').upper() #remove all N
    # seq2 = seq2.replace('-','').upper()
    # seq2 = seq2.replace('N','').upper() #remove all N
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    return 1 - lv.ratio(seq1, seq2)


def JCdistance(
    seq1, seq2
):  # https://github.com/hboutou/jukes_cantor/blob/master/jukes_cantor.py
    """
    distance = -b log(1 - p / b)
    where:
    b = 3/4
    and p = p-distance, i.e. uncorrected distance between seq1 and seq2
    """
    from math import log

    b = 0.75
    p = pdistance(seq1, seq2)
    try:
        d = -b * log(1 - p / b)
    except ValueError:
        print("Tried to take log of a negative number")
        return None
    return d


def TNdistance(seq1, seq2):  # this doesn't work
    """
    Tajima-Nei distance = -b log(1 - p / b)
    where:
    b = 0.5 * [ 1 - Sum i from A to T(Gi^2+p^2/h) ]
    h = Sum i from A to G( Sum j from C to T (Xij^2/2*Gi*Gj))
    p = p-distance, i.e. uncorrected distance between seq1 and seq2
    Xij = frequency of pair (i,j) in seq1 and seq2, with gaps removed
    Gi = frequency of base i over seq1 and seq2"""
    from math import log

    ns = ["A", "C", "G", "T"]
    G = estimate_nucleotide_frequencies(seq1 + seq2)
    p = pdistance(seq1, seq2)
    pairs = []
    h = 0

    # collect ungapped pairs
    for x in zip(seq1, seq2):
        pairs.append(x)

    # pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns) - 1):
        for j in range(i + 1, len(ns)):
            if i != j:
                paircount = pairs.count((ns[i], ns[j])) + pairs.count((ns[j], ns[i]))
            Xij_sq = (float(paircount) / len(pairs)) ** 2
            GiGj = G[i] * G[j]
            h += 0.5 * Xij_sq / GiGj  # h value used to calculate b

    b = 0.5 * (1 - sum([x**2 for x in G]) + p**2 / h)
    try:
        d = -b * log(1 - p / b)
    except ValueError:
        print("Tried to take log of a negative number")
        return None
    return d


def K2Pdistance(seq1, seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt

    pairs = []
    seq1, seq2 = seq1.upper(), seq2.upper()
    # collect ungapped pairs
    for x in zip(seq1, seq2):
        pairs.append(x)

    ts_count = 0
    tv_count = 0
    length = len(pairs)

    transitions = ["AG", "GA", "CT", "TC"]
    transversions = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]

    for x, y in pairs:
        if x + y in transitions:
            ts_count += 1
        elif x + y in transversions:
            tv_count += 1

    p = float(ts_count) / length
    q = float(tv_count) / length
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return None
    return d


def Tamuradistance(seq1, seq2):
    """
    Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
    where:
    P = transition frequency
    Q = transversion frequency
    C = GC1 + GC2 - 2 * GC1 * GC2
    GC1 = GC-content of sequence 1
    GC2 = GC-coontent of sequence 2
    """
    from math import log

    pairs = []
    seq1, seq2 = seq1.upper(), seq2.upper()
    # collect ungapped pairs
    for x in zip(seq1, seq2):
        pairs.append(x)

    ts_count = 0
    tv_count = 0
    length = len(pairs)

    transitions = ["AG", "GA", "CT", "TC"]
    transversions = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]

    for x, y in pairs:
        if x + y in transitions:
            ts_count += 1
        elif x + y in transversions:
            tv_count += 1

    p = float(ts_count) / length
    q = float(tv_count) / length
    gc1 = sum(estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(estimate_nucleotide_frequencies(seq2)[1:3])
    c = gc1 + gc2 - 2 * gc1 * gc2

    try:
        d = -c * log(1 - p / c - q) - 0.5 * (1 - c) * log(1 - 2 * q)
    except ValueError:
        print("Tried to take log of a negative number")
        return None
    return d


# refer: https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py


def distancearray(dictG, fn):
    out = {}
    dictKeys = list(dictG.keys())
    for i in dictG:
        out[i] = []
    for key in dictKeys:
        for j in range(len(dictKeys)):
            out[key].append(fn(dictG[key], dictG[dictKeys[j]]))
    print("dones")
    fnnames = {
        "Tamuradistance": "TamuraArray",
        "levendistance": "LevenArray",
        "pdistance": "RawArray",
        "JCdistance": "JCArray",
        "K2Pdistance": "K2PArray",
    }
    fnname = fnnames[str(fn.__name__)]
    os.makedirs(PATH + "/distancearrays", exist_ok=True)
    with open(PATH + "/distancearrays/" + str(fnname) + ".txt", "w") as f:
        f.write(str(out))
    return out


fns = [Tamuradistance, levendistance, pdistance, JCdistance, K2Pdistance]
for k in fns:
    # if you only want tamura distance
    if k == Tamuradistance:
        distancearray(dictSamples, k)
        
# generate output file for smk
open("res/1distanceMatrixFormation.done", "x")


#################################
#     Phylogenetic Distance     #
#################################
# finds branch distance between nodes a and b
def branchdistance(a, b):
    if a == b:
        return 0
    ancestora = a.ancestor
    ancestoralst = [ancestora]
    # print(ancestora, b)
    # if b in ancestora.descendants:
    #     return a.length + b.length
    while b not in ancestora.get_leaves():
        ancestora = ancestora.ancestor
        ancestoralst.append(ancestora)

    ancestorb = b.ancestor
    ancestorblst = [ancestorb]
    # print(ancestora, b)
    # if b in ancestora.descendants:
    #     return a.length + b.length
    while a not in ancestorb.get_leaves():
        ancestorb = ancestorb.ancestor
        ancestorblst.append(ancestorb)

    totaldistance = a.length + b.length
    if ancestoralst:
        ancestoralst.pop()
    for node in ancestoralst:
        totaldistance += node.length

    if ancestorblst:
        ancestorblst.pop()
    for node in ancestorblst:
        totaldistance += node.length
    return round(totaldistance, 10)


# takes in src, which is the ABSOLUTE path to the tree
# outputs phylogenetic distance matrix
def treeDistMatrix(src):
    import newick, os

    with open(src) as f1:
        tree = newick.load(f1)[0]

    keys = tree.get_leaves()
    out = {}
    for i in keys:
        temp = []
        out[i] = temp
        for j in keys:
            # branchdistance is distance by branch lengths between i and j
            temp.append(branchdistance(i, j))
        # print(i)
    # print(out)
    PATH = os.path.dirname(src)
    with open(PATH + "/topdownDistMatrix.txt", "w") as f:
        f.write(str(out))
    return out


# example:
# treepath = 'dengue/K1.raxml.bestTree'
# treeDistMatrix(treepath)
