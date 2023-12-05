from sklearn.metrics import (
    adjusted_rand_score,
    rand_score,
    adjusted_mutual_info_score,
    v_measure_score,
)
import csv
from collections import defaultdict
import matplotlib as plt
import numpy as np


def readColumn(f):  # reads the partitioning from csv file
    columns = defaultdict(list)
    with open(f) as file:
        reader = csv.reader(file)
        for row in reader:
            for i, v in enumerate(row):
                columns[i].append(v)
            # print(i,v)
    return list(columns[1])


def comparison(f1, f2):  # gives ARI of partitioning from f1 and f2
    return adjusted_rand_score(readColumn(f1), readColumn(f2))


# Please do this for phylo distance methods first before running ARI
# src is folder w Phylo array
def correctPhylo(src):
    correctrows = []
    # any other cluster label of correct genetic distance, example:
    with open("clusterlabels/Tamura/Kmeans2.csv", "r") as f1:
        csv_reader = csv.reader(f1)
        for row in csv_reader:
            correctrows.append([row[0]])
    # print(correctrows)
    import os

    PATH = os.path.dirname(src)
    # print(correctrows)
    with open(src, "r") as f:
        csv_reader1 = csv.reader(f)
        for row in csv_reader1:
            for j in range(len(correctrows)):
                if row[0] == correctrows[j][0]:
                    correctrows[j].append(row[1])
    out = open(PATH + "/_" + os.path.basename(src), "w")
    with out:
        write = csv.writer(out)
        for entry in correctrows:
            write.writerow(entry)


# src1 = "clusterlabels/Phylo/"
# import os
# for filename in os.listdir(src1):
#     if filename[0] == '.':
#         continue
#     file = src1 + filename
#     print(file)
#     correctPhylo(file)


# input ABSOLUTE Path of grd truth
grdtruth = "/clusterlabels/grdtruth2.csv"


# src is ABOLUTE PATH where distance matrices are, and where the ARI.csv output file will be written to.
def compTable(src):
    import os

    out = []
    for jj in range(4, 5):  # modify numbers according to which distances you want
        # eg. if you only want K2P, then put range(3,4)
        arraypath = ["Leven", "Raw", "JC", "K2P", "Tamura", "Phylo"][jj]
        src1 = src + arraypath + "/"
        tmp = []
        for filename in sorted(os.listdir(src1)):
            if filename[0] == "." or filename[-1] != "v":  # removes random stuff
                continue
            src2 = src1 + filename
            print(src2)
            tmp.append(comparison(grdtruth, src2))
        out.append(tmp)
    file = open(src + "/ARI.csv", "w+")
    with file:
        write = csv.writer(file)
        for entry in out:
            write.writerow(entry)
    # print(out)


compTable("/clusterlabels/")
# if you only want a single comparison between foo.csv and bar.csv:
# print(comparison("foo.csv","bar.csv"))
