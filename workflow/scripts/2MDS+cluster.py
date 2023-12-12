from ast import literal_eval
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import argparse
from sklearn import manifold  # multidimensional scaling
from sklearn.cluster import OPTICS, KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import fcluster, linkage

parser = argparse.ArgumentParser()
parser.add_argument("source", help="The path of the directory with the distance array")
parser.add_argument("--useMDS", help="Whether or not to use multi-dimentional scaling", default= "noMDS")
args = parser.parse_args()

os.makedirs(args.source + "/clusterlabels", exist_ok=True)
os.makedirs(args.source + "/silhouette", exist_ok=True)

# MDS from distance matrix
# src is ABSOLUTE path of the distance matrix (eg of Tamura array)
# please make sure that there is a clusterlabels and MDSplots in the PATH folder

for jj in range(4, 5):  # modify numbers according to which distances you want
    # eg. if you only want K2P, then put range(3,4)
    src = args.source
    mds = (args.useMDS == "MDS")
    arraypath = ["Leven", "Raw", "JC", "K2P", "Tamura", "Phylo"][jj]
    distanceMeasure = ["Leven", "Raw", "JC", "K2P", "Tamura", "Phylo"][jj]
    os.makedirs(src + "/clusterlabels/" + distanceMeasure, exist_ok=True)
    os.makedirs(src + "/silhouette/" + distanceMeasure, exist_ok=True)
    src = src + "/" + arraypath + "Array.txt"

    with open(src) as f:
        data = f.read()
    d = literal_eval(data)  # d is the distance matrix
    PATH = os.path.dirname(src)
    vals = np.array(list(d.values()))
    labels = d.keys()
    # Performs MDS onto n_components dimensions
    mds_model = manifold.MDS(
        n_components=2, random_state=0, dissimilarity="precomputed"
    )
    mds_fit = mds_model.fit(vals)
    mds_coords = mds_model.fit_transform(vals)

    # sets colours for the clusterings below
    # clist = ["#D81B60", "#1E88E5", "#FFC107", "#004D40"]
    clist = ["#D81B60", "#1E88E5", "#FFC107", "#8B4000", "#B19CD8", "#004D40"]
    mycmap = cm.ListedColormap(clist, name="custom_cmap")

    def filterMDScoords(lst):  # outputs MDS coordinates of all/filtered points
        coords_dict = {}
        for label, x, y in zip(labels, mds_coords[:, 0], mds_coords[:, 1]):
            if label in lst:
                coords_dict[label] = (x, y)
        outcoords = PATH + "/filteredMDScoords.txt"
        with open(outcoords, "w") as f1:
            f1.write(str(coords_dict))
        return

    # Function calculates distance between two points
    def dist(p1, p2):
        x0 = p1[0] - p2[0]
        y0 = p1[1] - p2[1]
        return x0 * x0 + y0 * y0

    # Function to find the maximum distance between any two points on MDS.
    # Takes in output from filterMDScoords
    def maxDist(src):
        with open(src) as f:
            data = f.read()
        p = literal_eval(data)  # p is the coordinates dictionary
        points = list(p.values())
        n = len(points)
        print(n)
        maxm = 0
        # Iterate over all possible pairs
        for i in range(n):
            for j in range(i + 1, n):
                # Update maxm
                maxm = max(maxm, dist(points[i], points[j]))
        from math import sqrt

        # Return actual distance
        return sqrt(maxm)

    # Function that plots MDS, colouring them based on pre-defined clusters
    # Takes in lst of labels. Eg. backplot(labels)
    def backplot(lst=labels):
        # lst is the list of all samples labels. Eg:
        # lst = ['OV099745', 'OV099735', 'OV099674', 'OV099560', 'OV099489', 'OV099461', 'OV099191', 'OV097299', 'OV089026', 'OV072861', 'OV072859', 'OV057511', 'OV057131', 'OV056904', 'OV056601', 'OV053610', 'OV049178', 'OV048808', 'OV070343', 'OV067676', 'OV028184', 'OV637479', 'OV635995', 'OV635135', 'OV635070', 'OV592060', 'OV587075', 'OV551194', 'OV551371', 'OV550356', 'OV550329', 'OV550656', 'OV549761', 'OV546677', 'OV536668', 'OV570666', 'OV570604', 'OV562129', 'OV572664', 'OV552333', 'OV552328', 'OV550991', 'OV550600', 'OV514006', 'OV508297', 'OV506314', 'OV504985', 'OV489961', 'OV486675', 'OV486442', 'OV485436', 'OV483284', 'OV477890', 'OV477359', 'OV505058', 'OV498414', 'OV489728', 'OV424790', 'OV467035', 'OV466951', 'OV462420', 'OV462361', 'OV462228', 'OV460403', 'OV413903', 'OV458339', 'OV457202', 'OV454176', 'OV425093', 'OV435145', 'OV435048', 'OV433362', 'OV466210', 'OV454858', 'OV454443', 'OV454022', 'OV449046', 'OV432616', 'OV426449', 'OV371283', 'OV370982', 'OV395815', 'OV093011', 'OV048775', 'OV048181', 'OV044774', 'OV065431', 'OV064961', 'OV011900', 'OV030712', 'OV029243', 'OV027233', 'OV027166', 'OV026607', 'OV025329', 'OV009389', 'OV008565', 'OV002576', 'OV002530', 'OU983369', 'OV343699', 'OV343562', 'OV343347', 'OV343321', 'OV341379', 'OV341339', 'OV337395', 'OV336288', 'OV330378', 'OV328759', 'OV328554', 'OV365663', 'OV364686', 'OV364625', 'OV404680', 'OV404500', 'OV404077', 'OV400243', 'OV346414', 'OV345069', 'OV344182', 'OV344031', 'OV300520', 'OV283801', 'OV283777', 'OV283090', 'OV281715', 'OV280088', 'OV276366', 'OV274308', 'OV273905', 'OV272478', 'OV284143', 'OV260161', 'OV259968', 'OV259217', 'OV258808', 'OV258716', 'OV258619', 'OV258384', 'OV257417', 'OV257229', 'OV257183', 'OV257132', 'OV256815', 'OV255909', 'OV254168', 'OV250885', 'OV165399', 'OV247551', 'OV247347', 'OV163560', 'OV246872', 'OV246811', 'OV246746', 'OV246648', 'OV163107', 'OV246390', 'OV162928', 'OV245954', 'OV162593', 'OV162205', 'OV161557', 'OV243894', 'OV236619', 'OV235596', 'OV235545', 'OV234825', 'OV184026', 'OV184022', 'OV266133', 'OV182212', 'OV264516', 'OV264503', 'OV180487', 'OV180047', 'OV260645', 'OV137486', 'OV160654', 'OV156845', 'OV156305', 'OV155979', 'OV150314', 'OV148830', 'OV146009', 'OV145890', 'OV135697', 'OV135648', 'OV133066', 'OV132856', 'OV132313', 'OV132156', 'OV100635', 'OV100321', 'OV246790', 'OV161504', 'OV161414', 'OV206054', 'OV183038', 'OV180816', 'OV137488', 'OV160644', 'OV160572', 'OV143628', 'OV132530', 'OV100278', 'OV100277', 'OV099897', 'OV099598', 'OV098835', 'OV072615', 'OV058614', 'OV056667', 'OV049121', 'OV067654', 'OV041821', 'OV637474', 'OV635517', 'OV634725', 'OV640981', 'OV597273', 'OV596364', 'OV593099', 'OV591436', 'OV591308', 'OV595083', 'OV543347', 'OV559737', 'OV551966', 'OV552266', 'OV537746', 'OV554550', 'OV572171', 'OV568332', 'OV570704', 'OV550001', 'OV545494', 'OV542731', 'OV514472', 'OV517987', 'OV517847', 'OV508496', 'OV508383', 'OV497113', 'OV486629', 'OV486432', 'OV432092', 'OV474337', 'OV474267', 'OV466338', 'OV465846', 'OV461460', 'OV461442', 'OV456140', 'OV454702', 'OV454468', 'OV424743', 'OV449428', 'OV448418', 'OV433115', 'OV437259', 'OV467671', 'OV462848', 'OV438263', 'OV432675', 'OV384477', 'OV384472', 'OV384238', 'OV383577', 'OV383558', 'OV383421', 'OV329574', 'OV370148', 'OV364248', 'OV322117', 'OV389360', 'OV315811', 'OU487169', 'OU474865', 'OU443123', 'OU439599', 'OU435062', 'OU434980', 'OU434946', 'OU658401', 'OU658209', 'OU654930', 'OU597700', 'OU638229', 'OU628594', 'OU627959', 'OU628105', 'OU626992', 'OU626101', 'OU626065', 'OU611034', 'OU607786', 'OU603737', 'OU588191', 'OU586114', 'OU565233', 'OU565061', 'OU563885', 'OU563698', 'OU562398', 'OU562184', 'OU562040', 'OU561979', 'OU561970', 'OU553242', 'OU836787', 'OU819425', 'OU818925', 'OU818599', 'OU781380', 'OU778109', 'OU778098', 'OU778027', 'OU774564', 'OU754561', 'OU750483', 'OU732726', 'OU732375', 'OU730556', 'OU730109', 'OU720925', 'OU720174', 'OU719279', 'OU711333', 'OU709191', 'OU679188', 'OU688270', 'OU683458', 'OV047330', 'OV040790', 'OV025779', 'OV009123', 'OV000075', 'OU997825', 'OU951706', 'OU951109', 'OU918950', 'OU914387', 'OU914058', 'OU913485', 'OU893905', 'OU892938', 'OU892554', 'OU892139', 'OU891580', 'OU856516', 'OU884932', 'OU884728', 'OU884233', 'OU884214', 'OU884129', 'OU883121', 'OU882799', 'OU881476', 'OU838574', 'OV331651', 'OV331630', 'OV331365', 'OV331047', 'OV405519', 'OV344024', 'OV275119', 'OV271852', 'OV257398', 'OV170058', 'OV168404', 'OV166264', 'OV165653', 'OV505024', 'OV496410', 'OV487754', 'OV487655', 'OV487449', 'OV487029', 'OV486970', 'OV485718', 'OV485132', 'OV484903', 'OV483994', 'OV483646', 'OV483588', 'OV482532', 'OV480623', 'OV489556', 'OV437573', 'OV433553', 'OV470047', 'OV467643', 'OV420592', 'OV463328', 'OV461338', 'OV414750', 'OV457413', 'OV453898', 'OV453848', 'OV453259', 'OV446481', 'OV440212', 'OV438871', 'OV437816', 'OV437141', 'OV437534', 'OV436157', 'OV432493', 'OV431166', 'OV429580', 'OV426893', 'OV466589', 'OV462336', 'OV449012', 'OV439234', 'OV439209', 'OV438231', 'OV434784', 'OV432516', 'OV384301', 'OV342824', 'OV342498', 'OV338810', 'OV329890', 'OV329839', 'OV407251', 'OV407225', 'OV394223', 'OV389904', 'OV389903', 'OV389273', 'OV345857', 'OV345186', 'OV308894', 'OV315819', 'OV315503', 'OV315470', 'OV639388', 'OV638818', 'OV638710', 'OV638520', 'OV638518', 'OV637692', 'OV637446', 'OV635800', 'OV635784', 'OV635782', 'OV635780', 'OV635146', 'OV635515', 'OV634410', 'OV634254', 'OV634521', 'OV633503', 'OV632836', 'OV669333', 'OV666151', 'OV666145', 'OV665552', 'OV665572', 'OV622436', 'OV619656', 'OV657247', 'OV657173', 'OV656438', 'OV656398', 'OV655643', 'OV652448', 'OV651794', 'OV650762', 'OV650719', 'OV650701', 'OV650429', 'OV649767', 'OV693781', 'OV649468', 'OV647432', 'OV645186', 'OV645006', 'OV643734', 'OV642700', 'OV641563', 'OV641176', 'OV640444', 'OV668561', 'OV667427', 'OV639745', 'OV639489', 'OV605149', 'OV603445', 'OV599633', 'OV599632', 'OV599354', 'OV598410', 'OV597680', 'OV596727', 'OV596705', 'OV596365', 'OV595876', 'OV594704', 'OV593447', 'OV593379', 'OV588044', 'OV583679', 'OV582019', 'OV600678', 'OV594759', 'OV560354', 'OV562933', 'OV560170', 'OV551240', 'OV559811', 'OV559702', 'OV552904', 'OV551598', 'OV550947', 'OV550921', 'OV551312', 'OV550375', 'OV550274', 'OV550225', 'OV547632', 'OV545614', 'OV536084', 'OV534486', 'OV566639', 'OV566420', 'OV561493', 'OV565350', 'OV563943', 'OV562090', 'OV561785', 'OV572191', 'OV571163', 'OV568300', 'OV566465', 'OV563752', 'OV552672', 'OV552056', 'OV551051', 'OV534684', 'OV534153', 'OV511185', 'OV517725', 'OV518512', 'OV497337', 'OV508726', 'OV508683', 'OV507744']
        with open("AY3.txt", "r") as f2:
            inp = f2.readlines()
            AY3 = [x.rstrip("\n") for x in inp]
        with open("AY122.txt", "r") as f3:
            inp = f3.readlines()
            AY122 = [x.rstrip("\n") for x in inp]
        # print(len(AY3), len(AY122))
        # print(lst)
        labels = []
        for i in lst:
            if i in AY3:
                labels.append(2)
            elif i in AY122:
                labels.append(1)
        with open("clusterlabels/grdtruth.txt", "w") as gf:
            gf.write(str(labels))
        # return labels
        plt.clf()
        plt.scatter(
            mds_coords[:, 0],
            mds_coords[:, 1],
            c=labels,
            cmap=mycmap,
            marker="o",
            picker=True,
        )
        plt.title("True clusters")
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        # plt.savefig(PATH+'/MDSplots/trueclusters.png', dpi=500)
        plt.show()

    # backplot()

    def OPTICScluster(size, useMDS=True):
        if useMDS:
            clustering = OPTICS(min_samples=size).fit(mds_coords)
        else:
            clustering = OPTICS(min_samples=size).fit(vals)
        lkeys, llabels = list(d.keys()), list(clustering.labels_)
        from csv import writer

        outcsv = (
            PATH + "/clusterlabels/" + distanceMeasure + "/OPTICS" + str(size) + ".csv"
        )
        if not useMDS:
            outcsv = (
                PATH
                + "/clusterlabels/"
                + distanceMeasure
                + "/noMDSOPTICS"
                + str(size)
                + ".csv"
            )
        with open(outcsv, "w") as f1:
            wr = writer(f1)
            for i in range(len(lkeys)):
                wr.writerow([lkeys[i], llabels[i]])
        if useMDS:
            plt.clf()
            plt.scatter(
                mds_coords[:, 0],
                mds_coords[:, 1],
                c=clustering.labels_,
                cmap=mycmap,
                marker="o",
                picker=True,
            )
            plt.title("OPTICS cluster, size = " + str(size))
            plt.savefig(
                PATH + "/MDSplots/" + distanceMeasure + "/OPTICS" + str(size) + ".png",
                dpi=500,
            )
        # plt.show()

    def Kmeanscluster(K, useMDS=True):
        if useMDS:
            clustering = KMeans(
                init="k-means++", n_clusters=K, n_init=10, random_state=42
            ).fit(mds_coords)
        else:
            clustering = KMeans(
                init="k-means++", n_clusters=K, n_init=10, random_state=42
            ).fit(vals)
        lkeys, llabels = list(d.keys()), list(clustering.labels_)
        from csv import writer

        outcsv = (
            PATH + "/clusterlabels/" + distanceMeasure + "/Kmeans" + str(K) + ".csv"
        )
        if not useMDS:
            outcsv = (
                PATH
                + "/clusterlabels/"
                + distanceMeasure
                + "/noMDSKmeans"
                + str(K)
                + ".csv"
            )
        with open(outcsv, "w") as f1:
            wr = writer(f1)
            for i in range(len(lkeys)):
                wr.writerow([lkeys[i], llabels[i]])
        if useMDS:
            plt.clf()
            plt.scatter(
                mds_coords[:, 0],
                mds_coords[:, 1],
                c=clustering.labels_,
                cmap=mycmap,
                marker="o",
                picker=True,
            )
            plt.title("K means cluster, K = " + str(K))
            plt.savefig(
                PATH + "/MDSplots/" + distanceMeasure + "/Kmeans" + str(K) + ".png",
                dpi=500,
            )
            # plt.show()

    def aggloCluster(distance, useMDS=True):
        if useMDS:
            clustering = AgglomerativeClustering(
                n_clusters=None, linkage="single", distance_threshold=distance
            ).fit(mds_coords)
        else:
            clustering = AgglomerativeClustering(
                n_clusters=None, linkage="single", distance_threshold=distance
            ).fit(vals)

        lkeys, llabels = list(d.keys()), list(clustering.labels_)
        from csv import writer

        outcsv = (
            PATH
            + "/clusterlabels/"
            + distanceMeasure
            + "/agglo"
            + str(distance)
            + ".csv"
        )
        if not useMDS:
            outcsv = (
                PATH
                + "/clusterlabels/"
                + distanceMeasure
                + "/noMDSagglo"
                + str(distance)
                + ".csv"
            )
        with open(outcsv, "w") as f1:
            wr = writer(f1)
            for i in range(len(lkeys)):
                wr.writerow([lkeys[i], llabels[i]])
        if useMDS:
            plt.clf()
            plt.scatter(
                mds_coords[:, 0],
                mds_coords[:, 1],
                c=clustering.labels_,
                cmap=mycmap,
                marker="o",
                picker=True,
            )
            plt.title("Agglomerative cluster, distance = " + str(distance))
            plt.savefig(
                PATH
                + "/MDSplots/"
                + distanceMeasure
                + "/agglomerative"
                + str(distance)
                + ".png",
                dpi=500,
            )

    """
    Example below. Note the clustering variables should be obtained from silhouette analysis code.

    # Calling for dengue dataset
    if jj == 3: #K2P
        print('K2P')

        Kmeanscluster(2,useMDS=False)
        Kmeanscluster(2,useMDS=True)
        aggloCluster(0.02875,useMDS=False)
        aggloCluster(0.0016,useMDS=True)
        OPTICScluster(0.07,useMDS=False)
        OPTICScluster(0.16,useMDS=True)

    if jj == 4: #Tamura
        print('Tamura')
        aggloCluster(0.0044,useMDS=False)
        aggloCluster(0.0026,useMDS=False)

        Kmeanscluster(2,useMDS=False)
        Kmeanscluster(2,useMDS=True)
        aggloCluster(0.02875,useMDS=False)
        aggloCluster(0.0016,useMDS=True)
        OPTICScluster(0.07,useMDS=False)
        OPTICScluster(0.16,useMDS=True)

    elif jj == 5: #phylo
        print('phylo')

        Kmeanscluster(2,useMDS=False)
        Kmeanscluster(3,useMDS=True)
        aggloCluster(0.20705,useMDS=False)
        aggloCluster(0.0412,useMDS=True)
        OPTICScluster(0.13,useMDS=False)
        OPTICScluster(0.12,useMDS=True)
    """

    # performs silhouette analysis
    # input which algorithm, and whether or not MDS is used
    def silhouette(alg, useMDS=True):
        if alg == "Kmeans":
            score, clusters = [], []
            for i in range(2, 31):
                model = KMeans(
                    init="k-means++", n_clusters=i, n_init=10, random_state=42
                )
                # visualizer = SilhouetteVisualizer(model, colors='yellowbrick')
                visualizer = model
                if useMDS:
                    clustering = visualizer.fit(mds_coords)
                    score.append(silhouette_score(mds_coords, clustering.labels_))
                    clusters.append(len(set(clustering.labels_)))
                else:
                    clustering = visualizer.fit(vals)
                    score.append(silhouette_score(vals, clustering.labels_))
                    clusters.append(len(set(clustering.labels_)))
            plt.clf()
            fig, ax = plt.subplots()
            ax.plot(list(range(2, 31)), score, label="Silhouette score", color="red")
            ax.set_xlabel("K value", fontsize=12)
            ax.set_ylabel("Silhouette score", color="red", fontsize=12)
            xmax = list(range(2, 31))[score.index(max(score))]
            ymax = max(score)
            ax.axvline(xmax, ls=":", c="k")
            plt.text(xmax, ymax * 1.0005, f"({xmax}, {round(ymax,3)})", fontsize=8)
            ax2 = ax.twinx()
            ax2.plot(
                list(range(2, 31)), clusters, label="Number of clusters", color="blue"
            )
            ax2.set_ylabel("Number of clusters", color="blue", fontsize=12)

            # plt.legend()
            if useMDS:
                plt.title("Silhouette scores for K means clustering with MDS")
                plt.savefig(
                    PATH + "/silhouette/" + distanceMeasure + "/Kmeanssilhouette.png",
                    dpi=500,
                )
            else:
                plt.title("Silhouette scores for K means clustering without MDS")
                plt.savefig(
                    PATH
                    + "/silhouette/"
                    + distanceMeasure
                    + "/noMDSKmeanssilhouette.png",
                    dpi=500,
                )
            best = list(range(2, 31))[score.index(max(score))]
            print(best)
            return best

        if alg == "OPTICS":
            score, clusters = [], []
            for i in [x / 100.0 for x in range(1, 100)]:
                model = OPTICS(min_samples=i)
                # visualizer = SilhouetteVisualizer(model, colors='yellowbrick')
                visualizer = model
                if useMDS:
                    clustering = visualizer.fit(mds_coords)
                    if len(set(clustering.labels_)) > 1:
                        score.append(silhouette_score(mds_coords, clustering.labels_))
                        clusters.append(len(set(clustering.labels_)))
                    else:
                        j = int(i * 100.0)
                        break
                else:
                    clustering = visualizer.fit(vals)
                    if len(set(clustering.labels_)) > 1:
                        score.append(silhouette_score(vals, clustering.labels_))
                        clusters.append(len(set(clustering.labels_)))
                    else:
                        j = int(i * 100.0)
                        break

            plt.clf()
            fig, ax = plt.subplots()
            ax.plot(
                [x / 100.0 for x in range(1, j)],
                score,
                label="Silhouette score",
                color="red",
            )
            ax.set_xlabel("OPTICS size", fontsize=12)
            ax.set_ylabel("Silhouette score", color="red", fontsize=12)
            xmax = [x / 100.0 for x in range(1, j)][score.index(max(score))]
            ymax = max(score)
            ax.axvline(xmax, ls=":", c="k")
            plt.text(xmax, ymax * 1.02, f"({xmax}, {round(ymax,3)})", fontsize=8)
            ax2 = ax.twinx()
            ax2.plot(
                [x / 100.0 for x in range(1, j)],
                clusters,
                label="Number of clusters",
                color="blue",
            )
            ax2.set_ylabel("Number of clusters", color="blue", fontsize=12)

            # plt.legend()
            if useMDS:
                plt.title("Silhouette scores for OPTICS clustering with MDS")
                plt.savefig(
                    PATH + "/silhouette/" + distanceMeasure + "/OPTICSsilhouette.png",
                    dpi=500,
                )
            else:
                plt.title("Silhouette scores for OPTICS clustering without MDS")
                plt.savefig(
                    PATH
                    + "/silhouette/"
                    + distanceMeasure
                    + "/noMDSOPTICSsilhouette.png",
                    dpi=500,
                )
            best = [x / 100.0 for x in range(1, j)][score.index(max(score))]
            print(best)
            return best

        if alg == "Agglomerative":
            score, clusters = [], []
            for i in [x / 20000.0 for x in range(1, 20001)]:
                model = AgglomerativeClustering(
                    n_clusters=None, linkage="single", distance_threshold=i
                )
                visualizer = model
                if useMDS:
                    clustering = visualizer.fit(mds_coords)
                    if len(set(clustering.labels_)) > 1:
                        score.append(silhouette_score(mds_coords, clustering.labels_))
                        clusters.append(len(set(clustering.labels_)))
                        # j = int(i * 20000)
                    else:
                        j = round(i * 20000.0)
                        break
                else:
                    clustering = visualizer.fit(vals)
                    if len(set(clustering.labels_)) > 1:
                        if False:
                            # if i > 0.00181:
                            j = int(i * 20000)
                            break
                        else:
                            score.append(silhouette_score(vals, clustering.labels_))
                            clusters.append(len(set(clustering.labels_)))
                            j = int(i * 20000)
                    else:
                        j = round(i * 20000.0)
                        break
            plt.clf()
            fig, ax = plt.subplots()
            plt.plot(
                [x / 20000.0 for x in range(1, j)],
                score,
                label="Silhouette score",
                color="red",
            )
            ax.set_xlabel("Agglomerative distance", fontsize=12)
            ax.set_ylabel("Silhouette score", color="red", fontsize=12)
            xmax = [x / 20000.0 for x in range(1, j)][score.index(max(score))]
            ymax = max(score)
            ax.axvline(xmax, ls=":", c="k")
            plt.text(xmax, ymax * 1.01, f"({xmax}, {round(ymax,3)})", fontsize=8)
            ax2 = ax.twinx()
            plt.plot(
                [x / 20000.0 for x in range(1, j)],
                clusters,
                label="Number of clusters",
                color="blue",
            )
            ax2.set_ylabel("Number of clusters", color="blue", fontsize=12)

            # plt.legend()
            if useMDS:
                plt.title("Silhouette scores for Agglomerative clustering with MDS")
                plt.savefig(
                    PATH
                    + "/silhouette/"
                    + distanceMeasure
                    + "/Agglomerativesilhouette.png",
                    dpi=500,
                )
            else:
                plt.title("Silhouette scores for Agglomerative clustering without MDS")
                plt.savefig(
                    PATH
                    + "/silhouette/"
                    + distanceMeasure
                    + "/noMDSAgglomerativesilhouette.png",
                    dpi=500,
                )
            best = [x / 20000.0 for x in range(1, j)][score.index(max(score))]
            print(best)
            return best

    # Calling of silhouette analysis
    for i in ["Kmeans", "OPTICS", "Agglomerative"]:
        if i == "Agglomerative":
            bestdist = silhouette(i, useMDS=mds)
            aggloCluster(bestdist, useMDS=mds)

    # generate output file for smk
    # open("res/2MDS+cluster.done", "x")

    # Generating cluster based on best dstance from silhouette score
