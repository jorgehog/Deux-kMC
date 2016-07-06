import sys
import os
import shutil
from os.path import join
import numpy as np
from scipy.stats import linregress
from numpy import where, empty, save, zeros
from skimage import measure

from matplotlib import pylab

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

min_area = 4

def periodic_image_analysis(pim, labels):
    L, W = pim.shape
    expanded = zeros(shape=[2*d for d in pim.shape], dtype=int)
    for x in range(2):
        for y in range(2):
            expanded[x*L:(x+1)*L, y*W:(y+1)*W] = pim

    n_labels = len(labels)

    idx_map = {}

    n = 0
    for l in labels:
        idx_map[l] = n
        n += 1

    winners = dict(zip(labels, [[] for _ in range(n_labels)]))
    winner_areas = [0 for _ in range(n_labels)]

    l = measure.label(expanded)
    l[np.where(expanded == 0)] = -1
    l += 1

    def find_fast(im, v, x0, y0, x1, y1):
        for _x in range(x0, x1):
            for _y in range(y0, y1):
                if im[_x, _y] == v:
                    return _x, _y

    orig_labels = []
    expanded_props = measure.regionprops(l)
    for prop in expanded_props:

        x, y = find_fast(l, prop.label, *prop.bbox)
        orig_label = expanded[x, y]
        orig_labels.append(orig_label)

        idx = idx_map[orig_label]

        if prop.area > winner_areas[idx]:
            winner_areas[idx] = prop.area

    for orig_label, prop in zip(orig_labels, expanded_props):

        idx = idx_map[orig_label]

        if prop.area == winner_areas[idx]:
            winners[orig_label].append([prop.perimeter, prop.centroid])

    return winners

def changefill(im, source, dest):
    if dest == -1 or source == -1:
        return

    if dest != source:
        I = where(im == dest)
        im[I] = source


def periodify(im):
    l, w = im.shape

    for i in range(l):
        source = im[i, 0]
        dest = im[i, w-1]

        changefill(im, source, dest)

    for j in range(w):
        source = im[0, j]
        dest = im[l-1, j]

        changefill(im, source, dest)

    return im


def analyze(data,
            do_cluster_analysis,
            i,
            covs,
            n_broken,
            n_gained,
            *cluster_measures):

    _s = data[i, :, :]

    covs[i] = _s.mean()

    if i == data.shape[0] - 1:
        return

    _s_next = data[i+1, :, :]

    delta = _s_next - _s

    n_broken[i] = len(where(delta == -1)[0])/float(delta.size)
    n_gained[i] = len(where(delta == 1)[0])/float(delta.size)

    if not do_cluster_analysis:
        return

    avg_cluster_sizes, n_clusters_list, avg_circs, centeroids = cluster_measures

    im_label = measure.label(_s)
    im_label[where(_s == 0)] = -1

    im_label_per = periodify(im_label) + 1

    props = measure.regionprops(im_label_per)

    all_areas = []
    all_circs = []
    centeroids.append([])

    if props:
        perimeter_data = periodic_image_analysis(im_label_per, [prop.label for prop in props])

        for m in props:
            if m.label != 0 and m.area > min_area:
                all_areas.append(m.area)
                all_circs.append(perimeter_data[m.label][0][0])

                all_images = []
                for image in perimeter_data[m.label]:
                    all_images.append(image[1])

                centeroids[-1].append([m.area, all_images])

    n_clusters = len(all_areas)

    if n_clusters == 0:
        avg_cluster_size = 0
        avg_circ = 0
    else:
        avg_cluster_size = sum(all_areas)/len(all_areas)
        avg_circ = sum(all_circs)/len(all_circs)

    avg_cluster_sizes[i] = avg_cluster_size
    n_clusters_list[i] = n_clusters
    avg_circs[i] = avg_circ

colors = ['r', 'g', 'b', 'k', 'm', 'y']

def calculate_cluster_trace(data, L, W):

    max_cluster_size = max(data, key=lambda x: 0 if x == [] else max(x, key=lambda y: y[0]))

    if not max_cluster_size:
        return
    else:
        max_cluster_size = max_cluster_size[0][0]

    i = len(data) - 1
    while not data[i]:
        i -= 1

        if i == -1:
            return

    last_centeroid = data[i][0][1]

    xc, yc = last_centeroid[0]

    x0 = xc - L/2
    x1 = xc + L/2
    y0 = yc - W/2
    y1 = yc + W/2

    all_points = []

    for zi, c in enumerate(data):
        for ci, cluster in enumerate(c):
            size = cluster[0]
            ri = float(size)/max_cluster_size

            for xi, yi in cluster[1]:

                if (xi < x0) or (xi >= x1):
                    pass
                elif (yi < y0) or (yi >= y1):
                    pass
                else:
                    all_points.append([xi-x0, yi-y0, zi, ri, ci])

    return all_points


def makeXYZ_single(data, xyz_dir, n):
    XYZ = ""
    L, W = data.shape

    h_min = data.min()

    count = 0
    for x in range(L):
        for y in range(W):
            for z in range(h_min, data[x, y]+1):
                XYZ += "0 %g %g %g\n" % (x, y, z)
                count += 1

    name = "%s/cluster%d.xyz" % (xyz_dir, n)
    with open(name, 'w') as f:
        f.write("%d\n-\n%s" % (count, XYZ))

    return name


def makeXYZ(heights, idx, dir):

    tot = 10
    N = len(heights)
    every = N/tot

    xyz_dir = "%s/XYZ/" % dir

    if not os.path.exists(xyz_dir):
        os.mkdir(xyz_dir)

    n = 0
    for step in idx[::every]:
        data = heights[step][()]
        name = makeXYZ_single(data, xyz_dir, n)
        n += 1

        sys.stdout.flush()
        print "\rStoring %d/%d" % (n, len(heights)/every),
    print

    return name

def filter_background(n):

    return n

    x = np.arange(len(n))
    a, b, _, _, _ = linregress(x, n)

    return n - (a*x + b)

def calcHeight(all_coverage, all_heights, idx, every):

    neq = (3*len(all_coverage))/4

    mean = 0
    for i in range(neq, len(all_coverage), every):
        coverage = all_coverage[i]
        heights = all_heights[idx[i]][()]

        if coverage.sum() == 0:
            mean += 0
        else:
           # print heights[np.where(coverage == 0)].mean(), heights[np.where(coverage == 1)].mean()
            mean += heights[np.where(coverage == 1)].mean() - heights[np.where(coverage == 0)].mean()

    return mean / ((len(all_coverage)-neq)/every)

def get_cov_std(covs, nbroken, ngained, time):

    neq = (3*len(covs))/4

    #pylab.plot(np.array(covs[neq:]))
    #pylab.show()

    cov = covs[neq:].mean()
    stdcov = np.array(covs[neq:]).std()
    stdbroken = np.array(nbroken[neq:]).std()
    stdgained = np.array(ngained[neq:]).std()

    return stdcov, stdbroken, stdgained, cov

def main():

    input_file = sys.argv[1]

    single = sys.argv[2] == "single"
    if single:
        s_alpha = float(sys.argv[3])
        s_F0 = float(sys.argv[4])
    else:
        path = sys.argv[2]
        every = int(sys.argv[3])
        every2 = int(sys.argv[4])

    do_cluster_analysis = len(sys.argv) <= 5
    print "cluster analysis:", do_cluster_analysis

    parser = ParseKMCHDF5(input_file)

    alphas = []
    F0s = []
    lmax = 0
    N = 0
    for data, L, W, run_id in parser:
        l = len(data["eq_storedEventValues"][0])

        if l > lmax:
            lmax = l

        alpha = round(data.attrs["alpha"], 5)
        F0 = round(data.attrs["F0"], 5)

        if alpha not in alphas:
            alphas.append(alpha)
        if F0 not in F0s:
            F0s.append(F0)

        if single:
            if alpha == s_alpha and F0 == s_F0:
                time = data["eq_storedEventValues"][1][()]

                covs = np.zeros(len(time))
                nbroken = np.zeros(len(time) - 1)
                ngained = np.zeros(len(time) - 1)

                coverage_matrix = data["eq_coverage_matrix"]

                for i in range(len(time)):
                    analyze(coverage_matrix, do_cluster_analysis, i, covs, nbroken, ngained)

                save("/tmp/extran_cluster_time.npy", time)
                save("/tmp/extran_cluster_covs.npy", covs)
                save("/tmp/extran_cluster_nbroken.npy", nbroken)
                save("/tmp/extran_cluster_ngained.npy", ngained)

                return
        N += 1

    print "lmax:", lmax

    alphas = sorted(alphas)
    F0s = sorted(F0s)

    if single:
        print "alpha=%g and F0=%g not found in set." % (s_alpha, s_F0)
        for alpha in alpha:
            print alpha
        for F0 in F0s:
            print F0
        return

    print "Found %d entries" % N

    stds = np.zeros(shape=(len(alphas), len(F0s), 3))
    all_covs = np.zeros(shape=(len(alphas), len(F0s)))
    heigts = np.zeros_like(all_covs)

    if not os.path.exists(path):
        os.mkdir(path)

    all_xyz_path = os.path.join(path, "all_xyz")
    if not os.path.exists(all_xyz_path):
        os.mkdir(all_xyz_path)

    std_path = os.path.join(path, "stds")
    if not os.path.exists(std_path):
        os.mkdir(std_path)

    all_counts = {}
    last_xyz_txt = ""
    stored_heights_keys = None
    for count, (data, L, W, run_id) in enumerate(parser):

        alpha = round(data.attrs["alpha"], 5)
        ai = alphas.index(alpha)

        F0 = round(data.attrs["F0"], 5)
        F0i = F0s.index(F0)

        dir_tag = "a%.3fF0%.3f" % (alpha, F0)

        if dir_tag in all_counts.keys():
            all_counts[dir_tag] += 1
            dir_tag += "_%d" % all_counts[dir_tag]
        else:
            all_counts[dir_tag] = 0

        dir = "%s/cluster_%s" % (path, dir_tag)
        xyz_name = os.path.join(all_xyz_path, "all_xyz%d.xyz" % count)

        print "%d/%d" % (count+1, N), L, W, run_id, alpha, F0

        coverage_matrix_h5 = data["eq_coverage_matrix"]

        event_values_h5 = data["eq_storedEventValues"]
        stored_heights = data["stored_heights"]

        coverage_matrix = np.zeros(shape=(len(coverage_matrix_h5)/every, L, W))
        time = np.zeros(shape=len(coverage_matrix_h5)/every)

        if len(event_values_h5[1]) != lmax:
            last_xyz = makeXYZ_single(coverage_matrix_h5[-1], "/tmp", 0)
            cov = 0
        else:
            if not stored_heights_keys:
                stored_heights_keys = sorted(stored_heights.keys(), key=lambda x: int(x))
            heigts[ai, F0i] = calcHeight(coverage_matrix_h5, stored_heights, stored_heights_keys, every)
            #continue

            if not os.path.exists(dir):
                os.mkdir(dir)
            for i, n in enumerate(range(0, len(coverage_matrix_h5), every)):
                coverage_matrix[i] = coverage_matrix_h5[n]
                time[i] = event_values_h5[1][n]

            n = coverage_matrix.shape[0]

            covs = empty(n)
            avg_cluster_sizes = empty(n)
            n_clusters_list = empty(n)
            avg_circs = empty(n)
            centeroids = []
            n_broken = empty(n - 1)
            n_gained = empty(n - 1)

            measures = [covs,
                        n_broken,
                        n_gained,
                        avg_cluster_sizes,
                        n_clusters_list,
                        avg_circs,
                        centeroids]


            for i in range(n):
                analyze(coverage_matrix, do_cluster_analysis, i, *measures)

            stdcov, stdbroken, stdgained, cov = get_cov_std(covs,
                                                            n_broken,
                                                            n_gained,
                                                            time)

            stds[ai, F0i, :] = [stdcov, stdbroken, stdgained]
            all_covs[ai, F0i] = cov

            if do_cluster_analysis:
                trace = calculate_cluster_trace(centeroids[::(every2/every)], L, W)
            else:
                trace = []

            last_xyz = makeXYZ(stored_heights, stored_heights_keys, dir)

            save("%s/extran_cluster_time.npy" % dir, time)
            save("%s/extran_cluster_covs.npy" % dir, covs)
            save("%s/extran_cluster_nbroken.npy" % dir, n_broken)
            save("%s/extran_cluster_ngained.npy" % dir, n_gained)
            if do_cluster_analysis:
                save("%s/extran_cluster_trace.npy" % dir, trace)
            save("%s/extran_cluster_size.npy" % dir, avg_cluster_sizes)
            save("%s/extran_cluster_n.npy" % dir, n_clusters_list)
            save("%s/extran_cluster_circs.npy" % dir, avg_circs)

        shutil.copy(last_xyz, xyz_name)
        last_xyz_txt += "%d %.3f %.3f %g\n" % (count, alpha, F0, cov)

    # for ia, alpha in enumerate(alphas):
    #     K = np.where(heigts[ia, :] != 0)
    #     pylab.plot(np.array(F0s)[K], heigts[ia, :][K], "s", label="%g" % alpha)
    # pylab.ylim(0, 10)
    # pylab.legend()
    # pylab.show()

    np.save("%s/cov_clusters_alphas.npy" % std_path, alphas)
    np.save("%s/cov_clusters_F0.npy" % std_path, F0s)
    np.save("%s/cov_clusters_stds.npy" % std_path, stds)
    np.save("%s/cov_clusters_covs.npy" % std_path, all_covs)
    np.save("%s/cov_clusters_heights.npy" % std_path, heigts)

    with open(os.path.join(all_xyz_path, "desc.txt"), 'w') as f:
        f.write(last_xyz_txt)

if __name__ == "__main__":
    main()