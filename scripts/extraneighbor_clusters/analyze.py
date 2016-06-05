import sys
import os
from os.path import join
import numpy as np
from numpy import where, empty, save, zeros
from skimage import measure
from threading import Thread
from matplotlib import pylab as plab
from mpl_toolkits.mplot3d import Axes3D

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


def analyze(data, i,
            covs,
            avg_cluster_sizes,
            n_clusters_list,
            avg_circs,
            n_broken,
            n_gained,
            centeroids):

    _s = data[i, :, :]

    covs[i] = _s.mean()

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

    if i == data.shape[0] - 1:
        return

    _s_next = data[i+1, :, :]

    delta = _s_next - _s

    n_broken[i] = len(where(delta == -1)[0])
    n_gained[i] = len(where(delta == 1)[0])


class AnalyzeThread(Thread):

    def __init__(self, data, thread_i, *m):

        self.data = data
        self.thread_i = thread_i

        self.m = m

        super(AnalyzeThread, self).__init__()

    def run(self):

        for i in self.thread_i:
            analyze(self.data, i, *self.m)

def calculate_cluster_trace(data, L, W):
    max_size = 20
    min_size = 5

    proj_props = {
        "s" : max_size,
        "marker" : "x",
        "c" : "k",
        "linewidth" : 1
    }

    max_cluster_size = max(data, key=lambda x: 0 if x == [] else max(x, key=lambda y: y[0]))[0][0]

    colors = ['r', 'g', 'b', 'k', 'm', 'y']

    fig = plab.figure()
    ax = fig.gca(projection='3d')

    xc, yc = data[-1][0][1][0]

    x0 = xc - L/2
    x1 = xc + L/2
    y0 = yc - W/2
    y1 = yc + W/2

    all_points = []

    for zi, c in enumerate(data):
        for ci, cluster in enumerate(c):
            size = cluster[0]
            ri = float(size)/max_cluster_size
            si = min_size + ri*(max_size-min_size)

            for xi, yi in cluster[1]:

                if (xi < x0) or (xi >= x1):
                    pass
                elif (yi < y0) or (yi >= y1):
                    pass
                else:
                    all_points.append([xi-x0, yi-y0, zi, ri, ci])

                ax.scatter(zi, xi, yi, s=si, c=colors[ci%len(colors)], edgecolors='none')

                if zi % 10 == 0:
                    #ax.scatter(len(data), xi, yi, s=sproj, c='k', edgecolors='none')
                    # ax.scatter(zi, x1, yi, **proj_props)
                    ax.scatter(zi, xi, y0, **proj_props)

    ax.set_xbound(0)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    pad = 1
    ax.xaxis._axinfo['label']['space_factor'] = pad
    ax.yaxis._axinfo['label']['space_factor'] = pad
    ax.zaxis._axinfo['label']['space_factor'] = pad

    fs = 30
    ax.set_xlabel(r"$\nu t$", size=fs)
    ax.set_ylabel(r"$x$", size=fs)
    ax.set_zlabel(r"$y$", size=fs)

    plab.show()

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

    with open("%s/cluster%d.xyz" % (xyz_dir, n), 'w') as f:
        f.write("%d\n-\n%s" % (count, XYZ))


def makeXYZ(heights, dir, every):

    xyz_dir = "%s/XYZ/" % dir

    if not os.path.exists(xyz_dir):
        os.mkdir(xyz_dir)

    n = 0
    for step in sorted(heights.keys(), key=lambda x: int(x))[::every]:
        data = heights[step][()]
        makeXYZ_single(data, xyz_dir, n)
        n += 1

        sys.stdout.flush()
        print "\rStoring %d/%d" % (n, len(heights)/every),
    print

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    nThreads = int(sys.argv[2])

    if len(sys.argv) > 3:
        every = int(sys.argv[3])
    else:
        every = 1

    for data, L, W, run_id in parser:

        dir = "/tmp/cluster_%s" % run_id

        if not os.path.exists(dir):
            os.mkdir(dir)

        print L, W, run_id, data.attrs["F0"]

        coverage_matrix_h5 = data["eq_coverage_matrix"]
        event_values_h5 = data["eq_storedEventValues"]
        coverage_matrix = np.zeros(shape=(len(coverage_matrix_h5)/every, L, W))
        time = np.zeros(shape=len(coverage_matrix_h5)/every)

        makeXYZ(data["stored_heights"], dir, 10*every)

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
                    avg_cluster_sizes,
                    n_clusters_list,
                    avg_circs,
                    n_broken,
                    n_gained,
                    centeroids]

        if nThreads != 1:
            all_threads = []
            all_i = range(n)
            all_given_i = []
            for thread in range(nThreads):

                if thread != nThreads - 1:
                    thread_i = all_i[thread::nThreads]
                    all_given_i += thread_i
                else:
                    thread_i = []
                    for i in all_i:
                        if i not in all_given_i:
                            thread_i.append(i)

                all_threads.append(AnalyzeThread(coverage_matrix, thread_i, *measures))

            for thread in all_threads:
                thread.start()

            for thread in all_threads:
                thread.join()
        else:
            for i in range(n):
                analyze(coverage_matrix, i, *measures)

        trace = calculate_cluster_trace(centeroids[::10], L, W)

        save("%s/extran_cluster_time.npy" % dir, time)
        save("%s/extran_cluster_trace.npy" % dir, trace)
        save("%s/extran_cluster_covs.npy" % dir, covs)
        save("%s/extran_cluster_size.npy" % dir, avg_cluster_sizes)
        save("%s/extran_cluster_n.npy" % dir, n_clusters_list)
        save("%s/extran_cluster_circs.npy" % dir, avg_circs)
        save("%s/extran_cluster_nbroken.npy" % dir, n_broken)
        save("%s/extran_cluster_ngained.npy" % dir, n_gained)


if __name__ == "__main__":
    main()