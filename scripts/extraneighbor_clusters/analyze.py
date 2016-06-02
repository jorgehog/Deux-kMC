import sys
import os
from os.path import join
import numpy as np
from numpy import where, empty, save, zeros, unique
from skimage import measure
from threading import Thread
from matplotlib import pylab as plab

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

min_area = 4

def find_perimeters(pim, labels):
    L, W = pim.shape
    expanded = zeros(shape=[3*d for d in pim.shape], dtype=int)
    for x in range(3):
        for y in range(3):
            expanded[x*L:(x+1)*L, y*W:(y+1)*W] = pim

    n_labels = len(labels)

    idx_map = {}

    n = 0
    for l in labels:
        idx_map[l] = n
        n += 1

    winners = [None for _ in range(n_labels)]
    winner_areas = [0 for _ in range(n_labels)]

    l = measure.label(expanded)
    l[np.where(expanded == 0)] = -1
    l += 1
    #figure()
    #pcolor(l)

    def find_fast(im, v, x0, y0, x1, y1):
        for x in range(x0, x1):
            for y in range(y0, y1):
                if im[x, y] == v:
                    return x, y

    expanded_props = measure.regionprops(l)
    for prop in expanded_props:

        x, y = find_fast(l, prop.label, *prop.bbox)
        orig_label = expanded[x, y]

        idx = idx_map[orig_label]

        if prop.area > winner_areas[idx]:
            winner_areas[idx] = prop.area
            winners[idx] = [orig_label, (prop.perimeter, prop.eccentricity)]

    return dict(winners)

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
            eccentricities,
            n_broken,
            n_gained):

    _s = data[i, :, :]

    covs[i] = _s.mean()

    im_label = measure.label(_s)
    im_label[where(_s == 0)] = -1

    im_label_per = periodify(im_label) + 1

    props = measure.regionprops(im_label_per)

    all_areas = []
    all_circs = []
    all_eccs = []

    if props:
        perimeter_data = find_perimeters(im_label_per, [prop.label for prop in props])

        for m in props:
            if m.label != 0 and m.area > min_area:
                all_areas.append(m.area)
                all_circs.append(perimeter_data[m.label][0])
                all_eccs.append(perimeter_data[m.label][1])

    n_clusters = len(all_areas)

    if n_clusters == 0:
        avg_cluster_size = 0
        avg_circ = 0
        avg_eccs = 0
    else:
        avg_cluster_size = sum(all_areas)/len(all_areas)
        avg_circ = sum(all_circs)/len(all_circs)
        avg_eccs = sum(all_eccs)/len(all_eccs)

    avg_cluster_sizes[i] = avg_cluster_size
    n_clusters_list[i] = n_clusters
    avg_circs[i] = avg_circ
    eccentricities[i] = avg_eccs

    if i == data.shape[0] - 1:

        plab.pcolor(_s)
        plab.show()

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

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    nThreads = int(sys.argv[2])

    if len(sys.argv) > 3:
        every = int(sys.argv[3])
    else:
        every = 1

    for data, L, W, run_id in parser:

        print L, W, run_id

        #coverage_matrix = data["eq_coverage_matrix"][()][::every]
        coverage_matrix = np.load("/tmp/eq_coverage_new.npy")[::every]

        n = coverage_matrix.shape[0]

        covs = empty(n)
        avg_cluster_sizes = empty(n)
        n_clusters_list = empty(n)
        avg_circs = empty(n)
        eccentricities = empty(n)
        n_broken = empty(n - 1)
        n_gained = empty(n - 1)

        measures = [covs,
                    avg_cluster_sizes,
                    n_clusters_list,
                    avg_circs,
                    eccentricities,
                    n_broken,
                    n_gained]

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

        save("/tmp/extran_cluster_covs.npy", covs)
        save("/tmp/extran_cluster_size.npy", avg_cluster_sizes)
        save("/tmp/extran_cluster_n.npy", n_clusters_list)
        save("/tmp/extran_cluster_circs.npy", avg_circs)
        save("/tmp/extran_cluster_eccs.npy", eccentricities)
        save("/tmp/extran_cluster_nbroken.npy", n_broken)
        save("/tmp/extran_cluster_ngained.npy", n_gained)

        break


if __name__ == "__main__":
    main()