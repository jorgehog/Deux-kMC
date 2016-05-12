import sys
import os
from os.path import join
from numpy import where, empty, save
from skimage import measure
from threading import Thread

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


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
            avg_cluster_sizes,
            n_clusters_list,
            avg_circs,
            n_broken,
            n_gained):

    _s = data[i, :, :]

    im_label = measure.label(_s)
    im_label[where(_s == 0)] = -1

    im_label_per = periodify(im_label)

    props = measure.regionprops(im_label_per)

    all_areas = []
    all_circs = []

    for m in props:
        if m.label != -1:
            all_areas.append(m.area)
            all_circs.append(m.perimeter)

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

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    nThreads = 8

    for data, L, W, run_id in parser:

        print L, W, run_id

        coverage_matrix = data["coverage_matrix_eq"][()]

        n = coverage_matrix.shape[0]

        avg_cluster_sizes = empty(n)
        n_clusters_list = empty(n)
        avg_circs = empty(n)
        n_broken = empty(n - 1)
        n_gained = empty(n - 1)

        measures = [avg_cluster_sizes,
                    n_clusters_list,
                    avg_circs,
                    n_broken,
                    n_gained]

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

        save("/tmp/extran_cluster_size.npy", avg_cluster_sizes)
        save("/tmp/extran_cluster_n.npy", n_clusters_list)

        break


if __name__ == "__main__":
    main()