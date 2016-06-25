import sys
import os
import numpy as np
from os.path import join
from matplotlib.pylab import *


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

from intercombinatorzor import ICZ

def find_front_pos(heights):
    return heights.mean()/float(heights.shape[1])

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    def skip(data):
        return data.attrs["flux"] != 2.40

    nbins = 10
    every = 1

    l = None
    n_entries = 0
    for data, L, W, run_id in parser:

        if skip(data):
            continue

        if not l:
            l = len(data["time"])
        n_entries += 1

    cmat = np.zeros(shape=(l/every, nbins))
    dy = W/float(nbins)

    combinator = ICZ("Time", "ys")

    entry_count = 0
    for data, L, W, run_id in parser:

        if skip(data):
            continue

        conf_height = data.attrs["height"]

        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]

        stored_heights_indices = sorted(stored_heights, key=lambda x: int(x))

        time = data["time"][()]

        ys_vec = np.zeros(len(time)/every)

        #Shift with 1 to translate from starting time till ending times
        t_prev = 0
        for hi, heights_id in enumerate(stored_heights_indices):

            if hi % every != 0:
                continue

            t_new = time[hi]
            dt = t_new - t_prev
            t_prev = t_new

            if heights_id in stored_particles:
                particles = stored_particles[heights_id][()]
            else:
                particles = []

            heights = stored_heights[heights_id][()].transpose()

            ys = find_front_pos(heights)
            ys_vec[hi/every] = ys

            #
            # for x, y, _ in particles:
            #     xl = round(x)
            #     yl = round(y)
            #     dh = conf_height - heights[xl, yl] - 1
            #
            #     cmat[hi/every, int((y+0.5)/dy)] += dt/dh
            #
            # for x in range(L):
            #     for y in range(W):
            #         height = heights[x, y]
            #         cmat[hi/every, y/dy] += dt*height

            if hi % 100 == 0:
                sys.stdout.flush()
                print "\r%d/%d" % (hi, len(stored_heights)),

        entry_count += 1
        combinator.feed(time[::every], ys_vec)

        sys.stdout.flush()
        print
        print "\rfin %d / %d" % (entry_count, n_entries)
        plot(time[::every], ys_vec)

        break
    print
    show()

    ti, ys_veci = combinator.intercombine("Time", "ys")

    np.save("/tmp/FelixSeqC_t.npy", ti)
    np.save("/tmp/FelixSeqC_ys.npy", ys_veci)

if __name__ == "__main__":
    main()
