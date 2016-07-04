import sys
import os
import numpy as np
from os.path import join
from __builtin__ import enumerate
from matplotlib.pylab import *


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

from intercombinatorzor import ICZ

def find_front_pos(heights):
    return heights.mean()

#9 14 17 39 43 56 59 61 62 64
def main():

    input_file = sys.argv[1]

    skiplist = []
    if len(sys.argv) > 2:
        skiplist = [int(x) for x in sys.argv[2:]]

    parser = ParseKMCHDF5(input_file)

    def skip(data):
        return data.attrs["flux"] != 2.40

    every = 1000
    thermRatio = 0.25

    l = None
    n_entries = 0
    for data, L, W, run_id in parser:

        if skip(data):
            continue

        if not l:
            l = len(data["time"])
        n_entries += 1

    therm = l*thermRatio
    nbins = 20
    cmat = np.zeros(shape=(n_entries, nbins))
    dy = W/float(nbins)

    combinator = ICZ("Time", "ys")

    Cmean = 0
    cmeancount = 0
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
        tot_weight = 0
        for hi, heights_id in enumerate(stored_heights_indices):

            if hi % every != 0:
                continue

            heights = stored_heights[heights_id][()].transpose()

            ys = find_front_pos(heights)
            ys_vec[hi/every] = ys

            if hi < len(time) - 1:
                t_new = time[hi+1]

                if hi >= therm:
                    if heights_id in stored_particles:
                        particles = stored_particles[heights_id][()]
                    else:
                        particles = []

                    dt = t_new - t_prev

                    for x, y, _ in particles:
                        xl = round(x)
                        yl = round(y)
                        dh = conf_height - heights[xl, yl] - 1

                        cmat[entry_count, int((y+0.5)/dy)] += dt/dh
                        tot_weight += dt

                t_prev = t_new

            if hi % 100 == 0:
                sys.stdout.flush()
                print "\r%d/%d" % (hi, len(stored_heights)),

        cmat[entry_count, :] /= tot_weight

        sys.stdout.flush()
        print
        print "\rfin %d / %d" % (entry_count+1, n_entries)

        if skiplist:
            if entry_count + 1 in skiplist:
                ans = "asd"
            else:
                ans = ""
        else:
            plot(time[::every], ys_vec)
            show()

            ans = raw_input("discard? (n)")

        if ans == "":
            combinator.feed(time[::every], ys_vec)
            Cmean += cmat[entry_count, :]
            cmeancount += 1

        entry_count += 1

    print

    Cmean /= cmeancount

    ti, ys_veci = combinator.intercombine("Time", "ys")

    np.save("/tmp/FelixSeqC_t.npy", ti)
    np.save("/tmp/FelixSeqC_ys.npy", ys_veci)
    np.save("/tmp/FelixSeqC_C.npy", Cmean)

if __name__ == "__main__":
    main()
