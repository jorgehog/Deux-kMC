import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), "..", ".."))

from parse_h5_output import ParseKMCHDF5
from intercombinatorzor import ICZ


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    supersaturations = []

    try:
        every = int(sys.argv[2])
    except:
        every = 1

    n = 0
    for data, _, _, _ in parser:

        supersaturation = data.attrs["supersaturation"]

        if not supersaturation in supersaturations:
            supersaturations.append(supersaturation)

        n += 1

    print "Parsed", n, "entries."

    N = n/len(supersaturations)

    print "Found", N, "repeats."

    combinators = [[ICZ("Time", "AverageHeight", "Concentration") for _ in range(len(supersaturations))] for _ in [0, 1]]

    C = np.zeros((2, len(supersaturations)))
    CMAX = 100

    for data, _, _, _ in parser:

        supersaturation = data.attrs["supersaturation"]

        heights = parser.get_ignis_data(data, "AverageHeight").copy()

        conc = parser.get_ignis_data(data, "Concentration").copy()

        time = parser.get_ignis_data(data, "Time").copy()

        dtype = data.attrs["diffuse"] - 3

        i = supersaturations.index(supersaturation)

        if C[dtype, i] >= CMAX:
            continue

        combinators[dtype][i].feed(time, heights, conc)
        


    for dtype, name in enumerate(["radial", "pathfind"]):

        n_steps = len(combinators[dtype][0]["Time"])
        print "found", n_steps, "sized data."

        all_times = np.zeros([len(supersaturations), n_steps])
        all_heights = np.zeros_like(all_times)
        all_conc = np.zeros_like(all_times)
        lengths = np.zeros(len(supersaturations))

        for i, supersaturation in enumerate(supersaturations):

            t, h = combinators[dtype][i].intercombine("Time", "AverageHeight")
            t, c = combinators[dtype][i].intercombine("Time", "Concentration")

            l = len(t)

            lengths[i] = l

            all_times[i, :l] = t
            all_heights[i, :l] = h
            all_conc[i, :l] = c

            print "\r%d/%d" % (i+1, len(supersaturations)),
            sys.stdout.flush()
        print

        if not os.path.exists("/tmp/%s" % name):
            os.mkdir("/tmp/%s" % name)

        np.save("/tmp/%s/confined_%s_supersaturations.npy" % (name, name), supersaturations)
        np.save("/tmp/%s/confined_%s_times.npy" % (name, name), all_times)
        np.save("/tmp/%s/confined_%s_heights.npy" % (name, name), all_heights)
        np.save("/tmp/%s/confined_%s_concentrations.npy" % (name, name), all_conc)
        np.save("/tmp/%s/confined_%s_lengths.npy" % (name, name), lengths)

if __name__ == "__main__":
    main()
