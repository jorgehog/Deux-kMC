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

    every = 250

    n = 0
    for data, _, _, _ in parser:

        supersaturation = data.attrs["supersaturation"]

        if not supersaturation in supersaturations:
            supersaturations.append(supersaturation)

        n += 1

    print "Parsed", n, "entries."

    N = n/len(supersaturations)

    print "Found", N, "repeats."

    combinators = [ICZ("Time", "AverageHeight", "Concentration") for _ in range(len(supersaturations))]

    for data, _, _, _ in parser:

        supersaturation = data.attrs["supersaturation"]\

        heights = parser.get_ignis_data(data, "AverageHeight")[::every]
        conc = parser.get_ignis_data(data, "Concentration")[::every]
        time = parser.get_ignis_data(data, "Time")[::every]

        i = supersaturations.index(supersaturation)

        combinators[i].feed(time, heights, conc)

    n_steps = len(combinators[0]["Time"])

    all_times = np.zeros([len(supersaturations), n_steps])
    all_heights = np.zeros_like(all_times)
    all_conc = np.zeros_like(all_times)
    lengths = np.zeros(len(supersaturations))

    supersaturations = sorted(supersaturations  )

    for i, supersaturation in enumerate(supersaturations):

        t, h = combinators[i].intercombine("Time", "AverageHeight")
        t, c = combinators[i].intercombine("Time", "Concentration")

        l = len(t)

        lengths[i] = l

        all_times[i, :l] = t
        all_heights[i, :l] = h
        all_conc[i, :l] = c

        print "\r%d/%d" % (i+1, len(supersaturations))
        sys.stdout.flush()

    print

    np.save("/tmp/cconc_supersaturations.npy", supersaturations)
    np.save("/tmp/cconc_times.npy", all_times)
    np.save("/tmp/cconc_heights.npy", all_heights)
    np.save("/tmp/cconc_conc.npy", all_conc)
    np.save("/tmp/cconc_lengths.npy", lengths)

if __name__ == "__main__":
    main()