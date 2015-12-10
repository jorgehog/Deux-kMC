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

    combinators = [ICZ("Time", "AverageHeight", "Concentration") for _ in range(len(supersaturations))]

    for data, _, _, _ in parser:

        supersaturation = data.attrs["supersaturation"]\

        heights = parser.get_ignis_data(data, "AverageHeight")[::every].copy()

        conc = parser.get_ignis_data(data, "Concentration")[::every].copy()

        time = parser.get_ignis_data(data, "Time")[::every].copy()

        i = supersaturations.index(supersaturation)

        try:
            combinators[i].feed(time, heights, conc)
        except:
            pass

    n_steps = len(combinators[0]["Time"])

    print "found", n_steps, "sized data."

    all_times = np.zeros([len(supersaturations), n_steps])
    all_heights = np.zeros_like(all_times)
    all_conc = np.zeros_like(all_times)
    lengths = np.zeros(len(supersaturations))

    for i, supersaturation in enumerate(supersaturations):

        t, h = combinators[i].mean("Time", "AverageHeight")
        t, c = combinators[i].mean("Time", "Concentration")

        l = len(t)

        lengths[i] = l

        all_times[i, :l] = t
        all_heights[i, :l] = h
        all_conc[i, :l] = c

        print "\r%d/%d" % (i+1, len(supersaturations)),
        sys.stdout.flush()

    print

    np.save("/tmp/confined_mfpt_supersaturations.npy", supersaturations)
    np.save("/tmp/confined_mfpt_times.npy", all_times)
    np.save("/tmp/confined_mfpt_heights.npy", all_heights)
    np.save("/tmp/confined_mfpt_concentrations.npy", all_conc)
    np.save("/tmp/confined_mfpt_lengths.npy", lengths)

if __name__ == "__main__":
    main()