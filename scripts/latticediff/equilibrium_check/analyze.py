import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), "..", ".."))

from parse_h5_output import ParseKMCHDF5
from intercombinatorzor import ICZ


def main():

    input_file = sys.argv[1]

    try:
        every = int(sys.argv[2])
    except:
        every = 1

    parser = ParseKMCHDF5(input_file)

    supersaturations = []

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
        concentrations = parser.get_ignis_data(data, "Concentration")[::every].copy()
        time = parser.get_ignis_data(data, "Time")
        time /= (1 + supersaturation)
        time -= time[0] #Normalize to real time and not state time
        time = time[::every].copy()

        i = supersaturations.index(supersaturation)

        combinators[i].feed(time, heights, concentrations)

    n_steps = len(combinators[0]["Time"])

    print "Found", n_steps, "size data."

    all_times = np.zeros([len(supersaturations), n_steps])
    all_heights = np.zeros_like(all_times)
    all_concentrations = np.zeros_like(all_times)
    lengths = np.zeros(len(supersaturations))

    for i, supersaturation in enumerate(supersaturations):

        t, h = combinators[i].intercombine("Time", "AverageHeight")
        t2, c = combinators[i].intercombine("Time", "Concentration")

        if not (t == t2).all():
            sys.exit("Failure")

        l = len(t)

        lengths[i] = l

        all_heights[i, :l] = h
        all_concentrations[i, :l] = c
        all_times[i, :l] = t


    np.save("/tmp/confined_lattice_supersaturations.npy", supersaturations)
    np.save("/tmp/confined_lattice_lengths.npy", lengths)
    np.save("/tmp/confined_lattice_heights.npy", all_heights)
    np.save("/tmp/confined_lattice_concentrations.npy", all_concentrations)
    np.save("/tmp/confined_lattice_times.npy", all_times)

if __name__ == "__main__":
    main()
