import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), "..", ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    supersaturations = []

    n = 0
    for data, _, _, _ in parser:

        supersaturation, speed = data.attrs["supersaturation"], data.attrs["GrowthSpeed"]

        if not supersaturation in supersaturations:
            supersaturations.append(supersaturation)

        n += 1

    print "Parsed", n, "entries."

    N = n/len(supersaturations)

    speeds = np.zeros(len(supersaturations))

    for data, _, _, _ in parser:

        supersaturation, speed = data.attrs["supersaturation"], data.attrs["GrowthSpeed"]

        i = supersaturations.index(supersaturation)

        speeds[i] += speed

    speeds /= N
    supersaturations = np.array(supersaturations)

    np.save("/tmp/lattice_supersaturations.npy", supersaturations)
    np.save("/tmp/lattice_growthspeeds.npy", speeds)

if __name__ == "__main__":
    main()