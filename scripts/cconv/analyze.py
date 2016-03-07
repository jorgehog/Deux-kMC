import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    radial = []
    pathfind = []
    heights = []

    n = 0
    for data, _, _, _ in parser:

        h0 = data.attrs["h0"]
        alpha = data.attrs["alpha"]
        type = data.attrs["type"]
        a = data.attrs["a"]
        b = data.attrs["b"]

        if h0 not in heights:
            heights.append(h0)

        i = heights.index(h0)

        if type == 0:
            radial.append(np.array([alpha, a, b, i]))
        else:
            pathfind.append(np.array([alpha, a, b, i]))

        n += 1

    print "Parsed", n, "entries."

    heights = np.array(heights)
    pathfind = np.array(pathfind)
    radial = np.array(radial)

    np.save("/tmp/cconv_heights.npy", heights)

    if len(pathfind) != 0:
        np.save("/tmp/cconv_pathfind.npy", pathfind)
    if len(radial) != 0:
        np.save("/tmp/cconv_radial.npy", radial)

if __name__ == "__main__":
    main()