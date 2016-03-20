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

    r_ids = []
    p_ids = []

    n = 0
    for data, _, _, run_id in parser:

        h0 = data.attrs["h0"]
        alpha = data.attrs["alpha"]
        type = data.attrs["type"]
        a = data.attrs["a"]
        b = data.attrs["b"]

        id = int(run_id)
        #
        # if id < 4E13:
        #     continue

        if h0 not in heights:
            heights.append(h0)

        i = heights.index(h0)

        arr = [alpha, a, b, i]

        add = True
        if type == 0:
            if arr in radial:
                j = radial.index(arr)

                if r_ids[j] < id:
                    radial[j] = arr
                    r_ids[j] = id
                    add = False

            if add:
                radial.append(arr)
                r_ids.append(id)

        else:
            if arr in pathfind:
                j = pathfind.index(arr)

                if p_ids[j] < id:
                    pathfind[j] = arr
                    p_ids[j] = id
                    add = False

            if add:
                pathfind.append(arr)
                p_ids.append(id)

        n += 1

    print "Parsed", n, "entries. Extracted", len(r_ids), len(p_ids), "values."

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