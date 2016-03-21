import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def dump_map(set, setname):

    s = "const maptype valueMap%s%s = {" % (setname[0].upper(), setname[1:])

    set = sorted(set, key=lambda x: x[2]*100 + x[3])

    for a, b, i, j in set:

        c = (a+b)/2.

        s += "{{%2d, %2d}, %.3f},\n" % (i, j, c)

    s = s.strip(",\n") + "};"

    print s

def dump_array(arr, name):

    s = "const vec %s = {" % name

    for i, v in enumerate(arr):
        s += "%g, " % v

        if (i+1) % 4 == 0:
            s += "\n"

    s = s.strip(", ").strip(", \n") + "};"

    print s

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    radial = []
    pathfind = []
    heights = []
    alphas = []

    r_ids = []
    p_ids = []

    n = 0
    for data, _, _, run_id in parser:

        h0 = data.attrs["h0"]
        alpha = data.attrs["alpha"]

        alpha = round(alpha, 3)

        if h0 not in heights:
            heights.append(h0)
        if alpha not in alphas:
            alphas.append(alpha)

    heights = sorted(heights)
    alphas = sorted(alphas)

    for data, _, _, run_id in parser:

        h0 = data.attrs["h0"]
        alpha = data.attrs["alpha"]
        type = data.attrs["type"]
        a = data.attrs["a"]
        b = data.attrs["b"]

        i = heights.index(h0)
        j = alphas.index(alpha)

        arr = [a, b, i, j]

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
    alphas = np.array(alphas)
    pathfind = np.array(pathfind)
    radial = np.array(radial)

    np.save("/tmp/cconv_heights.npy", heights)
    np.save("/tmp/cconv_alphas.npy", alphas)

    if len(pathfind) != 0:
        np.save("/tmp/cconv_pathfind.npy", pathfind)
    if len(radial) != 0:
        np.save("/tmp/cconv_radial.npy", radial)

    dump_map(radial, "radial")
    dump_map(pathfind, "pathfind")
    dump_array(heights, "knownHeights")
    dump_array(alphas, "knownAlphas")

if __name__ == "__main__":
    main()