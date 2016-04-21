import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), "..", ".."))

from parse_h5_output import ParseKMCHDF5

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    cs = []

    n = 0
    for data, _, _, _ in parser:

        c = data.attrs["depRateConstant"]

        if c not in cs:
            cs.append(c)

        n += 1

    print "Parsed", n, "entries."

    N = n/len(cs)

    print "Found", N, "repeats."

    speeds = np.zeros((2, len(cs)))

    for data, _, _, _ in parser:

        c = data.attrs["depRateConstant"]

        i = cs.index(c)

        type = data.attrs["diffuse"]

        speeds[type-3, i] += data.attrs["GrowthSpeed"]

    speeds /= N

    np.save("/tmp/mfptc_cs.npy", cs)
    np.save("/tmp/mfptc_growthspeeds.npy", speeds)

if __name__ == "__main__":
    main()