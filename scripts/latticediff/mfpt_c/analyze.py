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

    cs = []

    n = 0
    for data, _, _, _ in parser:

        c = data.attrs["depRateConstant"]

        if not c in cs:
            cs.append(c)

        n += 1

    print "Parsed", n, "entries."

    N = n/len(cs)

    print "Found", N, "repeats."

    speeds = np.zeros_like(cs)

    for data, _, _, _ in parser:

        c = data.attrs["depRateConstant"]

        i = cs.index(c)

        speeds[i] += data.attrs["GrowthSpeed"]

    speeds /= N

    np.save("/tmp/mfptc_cs.npy", cs)
    np.save("/tmp/mfptc_growthspeeds.npy", speeds)

if __name__ == "__main__":
    main()