import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    autocorr = 0.0

    n = 0
    for data, _, _, run_id in parser:

        autocorr += data["autocorrelation"].__array__()

        n += 1

    print "Parsed", n, "entries."

    np.save("/tmp/autocorr.npy", autocorr/n)


if __name__ == "__main__":
    main()