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

    system_sizes = []
    alpha_values = []

    n = 0
    for data, L, W, _ in parser:

        if L != W:
            raise RuntimeError("Incorrect system size.")

        if not L in system_sizes:
            system_sizes.append(L)

        alpha = data.attrs["alpha"]

        if alpha not in alpha_values:
            alpha_values.append(alpha)

        n += 1

    print "Parsed", n, "entries."

    N = n/(len(system_sizes)*len(alpha_values))

    print "Found", N, "repeats."

    sizes = np.zeros((len(system_sizes), len(alpha_values)))
    counts = np.zeros((len(system_sizes), len(alpha_values)))

    alpha_values = sorted(alpha_values)

    for data, L, W, _ in parser:

        alpha = data.attrs["alpha"]
        size = data.attrs["size"]

        i = system_sizes.index(L)
        j = alpha_values.index(alpha)

        sizes[i, j] += size
        counts[i, j] += 1

    counts[np.where(counts == 0)] = 1

    sizes /= counts

    np.save("/tmp/size_test_system_sizes.npy", np.array(system_sizes))
    np.save("/tmp/size_test_alpha_values.npy", np.array(alpha_values))
    np.save("/tmp/size_test_sizes.npy", sizes)

if __name__ == "__main__":
    main()