import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np
from math import floor

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    Pls = []
    s0s = []

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)

        Pl = data.attrs["Pl"]
        Pl = round(Pl, 3)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)

        if alpha not in alphas:
            alphas.append(alpha)

        if Pl not in Pls:
            Pls.append(Pl)

        if s0 not in s0s:
            s0s.append(s0)


    cmat = np.zeros(shape=(len(s0s), len(alphas), len(Pls)))
    ccounts = np.zeros_like(cmat)

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        ia = alphas.index(alpha)

        Pl = data.attrs["Pl"]
        Pl = round(Pl, 3)
        ipl = Pls.index(Pl)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        is0 = s0s.index(s0)

        coverage = data["coverage"]
        L = len(coverage)

        cval = np.array(coverage)[-1]

        cmat[is0, ia, ipl] += cval
        ccounts[is0, ia, ipl] += 1

    cmat /= ccounts

    for is0, s0 in enumerate(s0s):
        plab.figure()
        plab.title(s0)
        plab.imshow(cmat[0, :, :], interpolation='none')
        plab.show()

    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)
    np.save("/tmp/extraneighbor_cmat.npy", cmat)


if __name__ == "__main__":
    main()