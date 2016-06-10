from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np
from scipy.signal import argrelextrema

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def find_peaks(coverage):

    diff = coverage[1:] - coverage[:-1]

    m = coverage.max()
    return np.where(diff < -m/2)

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    F0s = []
    s0s = []
    lmax = 0

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        if alpha not in alphas:
            alphas.append(alpha)

        F0 = data.attrs["F0"]
        F0 = round(F0, 3)
        if F0 not in F0s:
            F0s.append(F0)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        if s0 not in s0s:
            s0s.append(s0)

        l = data["eq_storedEventValues"].shape[1]
        if l > lmax:
            lmax = l

    F0s = sorted(F0s)
    s0s = sorted(s0s)
    alphas = sorted(alphas)

    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_F0s.npy", F0s)

    cmat = np.zeros(shape=(len(s0s), len(alphas), len(F0s)))
    ccounts = np.zeros(shape=(len(s0s), len(alphas), len(F0s)))

    facSmall = 9/10.
    facBig = 9/10.
    plotevery = 100
    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        ia = alphas.index(alpha)

        F0 = data.attrs["F0"]
        F0 = round(F0, 3)
        iF0 = F0s.index(F0)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        is0 = s0s.index(s0)

        eq_coverage = data["eq_storedEventValues"][(0)]

        l = len(eq_coverage)

        if l != lmax:
            #The system clogged
            if eq_coverage[-1] > (L*W)/10:
                cval = L*W

            #Because there was no change in coverage
            else:
                cval = 0
        else:
            if eq_coverage[:l/2].mean() < L*W/10:
                start = int(facSmall*l)
                plab.plot(eq_coverage[::plotevery]/float(L*W))
            else:
                start = int(facBig*l)

            cval = eq_coverage[start:].mean()

                #
                # else:
                #     X = find_peaks(coverage)
                #
                #     if len(X[0]) != 0:
                #         cval = coverage[X].mean()
                #     else:
                #         cval = coverage[start:].mean()


        cmat[is0, ia, iF0] += cval/float(L*W)
        ccounts[is0, ia, iF0] += 1.

    print ccounts.max(), ccounts.min()

    I = np.where(ccounts != 0)
    J = np.where(ccounts == 0)

    cmat[I] /= ccounts[I]
    cmat[J] = -1
    np.save("/tmp/extraneighbor_cmat.npy", cmat)

    plab.plot([lmax*facSmall/plotevery, lmax*facSmall/plotevery], [0, 1])
    plab.show()

if __name__ == "__main__":
    main()