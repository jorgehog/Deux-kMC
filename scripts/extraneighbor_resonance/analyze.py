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
    Pls = []
    s0s = []
    lmax = 0
    omegashift = None

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        if alpha not in alphas:
            alphas.append(alpha)

        Pl = data.attrs["Pl"]
        Pl = round(Pl, 3)
        if Pl not in Pls:
            Pls.append(Pl)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        if s0 not in s0s:
            s0s.append(s0)

        if not omegashift:
            omegashift = data.attrs["omegaShift"]
        else:
            if data.attrs["omegaShift"] != omegashift:
                raise RuntimeError("invalid omegashift in datasets.")

        l = data["eq_storedEventValues"].shape[1]
        if l > lmax:
            lmax = l

    Pls = sorted(Pls)
    s0s = sorted(s0s)
    alphas = sorted(alphas)

    np.save("/tmp/extraneighbor_omegas.npy", [0, omegashift])
    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)

    cmat = np.zeros(shape=(2, len(s0s), len(alphas), len(Pls)))
    ccounts = np.zeros(shape=(2, len(s0s), len(alphas), len(Pls)))

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

        eq_coverage = data["eq_storedEventValues"][(0)]
        omega_coverage = data["omega_storedEventValues"][(0)]

        cval = None
        for i, coverage in enumerate([eq_coverage, omega_coverage]):

            l = len(eq_coverage)
            start = (9*l)/10

            if l != lmax:
                #The system clogged
                if coverage[-1] > (L*W)/10:
                    cval = L*W

                #Because there was no change in coverage
                else:
                    cval = 0
            else:
                #equilibrium
                if i == 0:
                    cval = coverage[start:].mean()

                    if Pl == 0.861 and s0 == 1.5 and alpha >= 1.5:
                        print alpha, cval
                        plab.plot(coverage)
                        plab.show()

                else:
                    X = find_peaks(coverage)

                    if len(X[0]) != 0:
                        cval = coverage[X].mean()
                    else:
                        cval = coverage[start:].mean()

            cmat[i, is0, ia, ipl] += cval/float(L*W)
            ccounts[i, is0, ia, ipl] += 1.

    print ccounts.max(), ccounts.min()


    I = np.where(ccounts != 0)
    J = np.where(ccounts == 0)

    cmat[I] /= ccounts[I]
    cmat[J] = -1
    np.save("/tmp/extraneighbor_cmat.npy", cmat)

if __name__ == "__main__":
    main()