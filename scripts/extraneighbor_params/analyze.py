from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np
from scipy.signal import argrelextrema

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


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

        try:
            s0 = data.attrs["s0"]
            s0 = round(s0, 3)
        except Exception as inst:
            print inst
            s0= 1.0

        if s0 not in s0s:
            s0s.append(s0)

        if not omegashift:
            omegashift = data.attrs["omegaShift"]
        else:
            if data.attrs["omegaShift"] != omegashift:
                raise RuntimeError("invalid omegaShift in datasets.")

        l = data["eq_coverage"].shape[1]
        if l > lmax:
            lmax = l

    omegas = [-omegashift, 0, omegashift]

    Pls = sorted(Pls)
    s0s = sorted(s0s)
    alphas = sorted(alphas)
    omegas = sorted(omegas)

    np.save("/tmp/extraneighbor_omegas.npy", omegas)
    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)

    all_cmat = {}
    all_ccounts = {}

    all_cmat["neg"] = np.zeros(shape=(len(s0s), len(alphas), len(Pls)))
    all_cmat["eq"] = np.zeros_like(all_cmat["neg"])
    all_cmat["pos"] = np.zeros_like(all_cmat["neg"])

    all_ccounts["neg"] = np.zeros_like(all_cmat["neg"])
    all_ccounts["eq"] = np.zeros_like(all_cmat["neg"])
    all_ccounts["pos"] = np.zeros_like(all_cmat["neg"])

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

        for io, name in enumerate(["neg", "eq", "pos"]):
            cmat = all_cmat[name]
            ccounts = all_ccounts[name]

            coverage = data[name + "_coverage"][(0)]

            l = len(coverage)
            start = -1000

            #Completely sealed
            if l != lmax:
                if coverage[-1] == 0:
                    cval = 0
                else:
                    cval = L*W
            else:
                #equilibrium
                if io == 1:
                    cval = coverage[start:].mean()

                    plab.plot(coverage)
                    #print alpha, Pl, s0, cval
                    #plab.show()

                else:
#                    print name, alpha, Pl, s0
                    #plab.plot(coverage)
                    #plab.show()

                    X = argrelextrema(coverage[start:], np.greater)

                    if len(X[0]) != 0:
                        cval = coverage[start:][X].mean()

                        # print name, s0, Pl, alpha
                        # plab.plot(coverage[start:])
                        # plab.hold('on')
                        # plab.plot(X[0], coverage[start:][X], 'ro')
                        # plab.plot([0, len(coverage[start:]) - 1], [cval, cval], "k-")
                        # plab.show()
                    else:
                        cval = coverage[start:].mean()

            cmat[is0, ia, ipl] += cval/float(L*W)
            ccounts[is0, ia, ipl] += 1.

            # plab.plot(coverage)
            # print alpha, Pl, s0, cval
    #plab.show()

    for io, name in enumerate(["neg", "eq", "pos"]):
        cmat = all_cmat[name]
        ccounts = all_ccounts[name]

        I = np.where(ccounts != 0)
        J = np.where(ccounts == 0)

        cmat[I] /= ccounts[I]
        cmat[J] = -1
        np.save("/tmp/extraneighbor_cmat_omega%d.npy" % (2-io), cmat)

if __name__ == "__main__":
    main()