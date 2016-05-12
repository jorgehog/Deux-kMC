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
    omegas = []
    lmax = 0

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

        try:
            omega = data.attrs["omega"]
            omega = round(omega, 3)
        except Exception as inst:
            print inst
            omega = 0

        if omega not in omegas:
            omegas.append(omega)

        l = len(data["coverage"])
        if l > lmax:
            lmax = l

    Pls = sorted(Pls)
    s0s = sorted(s0s)
    alphas = sorted(alphas)
    omegas = sorted(omegas)

    np.save("/tmp/extraneighbor_omegas.npy", omegas)
    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)

    for this_io, omega in enumerate(omegas):

        cmat = np.zeros(shape=(len(s0s), len(alphas), len(Pls)))
        ccounts = np.zeros_like(cmat)

        for data, L, W, run_id in parser:

            try:
                omega = data.attrs["omega"]
                omega = round(omega, 3)
            except Exception as inst:
                print inst
                omega = 0.0

            io = omegas.index(omega)

            if io != this_io:
                continue

            alpha = data.attrs["alpha"]
            alpha = round(alpha, 3)
            ia = alphas.index(alpha)

            Pl = data.attrs["Pl"]
            Pl = round(Pl, 3)
            ipl = Pls.index(Pl)

            try:
                s0 = data.attrs["s0"]
                s0 = round(s0, 3)
            except Exception as inst:
                print inst
                s0 = 1.0

            is0 = s0s.index(s0)

            coverage = data["coverage"][()]

            l = len(coverage)
            start = (9*l)/10

            #Completely sealed
            if l != lmax:
                cval = L*W
            else:
                if omega == 0:
                    cval = coverage[start:].mean()

                else:
                    X = argrelextrema(coverage[start:], np.greater)

                    if len(X[0]) != 0:
                        cval = coverage[start:][X].mean()

                        # if omega < 0 and cval > 100 and cval < 800:
                        #     print omega, s0, Pl, alpha
                        #     plab.plot(coverage[start:])
                        #     plab.hold('on')
                        #     plab.plot(X[0], coverage[start:][X], 'ro')
                        #     plab.plot([0, len(coverage[start:]) - 1], [cval, cval], "k-")
                        #     plab.show()
                    else:
                        cval = coverage[start:].mean()

            cmat[is0, ia, ipl] += cval/float(L*W)
            ccounts[is0, ia, ipl] += 1.

            # plab.plot(coverage)
            # print alpha, Pl, s0, cval
            # plab.show()

        I = np.where(ccounts != 0)
        J = np.where(ccounts == 0)

        cmat[I] /= ccounts[I]
        cmat[J] = -1

        np.save("/tmp/extraneighbor_cmat_omega%d.npy" % this_io, cmat)


if __name__ == "__main__":
    main()