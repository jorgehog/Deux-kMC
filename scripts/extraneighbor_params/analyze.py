from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np
from scipy.signal import argrelextrema

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def find_peaks(coverage, sign):

    diff = coverage[1:] - coverage[:-1]

    if sign == 1:
        m = coverage.max()
        return np.where(diff < -m/2)
    else:
        if coverage.mean() > coverage.max()/4:
            return find_peaks(coverage, 1)

        return [[]]

        X = np.arange(len(coverage))
        J = np.where(coverage == 1)[0][0]

        m = ((diff**2).mean() - diff.mean()**2)**0.5

        plab.plot(X[J:], coverage[J:])


        K = np.where(coverage[J:] > coverage[J:].mean())
        plab.plot(X[J:][K], coverage[J:][K])
        plab.show()

        I = list(argrelextrema(coverage, np.greater)[0])

        sheds = []
        for idx in range(1, len(I) - 1):
            elem = coverage[I[idx]]
            prev = coverage[I[idx-1]]
            next = coverage[I[idx+1]]

            if not ((prev < elem) and (next < elem)):
                sheds.append(idx)

        bigger = False
        this = None
        for idx in range(len(I)):
            if idx in sheds:
                continue

            if this:
                if coverage[I[idx]] > this:
                    bigger = True
                else:
                    this = coverage[I[idx]]
                    sheds.append(idx)
            else:
                this = coverage[I[idx]]
                sheds.append(idx)
        sheds = sorted(sheds)

        for idx in sheds[::-1]:
            I.pop(idx)



        return [np.array(I)]



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

    all_cmat = [np.zeros(shape=(len(s0s), len(alphas), len(Pls))),
                np.zeros(shape=(len(s0s), len(alphas), len(Pls))),
                np.zeros(shape=(len(s0s), len(alphas), len(Pls)))]
    all_ccounts = [np.zeros(shape=(len(s0s), len(alphas), len(Pls))),
                np.zeros(shape=(len(s0s), len(alphas), len(Pls))),
                np.zeros(shape=(len(s0s), len(alphas), len(Pls)))]

    names = ["neg", "eq", "pos"]

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

        tresh = 100

        for io, name in enumerate(names):

            coverage = data[name + "_coverage"][(0)]

            cmat = all_cmat[io]
            ccounts = all_ccounts[io]

            l = len(coverage)
            start = (9*l)/10

            #Exited premturely
            if l != lmax:

                #plab.plot(coverage)
                #plab.show()
                #Because there was no change in coverage
                if coverage[start:].mean() < tresh:
                    cval = 0
                #The system clogged
                else:
                    cval = L*W
            else:
                #equilibrium
                if io == 1:
                    cval = coverage[start:].mean()

                    #plab.plot(coverage)
                    #print alpha, Pl, s0, cval
                    #plab.show()

                else:

                    if io == 0:
                        X = find_peaks(coverage, 1)
                    else:
                        # print name, alpha, Pl, s0
                        # if not (alpha == 1.9 and Pl == 0.816 and s0 == 0.5 ):
                        #     continue
                        # plab.plot(coverage)
                        # plab.show()
                        X = find_peaks(coverage, -1)

                    if len(X[0]) != 0:
                        cval = coverage[X].mean()

                        # print name, s0, Pl, alpha
                        # plab.plot(coverage)
                        # plab.hold('on')
                        # plab.plot(X[0], coverage[X], 'ro')
                        # plab.plot([0, len(coverage) - 1], [cval, cval], "k-")
                        # plab.show()
                    else:
                        cval = coverage.mean()

            cmat[is0, ia, ipl] += cval/float(L*W)
            ccounts[is0, ia, ipl] += 1.

            # plab.plot(coverage)
            # print alpha, Pl, s0, cval
    #plab.show()

    print ccounts.max(), ccounts.min()

    for io in range(3):
        cmat = all_cmat[io]
        ccounts = all_ccounts[io]

        I = np.where(ccounts != 0)
        J = np.where(ccounts == 0)

        cmat[I] /= ccounts[I]
        cmat[J] = -1
        np.save("/tmp/extraneighbor_cmat_omega%d.npy" % io, cmat)

if __name__ == "__main__":
    main()