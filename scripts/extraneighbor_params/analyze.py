from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def Plify(x):
    return x
    return np.log(x)

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    Pls = []
    s0s = []

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        if alpha not in alphas:
            alphas.append(alpha)

        Pl = Plify(data.attrs["Pl"])
        Pl = round(Pl, 3)
        if Pl not in Pls:
            Pls.append(Pl)

        try:
            s0 = data.attrs["s0"]
            s0 = round(s0, 3)
            if s0 not in s0s:
                s0s.append(s0)
        except:
            s0s = [1.0]

    cmat = np.zeros(shape=(len(s0s), len(alphas), len(Pls)))
    ccounts = np.zeros_like(cmat)

    Pls = sorted(Pls)
    s0s = sorted(s0s)
    alphas = sorted(alphas)

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        ia = alphas.index(alpha)

        Pl = Plify(data.attrs["Pl"])
        Pl = round(Pl, 3)
        ipl = Pls.index(Pl)

        # try:
        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        is0 = s0s.index(s0)
        # except:
        #     s0 = 1
        #     is0 = 0

        coverage = data["coverage"]
        l = len(coverage)

        start = (9*l)/10

        k = 0
        cval = sum(coverage[start:])/float(L*W*(l-start))
        #
        # for i in range(start, L):
        #     cval += coverage[i]
        #     k += 1

        if cval != 0:
            cmat[is0, ia, ipl] += cval
            ccounts[is0, ia, ipl] += 1.

        # plab.plot(coverage)
        # print alpha, Pl, s0, cval
        # plab.show()

    I = np.where(ccounts != 0)

    cmat[I] /= ccounts[I]
    #
    # for is0, s0 in enumerate(s0s):
    #     f = plab.figure()
    #     ax = Axes3D(f)
    #     xpos, ypos = np.meshgrid(Pls, alphas)
    #
    #     print ypos.shape, xpos.shape, len(alphas), len(Pls), cmat[is0, :, :].shape
    #
    #     ax.plot_surface(xpos, ypos, cmat[is0, :, :], cstride=1, rstride=1)
    #
    #     print cmat[is0, :, :].max()
    #     plab.figure()
    #     plab.title(s0)
    #     plab.imshow(cmat[is0, :, :], interpolation='none', vmin=0, vmax=1, extent=[alphas[0], alphas[-1], Pls[0], Pls[-1]])
    #     plab.show()

    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)
    np.save("/tmp/extraneighbor_cmat.npy", cmat)


if __name__ == "__main__":
    main()