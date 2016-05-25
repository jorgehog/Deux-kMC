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

    F0s = []
    lmax = 0

    for data, L, W, run_id in parser:

        F0 = data.attrs["F0"]
        F0 = round(F0, 5)
        if F0 not in F0s:
            F0s.append(F0)

        l = data["eq_storedEventValues"].shape[1]
        if l > lmax:
            lmax = l

    F0s = sorted(F0s)

    np.save("/tmp/resonance_F0s.npy", F0s)

    cvec = np.zeros(len(F0s))
    counts = np.zeros_like(cvec)

    for data, L, W, run_id in parser:
        F0 = data.attrs["F0"]
        F0 = round(F0, 5)
        if0 = F0s.index(F0)

        eq_coverage = data["eq_storedEventValues"][(0)]

        l = len(eq_coverage)
        start = (9*l)/10

        if l != lmax:
            if sum(eq_coverage) != 0:
                print eq_coverage.mean(), l, lmax
            cval = 0
        else:
            cval = eq_coverage[start:].mean()

        #plab.plot(eq_coverage)

        cvec[if0] += cval/float(L*W)
        counts[if0] += 1.

    #plab.show()
    print counts.max(), counts.min()

    if counts.max() != counts.min():
        print np.array(F0s)[np.where(counts == counts.max())]
        print np.array(F0s)[np.where(counts == counts.min())]


    cvec /= counts

    np.save("/tmp/resonance_cvec.npy", cvec)

if __name__ == "__main__":
    main()