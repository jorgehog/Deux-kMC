import sys
import os
from os.path import join
import numpy as np
from scipy.optimize import curve_fit
from math import exp, sqrt

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def find_corrl(x, y):
    I = np.where(y > exp(-2))

    f = lambda _x, b: np.exp(-np.abs(_x)/b)

    # from matplotlib.pylab import plot, show
    # plot(x, y)
    # show()

    p0 = (1.)
    pl, _ = curve_fit(f, x[I], y[I], p0)

    return pl[0]

def analyze(input_file, typeint, typestr):

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    heights = []

    for data, L, W, run_id in parser:

        if data.attrs["diffuse"] != typeint:
            continue

        alpha = data.attrs["alpha"]
        height = data.attrs["confiningSurfaceHeight"]

        if alpha not in alphas:
            alphas.append(alpha)

        if height not in heights:
            heights.append(height)


    autocorrs = np.zeros(shape=(len(heights), len(alphas), L+1, W+1))
    RMSes = np.zeros(shape=(len(heights), len(alphas)))
    counts = np.zeros_like(RMSes)

    for data, L, W, run_id in parser:

        if data.attrs["diffuse"] != typeint:
            continue

        alpha = data.attrs["alpha"]
        height = data.attrs["confiningSurfaceHeight"]

        ia = alphas.index(alpha)
        ih = heights.index(height)

        xl = np.linspace(-L/2, L/2, L + 1)
        xw = np.linspace(-W/2, W/2, W + 1)

        x1 = sqrt(2)*xl
        x2 = sqrt(2)*xw

        autocorr = np.array(data["autocorrelation"])

        autocorrs[ih, ia, :, :] += autocorr

        labels = data["ignisEventDescriptions"]
        i = list(list(labels)[0]).index("HeightRMS@MainMesh")
        RMS = data["ignisData"][i, :]
        rms = (RMS[len(RMS)/2:]).mean()

        RMSes[ih, ia] += rms
        counts[ih, ia] += 1

    RMSes /= counts

    for j, height in enumerate(heights):

        Cs = []

        for i in range(len(alphas)):

            autocorr = autocorrs[j, i, :, :]

            autocorr /= autocorr.max()

            dl = autocorr[:, W/2]
            dw = autocorr[L/2, :]

            d1 = np.diag(autocorr)
            d2 = np.diag(np.flipud(autocorr))

            cl = find_corrl(xl, dl)
            cw = find_corrl(xw, dw)
            c1 = find_corrl(x1, d1)
            c2 = find_corrl(x2, d2)

            Cs.append([(cl + cw)/2., (c1 + c2)/2.])

        alphas = np.array(alphas)
        Cs = np.array(Cs)

        np.save("/tmp/acf_h%d_%s_alphas.npy" % (j, typestr), alphas)
        np.save("/tmp/acf_h%d_%s_RMSes.npy" % (j, typestr), RMSes[j, :])
        np.save("/tmp/acf_h%d_%s_Cs.npy" % (j, typestr), Cs)

    np.save("/tmp/acf_heights.npy", np.array(heights))

def main():

    input_file = sys.argv[1]

    analyze(input_file, 6, "uniform")
    analyze(input_file, 2, "lattice")
    analyze(input_file, 3, "radial")
    analyze(input_file, 4, "pathfind")


if __name__ == "__main__":
    main()