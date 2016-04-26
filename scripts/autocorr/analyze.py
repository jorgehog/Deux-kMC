from matplotlib.pylab import plot, show
import sys
import os
from os.path import join
import numpy as np
from scipy.optimize import curve_fit
from math import exp, sqrt

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

from scipy.stats import linregress


def find_match(x, y, start, end):
    slope, intercept, r, p, stderr = linregress(x[start:end], y[start:end])
    return abs(1-r**2), slope

def fit_best_line(x, y):

    all_err = []
    ends = range(3, len(x))

    for end in ends:
        err, slope = find_match(x, y, 0, end)

        all_err.append([end, err, slope])

    return min(all_err, key=lambda x: x[1])

def find_corrl(x, y):
    I = np.where(y > exp(-2))

    f = lambda _x, b: np.exp(-np.abs(_x)/b)

    # from matplotlib.pylab import plot, show
    # plot(x, y)
    # show()

    p0 = (1.)
    pl, _ = curve_fit(f, x[I], y[I], p0)

    return pl[0]

def find_corrl2(x, y, doPlot=False):

    I = np.where(y > 0)
    J = np.where(x[I] >= 0)

    X = x[I][J]
    Y = np.log((y[I][J]))

    # print "[",
    # for _y in y:
    #     print _y, ", ",
    # print "]"

    try:
        end, err, slope = fit_best_line(X, Y)
    except:
        print "ERROR", X, Y
        return 0

    if doPlot:
        plot(X, Y)
        plot(X[:end], Y[:end], 'rx')
        print end, err, slope
        print -1/slope
        show()

    return -1./slope

def analyze(input_file, typeint, typestr):

    print typeint

    parser = ParseKMCHDF5(input_file)

    alphas = []
    heights = []

    n = 0
    for data, L, W, run_id in parser:

        if data.attrs["diffuse"] != typeint:
            continue

        n+=1
     #   continue
        alpha = data.attrs["alpha"]
        height = data.attrs["confiningSurfaceHeight"]

        if alpha not in alphas:
            alphas.append(alpha)

        if height not in heights:
            heights.append(height)
    print "n:",n
    #return

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

        autocorr = np.array(data["autocorrelation"])

        autocorrs[ih, ia, :, :] += autocorr

        labels = data["ignisEventDescriptions"]
        i = list(list(labels)[0]).index("HeightRMS@MainMesh")
        RMS = data["ignisData"][i, :]
        rms = (RMS[len(RMS)/2:]).mean()

        RMSes[ih, ia] += rms
        counts[ih, ia] += 1
    show()

    RMSes /= counts

    xl = np.linspace(-L/2, L/2, L + 1)
    xw = np.linspace(-W/2, W/2, W + 1)

    x1 = sqrt(2)*xl
    x2 = sqrt(2)*xw

    heights = sorted(heights)

    for j, height in enumerate(heights):

        Cs = []

        for i in range(len(alphas)):

            autocorr = autocorrs[j, i, :, :]

            autocorr /= autocorr.max()

            dl = autocorr[:, W/2]
            dw = autocorr[L/2, :]

            d1 = np.diag(autocorr)
            d2 = np.diag(np.flipud(autocorr))

            print height, alphas[i]

            doPlot = (typeint == 2) and (alphas[i] == 0.6) and False

            cl = find_corrl2(xl, dl, doPlot)
            cw = find_corrl2(xw, dw, doPlot)
            c1 = find_corrl2(x1, d1, doPlot)
            c2 = find_corrl2(x2, d2, doPlot)

            Cs.append([(cl + cw)/2., (c1 + c2)/2.])

        alphas = np.array(alphas)
        Cs = np.array(Cs)

        print alphas.shape, RMSes[j, :].shape, Cs.shape
        np.save("/tmp/acf_h%d_%s_alphas.npy" % (j, typestr), alphas)
        np.save("/tmp/acf_h%d_%s_RMSes.npy" % (j, typestr), RMSes[j, :])
        np.save("/tmp/acf_h%d_%s_Cs.npy" % (j, typestr), Cs)

    if heights:
        np.save("/tmp/acf_heights.npy", np.array(heights))

def main():

    input_file = sys.argv[1]

    analyze(input_file, 1, "uniform")
    analyze(input_file, 2, "lattice")
    analyze(input_file, 3, "radial")
    analyze(input_file, 4, "pathfind")


if __name__ == "__main__":
    main()
