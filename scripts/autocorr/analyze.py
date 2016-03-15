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

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    RMSes = []
    Cs = []

    n = 0
    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]

        alphas.append(alpha)

        xl = np.linspace(-L/2, L/2, L + 1)
        xw = np.linspace(-W/2, W/2, W + 1)

        x1 = sqrt(2)*xl
        x2 = sqrt(2)*xw

        autocorr = np.array(data["autocorrelation"])
        autocorr /= autocorr.max()

        dl = autocorr[:, W/2]
        dw = autocorr[L/2, :]

        d1 = np.diag(autocorr)
        d2 = np.diag(np.flipud(autocorr))

        cl = find_corrl(xl, dl)
        cw = find_corrl(xw, dw)
        c1 = find_corrl(x1, d1)
        c2 = find_corrl(x2, d2)

        Cs.append([cl, cw, c1, c2])

        labels = data["ignisEventDescriptions"]
        i = list(list(labels)[0]).index("HeightRMS@MainMesh")
        RMS = data["ignisData"][i, :]
        rms = (RMS[len(RMS)/2:]).mean()

        RMSes.append(rms)

        n += 1

    print "Parsed", n, "entries."

    alphas = np.array(alphas)
    RMSes = np.array(RMSes)
    Cs = np.array(Cs)

    np.save("/tmp/acf_alphas.npy", alphas)
    np.save("/tmp/acf_RMSes.npy", RMSes)
    np.save("/tmp/acf_Cs.npy", Cs)


if __name__ == "__main__":
    main()