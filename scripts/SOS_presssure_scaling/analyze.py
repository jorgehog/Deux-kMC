import sys
import os
import re
from os.path import join
from scipy.stats import linregress

from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

mean_err1 = 0
mean_err2 = 0
c1 = 0
c2 = 0
max_err1 = 0
max_err2 = 0


def get_coefficient_from_data(data):

    global mean_err1
    global mean_err2
    global c1
    global c2
    global max_err1
    global max_err2

    d_mu_d_alpha_values = []
    E0_values = []

    for E0, raw_data in data.items():

        d_mu_d_alpha, a, r1, c, err1 = linregress(raw_data["alphas"], raw_data["muEqs"])

        d_mu_d_alpha_values.append(d_mu_d_alpha)
        E0_values.append(E0)

        mean_err1 += err1
        c1 += 1

        if err1 > max_err1:
            max_err1 = err1

    coefficient, a, r2, c, err2 = linregress(E0_values, d_mu_d_alpha_values)

    mean_err2 += err2
    c2 += 1

    if err2 > max_err2:
        max_err2 = err2

    return 1./coefficient


def main():
    parsed_data = {}

    input_file = sys.argv[1]
    dir_name = os.path.basename(os.path.dirname(input_file))

    parser = ParseKMCHDF5(input_file)

    n = 0
    for stuff in parser:

        data = stuff[-1]

        L, W = stuff[:2]
        area = L*W

        alpha, mu, E0, s0, r0, neighbors = stuff[3:-2]

        E0 /= area

        muEq = data.attrs["muEq"]
        muEqError = data.attrs["muEqError"]

        if muEq != muEq:
            print "NAN", muEq, alpha, r0, s0
            import time
            time.sleep(0.1)
            continue

        if s0 not in parsed_data.keys():
            parsed_data[s0] = {}
        if r0 not in parsed_data[s0].keys():
            parsed_data[s0][r0] = {}
        if E0 not in parsed_data[s0][r0].keys():
            parsed_data[s0][r0][E0] = {"alphas": [],
                                       "muEqs": [],
                                       "muEqErrors": [],
                                       "neighbors": []}


        parsed_data[s0][r0][E0]["alphas"].append(alpha)
        parsed_data[s0][r0][E0]["muEqs"].append(muEq)
        parsed_data[s0][r0][E0]["muEqErrors"].append(muEqError)
        parsed_data[s0][r0][E0]["neighbors"].append(neighbors)


    s0_values = np.array(sorted(parsed_data.keys()))
    r0_values = np.array(sorted(parsed_data[s0_values[0]].keys()))
    coefficients = np.zeros(shape=[len(s0_values), len(r0_values)])

    for i, s0 in enumerate(s0_values):

        for j, r0 in enumerate(r0_values):

            coefficients[i][j] = get_coefficient_from_data(parsed_data[s0][r0])

            print s0, r0, coefficients[i][j]

    print "linearization errors (alpha mu) ( E0 ) : ", mean_err1/c1, mean_err2/c2
    print "max: ", max_err1, max_err2

    np.save("/tmp/analyze_%s_s0_values.npy" % dir_name, s0_values)
    np.save("/tmp/analyze_%s_r0_values.npy" % dir_name, r0_values)
    np.save("/tmp/analyze_%s_C_values.npy" % dir_name, coefficients)

    X, Y = np.meshgrid(s0_values, r0_values)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, coefficients.transpose(), rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    plt.xlabel(r"$\sigma_0$")
    plt.ylabel(r"$r_0$")

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.figure()

    plt.contour(coefficients)

    plt.xlabel(r"$\sigma_0$")
    plt.ylabel(r"$r_0$")


    plt.show()



if __name__ == "__main__":
    main()