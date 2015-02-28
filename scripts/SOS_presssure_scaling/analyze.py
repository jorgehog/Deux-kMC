import sys
import os
from os.path import join
from scipy.stats import linregress

from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def get_coefficient_from_data(data):

    d_mu_d_alpha_values = []
    E0_values = []

    for E0, raw_data in data.items():

        d_mu_d_alpha, a, b, c, d = linregress(raw_data["alphas"],raw_data["muEqs"])

        d_mu_d_alpha_values.append(d_mu_d_alpha)
        E0_values.append(E0)

    coefficient, a, b, c, d = linregress(E0_values, d_mu_d_alpha_values)

    return 1./coefficient


def main():
    parsed_data = {}

    parser = ParseKMCHDF5(sys.argv[1])

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

    X, Y = np.meshgrid(s0_values, r0_values)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, coefficients, rstride=1, cstride=1, cmap=cm.coolwarm,
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