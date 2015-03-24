from mercurial.ignore import ignore
import sys
import os
import re
import random

from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5



def main():
    parsed_data = {}

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    n = 0
    for stuff in parser:
        n += 1

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data = stuff

        area = L*W

        E0 /= area

        muEq = data.attrs["muEq"]

        if muEq != muEq:
            print "NAN", muEq, alpha, r0, s0
            import time
            time.sleep(0.1)
            continue

        if E0 not in parsed_data.keys():
            parsed_data[E0] = {"alphas": [],
                               "mean_s": [],
                               "var_s": []}

        mean_s = data.attrs["size"]
        var_s = data.attrs["var"]

        parsed_data[E0]["alphas"].append(alpha)
        parsed_data[E0]["mean_s"].append(mean_s)
        parsed_data[E0]["var_s"].append(var_s)

    E0_array = []
    alpha_array = []
    mean_s_array = []
    var_s_array = []

    for E0, data in parsed_data.items():

        E0_array.append(E0)
        alpha_array.append(data["alphas"])
        mean_s_array.append(data["mean_s"])
        var_s_array.append(data["var_s"])

    print "Parsed", n, "entries."

    E0_array = np.asarray(E0_array)
    alpha_array = np.asarray(alpha_array)
    mean_s_array = np.asarray(mean_s_array)
    var_s_array = np.asarray(var_s_array)

    np.save("/tmp/pressure_plots_E0.npy", E0_array)
    np.save("/tmp/pressure_plots_alphas.npy", alpha_array)
    np.save("/tmp/pressure_plots_mean_s.npy", mean_s_array)
    np.save("/tmp/pressure_plots_var_s.npy", var_s_array)


if __name__ == "__main__":
    main()