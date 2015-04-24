import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def roundify(lim, *args):

    return (round(arg, lim) for arg in args)

def main():

    parsed_data = {}

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    E0_values = set([])
    alpha_values = set([])
    mu_shift_values = set([])
    r0_values = set([])
    s0_values = set([])

    lim = 5

    n = 0
    for stuff in parser:

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data, repeat = stuff

        if E0 == 0:
            mu_shift = mu
            s0 = 0
            r0 = 0
        else:
            mu_shift = data.attrs["muShift"]

        area = L*W
        E0 /= area

        E0, alpha, mu_shift, r0, s0 = roundify(lim, E0, alpha, mu_shift, r0, s0)

        if E0 not in parsed_data.keys():
            parsed_data[E0] = {}
            E0_values.add(E0)

        if alpha not in parsed_data[E0].keys():
            parsed_data[E0][alpha] = {}
            alpha_values.add(alpha)

        if mu_shift not in parsed_data[E0][alpha].keys():
            parsed_data[E0][alpha][mu_shift] = {}
            mu_shift_values.add(mu_shift)

        if r0 not in parsed_data[E0][alpha][mu_shift].keys():
            parsed_data[E0][alpha][mu_shift][r0] = {}
            r0_values.add(r0)

        if s0 not in parsed_data[E0][alpha][mu_shift][r0].keys():
            parsed_data[E0][alpha][mu_shift][r0][s0] = [0, 0, 0, 0]
            s0_values.add(s0)

        pdata = parsed_data[E0][alpha][mu_shift][r0][s0]

        pdata[0] += data.attrs["GrowthSpeed"]
        pdata[1] += mu
        pdata[2] += 1
        pdata[3] += data.attrs["nNeighbors"]

        n += 1

    print "Parsed", n, "entries."

    E0_values = sorted(list(E0_values))
    alpha_values = sorted(list(alpha_values))
    mu_shift_values = sorted(list(mu_shift_values))
    r0_values = sorted(list(r0_values))
    s0_values = sorted(list(s0_values))

    v_values = np.zeros(shape=(len(E0_values), len(alpha_values), len(mu_shift_values), len(r0_values), len(s0_values)))
    mu_values = np.zeros_like(v_values)
    n_values = np.zeros_like(mu_values)

    for i, E0 in enumerate(E0_values):
        for j, alpha in enumerate(alpha_values):
            for k, mu_shift in enumerate(mu_shift_values):
                for l, r0 in enumerate(r0_values):
                    for m, s0 in enumerate(s0_values):
                        E0, alpha, mu_shift, r0, s0 = roundify(lim, E0, alpha, mu_shift, r0, s0)

                        if E0 == 0:
                            s0 = 0
                            r0 = 0
                        elif s0 == 0 or r0 == 0:
                            continue

                        try:
                            pdata = parsed_data[E0][alpha][mu_shift][r0][s0]
                        except:
                            print "fail"
                            print E0, alpha, mu_shift, r0, s0
                            print parsed_data[E0][alpha][mu_shift].keys()
                            return
                        v_values[i][j][k][l][m] = pdata[0]/pdata[2]
                        mu_values[i][j][k][l][m] = pdata[1]/pdata[2]
                        n_values[i][j][k][l][m] = pdata[3]/pdata[2]

    E0_values = np.asarray(E0_values)
    alpha_values = np.asarray(alpha_values)
    mu_shift_values = np.asarray(mu_shift_values)
    r0_values = np.asarray(r0_values)
    s0_values = np.asarray(s0_values)
    v_values = np.asarray(v_values)
    n_values = np.asarray(n_values)
    mu_values = np.asarray(mu_values)

    np.save("/tmp/growthspeed_E0.npy", E0_values)
    np.save("/tmp/growthspeed_alpha.npy", alpha_values)
    np.save("/tmp/growthspeed_mu_shift.npy", mu_shift_values)
    np.save("/tmp/growthspeed_r0.npy", r0_values)
    np.save("/tmp/growthspeed_s0.npy", s0_values)
    np.save("/tmp/growthspeed_mu.npy", mu_values)
    np.save("/tmp/growthspeed_v.npy", v_values)
    np.save("/tmp/growthspeed_n.npy", n_values)

if __name__ == "__main__":
    main()