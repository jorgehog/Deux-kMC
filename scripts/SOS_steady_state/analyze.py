import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

from intercombinatorzor import ICZ


def main():

    parsed_data = {}

    input_file = sys.argv[1]
    # N = int(sys.argv[2])

    parser = ParseKMCHDF5(input_file)

    data_keys = ["Time", "HeightRMS", "SurfaceSize", "PressureWall"]



    E0_array = []
    alpha_array = []
    mu_shift_array = []

    # parser.skipWhen = lambda p_alpha, p_mu, p_E0, p_s0, p_r0, p_n: p_n > 9

    n = 0
    for stuff in parser:

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data, repeat = stuff

        area = L*W

        E0 /= area

        mu_shift = data.attrs["muShift"]

        if E0 not in E0_array:
            print "n=%d located E0=%g" % (n, E0)
            E0_array.append(E0)

        if E0 not in parsed_data:
            parsed_data[E0] = {}

        if alpha not in alpha_array:
            print "n=%d located alpha=%g" % (n, alpha)
            alpha_array.append(alpha)

        if alpha not in parsed_data[E0]:
            parsed_data[E0][alpha] = {}

        if mu_shift not in mu_shift_array:
            print "n=%d located muShift=%g" % (n, mu_shift)
            mu_shift_array.append(mu_shift)

        if mu_shift not in parsed_data[E0][alpha]:
            parsed_data[E0][alpha][mu_shift] = []

        parsed_data[E0][alpha][mu_shift].append([parser.file_number, data.name])

        n += 1

        if n % 1000 == 0:
            print "n =", n

    E0_array = sorted(E0_array)
    mu_shift_array = sorted(mu_shift_array)
    alpha_array = sorted(alpha_array)

    E0_array = np.asarray(E0_array)
    alpha_array = np.asarray(alpha_array)
    mu_shift_array = np.asarray(mu_shift_array)


    N = len(parsed_data[E0][alpha][mu_shift])
    print "N =", N

    np.save("/tmp/steadystate_E0.npy", E0_array)
    np.save("/tmp/steadystate_alpha.npy", alpha_array)
    np.save("/tmp/steadystate_mu_shift.npy", mu_shift_array)

    combined_data = []

    combinator = ICZ(*data_keys)
    count = 0

    eqMu = np.zeros(shape=(len(E0_array), len(alpha_array)))

    for i, E0_match in enumerate(E0_array):
        for j, alpha_match in enumerate(alpha_array):

            nEq = 0
            for mu_shift_match in mu_shift_array:

                for file_number, name in parsed_data[E0_match][alpha_match][mu_shift_match]:
                    data = parser.get_data(file_number, name)

                    try:
                        eqmu_local = data.attrs["muEq"]
                        eqMu[i, j] += eqmu_local
                    except KeyError:
                        print data.attrs.keys()
                        print "FAIL", E0_match, alpha_match, mu_shift_match
                        return

                    nEq += 1

                    for key in data_keys:
                        # if key == "Time":
                        #     combinator[key].append(data["ignisData"][ignis_map[key]]/np.exp(eqmu_local + mu_shift_match - 2*alpha))
                        # else:
                        combinator[key].append(data["ignisData"][ignis_map[key]])

                    n += 1

                    if n % 100 == 0:
                        sys.stdout.flush()
                        print "\rn=", n,

                for key in data_keys:
                    if key == "Time":
                        continue

                    t, measure = combinator.intercombine("Time", key)
                    combined_data.append(["/tmp/steadystate_%s_%d.npy" % (key, count), measure])

                combined_data.append(["/tmp/steadystate_Time_%d.npy" % count, t])

                combinator.clear()

                count += 1

                print
                print "progress: ", float(count) / (len(E0_array)*len(alpha_array)*len(mu_shift_array))*100, "%"

        eqMu[i, j] /= nEq

    np.save("/tmp/steadystate_eqMu.npy", eqMu)

    for name, data in combined_data:
        np.save(name, data)


if __name__ == "__main__":
    main()