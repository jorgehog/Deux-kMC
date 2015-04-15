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

    parser = ParseKMCHDF5(input_file)

    data_keys = ["Time", "HeightRMS", "SurfaceSize", "PressureWall"]

    n = 0
    for stuff in parser:
        n += 1

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data, repeat = stuff

        area = L*W

        E0 /= area

        mu_shift = data.attrs["muShift"]

        if E0 not in parsed_data.keys():
            parsed_data[E0] = {}
        if alpha not in parsed_data[E0].keys():
            parsed_data[E0][alpha] = {}
        if mu_shift not in parsed_data[E0][alpha].keys():
            parsed_data[E0][alpha][mu_shift] = ICZ(*data_keys)

        for key in data_keys:
            parsed_data[E0][alpha][mu_shift][key].append(data["ignisData"][ignis_map[key]])

    print "Parsed", n, "entries."

    E0_array = []
    alpha_array = []
    mu_shift_array = []

    count = 0

    for i, (E0, data) in enumerate(sorted(parsed_data.items(), key=lambda x: x[0])):

        E0_array.append(E0)

        for j, (alpha, data2) in enumerate(sorted(data.items(), key=lambda x: x[0])):

            if i == 0:
                alpha_array.append(alpha)

            for mu_shift, data3 in sorted(data2.items(), key=lambda x: x[0]):

                if i == 0 and j == 0:
                    mu_shift_array.append(mu_shift)

                combzor = parsed_data[E0][alpha][mu_shift]
                print E0, alpha, mu_shift, "%d/%d" % (count+1,  n/len(combzor))

                for key in data_keys:
                    if key == "Time":
                        continue

                    t, measure = combzor.intercombine("Time", key)
                    np.save("/tmp/steadystate_%s_%d.npy" % (key, count), measure)

                np.save("/tmp/steadystate_Time_%d.npy" % count, t)

                count += 1

    E0_array = np.asarray(E0_array)
    alpha_array = np.asarray(alpha_array)
    mu_shift_array = np.asarray(mu_shift_array)

    np.save("/tmp/steadystate_E0.npy", E0_array)
    np.save("/tmp/steadystate_alpha.npy", alpha_array)
    np.save("/tmp/steadystate_mu_shift.npy", mu_shift_array)

if __name__ == "__main__":
    main()