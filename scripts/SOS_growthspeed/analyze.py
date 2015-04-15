import sys
import os
from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def main():

    parsed_data = {}

    input_file = sys.argv[1]

    if len(sys.argv) > 2:
        unloaded_input_file = sys.argv[2]

        print "unloaded data found."

        unloaded_parser = ParseKMCHDF5(unloaded_input_file)

        _data = []

        n = 0
        for chunk in unloaded_parser:
            _, _, _, _, mu, _, _, _, _, ignis_map, data, _ = chunk

            L = len(data["ignisData"][ignis_map["GrowthSpeed"]])

            v = data["ignisData"][ignis_map["GrowthSpeed"]][L/2:].mean()

            _data.append([mu, v])

        v_unloaded_array = np.array(sorted(_data, key=lambda x: x[0]))[:, 1]

        np.save("/tmp/growthspeed_unloaded.npy", v_unloaded_array)


    parser = ParseKMCHDF5(input_file)

    n = 0
    for stuff in parser:
        n += 1

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data, repeat = stuff

        if E0 == 0:
            continue

        area = L*W

        E0 /= area

        mu_shift = data.attrs["muShift"]

        if E0 not in parsed_data.keys():
            parsed_data[E0] = {}
        if mu_shift not in parsed_data[E0].keys():
            parsed_data[E0][mu_shift] = [0, 0, 0]

        pdata = parsed_data[E0][mu_shift]
        # print data
        L = len(data["ignisData"][ignis_map["GrowthSpeed"]])

        pdata[0] += data["ignisData"][ignis_map["GrowthSpeed"]][L/2:].mean()
        pdata[1] += mu
        pdata[2] += 1

    print "Parsed", n, "entries."

    E0_array = sorted(parsed_data.keys())
    mu_shift_array = sorted(parsed_data.values()[0].keys())
    v_array = []
    mu_array = []

    for E0, data in sorted(parsed_data.items(), key=lambda x: x[0]):
        v_array.append([])
        mu_array.append([])
        for mu_shift, data2 in sorted(data.items(), key=lambda x: x[0]):

            combzor = parsed_data[E0][mu_shift]

            v_array[-1].append(combzor[0]/combzor[2])
            mu_array[-1].append(combzor[1]/combzor[2])

    E0_array = np.asarray(E0_array)
    mu_shift_array = np.asarray(mu_shift_array)
    v_array = np.asarray(v_array)
    mu_array = np.asarray(mu_array)

    np.save("/tmp/growthspeed_E0.npy", E0_array)
    np.save("/tmp/growthspeed_mu_shift.npy", mu_shift_array)
    np.save("/tmp/growthspeed_mu.npy", mu_array)
    np.save("/tmp/growthspeed_v.npy", v_array)

if __name__ == "__main__":
    main()