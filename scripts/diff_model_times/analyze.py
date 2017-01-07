import sys
import os

from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():
    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    parsed_data = [{}, {}, {}, {}]

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        height = data.attrs["height"]

        diff = data.attrs["diffusionType"]
        cpu = data.attrs["cpuTime"]

        entry = parsed_data[diff-1]

        if alpha not in entry.keys():
            entry[alpha] = {}
        if height not in entry[alpha].keys():
            entry[alpha][height] = {}

        entry[alpha][height]["cpu"] = cpu

    alphas = sorted(entry.keys())
    heights = sorted(entry.values()[0].keys())

    np.save("/tmp/diff_model_times_alphas.npy", alphas)
    np.save("/tmp/diff_model_times_heights.npy", heights)

    results = np.zeros(shape=(4, len(alphas), len(heights), 1))

    for i in range(4):
        for j, alpha in enumerate(alphas):
            for k, height in enumerate(heights):
                for l, field in enumerate(["cpu"]):
                    results[i, j, k, l] = (parsed_data[i][alpha][height][field])

    np.save("/tmp/diff_model_times_results.npy", results)


if __name__ == "__main__":
    main()
