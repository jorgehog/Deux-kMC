import sys
import os

from matplotlib.pylab import plot, show, figure

from os.path import join
import numpy as np

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def interpolate(x_start, x_end, y_start, y_end, x):
    slope = (y_end - y_start)/(x_end - x_start)

    return y_start + (x-x_start)*slope


def combine_results(all_x, all_y):

    n_runs = len(all_x)
    size = len(all_x[0])

    min_x = min(all_x, key=lambda z: z[0])[0]
    max_x = max(all_x, key=lambda z: z[-1])[-1]

    new_x = np.linspace(min_x, max_x, size)

    new_y = np.zeros_like(new_x)
    counts = np.zeros_like(new_x)

    for count, (x, y) in enumerate(zip(all_x, all_y)):

        print "%d / %d" % (count+1, n_runs)

        for i, xi in enumerate(x):

            closest_x_new = 0
            minimum_gap = 10000000
            loc = size+1
            for j, x_new_i in enumerate(new_x):
                gap = abs(x_new_i - xi)

                if gap < minimum_gap:
                    minimum_gap = gap
                    closest_x_new = x_new_i
                    loc = j

            if closest_x_new >= xi:

                start = i
                end = i + 1

                if end == size:
                    start = i - 1
                    end = i

            else:
                start = i-1
                end = i

                if start == -1:
                    start = 0
                    end = 1

            y_new_i = interpolate(x[start], x[end], y[start], y[end], closest_x_new)

            new_y[loc] += y_new_i
            counts[loc] += 1

    nonzero_slice = np.where(counts != 0)

    figure()
    plot(new_x, counts)

    new_x = new_x[nonzero_slice]
    new_y = new_y[nonzero_slice]/counts[nonzero_slice]


    return new_x, new_y

def testbed():

    all_x = []
    all_y = []

    n_runs = 100
    x_res = 1000
    x_max = 6.28

    figure()
    for i in range(n_runs):
        x = sorted(np.random.uniform(0, x_max, size=x_res))
        y = np.sin(x) + 0.1*(-1 + 2*np.random.uniform(size=x_res))

        all_x.append(x)
        all_y.append(y)

        plot(x, y)

    x, y = combine_results(all_x, all_y)

    figure()
    plot(x, y, 'rx', linewidth=3)
    show()

def main():

    #return testbed()

    parsed_data = {}

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    n = 0
    for stuff in parser:
        n += 1

        L, W, potential, alpha, mu, E0, s0, r0, neighbors, ignis_map, data, repeat = stuff

        print repeat

        area = L*W

        E0 /= area

        mu_shift = data.attrs["muShift"]

        if E0 not in parsed_data.keys():
            print ignis_map, n
            parsed_data[E0] = {}
        if alpha not in parsed_data[E0].keys():
            parsed_data[E0][alpha] = {}
        if mu_shift not in parsed_data[E0][alpha].keys():
            parsed_data[E0][alpha][mu_shift] = {"data": 0, "count": 0}

        parsed_data[E0][alpha][mu_shift]["data"] += np.array(data["ignisData"])
        parsed_data[E0][alpha][mu_shift]["count"] += 1

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

                data3["data"]/=data3["count"]

                np.save("/tmp/steadystate_data_%d.npy" % count, data3["data"])

                count += 1


    print "Parsed", n, "entries."

    E0_array = np.asarray(E0_array)
    alpha_array = np.asarray(alpha_array)
    mu_shift_array = np.asarray(mu_shift_array)

    np.save("/tmp/steadystate_E0.npy", E0_array)
    np.save("/tmp/steadystate_alpha.npy", alpha_array)
    np.save("/tmp/steadystate_mu_shift.npy", mu_shift_array)

if __name__ == "__main__":
    main()