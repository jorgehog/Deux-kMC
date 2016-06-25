import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def make_xyz(dirname, heights, idx, conf_height):
    n_xyz = 100
    l = len(heights)
    every = l/n_xyz

    for i, key in enumerate(idx):

        if i % every != 0:
            continue

        hi = heights[key][()].transpose()

        L, W = hi.shape
        bottom = hi.min()
        xyz = ""
        n = 0

        for x in range(L):
            for y in range(W):
                h = hi[x, y]

                for z in range(bottom, h+1):
                    xyz += "0 %d %d %d\n" % (x, y, z)
                    n += 1

                xyz += "1 %d %d %g\n" % (x, y, conf_height)

                n += 1

        xyz_file = "%d\n---\n%s" % (n, xyz)

        with open("%s/surfaces%d.xyz" % (dirname, i/every), 'w') as f:
            f.write(xyz_file)


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    def skip(data):
        return data.attrs["flux"] != 2.00

    nbins = 30
    every = 1

    l = None
    n_entries = 0
    for data, L, W, run_id in parser:

        if skip(data):
            continue

        if not l:
            l = len(data["time"])
        n_entries += 1

    cmat = np.zeros(shape=(l/every, nbins))
    hmat = np.zeros(shape=(l/every, nbins))
    dy = W/float(nbins)

    t_tot = 0
    entry_count = 0
    for data, L, W, run_id in parser:

        if skip(data):
            continue

        conf_height = data.attrs["height"]

        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]

        stored_heights_indices = sorted(stored_heights, key=lambda x: int(x))

        time = data["time"][()]

        #Shift with 1 to translate from starting time till ending times
        t_tot += time[::every] + time[1]
        t_prev = 0

        if entry_count == 0:
            xyz_dir = "/tmp/first_cav_front"
            if not os.path.exists(xyz_dir):
                os.mkdir(xyz_dir)
            make_xyz(xyz_dir, stored_heights, stored_heights_indices, conf_height)

        for hi, heights_id in enumerate(stored_heights):

            if hi % every != 0:
                continue

            t_new = time[hi]
            dt = t_new - t_prev
            t_prev = t_new

            if heights_id in stored_particles:
                particles = stored_particles[heights_id][()]
            else:
                particles = []

            heights = stored_heights[heights_id][()].transpose()

            for x, y, _ in particles:
                xl = round(x)
                yl = round(y)
                dh = conf_height - heights[xl, yl] - 1

                cmat[hi/every, int((y+0.5)/dy)] += dt/dh

            for x in range(L):
                for y in range(W):
                    height = heights[x, y]
                    hmat[hi/every, y/dy] += dt*height

            if hi % 100 == 0:
                sys.stdout.flush()
                print "\r%d/%d" % (hi, len(stored_heights)),

        entry_count += 1

        sys.stdout.flush()
        print "\rfin %d / %d" % (entry_count, n_entries)
    print

    for i, ti in enumerate(t_tot):
        cmat[i] /= ti
        hmat[i] /= ti

    t_avg = t_tot/entry_count

    np.save("/tmp/FelixSeqC_t.npy", t_avg)
    np.save("/tmp/FelixSeqC_c.npy", cmat)
    np.save("/tmp/FelixSeqC_h.npy", hmat)

if __name__ == "__main__":
    main()
