import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    def skip(data):
        return data.attrs["flux"] != 1.25

    nbins = 40
    every = 1000

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

        time = data["time"][()]

        #Shift with 1 to translate from starting time till ending times
        t_tot += time[::every] + time[1]
        t_prev = 0

        for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):

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

        xyz = ""
        n = 0

        bottom = heights.min()

        for x in range(L):
            for y in range(W):
                h = heights[x, y]

                for z in range(bottom, h+1):
                    xyz += "0 %d %d %d\n" % (x, y, z)
                    n += 1

                xyz += "1 %d %d %g\n" % (x, y, conf_height)

                n += 1

        xyz_file = "%d\n---\n%s" % (n, xyz)

        with open("/tmp/surfaces%d.xyz" % entry_count, 'w') as f:
            f.write(xyz_file)

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
