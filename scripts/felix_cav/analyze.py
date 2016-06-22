import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    if len(sys.argv) > 2:
        NXYZ = int(sys.argv[2])
    else:
        NXYZ = 1000

    parser = ParseKMCHDF5(input_file)

    if "cavity" in input_file:
        masterdir = "/tmp/cavity"
    else:
        masterdir = "/tmp/front"

    if not os.path.exists(masterdir):
        os.mkdir(masterdir)

    nbins = 20

    fluxes = []
    for data, L, W, run_id in parser:

        flux = data.attrs["flux"]

        if flux not in fluxes:
            fluxes.append(flux)

    fluxes = sorted(fluxes)

    hmat = np.zeros(shape=(len(fluxes), nbins, nbins))
    dx = L/nbins
    dy = W/nbins

    for data, L, W, run_id in parser:

        flux = data.attrs["flux"]

        io = fluxes.index(flux)

        conf_height = data.attrs["height"]

        pathname = "felix_o%.2f_h%d" % (flux, conf_height)

        dir = join(masterdir, pathname)
        hdir = join(dir, "x")

        if not os.path.exists(dir):
            os.mkdir(dir)
        if not os.path.exists(hdir):
            os.mkdir(hdir)

        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]

        n_file = 0
        every = max([1, len(stored_heights)/NXYZ])

        time = data["time"]
        t_prev = 0
        T = time[-1]

        for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):

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

                hmat[io, int((x+0.5)/dx), int((y+0.5)/dy)] += dt/dh

            if hi % every != 0:
                continue

            if hi != 0:
                np.save("%s/fcav_evo_%d.npy" % (hdir, int(hi/every) - 1), hmat[io]/time[hi])

            xyz = ""
            n = 0

            for x, y, z in particles:
                xyz += "2 %.3f %.3f %.3f\n" % (x, y, z)
                n += 1

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

            sys.stdout.flush()
            print "\rStored %s xyz %5d / %5d" % (dir, hi/every, len(stored_heights)/every)

            with open("%s/surfaces%d.xyz" % (dir, n_file), 'w') as f:
                f.write(xyz_file)

            n_file += 1

        del stored_heights

        hmat[io] /= T

        print "fin", dir

    # np.save("/tmp/felix_cav_phist_omegas.npy", omegas)
    # np.save("/tmp/felix_cav_phist_X.npy", X)
    # np.save("/tmp/felix_cav_phist_H.npy", H)

    print "fin"





if __name__ == "__main__":
    main()
