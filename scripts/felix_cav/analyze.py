import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def store_xyz(heights, conf_height, n_file, dir, particles=None):
    xyz = ""
    n = 0

    if particles is not None:
        for x, y, z in particles:
            xyz += "2 %.3f %.3f %.3f\n" % (x, y, z)
            n += 1

    L, W = heights.shape
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

    with open("%s/surfaces%d.xyz" % (dir, n_file), 'w') as f:
        f.write(xyz_file)


def main():

    input_file = sys.argv[1]

    if len(sys.argv) > 2:
        NXYZ = int(sys.argv[2])
    else:
        NXYZ = 1000

    whichDir = None
    whichFlux = None
    if len(sys.argv) > 3:
        whichFlux = float(sys.argv[3])
        whichDir = sys.argv[4]

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
    dx = L/float(nbins)
    dy = W/float(nbins)

    descs = ""
    these_flux = []

    for data, L, W, run_id in parser:

        flux = data.attrs["flux"]
        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]
        conf_height = data.attrs["height"]

        if flux in these_flux:
            continue

        these_flux.append(flux)

        every = max([1, len(stored_heights)/NXYZ])
        n_file = 0

        if whichFlux:
            if flux == whichFlux:

                os.system("mkdir -p %s" % whichDir)
                print parser.filename, run_id
                for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):
                    if hi % every != 0:
                        continue

                    if heights_id in stored_particles:
                        particles = stored_particles[heights_id][()]
                    else:
                        particles = None

                    heights = stored_heights[heights_id][()].transpose()
                    store_xyz(heights, conf_height, n_file, whichDir, particles)

                    n_file += 1
                return
            else:
                continue

        io = fluxes.index(flux)

        pathname = "felix_o%.2f_h%d" % (flux, conf_height)

        dir = join(masterdir, pathname)
        hdir = join(dir, "x")

        if not os.path.exists(dir):
            os.mkdir(dir)
        if not os.path.exists(hdir):
            os.mkdir(hdir)


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
                particles = None

            heights = stored_heights[heights_id][()].transpose()

            if particles is not None:
                for x, y, _ in particles:
                    xl = round(x)
                    yl = round(y)
                    dh = conf_height - heights[xl, yl] - 1

                    hmat[io, int((x+0.5)/dx), int((y+0.5)/dy)] += dt/dh

            if hi % every != 0:
                continue

            if hi != 0:
                np.save("%s/fcav_evo_%d.npy" % (hdir, int(hi/every) - 1), hmat[io]/time[hi])

            store_xyz(heights, conf_height, n_file, dir, particles)
            n_file += 1

            sys.stdout.flush()
            print "\rStored %s xyz %5d / %5d" % (dir, hi/every, len(stored_heights)/every)

        del stored_heights

        hmat[io] /= T

        print "fin", dir

        descs += "%s : %s %s\n" % (dir, parser.filename, run_id)

    desc_f = os.path.join(masterdir, "desc.txt")

    with open(desc_f, 'w') as f:
        f.write(descs)

    # np.save("/tmp/felix_cav_phist_omegas.npy", omegas)
    # np.save("/tmp/felix_cav_phist_X.npy", X)
    # np.save("/tmp/felix_cav_phist_H.npy", H)

    print "fin"





if __name__ == "__main__":
    main()
