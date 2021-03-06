import sys
import os
import numpy as np
from os.path import join
import matplotlib.pylab as plab

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
        NXYZ = 100

    parser = ParseKMCHDF5(input_file)

    masterdir = "/tmp/bunch"

    if not os.path.exists(masterdir):
        os.mkdir(masterdir)

    # nbins = 20

    repeats = 0
    nc = None
    conf_height = None
    for data, L, W, run_id in parser:

        # flux = data.attrs["flux"]

        repeats+=1

        if nc is None:
            nc = len(data["stored_heights"])
        if conf_height is None:
            conf_height = data.attrs["height"]

    # hmat = np.zeros(shape=(len(fluxes), nbins, nbins))
    # dx = L/float(nbins)
    # dy = W/float(nbins)

    # descs = ""
    these_flux = []

    frontposes = np.zeros(shape=(repeats, int(conf_height), nc))
    times = np.zeros(shape=(repeats, nc))

    repeater_count = 0

    for data, L, W, run_id in parser:

        flux = data.attrs["flux"]
        stored_heights = data["stored_heights"]
        # stored_particles = data["stored_particles"]
        conf_height = data.attrs["height"]

        every = max([1, len(stored_heights)/NXYZ])

        # io = fluxes.index(flux)

        if flux in these_flux:
            store=False
        else:
            store=True
            these_flux.append(flux)


        time = data["time"]
        times[repeater_count, :] = time
        t_prev = 0
        # T = time[-1]

        n_file = 0

        for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):

            # t_new = time[hi]
            # dt = t_new - t_prev
            # t_prev = t_new

            # if heights_id in stored_particles:
            #     particles = stored_particles[heights_id][()]
            # else:
            #     particles = None

            heights = stored_heights[heights_id][()].transpose()

            # if particles is not None:
            #     for x, y, _ in particles:
            #         xl = round(x)
            #         yl = round(y)
            #         dh = conf_height - heights[xl, yl] - 1
            #
            #         hmat[io, int((x+0.5)/dx), int((y+0.5)/dy)] += dt/dh
            #

            for hlevel in range(int(conf_height)):
                I = np.where(heights >= hlevel)

                frontposes[repeater_count, hlevel, hi] = I[0].size/float(L*W)

            if hi % every != 0 or not store:
                continue

            pathname = "felix_o%.2f_h%d" % (flux, conf_height)

            dir = join(masterdir, pathname)
            hdir = join(dir, "x")

            if not os.path.exists(dir):
                os.mkdir(dir)
            if not os.path.exists(hdir):
                os.mkdir(hdir)

            # if hi != 0:
            #     np.save("%s/fcav_evo_%d.npy" % (hdir, int(hi/every) - 1), hmat[io]/time[hi])

            particles = None
            store_xyz(heights, conf_height, n_file, dir, particles)
            n_file += 1

            sys.stdout.flush()
            print "\rStored %s xyz %5d / %5d" % (dir, hi/every, len(stored_heights)/every)
        repeater_count += 1

        del stored_heights

        # hmat[io] /= T

        print "fin", run_id

        # descs += "%s : %s %s\n" % (dir, parser.filename, run_id)

        # for hlevel in range(frontposes.shape[0]):
        #     plab.plot(time, frontposes[hlevel], label=str(hlevel))
        # plab.legend()
        # plab.show()
        # np.save("/tmp/flx_bunch_time.npy", time)
        # np.save("/tmp/flx_bunch_fpos.npy", frontposes)

    # desc_f = os.path.join(masterdir, "desc.txt")
    #
    # with open(desc_f, 'w') as f:
    #     f.write(descs)


    # np.save("/tmp/felix_cav_phist_omegas.npy", omegas)
    # np.save("/tmp/felix_cav_phist_X.npy", X)
    # np.save("/tmp/felix_cav_phist_H.npy", H)

    np.save("/tmp/flx_bunch_time.npy", times)
    np.save("/tmp/flx_bunch_fpos.npy", frontposes)

    print "fin"





if __name__ == "__main__":
    main()
