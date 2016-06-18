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

    omegas = []
    for data, L, W, run_id in parser:

        omega = data.attrs["omega"]

        if omega not in omegas:
            omegas.append(omega)

    omegas = sorted(omegas)

    X = np.linspace(-0.5, L - 0.5, nbins)
    Y = np.linspace(-0.5, W - 0.5, nbins)
    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    Hx = np.zeros(shape=(len(omegas), nbins))
    Hy = np.zeros_like(Hx)

    for data, L, W, run_id in parser:

        omega = data.attrs["omega"]
        alpha = data.attrs["alpha"]

        if alpha != 3:
            continue

        io = omegas.index(omega)
        conf_height = data.attrs["height"]

        pathname = "felix_o%.2f_h%d" % (omega, conf_height)

        dir = join(masterdir, pathname)
        xdir = join(dir, "x")
        ydir = join(dir, "y")

        if not os.path.exists(dir):
            os.mkdir(dir)
        if not os.path.exists(xdir):
            os.mkdir(xdir)
        if not os.path.exists(ydir):
            os.mkdir(ydir)

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

                if dh == 0:
                    print "error..."
                    continue

                Hx[io, int((x+0.5)/dx)] += dt/dh
                Hy[io, int((y+0.5)/dy)] += dt/dh

            if hi % every != 0:
                continue

            if hi != 0:
                np.save("%s/fcav_evo_%d.npy" % (xdir, int(hi/every) - 1), Hx[io]/time[hi])
                np.save("%s/fcav_evo_%d.npy" % (ydir, int(hi/every) - 1), Hy[io]/time[hi])

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

        Hx[io] /= T
        Hy[io] /= T

        print "fin", dir

    # np.save("/tmp/felix_cav_phist_omegas.npy", omegas)
    # np.save("/tmp/felix_cav_phist_X.npy", X)
    # np.save("/tmp/felix_cav_phist_H.npy", H)

    print "fin"





if __name__ == "__main__":
    main()
