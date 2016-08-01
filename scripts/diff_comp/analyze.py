import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    nbins = 10
    L = 30
    W = 60

    conf_heights = [5, 10, 20]
    H = np.zeros(shape=(2, 3, nbins, nbins))
    Ws = np.zeros(shape=(2, 3))
    cH = np.zeros(shape=(2, 3, L, W))

    therm = 0.75

    for data, L, W, run_id in parser:

        dx = L/float(nbins)
        dy = W/float(nbins)

        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]
        time = data["time"]

        conf_height = data.attrs["height"]
        refl = data.attrs["reflecting"]

        ih = conf_heights.index(conf_height)

        t_prev = 0
        T = time[-1]
        t0 = None
        hmat = np.zeros(shape=(nbins, nbins))
        heightmat = np.zeros(shape=(L, W))


        for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):

            t_new = time[hi]

            if (hi > therm*len(stored_heights)) and (t0 is None):
                t0 = t_new

            dt = t_new - t_prev
            t_prev = t_new

            if heights_id in stored_particles:
                particles = stored_particles[heights_id][()]
            else:
                particles = None

            heights = stored_heights[heights_id][()].transpose()

            heightmat += heights*dt

            if particles is not None:
                for x, y, _ in particles:
                    xl = round(x)
                    yl = round(y)
                    dh = conf_height - heights[xl, yl] - 1

                    hmat[int((x+0.5)/dx), int((y+0.5)/dy)] += dt/(dh*dx*dy)

            if hi % 100 == 0:
                sys.stdout.flush()
                print "\r%d/%d" % (hi, len(stored_heights)),

        del stored_heights

        H[refl, ih] += hmat/(T - t0)
        cH[refl, ih] += heightmat/(T - t0)
        Ws[refl, ih] += 1

        print "fin", run_id

    for refl in [0, 1]:
        for ih in range(len(conf_heights)):
            H[refl, ih, :, :] /= Ws[refl, ih]
            cH[refl, ih, :, :] /= Ws[refl, ih]

    np.save("/tmp/diff_comp_C.npy", H)
    np.save("/tmp/diff_comp_H.npy", cH)
    np.save("/tmp/diff_comp_heights.npy", conf_heights)

    print "fin"





if __name__ == "__main__":
    main()
