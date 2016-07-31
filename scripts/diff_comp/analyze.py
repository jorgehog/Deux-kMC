import sys
import os
import numpy as np
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    nbins = 20

    H = 0
    W = 0

    therm = 0.5

    for data, L, W, run_id in parser:

        dx = L/float(nbins)
        dy = W/float(nbins)

        stored_heights = data["stored_heights"]
        stored_particles = data["stored_particles"]
        time = data["time"]

        conf_height = data.attrs["height"]

        if conf_height != 5:
            continue

        t_prev = 0
        T = time[-1]
        t0 = None
        hmat = np.zeros(shape=(nbins, nbins))


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

            if particles is not None:
                for x, y, _ in particles:
                    xl = round(x)
                    yl = round(y)
                    dh = conf_height - heights[xl, yl] - 1

                    hmat[int((x+0.5)/dx), int((y+0.5)/dy)] += dt/dh

        del stored_heights

        H += hmat/(T - t0)
        W += 1

        print "fin", run_id

    H /= W

    np.save("/tmp/diff_comp_h.npy", H)

    print "fin"





if __name__ == "__main__":
    main()
