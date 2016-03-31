from matplotlib.pylab import imshow, show
import sys
import os
from os.path import join
import numpy as np
from math import floor

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)

    alphas = []
    Pls = []

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)

        Pl = data.attrs["Pl"]
        Pl = round(Pl, 3)

        if alpha not in alphas:
            alphas.append(alpha)

        if Pl not in Pls:
            Pls.append(Pl)

    cmat = np.zeros(shape=(len(alphas), len(Pls)))

    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        ia = alphas.index(alpha)

        Pl = data.attrs["Pl"]
        Pl = round(Pl, 3)
        ipl = Pls.index(Pl)

        height = data.attrs["h"]
        heights = data["heights"]

        # print np.array(heights).max(), floor(height)
        c = len(np.where(np.array(heights) > floor(height) - 2)[0])

        if c > 10:
            cmat[ia, ipl] = 1

    imshow(cmat)
    show()

    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_Pls.npy", Pls)
    np.save("/tmp/extraneighbor_cmat.npy", cmat)


if __name__ == "__main__":
    main()