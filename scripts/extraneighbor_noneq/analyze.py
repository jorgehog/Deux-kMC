import os
import sys
import numpy as np
from os.path import join

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def main():

    input_file = sys.argv[1]
    parser = ParseKMCHDF5(input_file)

    F0s = []
    omegas = []

    for data, L, W, run_id in parser:

        os = data.attrs["omegaSign"]
        omegaVal = data.attrs["omegaVal"]
        omega = os*omegaVal

        if omega not in omegas:
            omegas.append(omega)

        F0 = round(data.attrs["F0"], 5)
        if F0 not in F0s:
            F0s.append(F0)


    F0s = sorted(F0s)
    omegas = sorted(omegas)

    for data, L, W, run_id in parser:

        scale = L*W

        os = data.attrs["omegaSign"]
        omegaVal = data.attrs["omegaVal"]
        omega = os*omegaVal
        io = omegas.index(omega)

        F0 = round(data.attrs["F0"], 5)
        iF0 = F0s.index(F0)

        eq_cov = data["eq_storedEventValues"][0]/scale
        neq_cov = data["omega_storedEventValues"][0]/scale

        eq_time = data["eq_storedEventValues"][1]
        neq_time = data["omega_storedEventValues"][1]

        _tuple = (iF0, io)
        np.save("/tmp/noneq_neigz_eqcov%d%d.npy" % _tuple, eq_cov)
        np.save("/tmp/noneq_neigz_eqtime%d%d.npy" % _tuple, eq_time)
        np.save("/tmp/noneq_neigz_neqcov%d%d.npy" % _tuple, neq_cov)
        np.save("/tmp/noneq_neigz_neqtime%d%d.npy" % _tuple, neq_time)

    np.save("/tmp/noneq_neigz_F0s.npy", F0s)
    np.save("/tmp/noneq_neigz_omegas.npy", omegas)

if __name__ == "__main__":
    main()