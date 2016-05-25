import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    F0s = ParameterSet(cfg, ["F0\s*=\s*(.*)\;",
                             "nCycles\s*=\s*(.*)\;"])
    OS = ParameterSet(cfg, "omegaSign\s*=\s*(.*)\;")

    F0s.initialize_set([0.6, 1.0], [2000000, 30000])

    controller.register_parameter_set(F0s)
    controller.register_parameter_set(OS)

    # controller.set_repeats(10)

    print path
    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
