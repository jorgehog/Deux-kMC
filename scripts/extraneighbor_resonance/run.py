import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    F0s = ParameterSet(cfg, "F0\s*=\s*(.*)\;")
    F0s.initialize_set(np.linspace(0.5, 1.5, 21))

    controller.register_parameter_set(F0s)

    controller.set_repeats(10)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
