import sys
import os
from os.path import join, split

import time

sys.path.append(join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, ParameterSetController


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    E0_values = ParameterSet(cfg, "E0dA\s*\=\s*(.*)\;")
    E0_values.initialize_set_incr(0.05, 0.1, 0.05)

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set_incr(0.9, 1, 0.1)

    controller.register_parameter_set(E0_values)
    controller.register_parameter_set(alpha_values)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs)


if __name__ == "__main__":
    main()