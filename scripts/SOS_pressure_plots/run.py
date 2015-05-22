import sys
import os
import numpy as np
from os.path import join, split

import time

sys.path.append(join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, ParameterSetController


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    E0_values = ParameterSet(cfg, "E0dA\s*\=\s*(.*)\;")
    E0_values.initialize_set([0.01, 0.1, 1.0])

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")

    alpha_min = 0.1
    alpha_max = 2.0
    N = 20

    log_one_over_alpha = np.linspace(np.log(alpha_min), np.log(alpha_max), N)
    logspaces_alphas = np.exp(log_one_over_alpha)

    alpha_values.initialize_set(logspaces_alphas)

    controller.register_parameter_set(E0_values)
    controller.register_parameter_set(alpha_values)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs)


if __name__ == "__main__":
    main()