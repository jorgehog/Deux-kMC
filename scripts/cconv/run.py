import sys
import os

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, quick_replace

import numpy as np


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set(np.linspace(0.5, 2., 7))

    h_values = ParameterSet(cfg, "h0\s*=\s*(.*)\;")
    h_values.initialize_set([5., 10., 20., 30.])

    quick_replace(cfg, "type", 0)

    controller.register_parameter_set(alpha_values)
    controller.register_parameter_set(h_values)

    success = controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs)

    if success != 0:
        return

    controller.clear()

    quick_replace(cfg, "type", 1)
    controller.register_parameter_set(alpha_values)
    controller.register_parameter_set(h_values)

    return controller.run(run_kmc, path, app, cfg, ask=False, n_procs=n_procs)


if __name__ == "__main__":
    main()