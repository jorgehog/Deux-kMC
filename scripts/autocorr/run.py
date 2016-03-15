import sys
import os

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, quick_replace

import numpy as np


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set([1.0])

    controller.register_parameter_set(alpha_values)

    controller.set_repeats(100)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
