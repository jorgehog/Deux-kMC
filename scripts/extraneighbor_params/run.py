import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set(np.linspace(0.5, 2, 16))

    Pls = ParameterSet(cfg, "Pl\s*=\s*(.*)\;")
    Pls.initialize_set(np.linspace(0.01, 0.5, 20))

    controller.register_parameter_set(alpha_values)
    controller.register_parameter_set(Pls)

    # controller.set_repeats(10)

    print path
    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
