import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), "..", ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, quick_replace


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    size_values = ParameterSet(cfg, "[LW]\s*\=\s*(.*)\;")
    size_values.initialize_set([5, 10, 20, 30])

    alpha_values = ParameterSet(cfg, "alpha = (.*)\;")
    alpha_values.initialize_set(np.linspace(1, 1.5, 3))

    n = 20

    controller.register_parameter_set(size_values)
    controller.register_parameter_set(alpha_values)
    controller.set_repeats(n)

    success = controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

    if success != 0:
        print "Simulation failed."

    return success


if __name__ == "__main__":
    main()
