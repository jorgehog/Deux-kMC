import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), "..", ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    c_values = ParameterSet(cfg, "c\s*\=\s*(.*)\;")
    csinv = np.linspace(1, 100, 40)
    c_values.initialize_set(1/csinv)
    # supersaturation_values.initialize_set([0.])

    n = 10

    controller.register_parameter_set(c_values)
    controller.set_repeats(n)

    success = controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

    if success != 0:
        print "Simulation failed."

    return success


if __name__ == "__main__":
    main()
