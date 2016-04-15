import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    omegas = ParameterSet(cfg, "omega\s*=\s*(.*)\;")
    omegas.initialize_set(np.linspace(0, 2, 9))

    heights = ParameterSet(cfg, "height\s*=\s*(.*)\;")
    heights.initialize_set([5.0, 10.0, 15.0, 20.0, 25.0, 30.0])

    controller.register_parameter_set(omegas)
    controller.register_parameter_set(heights)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
