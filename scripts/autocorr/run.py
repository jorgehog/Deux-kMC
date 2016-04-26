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

    heights = ParameterSet(cfg, "confiningSurfaceHeight\s*=\s*(.*)\;")
    heights.initialize_set([20.])

    diffusions = ParameterSet(cfg, "diffuse\s*=\s*(.*)\;")
    diffusions.initialize_set([3])

    controller.register_parameter_set(alpha_values)
    controller.register_parameter_set(heights)
    controller.register_parameter_set(diffusions)

    controller.set_repeats(20)

    controller.run(run_kmc, path, app, cfg, ask=False, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
