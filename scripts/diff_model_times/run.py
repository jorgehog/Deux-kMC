import sys
import os

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    alphas = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alphas.initialize_set([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])

    heights = ParameterSet(cfg, "height\s*=\s*(.*)\;")
    heights.initialize_set([5., 10., 15., 20.])

    diff_types = ParameterSet(cfg, "diffusionType\s*=\s*(.*)\;")
    diff_types.initialize_set([1, 2, 3, 4])

    controller.register_parameter_set(alphas)
    controller.register_parameter_set(heights)
    controller.register_parameter_set(diff_types)

    controller.set_repeats(10)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
