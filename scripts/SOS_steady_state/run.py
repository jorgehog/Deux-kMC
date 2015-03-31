import sys
import os

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    E0_values = ParameterSet(cfg, "E0dA\s*\=\s*(.*)\;")
    E0_values.initialize_set([0.01, 0.1, 1.0])

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set([0.5, 1.0, 2.0])

    mu_shift_values = ParameterSet(cfg, "muShift\s*=\s*(.*)\;")
    mu_shift_values.initialize_set([-1., -0.5, 0., 0.5, 1.])

    controller.register_parameter_set(E0_values)
    controller.register_parameter_set(alpha_values)
    controller.register_parameter_set(mu_shift_values)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs)


if __name__ == "__main__":
    main()