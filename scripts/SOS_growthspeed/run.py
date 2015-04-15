import sys
import os
import numpy as np

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet, quick_replace


def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    quick_replace(cfg, "pressureWall", 1)

    E0_values = ParameterSet(cfg, "E0dA\s*\=\s*(.*)\;")
    E0_values.initialize_set(np.linspace(0.1, 1.0, 9))

    r0_values = ParameterSet(cfg, "r0\s*\=\s*(.*)\;")
    r0_values.initialize_set(np.linspace(0.5, 5, 10))

    s0_values = ParameterSet(cfg, "sigma0\s*\=\s*(.*)\;")
    s0_values.initialize_set(np.linspace(0.1, 2.0, 10))

    alpha_values = ParameterSet(cfg, "alpha\s*\=\s*(.*)\;")
    alpha_values.initialize_set(np.linspace(0.5, 3.0, 10))

    mu_shift_values = ParameterSet(cfg, "muShift\s*\=\s*(.*)\;")
    super_saturation = np.linspace(0.01, 4, 10)
    mu_shifts = np.log(super_saturation)
    mu_shift_values.initialize_set(mu_shifts)

    #For the unloaded system
    mu_values = ParameterSet(cfg, "mu\s*\=\s*(.*)\;")
    mu_values.initialize_set(mu_shifts)

    repeater_values = ParameterSet(cfg, "repeater\s*\=\s*(.*)\;")
    N = 1
    repeater_values.initialize_set_incr(0, N-1, 1)

    # controller.register_parameter_set(E0_values)
    # controller.register_parameter_set(r0_values)
    # controller.register_parameter_set(s0_values)
    # controller.register_parameter_set(alpha_values)
    # controller.register_parameter_set(mu_shift_values)
    controller.register_parameter_set(repeater_values)

    success = controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs)

    if success != 0:
        return

    controller.clear()

    quick_replace(cfg, "pressureWall", 0)

    # controller.register_parameter_set(alpha_values)
    # controller.register_parameter_set(mu_values)
    controller.register_parameter_set(repeater_values)

    print "Running unloaded system."
    controller.run(run_kmc, path, app, cfg, ask=False, n_procs=n_procs)


if __name__ == "__main__":
    main()