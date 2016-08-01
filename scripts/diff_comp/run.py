import sys
import os

sys.path.append(os.path.join(os.getcwd(), ".."))

from run_utils import run_kmc, parse_input
from ParameterJuggler import ParameterSet

def main():

    controller, path, app, cfg, n_procs = parse_input(sys.argv)

    refl = ParameterSet(cfg, "reflecting\s*=\s*(.*)\;")
    refl.initialize_set([0, 1])

    heights = ParameterSet(cfg, "height\s*=\s*(.*)\;")
    heights.initialize_set([5.0, 10.0, 20.0])

    controller.register_parameter_set(refl)
    controller.register_parameter_set(heights)

    controller.set_repeats(10)

    controller.run(run_kmc, path, app, cfg, ask=not controller.use_mpi, n_procs=n_procs, shuffle=True)

if __name__ == "__main__":
    main()
