import os
import sys
import time

from ParameterJuggler import ParameterSetController, quick_replace

from mpi4py import MPI

def run_kmc(proc, combination, path, app, cfg):

    print "Running ",
    for value in combination:
        print "%.3f" % value,
    print

    this_dir = os.getcwd()

    os.chdir(path)
    success = os.system("./%s %d %s_%d.cfg >> /tmp/kmc_dump_%d.txt" % (app, proc, cfg.split(".cfg")[0], proc, proc))
    os.chdir(this_dir)

    return success

def parse_input(argv):

    if len(argv) < 2:
        print "Error: Usage: path_to_app -n n_procs -path output_path"
        sys.exit(1)

    path, app = os.path.split(argv.pop(1))

    if os.path.exists("/tmp/kmc_dump"):
        os.remove("/tmp/kmc_dump")


    cfg = os.path.join(path, "infiles", app + ".cfg")
    n_procs = 1
    out_path = None
    use_mpi = False

    if argv:
        for i, input in enumerate(argv):

            if "-n" == input:
                try:
                    n_procs = int(argv[i+1])
                except:
                    raise ValueError("invalid value succeeding -n")
            elif "-path" == input:
                try:
                    out_path = argv[i+1]
                except:
                    raise ValueError("invalid value succeeding -path")

            elif "-mpi" == input:
                use_mpi = True

                if MPI.COMM_WORLD.rank == 0:
                    print "Running", MPI.COMM_WORLD.size, "processes."

            elif "-cfg" == input:
                cfg = argv[i+1]

    controller = ParameterSetController(use_mpi=use_mpi)

    if out_path:
        quick_replace(cfg, "path", out_path)
    MPI.COMM_WORLD.Barrier()

    return controller, path, app, cfg, n_procs
