import os
import sys
import inspect
import time

this_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

sys.path.append(os.path.join(this_dir, "..", "utils", "ParameterJuggler"))
print os.path.join(os.getcwd(), "..", "utils", "ParameterJuggler")
from ParameterJuggler import quick_replace


def run_kmc(proc, combination, path, app):
    time.sleep(proc/10.)

    print "Running ",
    for value in combination:
        print "%.3f" % value,
    print

    this_dir = os.getcwd()

    os.chdir(path)
    success = os.system("./%s %d >> /tmp/kmc_dump_%d.txt" % (app, proc, proc))
    os.chdir(this_dir)
 
    return success

def parse_input(argv):

    if len(argv) < 2:
        print "Error: Usage: path_to_app -n n_procs -path output_path"
        sys.exit(1)

    path, app = os.path.split(argv.pop(1))

    if os.path.exists("/tmp/kmc_dump"):
        remove("/tmp/kmc_dump")

    cfg = os.path.join(path, "infiles", app + ".cfg")

    n_procs = 1
    out_path = None

    if argv:
        for i, input in enumerate(argv):

            if "-n" in input:
                try:
                    n_procs = int(argv[i+1])
                except:
                    raise ValueError("invalid value succeeding -n")
            elif "-path" in input:
                if "-n" in input:
                    try:
                        out_path = int(argv[i+1])
                    except:
                        raise ValueError("invalid value succeeding -path")

    if out_path:
        quick_replace(cfg, "path", out_path)

    return path, app, cfg, n_procs
