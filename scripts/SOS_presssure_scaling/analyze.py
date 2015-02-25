import sys
import os
from os.path import join

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def main():
    parser = ParseKMCHDF5(sys.argv[1])

    for stuff in parser:
        print stuff

if __name__ == "__main__":
    main()