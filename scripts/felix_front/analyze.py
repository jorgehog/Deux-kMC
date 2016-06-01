import sys
import os
from os.path import join


sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5


def main():

    input_file = sys.argv[1]

    parser = ParseKMCHDF5(input_file)


    for data, L, W, run_id in parser:

        omega = data.attrs["omega"]
        conf_height = data.attrs["height"]

        pathname = "felix_o%.2f_h%d" % (omega, conf_height)

        dir = "/tmp/%s" % pathname
        if not os.path.exists(dir):
            os.mkdir(dir)

        stored_heights = data["stored_heights"]

        n_file = 0
        every = 100
        for hi, heights_id in enumerate(sorted(stored_heights, key=lambda x: int(x))):

            if hi % every != 0:
                continue

            heights = stored_heights[heights_id][()].transpose()

            bottom = heights.min()

            xyz = ""
            n = 0
            for x in range(L):
                for y in range(W):
                    h = heights[x, y]

                    for z in range(bottom, h+1):
                        xyz += "0 %d %d %d\n" % (x, y, z)
                        n += 1

                    xyz += "1 %d %d %g\n" % (x, y, conf_height)

                    n += 1

            xyz_file = "%d\n---\n%s" % (n, xyz)

            sys.stdout.flush()
            print "\rStored %s xyz %5d / %5d" % (dir, int(hi), len(stored_heights))

            with open("%s/surfaces%d.xyz" % (dir, n_file), 'w') as f:
                f.write(xyz_file)

            n_file += 1

        del stored_heights
        print "fin", dir

    print "fin"





if __name__ == "__main__":
    main()