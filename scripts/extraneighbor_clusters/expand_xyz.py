import sys
import os

def main():
    L, W = [int(x) for x in sys.argv[1:3]]
    in_files = sys.argv[3:]

    dirs = [[0, 0], [1, 0], [0, 1], [1, 1]]

    for in_file in in_files:
        path, name = os.path.split(in_file)

        out_path = "%s/%s" % (path, "expanded")

        if not os.path.exists(out_path):
            os.mkdir(out_path)

        expanded_xyz = ""

        with open(in_file, 'r') as f:
            N = int(f.readline())
            f.readline()

            expanded_xyz += "%d\n-\n" % (4*N)

            for line in f:
                i, x, y, z = [float(x) for x in line.split()]

                xp = x + L
                yp = y + W

                xm = x - L
                ym = y - W

                for _dir in dirs:
                    expanded_xyz += "%d %g %g %g\n" % (i, x + _dir[0]*L, y + _dir[1]*W, z)

        out_file = "%s/expanded_%s" % (out_path, name)

        with open(out_file, 'w') as f:
            f.write(expanded_xyz)

if __name__ == "__main__":
    main()

