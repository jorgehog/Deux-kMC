import ovito
from ovito.io import *
from ovito.vis import *
import os
import sys

print("Hello, this is OVITO %i.%i.%i" % ovito.version)

rs = RenderSettings(
    filename = '/tmp/image.png',
    size = (640,480),
    background_color = (1.0, 1.0, 1.0)
)
rs.renderer.antialiasing = True

def main():

    in_dir = sys.argv[1]
    out_dir = sys.argv[2]

    alphas = [0.5, 1.0, 2.0]

    frames = [[],
              [],
              []]

    with open(os.path.join(in_dir, "desc.txt"), 'r') as f:
        for line in f:
            n, alpha, F0, cov = line.split()

            frames[alphas.index(float(alpha))].append([float(F0), float(cov), int(n)])

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    node = ovito.dataset.selected_node

    for i in range(len(alphas)):

        conv = ""

        a_dir = os.path.join(out_dir, "a%1.f" % alphas[i])

        if not os.path.exists(a_dir):
            os.mkdir(a_dir)

        #sort by coverage
        a_frames = sorted(frames[i], key=lambda x: x[1])

        for frame, (F0, cov, n) in enumerate(a_frames):

            if n:
                xyz = os.path.join(in_dir, "all_xyz%d.xyz" % n)

                if not node:
                    node = import_file(xyz,
                                       columns=["Particle Type",
                                                "Position.X",
                                                "Position.Y",
                                                "Position.Z"])
                else:
                    node.source.load(xyz)
            else:
                node.source.load("null.xyz")

            rs.filename = os.path.join(a_dir, "cluster_grid_img_%d.png" % frame)
            ovito.dataset.viewports.active_vp.render(rs)

            conv += "%d %d %g %g\n" % (frame, n, F0, cov)

        with open(os.path.join(a_dir, "conv.txt"), 'w') as f:
            f.write(conv)

if __name__ == "__main__":
    main()
