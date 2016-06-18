import ovito
from ovito.io import *
from ovito.vis import *
import os
import glob
import sys
import re
from os.path import join

print("Hello, this is OVITO %i.%i.%i" % ovito.version)

rs = RenderSettings(
    filename = '/tmp/image.png',
    size = (320,240),
    background_color = (1.0, 1.0, 1.0)
)
rs.renderer.antialiasing = False

def makemovie(moviedir, path):
    files = sorted(glob.glob(path + "/*.xyz"), key= lambda x: int(re.findall("surfaces(\d+)\.xyz", x)[0]))
    name = path.split("/")[-1] + ".gif"

    imgs = []
    node = ovito.dataset.selected_node
    for frame, f in enumerate(files):
        if not node:
            node = import_file(f, columns=["Particle Type", "Position.X", "Position.Y", "Position.Z"])
        else:
            node.source.load(f)

        rs.filename = "%s/image%s.png" % (path, str(frame).rjust(4, "0"))
        ovito.dataset.viewports.active_vp.render(rs)

        imgs.append(rs.filename)

        sys.stdout.flush()
        print("\rRendered img %4d/%d : %s from %s" % (frame+1, len(files),
                                                      rs.filename.split("/")[-1],
                                                      f.split("/")[-1]))

    os.system("convert -delay 3 %s -loop 1 %s" % (" ".join(imgs), join(moviedir, name)))

    for img in imgs:
        os.remove(img)


def main():
    superdir = sys.argv[1]
    dirs = glob.glob(join(superdir, "felix_*"))
    moviedir = join(superdir, "movies")

    if not os.path.exists(moviedir):
        os.mkdir(moviedir)

    for i, dir in enumerate(dirs):
        print("Entering %s : %d/%d" % (dir, i+1, len(dirs)))
        makemovie(moviedir, dir)


if __name__ == "__main__":
    main()