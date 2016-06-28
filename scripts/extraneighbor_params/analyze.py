import matplotlib.pylab as plab
import sys
import os
from os.path import join
import numpy as np
from skimage import measure

sys.path.append(join(os.getcwd(), ".."))

from parse_h5_output import ParseKMCHDF5

def periodic_image_analysis(pim, labels):
    L, W = pim.shape
    expanded = np.zeros(shape=[2*d for d in pim.shape], dtype=int)
    for x in range(2):
        for y in range(2):
            expanded[x*L:(x+1)*L, y*W:(y+1)*W] = pim

    n_labels = len(labels)

    idx_map = {}

    n = 0
    for l in labels:
        idx_map[l] = n
        n += 1

    winners = dict(zip(labels, [[] for _ in range(n_labels)]))
    winner_areas = [0 for _ in range(n_labels)]

    l = measure.label(expanded)
    l[np.where(expanded == 0)] = -1
    l += 1

    def find_fast(im, v, x0, y0, x1, y1):
        for _x in range(x0, x1):
            for _y in range(y0, y1):
                if im[_x, _y] == v:
                    return _x, _y

    orig_labels = []
    expanded_props = measure.regionprops(l)
    for prop in expanded_props:

        x, y = find_fast(l, prop.label, *prop.bbox)
        orig_label = expanded[x, y]
        orig_labels.append(orig_label)

        idx = idx_map[orig_label]

        if prop.area > winner_areas[idx]:
            winner_areas[idx] = prop.area

    for orig_label, prop in zip(orig_labels, expanded_props):

        idx = idx_map[orig_label]

        if prop.area == winner_areas[idx]:
            winners[orig_label].append([prop.perimeter, prop.centroid])

    return winners

def changefill(im, source, dest):
    if dest == -1 or source == -1:
        return

    if dest != source:
        I = np.where(im == dest)
        im[I] = source


def periodify(im):
    l, w = im.shape

    for i in range(l):
        source = im[i, 0]
        dest = im[i, w-1]

        changefill(im, source, dest)

    for j in range(w):
        source = im[0, j]
        dest = im[l-1, j]

        changefill(im, source, dest)

    return im


def legal(coverage, x, y):
    if not coverage[x, y] == 1:
        return False

    L, W = coverage.shape

    for dx in [-1, 0, 1]:
        xn = (x + dx + L) % L
        for dy in [-1, 0, 1]:
            yn = (y + dy + W) % W

            if coverage[xn, yn] == 0:
                return True

    return False


def get_p_value(heights, id):

    # print id

    #if id != '151465468247664':
    #    return 0

   # plab.figure(0)
   # plab.pcolor(heights.transpose())

    L, W = heights.shape
    coverage = (heights == heights.max())

    l = measure.label(coverage, connectivity=1)
    l[np.where(coverage == 0)] = -1

    l_per = periodify(l) + 1
    props = measure.regionprops(l_per)

    max_area = 0
    winning_label = None
    for prop in props:
        if prop.area > max_area:
            max_area = prop.area
            winning_label = prop.label

    for prop in props:
        if prop.label != winning_label:
            for x, y in prop.coords:
                coverage[x, y] = 0
       # else:
        #    for x, y in prop.coords:
          #      plab.scatter(x + 0.5, y + 0.5, color='w', marker='x')
    area = max_area
    start_x = None
    start_y = None

    for prop in props:
        if prop.label == winning_label:
            sheds = None
            while sheds != 0:
                sheds = 0
                for x, y in prop.coords:
                    if coverage[x, y] == 0:
                        continue

                    xp = (x + 1) % L
                    xm = (x - 1 + L) % L
                    yp = (y + 1) % W
                    ym = (y - 1 + W) % W

                    n = 0
                    if coverage[xp, y] == 1:
                        n += 1
                    if coverage[xm, y] == 1:
                        n += 1
                    if coverage[x, yp] == 1:
                        n += 1
                    if coverage[x, ym] == 1:
                        n += 1

                    if n == 1:
                        coverage[x, y] = 0
                        sheds += 1
                        area -= 1
                       # plab.scatter(x + 0.5, y + 0.5, color='w', marker='s')

    if area < 9:
        return 0.

    pool = np.zeros_like(coverage) - 1
    pool[np.where(coverage == 0)] = 1
    lz = measure.label(pool, connectivity=1)
    lz_per = periodify(lz) + 1
    lz_props = measure.regionprops(lz_per)

    max_area = 0
    winning_label_2 = None
    for prop in lz_props:
        if prop.area > max_area:
            max_area = prop.area
            winning_label_2 = prop.label

    for prop in lz_props:
        if prop.label != winning_label_2:
            for x, y in prop.coords:
                coverage[x, y] = 1
              #  plab.scatter(x + 0.5, y + 0.5, color='w', marker='o')

    for prop in props:
        if prop.label == winning_label:
            for x, y in prop.coords:
                if legal(coverage, x, y):
                    start_x, start_y = x, y
                   #plab.scatter(x + 0.5, y + 0.5, color='w', marker='d', s=20)
                    break

    if start_x is None:
        print "ERROR: No start point.."
        #plab.show()
        exit(1)
    ori_to_idx = {
         0: {
             1: 1,
            -1: 3
         },
         1: {
            0:  0
         },
        -1: {
            0:  2
         }
    }

    idx_to_ori = {
        0: [ 1,  0],
        1: [ 0,  1],
        2: [-1,  0],
        3: [ 0, -1]
    }


    #0: right
    #1: up
    #2: left
    #3: down
    right_trans = [
        [ 0, -1],
        [ 1,  0],
        [ 0,  1],
        [-1,  0]
    ]

    forward_trans = [
        [ 1,  0],
        [ 0,  1],
        [-1,  0],
        [ 0, -1]
    ]

    left_trans = [
        [ 0,  1],
        [-1,  0],
        [ 0, -1],
        [ 1,  0]
    ]

    down_trans = [
        [-1,  0],
        [ 0, -1],
        [ 1,  0],
        [ 0,  1]
    ]

    paths = [right_trans,
             forward_trans,
             left_trans,
             down_trans]

    ori_int = 1
    for new_ori, path in enumerate(paths):
        dx, dy = path[ori_int]
        x = (start_x + dx + L) % L
        y = (start_y + dy + W) % W

        if legal(coverage, x, y):
            break

    start_x = x
    start_y = y
    ori_int = ori_to_idx[dx][dy]

    #print "found initial ori:", ori_int


    x = None
    y = None
    x_prev = start_x
    y_prev = start_y
    prev_turn = None

    r = 0
    l = 0

    #plab.figure(1)
    #plab.pcolor(coverage.transpose())
    #plab.scatter(start_x + 0.5, start_y + 0.5, c='w', marker='o')
    while not (start_x == x and start_y == y):

        #print "currently at", ori_int, x_prev, y_prev

        for new_ori, path in enumerate(paths):
            dx, dy = path[ori_int]
            x = (x_prev + dx + L) % L
            y = (y_prev + dy + W) % W

            #print "trying ", new_ori, x, y, coverage[x, y]

            if legal(coverage, x, y):
                #plab.scatter(x+0.5, y+0.5, c="w", marker='x')
                break



        if new_ori == 0:
            r += 1
           # print "went right"
        elif new_ori == 2:
          #  print "went left"
            l += 1
        elif new_ori == 3:
            print "went back"
            print "ERRRORR"
            #plab.show()
            exit(1)
        # else:
        #     print "went forwards"

        #raw_input()
        x_prev = x
        y_prev = y

        ori_int = ori_to_idx[dx][dy]

    if r == l:
        #print r, l
        #plab.show()
        return 2.
    else:
        #plab.clf()
        if coverage.mean() > 0.5:
            return 3.0
        else:
            return 1.0

    #print "fin"


def find_peaks(coverage):

    diff = coverage[1:] - coverage[:-1]

    m = coverage.max()
    return np.where(diff < -m/2)

def main():

    input_file = sys.argv[1]

    plot = False
    if len(sys.argv) > 2:
        plot = True

    parser = ParseKMCHDF5(input_file)

    alphas = []
    F0s = []
    s0s = []
    lmax = 0

    count = 0
    total = 0
    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        if alpha not in alphas:
            alphas.append(alpha)

        F0 = data.attrs["F0"]
        F0 = round(F0, 3)
        if F0 not in F0s:
            F0s.append(F0)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        if s0 not in s0s:
            s0s.append(s0)

        l = data["eq_storedEventValues"].shape[1]
        if l > lmax:
            lmax = l

        total += 1

    F0s = sorted(F0s)
    s0s = sorted(s0s)
    alphas = sorted(alphas)

    np.save("/tmp/extraneighbor_s0s.npy", s0s)
    np.save("/tmp/extraneighbor_alphas.npy", alphas)
    np.save("/tmp/extraneighbor_F0s.npy", F0s)

    cmat = np.zeros(shape=(len(s0s), len(alphas), len(F0s)))
    pmat = np.zeros(shape=(len(s0s), len(alphas), len(F0s)))  #0: no phase | 1: island | 2: band | 3: hole
    ccounts = np.zeros(shape=(len(s0s), len(alphas), len(F0s)))

    facSmall = 9/10.
    facBig = 9/10.
    plotevery = 100
    for data, L, W, run_id in parser:

        alpha = data.attrs["alpha"]
        alpha = round(alpha, 3)
        ia = alphas.index(alpha)

        F0 = data.attrs["F0"]
        F0 = round(F0, 3)
        iF0 = F0s.index(F0)

        s0 = data.attrs["s0"]
        s0 = round(s0, 3)
        is0 = s0s.index(s0)

        eq_coverage = data["eq_storedEventValues"][(0)]

        l = len(eq_coverage)

        if l != lmax:
            p_value = 0

            #The system clogged
            if eq_coverage[-1] > (L*W)/10:
                cval = L*W

            #Because there was no change in coverage
            else:
                cval = 0
        else:
            p_value = get_p_value(data["eq_heights"][()], run_id)

            if eq_coverage[:l/2].mean() < L*W/10:
                start = int(facSmall*l)
                if plot:
                    plab.plot(eq_coverage[::plotevery]/float(L*W))
            else:
                start = int(facBig*l)

            cval = eq_coverage[start:].mean()

        pmat[is0, ia, iF0] += p_value
        cmat[is0, ia, iF0] += cval/float(L*W)
        ccounts[is0, ia, iF0] += 1.

        count += 1

        sys.stdout.flush()
        print "\r%d/%d" % (count, total),
    print

    print ccounts.max(), ccounts.min()

    I = np.where(ccounts != 0)
    J = np.where(ccounts == 0)

    cmat[I] /= ccounts[I]
    cmat[J] = -1

    pmat[I] /= ccounts[I]
    pmat[J] = -1

    np.save("/tmp/extraneighbor_cmat.npy", cmat)
    np.save("/tmp/extraneighbor_pmat.npy", pmat)

    if plot:
        plab.plot([lmax*facSmall/plotevery, lmax*facSmall/plotevery], [0, 1])
        plab.show()

if __name__ == "__main__":
    main()