#DCVIZ

from DCViz_sup import DCVizPlotter

from matplotlib.ticker import FuncFormatter

from intercombinatorzor import ICZ

import numpy as np

def g_format_func(v, _):
    if int(v) == v:
        return r"$%d$" % v
    elif abs(int(10*v) - 10*v) < 1E-3:
        return r"$%.1f$" % v
    else:
        return ""

GFORMATTER = FuncFormatter(g_format_func)

class flx_bunch(DCVizPlotter):

    nametag = "flx_bunch_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    fig_size = [6, 8]

    specific_fig_size = {"figure": [6, 6],
                         "vfigure": [6, 4]}

    figMap = {"figure": ["subfigure", "subfigure_combo"],
              "vfigure": "vfig",
              "speedfig": "sfig",
              "allfig": "afig",
              "allfig2": "afig2",
              "allfig3": "afig3"}

    def adjust(self):
        fig = "figure"
        l = 0.16
        self.adjust_maps[fig]['top'] = 0.97
        self.adjust_maps[fig]['bottom'] = 0.13
        self.adjust_maps[fig]['left'] = l
        self.adjust_maps[fig]['right'] = 1-l
        self.adjust_maps[fig]['hspace'] = 0.05

    plotOnly = "allfig3"

    def find_half(self, t):

        i = 0
        while t[i] < 0.5:
            i+=1

        return i

    def plot(self, data):

        if self.argv:
            j = int(self.argv[0])
        else:
            j = 0


        combo = len(self.argv) > 1
        if combo:
            every = int(self.argv[1])

        all_t = self.get_family_member_data(data, "time")
        all_x = self.get_family_member_data(data, "fpos")


        nr, nl, nc = all_x.shape

        T = np.zeros(shape=(nc))
        X = np.zeros(shape=(nl, nc))
        nx = np.zeros_like(X)

        combinator = ICZ("Time", *range(1, nl))

        skip = set()

        for i in range(nr):
            t = all_t[i]
            x = all_x[i]

            if combo:
                combinator.feed(t[::every], *x[1:, ::every])

            for hlevel in range(1, nl):
                xf = x[hlevel]

                if xf.mean() < 0.01:
                    skip.add(hlevel)
                    continue

                if i == j:
                    self.subfigure.plot(t/t[-1], xf, "r-")

                X[hlevel] += xf
                nx[hlevel] += 1

            T += t

        assume_twos = True
        lim = 0.03
        back = 5
        get_s = lambda xvals: np.where(xvals < lim)[0][-1] - back
        if assume_twos:
            combinator = ICZ("Time", "h1")
            combinator2 = ICZ("Time", "h2")

            for i in range(nr):

                t = all_t[i]

                x1 = all_x[i][1]
                x2 = all_x[i][2]

                s1 = get_s(x1)
                s2 = get_s(x2)

                t_trans = (t - t[s1])
                t_trans2 = (t - t[s2])

                self.afig.plot(t_trans, x1)
                self.afig2.plot(t_trans2, x2)

                combinator.feed(t_trans/t[s2], x1)
                combinator2.feed(t_trans2/t[s2], x2)
            t1, h1 = combinator.intercombine("Time", "h1")
            t2, h2 = combinator2.intercombine("Time", "h2")

            self.afig3.plot(t1, h1)
            self.afig3.plot(t2+1, h2)
            self.afig3.plot(t2, h2*h1[self.find_half(t1)]/h2[self.find_half(t2)])

        seq_levels = []
        for hlevel in range(1, nl):
            if hlevel in skip:
                continue

            if combo:
                tf, xl = combinator.intercombine("Time", hlevel)
            else:
                xl = X[hlevel]/nx[hlevel]
                tf = T

            print hlevel
            seq_levels.append(xl)
            self.subfigure_combo.plot(tf/tf[-1], xl, "r-")

        #for xf in np.linspace(min(xfs), max(xfs), len(xfs)):
        # for i, xf in enumerate(xfs):
        #     hlevel = i+1
        #     self.subfigure.text(1.02, xf, r"$s_{%d}$" % hlevel,
        #                         verticalalignment="center",
        #                         horizontalalignment="left")


        # self.subfigure.set_xlim(0, 1.07)
        # self.subfigure.set_ylim(0, x[1:, :].max()+0.02)
        self.subfigure.set_xlim(0, 1)
        self.subfigure.set_ylim(0, all_x[:, 1:, :].max()+0.01)
        self.subfigure_combo.set_xlim(0, 1)
        self.subfigure_combo.set_ylim(0, all_x[:, 1:, :].max()+0.01)


        self.subfigure_combo.set_ylabel(r"$y/L_y$")
        self.subfigure_combo.set_xlabel(r"$t/t_\mathrm{end}$")

        self.subfigure.set_ylabel(r"$y/L_y$")

        self.subfigure.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure.yaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.yaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)

        self.subfigure.set_xticklabels([])

        for i in range(len(seq_levels)-1):
            diff = seq_levels[i+1]-seq_levels[i]
            self.vfig.plot(tf/tf[-1], diff)

        every = 100
        seq_levels[0] = seq_levels[0][::every]
        diff = np.zeros(len(seq_levels[0])-2)
        for i in range(1, len(seq_levels[0])-1):
            n = i+1
            p = i
            prev = seq_levels[0][p]
            next = seq_levels[0][n]

            tprev = tf[p]
            tnext = tf[n]

            diff[i-1] = tf[-1]*(next - prev)/(tnext-tprev)

        self.sfig.plot(tf[::every][1:-1]/tf[-1], diff)