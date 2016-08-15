#DCVIZ

from DCViz_sup import DCVizPlotter

from matplotlib.ticker import FuncFormatter

from intercombinatorzor import ICZ

import numpy as np

from scipy.optimize import curve_fit


my_props = {

    "fmt" : {
            "markersize" : 7,
            "markeredgewidth" : 1.5,
            "linewidth" : 1,
            "fillstyle" : "none"
    }

}

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
                         "vfigure": [6, 4],
                         "twostepfig": [6, 6],
                         "twosteplogfig": [6, 4],
                         "speedfig": [6, 3]}

    figMap = {"figure": ["subfigure", "subfigure_combo"],
              "vfigure": "vfig",
              "speedfig": "sfig",
              "allfig": "afig",
              "allfig2": "afig2",
              "twostepfig": ["afig3", "lafig3"]}

    def adjust(self):
        l = 0.18
        for fig in ["figure", "twostepfig", "vfigure", "speedfig"]:
            self.adjust_maps[fig]['left'] = l
            self.adjust_maps[fig]['right'] = 1-l
            self.adjust_maps[fig]['hspace'] = 0.125

        self.adjust_maps["figure"]['top'] = 0.97
        self.adjust_maps["figure"]['bottom'] = 0.13

        self.adjust_maps["twostepfig"]['top'] = 0.97
        self.adjust_maps["twostepfig"]['bottom'] = 0.13

        self.adjust_maps["speedfig"]['top'] = 0.97
        self.adjust_maps["speedfig"]['bottom'] = 0.24

    #plotOnly = "twostepfig"

    def find_half(self, t):

        i = 0
        while t[i] < 0.9:
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

        assume_twos = 3 in skip
        lim = 0.025
        back = 0
        get_s = lambda xvals: np.where(xvals < lim)[0][-1] - back
        if assume_twos:
            combinator = ICZ("Time", "h1", "h2")

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

                combinator.feed(t_trans/(t[s2]-t[s1]), x1, x2)
            t1, h1 = combinator.intercombine("Time", "h1")
            _, h2 = combinator.intercombine("Time", "h2")

            h2eq = 0.57
            h1eq = 0.7
            xmax = 12

            self.afig3.plot(t1, h1, "r-", label=r"$y_s^{(1)}(t')$")
            scale = h1[self.find_half(t1)]/h2[self.find_half(t1-1)]
            self.afig3.plot(t1-1, h2*scale, "k--", label=r"$%.2fy_s^{(2)}(t'-1)$" % scale)
            self.afig3.plot(t1, h2, "k-", linewidth=2, label=r"$y_s^{(2)}(t')$")

            def g(t, k, ys):
                return ys*(1-np.exp(-k*t))

            I = np.where(t1 > 1.5)
            J = np.where(t1[I] < xmax)
            t2 = t1[I][J]
            h2p = h2[I][J]

            p0 = (1,1)
            p, _ = curve_fit(g, t2, h2p, p0)
            print p
            #self.afig3.plot(t2, g(t2, *p), "k--", linewidth=2)

            self.lafig3.semilogy(t1, abs(1-h1/h1eq), "r-")
            self.lafig3.semilogy(t1, abs(1-h2/h2eq), "k-")
            self.lafig3.set_xlim(0, 9)
            self.lafig3.set_ylim(5E-3, 1.1)

            # self.lafig3.set_ylim(0, 0.8)
            # self.lafig3.set_xlim(0.4, 0.8)

            l = self.afig3.legend(loc="center",
                   numpoints=1,
                   ncol=2,
                   handlelength=1.25,
                   borderpad=0.2,
                   labelspacing=0.2,
                   columnspacing=-2.5,
                   handletextpad=0.25,
                   borderaxespad=0.0,
                   frameon=False,
                   fontsize=20,
                   bbox_to_anchor=(0.625, 0.25))

            l.get_frame().set_fill(not (self.toFile and self.transparent))

            self.afig3.set_xlim(-0.5, xmax)
            self.lafig3.set_xlim(-0.5, xmax)
            #self.afig3.set_xlabel(r"$t'=(t-t_1)/(t_2-t_1)$")
            self.afig3.set_xticklabels([])
            self.lafig3.set_xlabel(r"$t'=(t-t_1)/(t_2-t_1)$")
            self.afig3.set_xticks(range(xmax+1))
            self.lafig3.set_xticks(range(xmax+1))
            self.afig3.set_ylabel(r"$y_s/L_y$",labelpad=20)
            self.lafig3.set_ylabel(r"$|1-y_s/\langle y_s\rangle|$")
            self.afig3.set_ylim(0, 0.8)

            # def haxxorformatter(v, _):
            #     if v == 0:
            #         return r"$0$"
            #     else:
            #         if v % 2 == 0:
            #             return ""
            #         else:
            #             return r"$%g$" % v

            # self.afig3.xaxis.set_major_formatter(FuncFormatter(haxxorformatter))

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


        self.subfigure_combo.set_ylabel(r"$y_s/L_y$")
        self.subfigure_combo.set_xlabel(r"$t/t_\mathrm{end}$")

        self.subfigure.set_ylabel(r"$y_s/L_y$")

        self.subfigure.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure.yaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.yaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)
        self.subfigure_combo.xaxis.set_major_formatter(GFORMATTER)

        self.subfigure.set_xticklabels([])

        vs = []
        for i in range(len(seq_levels)-1):
            diff = seq_levels[i]-seq_levels[i+1]
            self.vfig.plot(tf/tf[-1], diff)
            vs.append(100*diff[len(diff)/2:].mean())
        self.sfig.plot(range(1, len(vs)+1), vs, "r--", **my_props['fmt'])
        self.sfig.plot(range(1, len(vs)+1), vs, "ks", **my_props['fmt'])
        self.sfig.set_xlabel(r"$n$")
        self.sfig.set_ylabel(r"$\langle y_s^{(n)}-y_s^{(n+1)}\rangle\,\,[l_0]$")
        self.sfig.set_xlim(0.9, len(vs)+0.1)
        self.sfig.set_ylim(np.floor(min(vs)), max(vs)+0.1)


        # def magicomg(v, _):
        #     return r"$\Delta \langle x\rangle_{%d,%d}$" % (v, v+1)
        #
        # self.sfig.xaxis.set_major_formatter(FuncFormatter(magicomg))

