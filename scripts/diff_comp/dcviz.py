#DCVIZ

from DCViz_sup import DCVizPlotter

from matplotlib.ticker import FuncFormatter
import numpy as np

my_props = {

    "fmt" : {
            "markersize" : 7,
            "markeredgewidth" : 1.5,
            "linewidth" : 1,
            "fillstyle" : "none"
    }

}
def g_format_func(v, _):
    if v == 0:
        return r"$0$"
    elif int(v) == v:
        return r"$%.1f$" % v
    else:
        return r"$%g$" % v

GFORMATTER = FuncFormatter(g_format_func)

class FelixParticleHDynCav(DCVizPlotter):

    nametag = "diff_comp_(.*)\.npy"
    isFamilyMember = True

    hugifyFonts = True

    figMap = {"c_open_figure": ['c00', 'c01', 'c02'],
              "c_refl_figure": ['c10', 'c11', 'c12'],
              "h_refl_figure": 'h_refl',
              "h_open_figure": 'h_open',
              "Ca_refl_figure": 'C_refl',
              "Ca_open_figure": 'C_open'}

    #plotOnly = ["c_open_figure", "c_refl_figure", "h_refl_figure", "h_open_figure"]

    col_fig_size = [5, 7]
    col_fig_size2 = [5, 4]

    specific_fig_size = {
         "c_refl_figure": col_fig_size,
         "c_open_figure": col_fig_size,
         "h_refl_figure": col_fig_size2,
         "h_open_figure": col_fig_size2,
         "Ca_refl_figure": col_fig_size2,
         "Ca_open_figure": col_fig_size2
    }

    styles = ["-.", "-", "--"]
    colors = ["k", "b", "r"]
    mstyles =["s", "^", "o"]

    def adjust(self):
        for figure in ["c_refl_figure", "c_open_figure"]:
            self.adjust_maps[figure]["bottom"] = 0.11
            self.adjust_maps[figure]["top"] = 0.94
            self.adjust_maps[figure]["hspace"] = 0.13

        self.adjust_maps["c_refl_figure"]["right"] = 1-0.24
        self.adjust_maps["c_refl_figure"]["left"] = 1-0.97

        self.adjust_maps["c_open_figure"]["right"] = 0.97
        self.adjust_maps["c_open_figure"]["left"] = 0.24

        for figure in ["h_refl_figure", "h_open_figure"]:
            self.adjust_maps[figure]["bottom"] = 0.19
            self.adjust_maps[figure]["top"] = 0.96

        for figure in ["Ca_refl_figure", "Ca_open_figure"]:
            self.adjust_maps[figure]["bottom"] = 0.13
            self.adjust_maps[figure]["top"] = 0.9

        self.adjust_maps["h_open_figure"]["left"] = 0.19
        self.adjust_maps["h_open_figure"]["right"] = 0.95
        self.adjust_maps["Ca_open_figure"]["left"] = 0.19
        self.adjust_maps["Ca_open_figure"]["right"] = 0.95

        self.adjust_maps["h_refl_figure"]["right"] = 1-0.19
        self.adjust_maps["h_refl_figure"]["left"] = 1-0.95
        self.adjust_maps["Ca_refl_figure"]["right"] = 1-0.19
        self.adjust_maps["Ca_refl_figure"]["left"] = 1-0.95

    def plot(self, data):

        C = self.get_family_member_data(data, "C")
        H = self.get_family_member_data(data, "H")
        heights = self.get_family_member_data(data, "heights")

        _, _, n, _ = C.shape
        _, _, W, L = H.shape
        X = np.arange(n, dtype=float)/(n-1)
        Xh = np.arange(L, dtype=float)/(L-1)
        titles = [r"$\mathrm{Open}$", r"$\mathrm{Reflecting}$"]

        hfigs = [self.h_open, self.h_refl]
        cfigs = [self.C_open, self.C_refl]

        self.h_open.set_ylabel("$h(x/L_x)/l_0$")
        self.C_open.set_ylabel("$10^4c(x/L_x)l_0^{3}$")

        ymaxes = [0,0,0]
        for ir in range(2):
            #hfigs[ir].set_title(titles[ir])
            hfigs[ir].set_xlabel("$x/L_x$")
            hfigs[ir].set_ylim(-2, 4)
            hfigs[ir].xaxis.set_major_formatter(GFORMATTER)

            cfigs[ir].set_title(titles[ir])
            #cfigs[ir].set_xlabel("$x/L_x$")
            cfigs[ir].xaxis.set_major_formatter(GFORMATTER)


            for ih in range(len(heights)):
                hlabel = r"$\Delta h/l_0=%d$" % heights[2-ih]

                hfigs[ir].plot(Xh, H[ir, 2-ih].mean(axis=0),
                               self.colors[ih] + self.styles[ih],
                               label=hlabel,
                               linewidth=2)

                Y = C[ir, 2-ih].mean(axis=0)*10**4

                cfigs[ir].plot(X, Y, self.colors[ih]+self.mstyles[ih]+"--", label=hlabel, **my_props["fmt"])

                ym = Y.max()
                if ym > ymaxes[ih]:
                    ymaxes[ih] = ym

                cfig = eval("self.c%d%d" % (ir, ih))

                cfig.plot(X, Y, "r--", linewidth=2)
                cfig.plot(X, Y, "ks", **my_props["fmt"])

                cfig.set_xlim(-0.02, 1.02)


                if ih == 0:
                    cfig.set_title(titles[ir])

                if ih != 2:
                    cfig.set_xticklabels([])
                else:
                    cfig.set_xlabel("$x/L_x$")

                    cfig.xaxis.set_major_formatter(GFORMATTER)
                cfig.yaxis.set_major_formatter(GFORMATTER)

                if ir == 1:
                    ax = cfig.axes.twinx()
                    ax.set_ylabel(hlabel, labelpad=15)
                    ax.yaxis.set_ticks([])
                    ax.yaxis.set_ticklabels([])
                    cfig.set_yticklabels([])
                else:
                    cfig.set_ylabel("$10^4c(x/L_x)l_0^{3}$")

        l = self.h_refl.axes.legend(loc="center",
                   numpoints=1,
                   ncol=1,
                   handlelength=1.5,
                   borderpad=0.2,
                   labelspacing=0.2,
                   columnspacing=0.3,
                   handletextpad=0.25,
                   borderaxespad=0.0,
                   frameon=False,
                   fontsize=25,
                   bbox_to_anchor=(0.6, 0.75))

        l.get_frame().set_fill(not (self.toFile and self.transparent))

        self.h_refl.set_yticklabels([])


        l = self.C_refl.axes.legend(loc="center",
                   numpoints=1,
                   ncol=1,
                   handlelength=1.5,
                   borderpad=0.2,
                   labelspacing=0.2,
                   columnspacing=0.3,
                   handletextpad=0.25,
                   borderaxespad=0.0,
                   frameon=False,
                   fontsize=25,
                   bbox_to_anchor=(0.4, 0.225))

        l.get_frame().set_fill(not (self.toFile and self.transparent))
        #
        # c0 = np.exp(-6)*100
        # for cfig in cfigs:
        #     cfig.plot([0, 1], [c0, c0], "k:")

        cfigs[0].set_ylim(0, max(ymaxes)*1.1)
        cfigs[1].set_ylim(0, max(ymaxes)*1.1)
        cfigs[1].set_yticklabels([])
        cfigs[0].set_xticklabels([])
        cfigs[1].set_xticklabels([])

        for ir in range(2):
            for ih in range(3):
                cfig = eval("self.c%d%d" % (ir, ih))
                cfig.set_ylim(0, ymaxes[ih]*1.1)