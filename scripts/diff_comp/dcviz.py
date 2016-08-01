#DCVIZ

from DCViz_sup import DCVizPlotter

from matplotlib.ticker import FuncFormatter

my_props = {

    "fmt" : {
            "markersize" : 7,
            "markeredgewidth" : 1.5,
            "linewidth" : 1,
            "fillstyle" : "none"
    }

}
def g_format_func(v, _):
    return r"$%g$" % v

GFORMATTER = FuncFormatter(g_format_func)

class FelixParticleHDynCav(DCVizPlotter):

    nametag = "diff_comp_(.*)\.npy"
    isFamilyMember = True

    hugifyFonts = True

    figMap = {"figure": ["subfigure_x",
                         "subfigure_xb",
                         "subfigure_y"],
              "allfigure": "planefigure",
              "c_refl_figure": ['c00', 'c01', 'c02'],
              "c_open_figure": ['c10', 'c11', 'c12'],
              "h_refl_figure": 'h_refl',
              "h_open_figure": 'h_open'}

    col_fig_size = [5, 7]
    col_fig_size2 = [6, 5]

    specific_fig_size = {
         "c_refl_figure": col_fig_size,
         "c_open_figure": col_fig_size,
         "h_refl_figure": col_fig_size2,
         "h_open_figure": col_fig_size2
    }

    styles = ["-.", "-", "--"]
    colors = ["k", "b", "r"]

    def adjust(self):
        for figure in ["c_refl_figure", "c_open_figure"]:
            self.adjust_maps[figure]["bottom"] = 0.11
            self.adjust_maps[figure]["top"] = 0.94
            self.adjust_maps[figure]["hspace"] = 0.13

        self.adjust_maps["c_open_figure"]["right"] = 0.87
        self.adjust_maps["c_open_figure"]["left"] = 0.14

        self.adjust_maps["c_refl_figure"]["right"] = 0.93
        self.adjust_maps["c_refl_figure"]["left"] = 0.20


        for figure in ["h_refl_figure", "h_open_figure"]:
            self.adjust_maps[figure]["bottom"] = 0.15
            self.adjust_maps[figure]["top"] = 0.92

        self.adjust_maps["h_open_figure"]["left"] = 0.16
        self.adjust_maps["h_open_figure"]["right"] = 0.96

        self.adjust_maps["h_refl_figure"]["right"] = 1-0.16
        self.adjust_maps["h_refl_figure"]["left"] = 1-0.96


    def plot(self, data):

        C = self.get_family_member_data(data, "C")
        H = self.get_family_member_data(data, "H")
        heights = self.get_family_member_data(data, "heights")

        if self.argv:
            ir = int(self.argv[0])
            ih = int(self.argv[1])

            h = C[ir, ih]
            print heights[ih]
            self.planefigure.pcolor(h)
            self.planefigure.set_xlabel(r"$x/L_x$")
            self.planefigure.set_ylabel(r"$y/L_y$")

            self.subfigure_x.plot(h.mean(axis=0), "k-s")
            self.subfigure_x.set_ylabel(r"$P(x/L_x)$")
            self.subfigure_x.set_ybound(0)

            self.subfigure_xb.plot(h.mean(axis=0), "k-s")
            self.subfigure_xb.set_ylabel(r"$P(x/L_x)$")
            self.subfigure_xb.set_xlabel(r"$x/L_x$")

            self.subfigure_y.plot(h.mean(axis=1), "k-s")
            self.subfigure_y.set_ylabel(r"$P(y/L_y)$")
            self.subfigure_y.set_xlabel(r"$y/L_y$")
            self.subfigure_y.set_ybound(0)

        else:
            _, _, n, _ = C.shape
            _, _, W, L = H.shape
            X = np.arange(n, dtype=float)/(n-1)
            Xh = np.arange(L, dtype=float)/(L-1)
            titles = [r"$\mathrm{Open}$", r"$\mathrm{Reflecting}$"]

            hfigs = [self.h_open, self.h_refl]

            self.h_open.set_ylabel("$h(x/L_x)$")
            self.h_open.yaxis.set_major_formatter(GFORMATTER)

            ymaxes = [0,0,0]
            for ir in range(2):
                hfigs[ir].set_title(titles[ir])
                hfigs[ir].set_xlabel("$x/L_x$")
                hfigs[ir].set_ylim(-5, 15)
                hfigs[ir].xaxis.set_major_formatter(GFORMATTER)

                for ih in range(len(heights)):
                    hlabel = r"$\Delta h=%d$" % heights[2-ih]

                    hfigs[ir].plot(Xh, H[ir, 2-ih].mean(axis=0),
                                   self.colors[ih] + self.styles[ih],
                                   label=hlabel,
                                   linewidth=2)

                    Y = C[ir, 2-ih].mean(axis=0)*100

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
                        cfig.set_ylabel("$c(x/L_x)\,\,[\%]$")

                    l = hfigs[ir].axes.legend(loc="center",
                               numpoints=1,
                               ncol=1,
                               handlelength=1.0,
                               borderpad=0.2,
                               labelspacing=0.2,
                               columnspacing=0.3,
                               handletextpad=0.25,
                               borderaxespad=0.0,
                               frameon=False,
                               fontsize=25,
                               bbox_to_anchor=(0.8, 0.8))

                    l.get_frame().set_fill(not (self.toFile and self.transparent))

            self.h_refl.set_yticklabels([])

            for ir in range(2):
                for ih in range(3):
                    cfig = eval("self.c%d%d" % (ir, ih))
                    cfig.set_ylim(0, ymaxes[ih]*1.1)