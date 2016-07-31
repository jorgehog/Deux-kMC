#DCVIZ

from DCViz_sup import DCVizPlotter


class FelixParticleHDynCav(DCVizPlotter):

    nametag = "diff_comp_(.*)\.npy"

    hugifyFonts = True

    figMap = {"figure": ["subfigure_x",
                         "subfigure_xb",
                         "subfigure_y"],
              "allfigure": "planefigure"}


    def plot(self, data):

        self.planefigure.pcolor(data)
        self.planefigure.set_xlabel(r"$x/L_x$")
        self.planefigure.set_ylabel(r"$y/L_y$")

        self.subfigure_x.plot(data.mean(axis=0), "k-s")
        self.subfigure_x.set_ylabel(r"$P(x/L_x)$")
        self.subfigure_x.set_ybound(0)

        self.subfigure_xb.plot(data.mean(axis=0), "k-s")
        self.subfigure_xb.set_ylabel(r"$P(x/L_x)$")
        self.subfigure_xb.set_xlabel(r"$x/L_x$")
        
        self.subfigure_y.plot(data.mean(axis=1), "k-s")
        self.subfigure_y.set_ylabel(r"$P(y/L_y)$")
        self.subfigure_y.set_xlabel(r"$y/L_y$")
        self.subfigure_y.set_ybound(0)

