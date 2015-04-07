from DCViz_sup import DCVizPlotter

class SteadyState(DCVizPlotter):

    nametag = "steadystate\_(.*)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {}

    def plot(self, data):

        clip = 1.1
        start = 1
        transform = True

        I, J, K = [int(x) for x in self.argv]

        E0_values = data[self.get_family_index_from_name("steadystate_E0.npy")].data
        alpha_values = data[self.get_family_index_from_name("steadystate_alpha.npy")].data
        mu_shift_values = data[self.get_family_index_from_name("steadystate_mu_shift.npy")].data

        count = 0
        for i, E0 in enumerate(E0_values):
            for j, alpha in enumerate(alpha_values):
                for k, mu_shift in enumerate(mu_shift_values):
                    if i == I and j == J and k == K:
                        print E0, alpha, mu_shift

                        time = data[self.get_family_index_from_name("steadystate_Time_%d.npy" % count)].data
                        L = int(len(time)/clip)
                        time = time[start:L]

                        rms = data[self.get_family_index_from_name("steadystate_HeightRMS_%d.npy" % count)].data[start:L]

                        ss = data[self.get_family_index_from_name("steadystate_SurfaceSize_%d.npy" % count)].data[start:L]

                        wh = data[self.get_family_index_from_name("steadystate_PressureWall_%d.npy" % count)].data[start:L]

                        rmslabel=r"$\sigma(h)$"
                        sslabel=r"$\sigma(s)$"
                        whlabel=r"$h_l(t) - \langle h(t) \rangle$"

                        if transform:
                            self.subfigure.plot(time[start:], (rms[start:] - rms[start])/(rms.max() - rms[start]), label=rmslabel)
                            self.subfigure.plot(time[start:], (ss[start:] - ss[start])/(ss.max() - ss[start]), label=sslabel)
                            self.subfigure.plot(time[start:], (wh[start:] - wh[start])/(wh.max() - wh[start]), label=whlabel)
                        else:
                            self.subfigure.plot(time, rms/rms.max(), label=rmslabel)
                            self.subfigure.plot(time, ss/ss.max(), label=sslabel)
                            self.subfigure.plot(time, wh/wh.max(), label=whlabel)



                        # self.subfigure.text(time[-1]/clip*xshift, 1./3 - yshift, whlabel,size=self.labelSize)


                    count += 1

        self.subfigure.set_xlabel(r"$t$")
        self.subfigure.legend(loc="lower right")
        self.subfigure.set_xlim(0, self.subfigure.get_xlim()[1])
        self.subfigure.set_ylim(0, 1.1)
