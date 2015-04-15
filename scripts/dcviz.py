#DCVIZ

from DCViz_sup import DCVizPlotter

import numpy as np
import os
import re
from numpy import exp
from scipy.stats import linregress
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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

        E0_values = data[self.get_family_index_from_name("steadystate_E0.npy")]
        alpha_values = data[self.get_family_index_from_name("steadystate_alpha.npy")]
        mu_shift_values = data[self.get_family_index_from_name("steadystate_mu_shift.npy")]

        count = 0
        for i, E0 in enumerate(E0_values):
            for j, alpha in enumerate(alpha_values):
                for k, mu_shift in enumerate(mu_shift_values):
                    if i == I and j == J and k == K:
                        print E0, alpha, mu_shift

                        time = data[self.get_family_index_from_name("steadystate_Time_%d.npy" % count)]
                        L = int(len(time)/clip)
                        time = time[start:L]

                        rms = data[self.get_family_index_from_name("steadystate_HeightRMS_%d.npy" % count)][start:L]

                        ss = data[self.get_family_index_from_name("steadystate_SurfaceSize_%d.npy" % count)][start:L]

                        wh = data[self.get_family_index_from_name("steadystate_PressureWall_%d.npy" % count)][start:L]

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

                    count += 1

        self.subfigure.set_xlabel(r"$t$")
        self.subfigure.legend(loc="lower right")
        self.subfigure.set_xlim(0, self.subfigure.get_xlim()[1])
        self.subfigure.set_ylim(0, 1.1)


class GrowthSpeed(DCVizPlotter):

    nametag = "growthspeed\_(.*)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"omega_vs_v": "subfigure", "slopes": "subfigure2"}

    def plot(self, data):

        plot_values = [0, 0.5, 1.0]
        shapes = ["s", "^", "v"]

        E0_values = data[self.get_family_index_from_name("growthspeed_E0.npy")]
        mu_shift_values = data[self.get_family_index_from_name("growthspeed_mu_shift.npy")]
        mu_values = data[self.get_family_index_from_name("growthspeed_mu.npy")]
        v_values = data[self.get_family_index_from_name("growthspeed_v.npy")]
        unloaded = self.get_family_member_data(data, "unloaded")

        X = np.array([mu_shift_values]).T
        V = np.array([unloaded]).T


        E0_values = np.insert(E0_values, 0, 0)
        mu_values = np.concatenate((X, mu_values), axis=1)
        v_values = np.concatenate((V, v_values), axis=1)
        print mu_shift_values
        print
        print mu_values

        omega = exp(mu_shift_values) - 1

        slopes = np.zeros_like(E0_values)
        mu_zero = np.zeros_like(slopes)


        print unloaded.shape, omega.shape, v_values.shape

        k = 0
        for i, E0 in enumerate(E0_values):

            mu = mu_values[:, i]
            mu0 = (mu - mu_shift_values).mean()

            c_over_c0 = exp(mu0 + mu_shift_values)

            v = v_values[:, i]*c_over_c0

            slope, intercept, _, _, _ = linregress(omega, v)

            if E0 in plot_values:
                # self.subfigure.plot(omega, intercept + slope*omega, "r--")
                self.subfigure.plot(omega, v, shapes[k], label="E0=%g" % E0,
                                    linewidth=1,
                                    fillstyle='none',
                                    markersize=7,
                                    markeredgewidth=1.5,
                                    color="black")
                k += 1

            print intercept

            slopes[i] = slope
            mu_zero[i] = mu0

        self.subfigure2.plot(E0_values, slopes, "ks--",
                             linewidth=1,
                             fillstyle='none',
                             markersize=7,
                             markeredgewidth=1.5)

        self.subfigure.set_xlabel(r"$\Omega$")
        self.subfigure.set_ylabel(r"$\dot{H}$")
        self.subfigure.legend(loc="upper left")

        self.subfigure2.set_xlabel(r"$E_0/L$")
        self.subfigure2.set_ylabel(r"$\partial \dot{H} / \partial \Omega$")



class Quasi2D_slopes_and_stuff(DCVizPlotter):

    nametag = "linearplots\_(.*)\_?\d*\.npy"

    numpyBin = True

    isFamilyMember = True

    figMap = {"fig" : "gammaslopes", "fig2" : "E0slopes"}

    hugifyFonts = True

    def n_runs(self):

        n_max = -1

        for name in self.familyFileNames:
            if not "alpha" in name:
                continue

            n = int(re.findall("linearplots_alpha_(\d+)\.npy", name)[0])

            if n > n_max:
                n_max = n

        return n_max + 1


    def plot(self, data):

        meta_data_file = os.path.join(self.familyHome, "linearplots_metadata.dat")

        with open(meta_data_file, 'r') as meta_data:
            print meta_data.read()

        N = 3
        shapes = ["s", "^", "v"]

        n_runs = self.n_runs()

        E0s = data[self.get_family_index_from_name("linearplots_E0.npy")]
        idxmap = range(len(E0s))

        E0s, idxmap = [np.array(x) for x in zip(*sorted(zip(E0s, idxmap), key=lambda x:x[0]))]

        slopes = data[self.get_family_index_from_name("linearplots_slopes.npy")][idxmap]
        print E0s
        n_plots = 0
        errors = []
        for n in range(n_runs):

            alpha_name = "linearplots_alpha_%d.npy" % idxmap[n]
            mu_name = "linearplots_muEqs_%d.npy" % idxmap[n]

            alphas = data[self.get_family_index_from_name(alpha_name)]
            mus = data[self.get_family_index_from_name(mu_name)]

            a, b, c, d, err = linregress(alphas, mus)

            errors.append(err)

            if n%(n_runs/(N-1)) != 0:
                continue
            print n, n_plots

            mu_error_name = "linearplots_muEqErrors_%d.npy" % n

            mu_errors = data[self.get_family_index_from_name(mu_error_name)]

            self.gammaslopes.errorbar(alphas, mus,
                                      yerr=mu_errors,
                                      fmt=shapes[n_plots],
                                      fillstyle='none',
                                      label=r"$E_0/L=%1.2f$" % (E0s[n]),
                                      markersize=7,
                                      markeredgewidth=1.5,
                                      linewidth=1,
                                      color="black")
            self.gammaslopes.plot([0, alphas.max()], [0, slopes[n]*alphas.max()], "r--")

            n_plots += 1

        self.gammaslopes.set_xbound(0)
        self.gammaslopes.set_ybound(0)
        self.gammaslopes.legend(loc="upper left")

        self.gammaslopes.set_xlabel(r"$\alpha$")
        self.gammaslopes.set_ylabel(r"$\gamma_\mathrm{eq} - \gamma_0$")

        self.E0slopes.errorbar(E0s, slopes,
                               yerr=errors,
                               fmt="s",
                               fillstyle='none',
                               markersize=7,
                               markeredgewidth=1.5,
                               linewidth=1,
                               color="black")

        sslope, a, b, c, d = linregress(E0s, slopes)

        self.E0slopes.plot([0, E0s.max()], [0, sslope*E0s.max()], 'r--')
        self.E0slopes.set_xbound(0)
        self.E0slopes.set_ylim(0, sslope*E0s.max()*1.05)

        self.E0slopes.set_xlabel(r"$E_0/L$")
        self.E0slopes.set_ylabel(r"$K(E_0/L)$")

        print sslope, d


class SOS_pressure_sizes(DCVizPlotter):

    nametag = "pressure_plots_.*\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"sfig" : "subfigure", "varfig" : "subfigure2", "cvfig" : "subfigure3"}

    def plot(self, data):

        p = 2

        shapes = ["s", "^", "o"]


        E0_array = data[self.get_family_index_from_name("pressure_plots_E0.npy")]
        alphas = data[self.get_family_index_from_name("pressure_plots_alphas.npy")]
        mean_s = data[self.get_family_index_from_name("pressure_plots_mean_s.npy")]
        var_s = data[self.get_family_index_from_name("pressure_plots_var_s.npy")]

        print E0_array.shape, alphas.shape, mean_s.shape, var_s.shape

        analytical_path = os.path.join(self.familyHome, "boltzmann_ascii_full256.arma")
        if os.path.exists(analytical_path):
            analytical = np.loadtxt(analytical_path)
            self.subfigure.plot(analytical[:, 0], analytical[:, 1], 'r-',
                                linewidth=3,
                                label="$E_0 = 0$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure2.plot(analytical[:, 0], analytical[:, 2], 'r-',
                                linewidth=3,
                                label="$E_0 = 0$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure3.plot(analytical[:, 0], analytical[:, 0]**p*analytical[:, 2], 'r-',
                                linewidth=3,
                                label="$E_0 = 0$",
                                fillstyle='none',
                                markersize=5)
        else:
            print("KEINE ANALYTISCH")


        nplots = 0
        for i, E0_value in sorted(enumerate(E0_array), key=lambda x: x[1]):
            #
            #
            # if i%(len(E0_array)/(N-1)) != 0:
            #     continue

            alpha_array = alphas[i, :]
            mean_s_array = mean_s[i, :]
            var_s_array = var_s[i, :]

            alpha_array, mean_s_array, var_s_array = [np.array(x) for x in zip(*sorted(zip(alpha_array, mean_s_array, var_s_array), key=lambda x: x[0]))]

            self.subfigure.plot(alpha_array, mean_s_array, 'k%s' % shapes[nplots],
                                 fillstyle='none',
                                 label="$E_0 = %g$" % E0_value,
                                 markersize=7,
                                 markeredgewidth=1.5,
                                 linewidth=1)

            self.subfigure2.plot(alpha_array, var_s_array, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$E_0 = %g$" % E0_value,
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            self.subfigure3.plot(alpha_array, var_s_array*alpha_array**p, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$E_0 = %g$" % E0_value,
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            nplots += 1

        self.subfigure.set_xlabel(r"$\alpha$")
        self.subfigure.set_ylabel(r"$\langle s \rangle / L$")
        ax = self.subfigure.axes.twinx()
        ax.set_ylabel(r"$\langle E \rangle /E_bL$")
        ax.yaxis.set_ticklabels([])
        self.subfigure.set_xbound(0)
        self.subfigure.legend(numpoints=1, handlelength=1)


        self.subfigure2.set_xlabel(r"$\alpha$")
        self.subfigure2.set_ylabel(r"$\sigma (s) / L$")
        ax2 = self.subfigure2.axes.twinx()
        ax2.set_ylabel(r"$\sigma (E) / E_bL$")
        ax2.yaxis.set_ticklabels([])
        self.subfigure2.set_xbound(0)
        self.subfigure2.legend(numpoints=1, handlelength=1)

        self.subfigure3.set_xlabel(r"$\alpha$")
        self.subfigure3.set_ylabel(r"$\alpha^2\sigma(s)^2$")
        ax3 = self.subfigure3.axes.twinx()
        ax3.set_ylabel(r"$C_V/k$")
        ax3.yaxis.set_ticklabels([])
        self.subfigure3.set_xbound(0)
        self.subfigure3.legend(numpoints=1, handlelength=1, loc="upper left")


class SOSanalyze(DCVizPlotter):

    nametag = "analyze_(.+)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"s0_figure" : "ss", "figure_2" : "mean_figure"}

    def scale(self, r0):
        return r0*(1-exp(-1./r0))

    def plot(self, data):

        dirname = re.findall("analyze\_(.+)\_.+\_values\.npy", self.familyHead)[0]

        print dirname

        C =  data[self.get_family_index_from_name("analyze_%s_C_values.npy" % dirname)]
        s0 =  data[self.get_family_index_from_name("analyze_%s_s0_values.npy" % dirname)]
        r0 =  data[self.get_family_index_from_name("analyze_%s_r0_values.npy" % dirname)]

        s0_cut = np.where(s0 >= 0.)
        r0_cut = np.where(r0 > 0.)
        r0 = r0[r0_cut]
        s0 = s0[s0_cut]

        #
        r0_idx, s0_idx = np.meshgrid(s0_cut, r0_cut)
        C = 1./C[s0_idx, r0_idx]

        # print s0_idx.shape, r0_idx.shape
        #
        # ax = Axes3D(self.figure)
        #
        # X, Y = np.meshgrid(s0[s0_cut], r0[r0_cut])
        #
        # Z = C[s0_idx, r0_idx]
        # Z *= self.scale(X)
        #
        # colormap = cm.hot
        #
        # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=colormap, linewidth=0)
        #
        # zdir = "x"
        # offset = -0.5
        #
        # cset = ax.contour(X, Y, Z, zdir=zdir, offset=offset*1.05, color='#008000')
        #
        # ax.set_xlim(-0.5, s0[s0_cut].max()*1.1)
        # ax.set_zlim(0, Z.max())
        #
        # ax.set_xlabel(r"$\sigma_0$")
        # ax.set_ylabel(r"$\lambda_D$")
        # ax.view_init(15, -50)
        # ax.set_zticks([0, 0.2, 0.4, 0.6, 0.8])
        #
        # self.subfigure.axes.contourf(X, Y, Z, zdir='z', cmap=colormap)
        # self.subfigure.set_xlabel(r"$\sigma_0$")
        # self.subfigure.set_ylabel(r"$\lambda_D$", rotation=0)
        # self.subfigure.set_ybound(0)
        # self.subfigure.axes.get_yaxis().get_label().set_fontsize(20)
        # self.subfigure.axes.get_xaxis().get_label().set_fontsize(20)

        r0_mean = C.sum(axis=1)/len(s0)*self.scale(r0)

        diff = np.zeros_like(r0)
        for i in range(len(r0)):
            avg = r0_mean[i]
            diff[i] = sum((C[i, :]*self.scale(r0[i]) - avg)**2)

        self.ss.plot(r0, 0.5*np.log(diff) - 0.5*np.log(len(s0) - 1) - np.log(r0_mean), 'ks',
                      linewidth=1,
                      fillstyle='none',
                      markersize=7,
                      markeredgewidth=1.5)

        # s0_min = s0.min()
        # s0_max = s0.max()
        # width = 2
        # height = 2

        # bbox_props = dict(boxstyle="square", fc="white", ec="k", lw=3)
        # self.ss.text(s0_max-width, height, r"$\sigma_0 \in [%g, %g]$" %(s0_min, s0_max),size=self.labelSize, bbox=bbox_props)

        self.ss.set_xlabel(r"$\lambda_D$")
        self.ss.set_ylabel(r"$\log \left[\mathrm{\sigma(K_3; \sigma_0})/K_3(\lambda_d)\right]$")

        self.mean_figure.plot(r0, 0*r0, "r-", label="Analytical", linewidth=3)

        self.mean_figure.plot(r0, -np.log(r0_mean), 'ks',
                              label="KMC",
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5)

        self.mean_figure.legend(loc="upper right", numpoints=1, handlelength=1)

        self.mean_figure.set_xlabel(r"$\lambda_D$")
        self.mean_figure.set_ylabel(r"$-\log K_3(\lambda_D)$")
        self.mean_figure.set_xbound(0)
        # self.mean_figure.set_ybound(-0.45)
        # self.mean_figure.axes.set_yticks([0, 1, 2, 3, 4])



class shifts(DCVizPlotter):

    nametag = "shifts\.arma|values\.arma"

    isFamilyMember = True

    armaBin = True

    hugifyFonts = True

    def plot(self, data):

        cut = 10

        # shifts = data[self.get_family_index_from_name("shifts.arma")][0]
        values = data[self.get_family_index_from_name("values.arma")][0]

        print sum(values[cut:])/(len(values) - cut)

        n = avg1 = avg2 = 0
        for v in values[cut:]:
            avg1 += v
            avg2 += v*v
            n+=1
        std = np.sqrt(1./(n-1)*avg2 - avg1*avg1/(n*(n-1)))
        print std

        # self.subfigure.plot(-shifts,'^', label="shifts")
        self.subfigure.plot([0 for i in range(len(values))], 'r--', linewidth=3, label=r"$\gamma_\mathrm{eq}$")

        self.subfigure.plot(values, "ks",
                            markersize=7,
                            markeredgewidth=1.5,
                            linewidth=1,
                            label=r"$\gamma_\mathrm{KMC}$",
                            fillstyle="none")

        self.subfigure.axes.set_xticks(range(cut + 1))

        self.subfigure.set_xlabel(r"$k$")
        self.subfigure.set_ylabel(r"$\gamma$")
        self.subfigure.axes.set_xlim(-0.1, cut)
        self.subfigure.axes.set_ylim(-0.05, 1.1)
        self.subfigure.legend(numpoints=1, handlelength=0.9)

        # zoomfac = len(values)/2
        #
        # if len(shifts) <= zoomfac + 3:
        #     return
        #
        # X = range(zoomfac, len(shifts))
        #
        # self.zoomed.plot(X, shifts[zoomfac:], marker='x', label="shifts")
        # self.zoomed.plot(X, values[zoomfac:], marker='o', label="values")
        # self.zoomed.plot(X, [0 for i in range(len(values) - zoomfac)], 'r--', label="exact")
        # self.zoomed.set_xlabel("n")
        #
        # avg = 0
        # for start in X:
        #     shift_sum = 0
        #     avg2 = 0
        #     for i, value in enumerate(values[start:]):
        #         avg2 += value
        #         shift_sum += 1
        #     avg2 /= shift_sum
        #
        #     print sum(values[start:])/(len(values) - start), avg2
        #
        #     avg += avg2
        #
        #
        # print avg/len(X), values[-1], sum(values)/len(values)



        # self.zoomed.legend()