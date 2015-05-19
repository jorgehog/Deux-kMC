#DCVIZ

from DCViz_sup import DCVizPlotter

import numpy as np
import os
import re
from numpy import exp
from scipy.stats import linregress
from matplotlib.ticker import FormatStrFormatter, MultipleLocator, FuncFormatter

class SteadyState(DCVizPlotter):

    nametag = "steadystate\_(.*)"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"converge_figure": ["rms_fig", "srms_fig", "dh_fig"], "value_figure": ["subfigure2", "subfigure3", "subfigure4"]}

    stack = "V"

    def adjust(self):
        self.adjust_maps["convergence_figure"]["hspace"] = self.adjust_maps["value_figure"]["hspace"] = 0.2
        self.adjust_maps["value_figure"]["top"] = 0.83

    def trans(self, v):
        return v
        return (v/v.min())**10

    shapes = ["^", "v", "o"]
    plot_values = [0.01, 0.1, 0.5]
    alpha_values = [1, 2, 3]
    om = -0.75
    E0 = 0.01
    alpha = 3

    unloaded = False

    fig_size = [8, 8]

    def plot(self, data):

        clip = 1.1
        start = 1
        transform = False

        rmslabel=r"\sigma(h)"
        sslabel=r"\sigma(s)"
        whlabel=r"\langle \delta h_l \rangle"

        sfigs = [self.subfigure2, self.subfigure3, self.subfigure4]

        E0_values = data[self.get_family_index_from_name("steadystate_E0.npy")]
        alpha_values = data[self.get_family_index_from_name("steadystate_alpha.npy")]
        mu_shift_values = data[self.get_family_index_from_name("steadystate_mu_shift.npy")]

        I = list(E0_values).index(self.E0)
        J = list(alpha_values).index(self.alpha)
        K = [round(x, 2) for x in list(exp(mu_shift_values))].index(1 + self.om)

        if self.unloaded:
            unloaded = data[self.get_family_index_from_name("steadystate_unloaded.dat")]

            the_alpha = alpha_values[J]

            unloaded_rms = []
            for alpha, mu, rms in unloaded:

                if alpha == the_alpha:
                    unloaded_rms.append([mu, rms])

            _, unloaded_rms = zip(*(sorted(unloaded_rms, key=lambda x: x[0])))

            self.subfigure2.plot(exp(mu_shift_values) - 1, unloaded_rms,
                                 label=r"$E_0=0.00$",
                                 marker="s",
                                 markeredgecolor="k",
                                 linestyle="--",
                                 color="r",
                                 linewidth=1,
                                 fillstyle='none',
                                 markersize=7,
                                 markeredgewidth=1.5)
        else:
            self.shapes = ["s"] + self.shapes

        count = 0
        slopes = []
        FFS = [0,0,0]
        for i, E0 in enumerate(E0_values):

            for j, alpha in enumerate(alpha_values):

                rms_values = []

                for k, mu_shift in enumerate(mu_shift_values):

                    if i == I and j == J and k == K:
                        print "E0 = %g, alpha= %g, mus = %g (om = %g)" % (E0, alpha, mu_shift, exp(mu_shift))

                    time = data[self.get_family_index_from_name("steadystate_Time_%d.npy" % count)]
                    # n = 1 + np.arange(len(time))*1000

                    L = int(len(time)/clip)
                    time = time[start:L]

                    n = time

                    rms = data[self.get_family_index_from_name("steadystate_HeightRMS_%d.npy" % count)][start:L]

                    ss = data[self.get_family_index_from_name("steadystate_SurfaceSize_%d.npy" % count)][start:L]

                    wh = data[self.get_family_index_from_name("steadystate_PressureWall_%d.npy" % count)][start:L]

                    if i == I and j == J and k == K:

                        if transform:
                            self.subfigure.loglog(n, rms, label=rmslabel)
                            self.subfigure.loglog(n, (ss - ss[0])/(ss.max() - ss[0]), label=sslabel)
                            self.subfigure.loglog(n, (wh - wh[0])/(wh.max() - wh[0]), label=whlabel)
                        else:
                            self.rms_fig.plot(n, np.log(self.trans(rms)), "--", label=rmslabel,
                                              linewidth=1,
                                              fillstyle='none',
                                              markersize=7,
                                              markeredgewidth=1.5,
                                              color="red")

                            self.dh_fig.plot(n, np.log(self.trans(wh)), "--", label=whlabel,
                                             linewidth=1,
                                             fillstyle='none',
                                             markersize=7,
                                             markeredgewidth=1.5,
                                             color="red")

                            self.srms_fig.plot(n, np.log(self.trans(ss)), "--", label=sslabel,
                                               linewidth=1,
                                               fillstyle='none',
                                               markersize=7,
                                               markeredgewidth=1.5,
                                               color="red")




                    L = len(rms)
                    L2 = L/3

                    rms_value = rms[L-L2:].mean()
                    ss_value = ss[L-L2:].mean()
                    wh_value = wh[L-L2:].mean()


                    rms_values.append(rms_value)

                    count += 1

                    if i == I and j == J and k == K:

                        rms_value_chosen = rms_value
                        srms_value_chosen = ss_value
                        dh_value_chosen = wh_value

                if E0 in self.plot_values and alpha in self.alpha_values:

                    n_sfig = self.alpha_values.index(alpha)
                    sfigs[n_sfig].plot(exp(mu_shift_values) - 1, rms_values,
                                                               label=r"$E_0=%.2f$" % E0,
                                                               marker=self.shapes[FFS[n_sfig]],
                                                               markeredgecolor="k",
                                                               linestyle="--",
                                                               color="r",
                                                               linewidth=1,
                                                               fillstyle='none',
                                                               markersize=7,
                                                               markeredgewidth=1.5)
                    FFS[n_sfig] += 1
            print E0

            slope, intercept, _, _, _ = linregress(exp(mu_shift_values), rms_values)
            slopes.append(slope)


        self.dh_fig.set_xlabel(r"$t\nu c$")

        self.srms_fig.axes.set_xscale('log')
        self.rms_fig.axes.set_xscale('log')
        self.dh_fig.axes.set_xscale('log')

        self.srms_fig.axes.xaxis.set_ticklabels([])
        self.rms_fig.axes.xaxis.set_ticklabels([])

        majorFormatter = FuncFormatter(self.formatter)
        self.rms_fig.axes.yaxis.set_visible(False)
        self.srms_fig.axes.yaxis.set_visible(False)
        self.dh_fig.axes.yaxis.set_visible(False)


        # majorFormatter = FormatStrFormatter('%.1f')
        # self.rms_fig.axes.yaxis.set_major_formatter(majorFormatter)
        # self.srms_fig.axes.yaxis.set_major_formatter(majorFormatter)
        # self.dh_fig.axes.yaxis.set_major_formatter(majorFormatter)
        #
        # majorLocator = MultipleLocator(1)
        # self.rms_fig.axes.yaxis.set_major_locator(majorLocator)
        # self.srms_fig.axes.yaxis.set_major_locator(majorLocator)
        # self.dh_fig.axes.yaxis.set_major_locator(majorLocator)

        xmax = 3E3
        self.rms_fig.set_xlim(0, xmax)
        self.srms_fig.set_xlim(0, xmax)
        self.dh_fig.set_xlim(0, xmax)


        rms_span = self.rms_fig.get_ylim()[1] + self.rms_fig.get_ylim()[0]
        srms_span = self.srms_fig.get_ylim()[1] + self.srms_fig.get_ylim()[0]
        dh_span = self.dh_fig.get_ylim()[1] + self.dh_fig.get_ylim()[0]

        xtext = 0.25E3

        self.rms_fig.text(xtext, rms_span/2, r"$%s = %.3f$" % (rmslabel, rms_value_chosen), fontsize=self.fontSize)
        self.srms_fig.text(xtext, srms_span/2, r"$%s = %.3f$" % (sslabel, srms_value_chosen), fontsize=self.fontSize)
        self.dh_fig.text(xtext, dh_span/2, r"$%s = %.3f$" % (whlabel, dh_value_chosen), fontsize=self.fontSize)

        r = [0.6, 0.65, 0.65]
        for i, sfig in enumerate(sfigs):
            sfig.set_xlim([-1, 2.1])
            sfig.axes.get_xaxis().set_ticks([-1, 0, 1, 2])
            sfig.axes.yaxis.set_major_formatter(majorFormatter)
            sfig.set_ylabel(rmslabel)
            sfig.set_ybound(0)

            if i != 2:
                sfig.axes.xaxis.set_ticklabels([])

            ax = sfig.axes.twinx()
            ax.set_ylabel(r"$\alpha=%d$" % self.alpha_values[i], labelpad=15)
            ax.yaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])

            # end = sfig.axes.get_ylim()[1]
            # start = sfig.axes.get_ylim()[0]
            #
            # loc_y = start + r[i]*(end - start)
            # loc_x = 1.5
            #
            # sfig.text(loc_x, loc_y, r"$\alpha=%d$" % self.alpha_values[i], fontsize=self.fontSize)


        self.subfigure4.set_xlabel(r"$\Omega = c/c_\mathrm{eq} - 1$")
        self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1.2, ncol=3, columnspacing=0.3, handletextpad=0.5, borderaxespad=0.3,  bbox_to_anchor=(-0.01, 1.55))

    @staticmethod
    def formatter(value, index):

        if index is None:
            return ""

        if index % 2 == 0:
            return r"$%g$" % value
        else:
            return ""

class UnloadedGrowthSpeed(DCVizPlotter):

    nametag = "unloaded_growthspeed\_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"omega_vs_v" : "subfigure",
              "slopes_vs_alpha": "subfigure2",
              "shifts_vs_alpha": "subfigure3",
              "omega_vs_v_param": "subfigure4",
              "n_Fig": "nfig"}

    def plot(self, data):

        alpha_values = self.get_family_member_data(data, "alpha")
        mu_values = self.get_family_member_data(data, "mu_shift")
        v_values = self.get_family_member_data(data, "v")

        c_over_c0 = exp(mu_values)
        omega = c_over_c0 - 1

        idx = np.where(omega <= 0)
        omega_growth = (1-omega[idx])[::-1]
        print omega_growth

        omega_powers = []
        log_ks = []

        for i, alpha in enumerate(alpha_values):

            v = v_values[i, :]*c_over_c0
            v_growth = -v[idx][::-1]
            print v_growth

            self.subfigure.loglog(omega_growth, v_growth, "--o", label=r"$\alpha=%.2f$" % alpha)
            slope, shift, _, _, _ = linregress(np.log(omega_growth), np.log(v_growth))

            omega_powers.append(slope)
            log_ks.append(shift)

        log_ks = np.asarray(log_ks)

        self.subfigure2.plot(alpha_values, omega_powers, "--o")
        self.subfigure3.loglog(alpha_values, -log_ks, "--o")

        omega_power_slope, power_zero_alpha, _, _, r1 = linregress(alpha_values, omega_powers)
        slope, shift, _, _, r2 = linregress(np.log(alpha_values), np.log(-log_ks))

        print omega_power_slope, power_zero_alpha, r1
        print slope, exp(shift), r2

        log_k = lambda alpha: -exp(shift)*alpha**slope

        omega_full = np.linspace(omega_growth[0], omega_growth[-1])
        for i, alpha in enumerate(alpha_values):

            v = v_values[i, :]*c_over_c0
            v_growth = v[idx]

            self.subfigure4.plot(omega_growth, v_growth, "o", label=r"$\alpha=%.2f$ n=%.2f" % (alpha, omega_powers[i]))
            self.subfigure4.plot(omega_full, exp(log_k(alpha))*omega_full**omega_powers[i], "k--")
            self.subfigure4.plot(omega_full, exp(log_k(alpha))*omega_full, "k-.")



        self.subfigure.legend(loc="upper left")
        self.subfigure.set_xlabel(r"$\Omega = c/c_\mathrm{eq} - 1$")
        self.subfigure.set_ylabel(r"v")

        self.subfigure2.set_xlabel(r"$\alpha")
        self.subfigure2.set_ylabel(r"n")

        self.subfigure3.set_xlabel(r"$\alpha$")
        self.subfigure3.set_ylabel(r"$-\log k_g$")

        self.subfigure4.legend(loc="upper left")
        self.subfigure4.set_xlabel(r"$\Omega$")
        self.subfigure4.set_ylabel(r"v")

        n_values = self.get_family_member_data(data, "n") - 2
        print n_values
        for i, alpha in enumerate(alpha_values):
            self.nfig.plot(omega, n_values[i, :], "--o", label=r"$\alpha=%.2f$" % alpha)
        self.nfig.set_xlabel(r"$\Omega$")
        self.nfig.set_ylabel(r"$\langle n \rangle$")
        self.nfig.legend()

class GrowthSpeed(DCVizPlotter):

    nametag = "^growthspeed\_(.*)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    # figMap = {"omega_vs_v": "subfigure", "slopes": "subfigure2", "mu0": "subfigure3"}

    figMap = {"omega_vs_v": "subfigure",
              "slopes": "subfigure2",
              "logk_slopes": "subfigure4",
              "logk_cutz": "subfigure5",
              "alpha_cuts_all": "subfigure6",
              "alpha_cutz_r0": "subfigure7",
              "alpha_slopes_full": "subfigure8",
              "alpha_slopes_comb": "subfigure9",
              "neighborstuff": "subfigure10",
              "abs_stuff": "subfigure11"}

    # figMap = {"omega_vs_v": ["subfigure",
    #           "subfigure2",
    #           "subfigure4",
    #           "subfigure5",
    #           "subfigure6",
    #           "subfigure7",
    #           "subfigure8",
    #           "subfigure9"]}

    # figMap = {"asd": "subfigure"}

    plot_values = [0.5, 1.0]
    shapes = ["s", "^", "v"]

    def plot_and_slopify(self, E0, alpha, omega, mu_shifts, mu, v):
        mu0 = (mu - mu_shifts).mean()
        c_over_c0 = exp(mu0 + mu_shifts)
        v *= c_over_c0*exp(-2*alpha)

        if E0 == 0:
            if not (omega == (c_over_c0-1)).all():

                print "OMEGA FAIL"
                print omega
                print c_over_c0
                self.Exit()

        idx_high = np.where(omega >= 0)
        idx_vhigh = np.where(omega >= 2)
        idx_low = np.where(omega < 0)
        idx_all = np.where(omega > -2)

        idx_chosen = idx_vhigh

        slope, intercept, _, _, err = linregress(omega[idx_chosen], v[idx_chosen])
        # slope_low, intercept, _, _, _ = linregress(omega[idx_low], v[idx_low])
        #
        # print E0, slope/slope_low - 1

        return mu0, slope, v, err

    def uberplot(self, omega, v, k, E0):

        idx_undersat = np.where(omega <= 0)
        idx_oversat = np.where(omega >= 0)

        omega_undersat_list = list(omega[idx_undersat])
        omega_oversat_list = list(omega[idx_oversat])

        v_undersat_list = list(v[idx_undersat])
        v_oversat_list = list(v[idx_oversat])

        omega_undersat = abs(omega[idx_undersat])
        omega_oversat = omega[idx_oversat]

        v_undersat = abs(v[idx_undersat])
        v_oversat = v[idx_oversat]

        v = v_undersat
        omega = omega_undersat

        v_log_under = np.sign(v)*np.log(abs(v))
        o_log_under = np.sign(omega)*np.log(abs(omega))
        self.subfigure11.plot(o_log_under, v_log_under, "k--" + self.shapes[k],
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="black")

        v = v_oversat
        omega = omega_oversat

        v_log_over = np.sign(v)*np.log(abs(v))
        o_log_over = np.sign(omega)*np.log(abs(omega))

        xshift = o_log_over[-1] - o_log_under[0]
        yshift = v_log_over[-1] - v_log_under[0]

        xshift = 0
        yshift = 0

        self.subfigure11.plot(o_log_over - xshift, v_log_over - yshift, "k-." + self.shapes[k],
                            label=r"$E_0=%.1f$" % E0,
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="black")

        self.subfigure.plot(omega_undersat_list[::3] + omega_oversat_list,
                            v_undersat_list[::3] + v_oversat_list,
                            label=r"$E_0=%.1f$" % E0,
                            marker=self.shapes[k],
                            markeredgecolor="k",
                            linestyle="--",
                            color="r",
                            linewidth=1,
                            fillstyle='none',
                            markersize=7,
                            markeredgewidth=1.5)

    def adjust(self):

        for map in self.adjust_maps.values():
            map["bottom"] = 0.2
            map["top"] = 0.94

    l_min = 2

    def plot(self, data):

        E0_values = self.get_family_member_data(data, "E0")
        alpha_values = self.get_family_member_data(data, "alpha")
        mu_shift_values = self.get_family_member_data(data, "mu_shift")
        r0_values = self.get_family_member_data(data, "r0")
        s0_values = self.get_family_member_data(data, "s0")
        mu_values = self.get_family_member_data(data, "mu")
        v_values = self.get_family_member_data(data, "v")
        n_values = self.get_family_member_data(data, "n")

        print E0_values.shape
        print alpha_values.shape
        print mu_shift_values.shape
        print r0_values.shape
        print s0_values.shape
        print mu_values.shape
        print v_values.shape
        print n_values.shape

        omega = exp(mu_shift_values) - 1

        k_values = np.zeros(shape=(len(E0_values), len(alpha_values), len(r0_values), len(s0_values)))
        mu_zero = np.zeros_like(k_values)
        errors = np.zeros_like(k_values)

        J, L, M = [int(x) for x in self.argv]

        print "alpha=%g, s0 = %g, r0 = %g" % (alpha_values[J], s0_values[M+1], r0_values[L+1])

        for j in range(len(alpha_values)):
            _, slope0, v0, err = self.plot_and_slopify(0, alpha_values[j], omega, mu_shift_values, mu_shift_values, v_values[0, j, :, 0, 0])

            k_values[0, j, :, :] = slope0
            mu_zero[0, j, :, :] = 0
            errors[0, j, :, :] = err

            if j == J:
                self.uberplot(exp(mu_shift_values) - 1, v0, 0, 0)
                # self.subfigure.loglog(abs(omega), abs(v0), "k--" + shapes[0], label="E0=0",
                #                     linewidth=1,
                #                     fillstyle='none',
                #                     markersize=7,
                #                     markeredgewidth=1.5,
                #                     color="black")

        k = 1
        for i, E0 in enumerate(E0_values[1:]):
            for j, alpha in enumerate(alpha_values):
                # for k, mu_shift in enumerate(mu_shift_values):
                    for l, r0 in enumerate(r0_values[1:]):
                        for m, s0 in enumerate(s0_values[1:]):

                            mu0, slope, v, error = self.plot_and_slopify(E0,
                                                                         alpha,
                                                                  omega,
                                                                  mu_shift_values,
                                                                  mu_values[i+1, j, :, l+1, m+1],
                                                                  v_values[i+1, j, :, l+1, m+1])

                            if slope < 0:
                                print "ERROR slope=%g, E0(%d)=%g, alpha(%d)=%g, r0(%d)=%g, s0(%d)=%g" % (slope, i, E0, j, alpha, l, r0, m, s0)

                            k_values[i+1][j][l+1][m+1] = slope
                            errors[i+1][j][l+1][m+1] = error
                            mu_zero[i+1][j][l+1][m+1] = mu0

                            if E0 in self.plot_values and j == J and l == L and m == M:
                                # self.subfigure.plot(omega, slope*omega, "r--")
                                self.uberplot(omega, v, k, E0)
                                # self.subfigure.loglog(abs(omega), abs(v), "k--" + shapes[k], label="E0=%g" % E0,
                                #                     linewidth=1,
                                #                     fillstyle='none',
                                #                     markersize=7,
                                #                     markeredgewidth=1.5,
                                #                     color="black")
                                k += 1

                            if i == 3 and l == L and m == M:
                                self.subfigure10.plot(omega, n_values[i+1, j, :, l+1, m+1], label="a=%.2f" % alpha)

        log_k_E0_slopes = np.zeros(shape=(len(alpha_values), len(r0_values), len(r0_values)))
        log_k_cutz = np.zeros_like(log_k_E0_slopes)
        log_k_slope_error = np.zeros_like(log_k_E0_slopes)

        ka = 0
        na = 4
        nm = na - 2
        d = len(alpha_values)/nm
        for j, alpha in enumerate(alpha_values):
            for l, r0 in enumerate(r0_values[1:]):
                F0_values = E0_values/(r0*(1 - np.exp(-1./r0)))

                for m, s0 in enumerate(s0_values[1:]):

                    if l == L and m == M:
                        if j == 0 or j == len(alpha_values) - 1 or j % d == 0:

                            self.subfigure2.plot(E0_values, k_values[:, j, l+1, m+1],
                                                 label=r"$\alpha=%g$" % round(alpha, 1),
                                                 marker=self.shapes[ka],
                                                 markeredgecolor="k",
                                                 linestyle="--",
                                                 color="r",
                                                 linewidth=1,
                                                 fillstyle='none',
                                                 markersize=7,
                                                 markeredgewidth=1.5)

                            ka += 1

                    log_k_slope, intercept, _, _, err = linregress(E0_values, np.log(k_values[:, j, l+1, m+1]))

                    log_k_E0_slopes[j, l+1, m+1] = log_k_slope
                    log_k_cutz[j, l+1, m+1] = intercept
                    log_k_slope_error[j, l+1, m+1] = err

        alpha_slopes = np.zeros(shape=(len(r0_values), len(s0_values)))
        alpha_cutz = np.zeros_like(alpha_slopes)

        for l, r0 in enumerate(r0_values[1:]):
                for m, s0 in enumerate(s0_values[1:]):
                    alpha_slope, intercept, _, _, err = linregress(alpha_values, log_k_E0_slopes[:, l+1, m+1])
                    print "slope", alpha_slope

                    alpha_slopes[l+1, m+1] = alpha_slope
                    alpha_cutz[l+1, m+1] = intercept


        self.subfigure4.plot(alpha_values, log_k_E0_slopes[:, L+1, M+1],
                             marker="s",
                             markeredgecolor="k",
                             linestyle="--",
                             color="r",
                             linewidth=1,
                             fillstyle='none',
                             markersize=7,
                             markeredgewidth=1.5)
        self.subfigure4.plot([0, alpha_values[0]],
                             [alpha_cutz[L+1, M+1], alpha_cutz[L+1, M+1] + alpha_slopes[L+1, M+1]*alpha_values[0]],
                             "r--")
        # X = 1
        # self.subfigure4.plot([0, X], [X, X], "k:")
        # self.subfigure4.plot([X, X], [0, X], "k:")


        S_tot = 0
        count2 = 0
        for l, r0 in enumerate(r0_values[1:]):
            S = 0
            count = 0

            if r0 < self.l_min:
                continue

            for m, s0 in enumerate(s0_values[1:]):
                S += log_k_cutz[:, l+1, m+1]
                count += 1

                self.subfigure5.loglog(alpha_values, -log_k_cutz[:, l+1, m+1], "kx",
                                     linewidth=1,
                                     fillstyle='none',
                                     markersize=7,
                                     markeredgewidth=1.5)


            S_tot += S
            count2 += count

        S_tot /= count2

        print alpha_values, S_tot

        power, log_constant, _, _, err = linregress(np.log(alpha_values), np.log(-S_tot))
        print "Power=%g, constant=%g, err=%g" % (power, exp(log_constant), err)

        self.subfigure5.loglog(alpha_values, exp(log_constant)*alpha_values**power, "g-^",
                               label=r"$%.2f\alpha^{%.2f}$" % (exp(log_constant), power))
        self.subfigure5.loglog(alpha_values, -S_tot, "r-^", label="Avg all param")

        self.subfigure5.legend(numpoints=1, handlelength=1.2, borderaxespad=0.3, )

        for l, r0 in enumerate(r0_values[1:]):
            self.subfigure6.plot(s0_values[1:], alpha_cutz[l+1, 1:], label="r0 = %g" % r0)

        self.subfigure6.set_xlabel(r"$\sigma_0$")
        self.subfigure6.set_ylabel(r"$\mathrm{alpha shift}$")
        # self.subfigure6.legend()

        cuts = alpha_cutz[1:, 1:].mean(axis=1)

        print r0_values[1:], cuts
        self.subfigure7.plot(r0_values[1:], cuts)
        self.subfigure7.set_xlabel(r"$\lambda_D$")
        self.subfigure7.set_ylabel(r"$\mathrm{avg alpha shift}$")

        for l, r0 in enumerate(r0_values[1:]):
            self.subfigure8.plot(s0_values[1:], alpha_slopes[l+1, 1:], label="r0=%g" % r0)
        # self.subfigure8.legend()
        self.subfigure8.set_xlabel(r"$\sigma_0$")
        self.subfigure8.set_ylabel(r"$\mathrm{k_values}$")

        comb = alpha_slopes[1:, 1:].mean(axis=1)

        self.subfigure9.plot(r0_values[1:], comb*(r0_values[1:]*(1 - np.exp(-1./r0_values[1:]))))
        self.subfigure9.set_xlabel(r"$\lambda_D$")
        self.subfigure9.set_ylabel(r"$\mathrm{avg slope}$")

        self.subfigure.set_xlabel(r"$\Omega = c/c_\mathrm{eq} - 1$")
        self.subfigure.set_ylabel(r"$\dot{H} = \Delta \langle h\rangle / \Delta t$")
        self.subfigure.legend(loc="upper left", numpoints=1, handlelength=1.2, borderaxespad=0.3 )
        self.subfigure.set_xlim(-1, 3)
        self.subfigure.axes.xaxis.set_ticks([-1, 0, 1, 2, 3])
        self.subfigure.axes.xaxis.set_major_formatter(FormatStrFormatter(r"$%d$"))
        self.subfigure.set_ylim(-0.2, 1)
        # self.subfigure.axes.set_xscale('log')
        # self.subfigure.axes.set_yscale('log')

        self.subfigure2.legend(loc="lower right", numpoints=1, handlelength=1.2, ncol=3, columnspacing=0.3, handletextpad=0.5, borderaxespad=0.3)
        self.subfigure2.axes.set_yscale('log')
        self.subfigure2.set_xlim(-0.05, 1.05)
        self.subfigure2.set_ybound(1E-4)

        self.subfigure2.set_xlabel(r"$E_0 \propto F_0/L$")
        self.subfigure2.set_ylabel(r"$k_g = \dot{H} / \Omega$")

        self.subfigure4.set_xlabel(r"$\alpha = E_b/kT$")
        self.subfigure4.set_ylabel(r"$\log (k_g/k_0) / E_0$")
        self.subfigure4.set_ybound(0)
        self.subfigure4.set_xlim(0, alpha_values.max()*1.1)


        self.subfigure5.set_xlabel(r"$\alpha = E_b/kT$")
        self.subfigure5.set_ylabel(r"$\mathrm{log k shift}$")

        self.subfigure11.set_xlabel(r"$\log (\Omega + 1)$")
        self.subfigure11.set_ylabel(r"$\log(|\dot{H}|)$")
        self.subfigure11.legend(loc="upper left")



        # self.subfigure3.set_xlabel(r"$E_0$")
        # self.subfigure3.set_ylabel(r"$\gamma_\mathrm{eq}$")


class Quasi2D_slopes_and_stuff(DCVizPlotter):

    nametag = "linearplots\_(.*)\_?\d*\.npy"

    numpyBin = True

    isFamilyMember = True

    figMap = {"fig" : "gammaslopes", "fig2" : "E0slopes"}

    hugifyFonts = True

    def adjust(self):

        for map in self.adjust_maps.values():
            map["top"] = 0.95
            map["bottom"] = 0.18

    fig_size = [8, 4]

    def n_runs(self):

        n_max = -1

        for name in self.familyFileNames:
            if not "alpha" in name:
                continue

            n = int(re.findall("linearplots_alpha_(\d+)\.npy", name)[0])

            if n > n_max:
                n_max = n

        return n_max + 1

    plot_values = [0.01, 0.1, 0.19]

    def plot(self, data):

        meta_data_file = os.path.join(self.familyHome, "linearplots_metadata.dat")

        with open(meta_data_file, 'r') as meta_data:
            print meta_data.read()

        N = 3
        shapes = ["^", "v", 'o']

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

            if E0s[n] not in self.plot_values:
                continue
            print n, n_plots

            mu_error_name = "linearplots_muEqErrors_%d.npy" % n

            mu_errors = data[self.get_family_index_from_name(mu_error_name)]

            self.gammaslopes.errorbar(alphas, mus,
                                      yerr=mu_errors,
                                      fmt=shapes[n_plots],
                                      fillstyle='none',
                                      label=r"$E_0=%1.2f$" % (E0s[n]),
                                      markersize=7,
                                      markeredgewidth=1.5,
                                      linewidth=1,
                                      color="black")

            kws = {}
            if n_plots == 0:
                kws["label"] = r"$\mathrm{Linear\,\,fit}$"

            self.gammaslopes.plot([0, alphas.max()], [0, slopes[n]*alphas.max()], "r--", **kws)

            n_plots += 1


        self.gammaslopes.set_xbound(0)
        self.gammaslopes.set_ybound(0)
        self.gammaslopes.legend(loc="upper left", numpoints=1, handlelength=0.8, ncol=2, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3)



        self.gammaslopes.set_xlabel(r"$\alpha = E_b/kT$")
        self.gammaslopes.set_ylabel(r"$\gamma_\mathrm{eq} = \ln c_\mathrm{eq}/c_0$")

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

        self.E0slopes.axes.yaxis.set_ticks([0, 0.1, 0.2, 0.3])
        majorFormatter = FormatStrFormatter(r'$%.1f$')
        self.E0slopes.axes.yaxis.set_major_formatter(majorFormatter)

        self.E0slopes.set_xlabel(r"$E_0 \propto F_0/L$")
        self.E0slopes.set_ylabel(r"$\gamma_\mathrm{eq}/\alpha$")

        print sslope, d


class SOS_pressure_sizes(DCVizPlotter):

    nametag = "pressure_plots_.*\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"sfig" : ["subfigure", "subfigure2", "subfigure3"]}

    # share_axis = True

    def adjust(self):
        for map in self.adjust_maps.values():
            map["hspace"] = 0.15
            map["top"] = 0.83

    fig_size = [8, 8]

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
            self.subfigure.loglog(1./analytical[:, 0], analytical[:, 1], 'r-',
                                linewidth=3,
                                # label="$E_0 = 0.00$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure2.loglog(1./analytical[:, 0], analytical[:, 2], 'r-',
                                linewidth=3,
                                # label="$E_0 = 0.00$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure3.loglog(1./analytical[:, 0], (analytical[:, 0]*analytical[:, 2])**p, 'r-',
                                linewidth=3,
                                # label="$E_0 = 0.00$",
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

            self.subfigure.loglog(1./alpha_array, mean_s_array, 'k%s' % shapes[nplots],
                                 fillstyle='none',
                                 label="$E_0 = %.2f$" % E0_value,
                                 markersize=7,
                                 markeredgewidth=1.5,
                                 linewidth=1)

            self.subfigure2.loglog(1./alpha_array, var_s_array, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$E_0 = %.2f$" % E0_value,
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            self.subfigure3.loglog(1./alpha_array, (var_s_array*alpha_array)**p, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$E_0 = %.2f$" % E0_value,
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            nplots += 1

        xmin = 0.5
        xmax = 6

        # self.subfigure.set_xlabel(r"$1/\alpha$")
        self.subfigure.set_xlim(xmin, xmax)
        # self.subfigure.axes.xaxis.set_visible(False)
        self.subfigure.axes.xaxis.set_ticklabels([])
        self.subfigure.set_ylabel(r"$\langle s_{\uparrow\downarrow} \rangle / L$")
        ax = self.subfigure.axes.twinx()
        ax.set_ylabel(r"$\langle E \rangle /E_bL$", labelpad=15)
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])
        self.subfigure.legend(loc="upper left", numpoints=1, handlelength=0.5, ncol=3, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3,  bbox_to_anchor=(0, 1.55))

        # self.subfigure2.set_xlabel(r"$1/\alpha$")
        self.subfigure2.set_xlim(xmin, xmax)
        self.subfigure2.set_ylabel(r"$\sigma (s_{\uparrow\downarrow}) / L$")
        # self.subfigure2.axes.xaxis.set_visible(False)
        self.subfigure2.axes.xaxis.set_ticklabels([])
        ax2 = self.subfigure2.axes.twinx()
        ax2.set_ylabel(r"$\sigma (E) / E_bL$", labelpad=15)
        ax2.yaxis.set_ticks([])
        ax2.yaxis.set_ticklabels([])
        # self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1)

        self.subfigure3.set_xlabel(r"$1/\alpha = kT/E_b$")
        self.subfigure3.set_xlim(xmin, xmax)
        self.subfigure3.axes.xaxis.set_ticklabels([r"$0.01$", r"$0.1$", r"$1$", r'$10$'])

        self.subfigure3.axes.xaxis.set_minor_formatter(FuncFormatter(self.axisFormat))
        self.subfigure3.set_ylabel(r"$\alpha^2\sigma(s_{\uparrow\downarrow})^2/L^2$")
        ax3 = self.subfigure3.axes.twinx()
        ax3.set_ylabel(r"$C_V/k$", labelpad=15)
        ax3.yaxis.set_ticks([])
        ax3.yaxis.set_ticklabels([])
        # self.subfigure3.legend(numpoints=1, handlelength=1, loc="upper right")

    @staticmethod
    def axisFormat(value, index):

        if index % 2 == 0:
            return r'$%g$' % value
        else:
            return ''

class SOSanalyze(DCVizPlotter):

    nametag = "analyze_(.+)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"figure_2" : "mean_figure"}

    fig_size = [8, 4]

    def adjust(self):
        for map in self.adjust_maps.values():
            map["top"] = 0.95
            map["bottom"] = 0.18

    def scale(self, r0):
        return r0*(1-exp(-1./r0))

    def plot(self, data):

        dirname = re.findall("analyze\_(.+)\_.+\_values\.npy", self.familyHead)[0]

        print dirname

        C =  data[self.get_family_index_from_name("analyze_%s_C_values.npy" % dirname)]
        s0 =  data[self.get_family_index_from_name("analyze_%s_s0_values.npy" % dirname)]
        r0 =  data[self.get_family_index_from_name("analyze_%s_r0_values.npy" % dirname)]

        # C = 1./C

        print r0.shape, s0.shape, C.shape
        r0_mean = C.sum(axis=0)/len(s0)

        # diff = np.zeros_like(r0)
        # for i in range(len(r0)):
        #     avg = r0_mean[i]
        #     diff[i] = sum((C[i, :]*self.scale(r0[i]) - avg)**2)
        #
        # self.ss.plot(r0, 0.5*np.log(diff) - 0.5*np.log(len(s0) - 1) - np.log(r0_mean), 'ks',
        #               linewidth=1,
        #               fillstyle='none',
        #               markersize=7,
        #               markeredgewidth=1.5)

        # s0_min = s0.min()
        # s0_max = s0.max()
        # width = 2
        # height = 2

        # bbox_props = dict(boxstyle="square", fc="white", ec="k", lw=3)
        # self.ss.text(s0_max-width, height, r"$\sigma_0 \in [%g, %g]$" %(s0_min, s0_max),size=self.labelSize, bbox=bbox_props)
        #
        # self.ss.set_xlabel(r"$\lambda_D$")
        # self.ss.set_ylabel(r"$\log \left[\mathrm{\sigma(K_3; \sigma_0})/K_3(\lambda_d)\right]$")

        self.mean_figure.plot(r0, self.scale(r0), "r-", label=r"$\mathrm{Analytical}$", linewidth=3)

        self.mean_figure.plot(r0, r0_mean, 'ks',
                              label=r"$\mathrm{KMC}$",
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5)

        self.mean_figure.legend(loc="lower right", numpoints=1, handlelength=0.8, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3)

        self.mean_figure.set_xlabel(r"$\lambda_D/l_0$")
        self.mean_figure.set_ylabel(r"$\gamma_\mathrm{eq}/\alpha E_0$")
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