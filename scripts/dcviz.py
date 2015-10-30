#DCVIZ

from DCViz_sup import DCVizPlotter

import numpy as np
import os
import re
import sys
from numpy import exp
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatter, MultipleLocator, FuncFormatter

from mpl_toolkits.mplot3d import Axes3D

E0_tex = r"P_\lambda"
fig_4_x = 0.05
fig_4_y = 0.85
fig_4_fs = 30

class SteadyState(DCVizPlotter):

    nametag = "steadystate\_(.*)"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"converge_figure": ["rms_fig", "srms_fig", "dh_fig"],
              "value_figure": ["subfigure1", "subfigure2", "subfigure3", "subfigure4"]}

    stack = "V"

    fig_size = [8, 8]

    def adjust(self):
        self.adjust_maps["value_figure"]["hspace"] = 0.15
        self.adjust_maps["converge_figure"]["hspace"] = 0.15
        self.adjust_maps["value_figure"]["top"] = 0.88

    def trans(self, v):
        return v
        return (v/v.min())**10

    shapes = ["^", "v", "o"]
    plot_values = [0.01, 0.1, 0.5]
    alpha_values = [0.5, 1, 2, 3]
    om = -0.75
    E0 = 0.01
    alpha = 3

    unloaded = False

    transparent = True

    def plot(self, data):

        clip = 1.1
        start = 1
        transform = False

        rmslabel=r"\sigma(h)"
        sslabel=r"\sigma(s)"
        whlabel=r"\langle \delta h_l \rangle"

        sfigs = [self.subfigure1, self.subfigure2, self.subfigure3, self.subfigure4]

        E0_values = data[self.get_family_index_from_name("steadystate_E0.npy")]
        alpha_values = data[self.get_family_index_from_name("steadystate_alpha.npy")]
        mu_shift_values = data[self.get_family_index_from_name("steadystate_mu_shift.npy")]
        eqMu_values = data[self.get_family_index_from_name("steadystate_eqMu.npy")]

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
                                 label=r"$%s=0.00$" % E0_tex,
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
        FFS = [0,0,0,0]
        for i, E0 in enumerate(E0_values):

            for j, alpha in enumerate(alpha_values):

                muEq = eqMu_values[i, j]

                rms_values = []

                for k, mu_shift in enumerate(mu_shift_values):

                    log_c = mu_shift+muEq -2*alpha

                    if i == I and j == J and k == K:
                        print "E0 = %g, alpha= %g, mus = %g (om = %g)" % (E0, alpha, mu_shift, exp(mu_shift))

                    log_time = np.log(data[self.get_family_index_from_name("steadystate_Time_%d.npy" % count)])
                    # n = 1 + np.arange(len(time))*1000

                    L = int(len(log_time)/clip)
                    log_time = log_time[start:L]

                    n = exp(log_time - log_c)

                    rms = data[self.get_family_index_from_name("steadystate_HeightRMS_%d.npy" % count)][start:L]


                    #multiply by two to convert from s_ud to s
                    ss = data[self.get_family_index_from_name("steadystate_SurfaceSize_%d.npy" % count)][start:L]*2

                    wh = data[self.get_family_index_from_name("steadystate_PressureWall_%d.npy" % count)][start:L]

                    if i == I and j == J and k == K:

                        if transform:
                            self.subfigure.loglog(n, rms, label=r"$%s$" % rmslabel)
                            self.subfigure.loglog(n, (ss - ss[0])/(ss.max() - ss[0]), label=r"$%s$" % sslabel)
                            self.subfigure.loglog(n, (wh - wh[0])/(wh.max() - wh[0]), label=r"$%s$" % whlabel)
                        else:
                            self.rms_fig.plot(n, np.log(self.trans(rms)), "-", label=r"$%s$" % rmslabel,
                                              linewidth=2,
                                              fillstyle='none',
                                              markersize=7,
                                              markeredgewidth=1.5,
                                              color="red")

                            self.dh_fig.plot(n, np.log(self.trans(wh)), "-", label=r"$%s$" % whlabel,
                                             linewidth=2,
                                             fillstyle='none',
                                             markersize=7,
                                             markeredgewidth=1.5,
                                             color="red")

                            self.srms_fig.plot(n, np.log(self.trans(ss)), "-", label=r"$%s$" % sslabel,
                                               linewidth=2,
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
                                                               label=r"$%s=%.2f$" % (E0_tex, E0),
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


        self.dh_fig.set_xlabel(r"$t\, [\nu]$")

        self.srms_fig.axes.set_xscale('log')
        self.rms_fig.axes.set_xscale('log')
        self.dh_fig.axes.set_xscale('log')

        majorFormatter = FuncFormatter(self.formatter)

        for sfig in [self.rms_fig, self.srms_fig, self.dh_fig]:
            sfig.axes.set_xscale('log')


        self.srms_fig.axes.xaxis.set_ticklabels([])
        self.rms_fig.axes.xaxis.set_ticklabels([])

        def srms_format_func(v, i):
            if i is None:
                return ""
            elif i % 2 == 0:
                return ""
            else:
                return r"$%.2f$" % v

        srms_formatter = FuncFormatter(srms_format_func)
        self.srms_fig.axes.yaxis.set_major_formatter(srms_formatter)

        def dh_format_func(v, i):
            if i is None:
                return ""

            if round(v, 2) == v:
                return r"$%.2f$" % v
            else:
                return ""

        dh_formatter = FuncFormatter(dh_format_func)
        self.dh_fig.axes.yaxis.set_major_formatter(dh_formatter)

        # majorFormatter = FormatStrFormatter('%.1f')
        # self.rms_fig.axes.yaxis.set_major_formatter(majorFormatter)
        # self.srms_fig.axes.yaxis.set_major_formatter(majorFormatter)
        # self.dh_fig.axes.yaxis.set_major_formatter(majorFormatter)
        #
        # majorLocator = MultipleLocator(1)
        # self.rms_fig.axes.yaxis.set_major_locator(majorLocator)
        # self.srms_fig.axes.yaxis.set_major_locator(majorLocator)
        # self.dh_fig.axes.yaxis.set_major_locator(majorLocator)
        #
        xmax = 4E6
        self.rms_fig.set_xlim(0, xmax)
        self.srms_fig.set_xlim(0, xmax)
        self.dh_fig.set_xlim(0, xmax)

        self.rms_fig.set_ylabel(r"$\log %s$" % rmslabel, labelpad=30)
        self.srms_fig.set_ylabel(r"$\log %s$" % sslabel)
        self.dh_fig.set_ylabel(r"$\log %s$" % whlabel, labelpad=20)



        rms_span = self.rms_fig.get_ylim()[1] + self.rms_fig.get_ylim()[0]
        srms_span = self.srms_fig.get_ylim()[1] + self.srms_fig.get_ylim()[0]
        dh_span = self.dh_fig.get_ylim()[1] + self.dh_fig.get_ylim()[0]

        xtext = 4E5

        self.rms_fig.text(xtext, rms_span/2, r"$%s = %.3f$" % (rmslabel, rms_value_chosen), fontsize=self.fontSize)
        self.srms_fig.text(xtext, srms_span/2, r"$%s = %.3f$" % (sslabel, srms_value_chosen), fontsize=self.fontSize)
        self.dh_fig.text(xtext, dh_span/2, r"$%s = %.3f$" % (whlabel, dh_value_chosen), fontsize=self.fontSize)

        labels = ["a", "b", "c"]
        sfigs_first = [self.rms_fig, self.srms_fig, self.dh_fig]
        for i, sfig in enumerate(sfigs_first):
            sfig.text(0.05, 0.8, r"$\mathrm{(%s)}$" % labels[i], verticalalignment="center", horizontalalignment="left", transform=sfig.axes.transAxes, fontsize=30)



        labels2 = labels + ["d"]

        self.subfigure1.set_ylim(0, 17)
        pad = 13
        labelpads = [None, pad, pad, pad]
        for i, sfig in enumerate(sfigs):
            sfig.set_xlim([-1.1, 2.5])
            sfig.axes.get_xaxis().set_ticks([-1, 0, 1, 2])
            sfig.axes.yaxis.set_major_formatter(majorFormatter)
            sfig.set_ylabel(r"$%s$" % rmslabel, labelpad=labelpads[i])
            sfig.set_ybound(0)
            sfig.text(1-0.1, 0.5, r"$\mathrm{(%s)}$" % labels2[i], horizontalalignment="left", verticalalignment="center", fontsize=30, transform=sfig.axes.transAxes)

            if i != len(sfigs) - 1:
                sfig.axes.xaxis.set_ticklabels([])

            ax = sfig.axes.twinx()
            ax.set_ylabel(r"$\alpha=%g$" % self.alpha_values[i], labelpad=15)
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
        lg = self.subfigure1.legend(loc="upper left", numpoints=1, handlelength=1.2, ncol=3, columnspacing=0.3, handletextpad=0.5, borderaxespad=0.3,  bbox_to_anchor=(-0.01, 1.6))
        lg.get_frame().set_fill(not (self.toFile and self.transparent))




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


class GrowthSpeed2(DCVizPlotter):

    nametag = "growthspeed_depr\_(.*)\.npy"

    hugifyFonts = True

    figMap = {"omega_vs_v" : "subfigure",
              "k_figure": "subfigure2",
              "asd" : "subfigure3",
              "slope_slopes": "subfigure4",
              "slope_ints": "subfigure5",
              "int_slopes": "subfigure6",
              "int_ints": "subfigure7",
              "asdASd": "subfigure8",
              "asssasdASd": "subfigure9"}

    isFamilyMember = True

    def plot(self, data):

        E0_values = self.get_family_member_data(data, "E0")
        alpha_values = self.get_family_member_data(data, "alpha")
        mu_shift_values = self.get_family_member_data(data, "mu_shift")
        r0_values = self.get_family_member_data(data, "r0")
        s0_values = self.get_family_member_data(data, "s0")
        mu_values = self.get_family_member_data(data, "mu")
        v_values = self.get_family_member_data(data, "v")

        # r0_values = r0_values[1:]             #strip zero at [0]
        # s0_values = s0_values[1:]             #strip zero at [0]
        # mu_values = mu_values[1:, :, :, :, :] #strip for E0=0 where mu = mu_shifts since mu_eq = 0
        # alpha_values = alpha_values[:]       #Skip < 1 because of messy behavior

        omega = exp(mu_shift_values) - 1
        # log_omega = np.sign(omega)*np.log(abs(omega))
        sat_idx = np.where(omega > 0)

        f = lambda x, a, b: a*(exp(-b*x) - 1)

        slopes = []
        ints = []
        for j, alpha in enumerate(alpha_values):

            remainder = self.find_remainder(E0=0,
                                            alpha=alpha,
                                            l_d=None,
                                            v=v_values[0, j, :, 0, 0],
                                            mu=mu_shift_values,
                                            d_mu=mu_shift_values,
                                            omega=omega,
                                            plot=j==2)
            #
            # params, _ = curve_fit(f, omega, remainder, (1, 1))
            #
            # print params
            # self.subfigure.plot(omega, f(omega, *params), "s", fillstyle="none")


            slope, intersection, _, _, err = linregress(np.log(omega+1), np.log(remainder))
            print intersection, err

            slopes.append(slope)
            ints.append(intersection)

        self.subfigure2.plot(alpha_values, slopes, 'r--o')

        self.subfigure.set_xlabel("rofldofl")
        self.subfigure.set_xlim([-1, 4])

        slope_slope0, slope_int0, _, _, err = linregress(alpha_values, slopes)

        # print slope_slope, slope_intersection, err

        self.subfigure3.plot(alpha_values, ints, "b--o")

        int_slope0, int_int0, _, _, err = linregress(alpha_values, ints)
        #
        # print int_slope0, int_int0, err
        #
        # print "log(rem) = s(a)*log(1 + omega) + b(a)"
        # print "s(a) = %g*a + %g" % (slope_slope0, slope_int0)
        # print "b(a) = %g*a + %g" % (int_slope0, int_int0)
        # self.Exit()

        int_ints = np.zeros(shape=(len(E0_values)-1, len(r0_values)-1))
        int_slopes = np.zeros_like(int_ints)
        slope_ints = np.zeros_like(int_ints)
        slope_slopes = np.zeros_like(int_ints)

        errmax = 0
        for i, E0 in enumerate(E0_values[1:]):

            slopes = []
            ints = []

            for j, alpha in enumerate(alpha_values):

                pressure_slopes = np.zeros(shape=(len(r0_values)-1, len(s0_values)-1))
                pressure_ints = np.zeros_like(pressure_slopes)

                for l, l_d in enumerate(r0_values[1:]):
                    for m, s0 in enumerate(s0_values[1:]):

                        doplot = j == 2 and l == 3 and m == 3

                        remainder = self.find_remainder(E0,
                                                        alpha,
                                                        l_d,
                                                        v_values[i+1, j, :, l+1, m+1],
                                                        mu_values[i+1, j, :, l+1, m+1],
                                                        mu_shift_values,
                                                        omega,
                                                        doplot)

                        #claim that log r = intersection + log(om + 1)*slope
                        slope, intersection, _, _, err = linregress(np.log(omega+1), np.log(remainder))

                        if err > errmax:
                            errmax = err
                            # print E0, alpha, l_d, s0, errmax

                        pressure_slopes[l, m] = slope
                        pressure_ints[l, m] = intersection

                #claim that slope and intersection is independent of s_0
                comb_slope = pressure_slopes.mean(axis=1)
                comb_ints = pressure_ints.mean(axis=1)

                slopes.append(comb_slope)
                ints.append(comb_ints)

            L = len(r0_values) - 2
            for l in range(len(r0_values)-1):
                l_slopes = [slope[l] for slope in slopes]
                l_ints = [_int[l] for _int in ints]

                if l == L:
                    self.subfigure2.plot(alpha_values, l_slopes, self.markers.next(), label="E_0=%g" % E0)
                    self.subfigure3.plot(alpha_values, l_ints, self.markers.next(), label="E_0=%g" % E0)

                # claim that slope is linear with alpha: slope = slope_slope*alpha + slope_int
                slope_slope, slope_int, _, _, err = linregress(alpha_values, l_slopes)
                # print slope_slope, slope_int, err

                #claim the same for intersection = int_slope*alpha + int_int
                int_slope, int_int, _, _, err = linregress(alpha_values, l_ints)

                slope_slopes[i, l] = slope_slope
                slope_ints[i, l] = slope_int
                int_slopes[i, l] = int_slope
                int_ints[i, l] = int_int
                # print int_slope, int_int, err
                #
                #
                # print "log(rem) = s(a)*log(1 + omega) + b(a)"
                # print "s(a) = %g*a + %g" % (slope_slope, slope_intersection)
                # print "b(a) = %g*a + %g" % (int_slope, int_int)
                #
                # raw_input()



                # print pressure_slopes, pressure_slopes.flatten().std()/pressure_slopes.flatten().mean(), comb_slope
                # print "----"
                # print pressure_ints, pressure_ints.flatten().std()/pressure_ints.flatten().mean(), comb_ints
                #
                # print errmax
                #
                # raw_input()

        #claim that int_ints and int_slopes are small and can be neglected (some dependency on lambda is seen)
        #we get: log r = log(om + 1)*(slope_slope*a + slope_int)
        # print int_int0, int_ints
        # self.Exit()
        mean_int_ints = (int_int0 + int_ints.flatten().sum())/((len(E0_values)-1)*(len(r0_values)-1) + 1)
        mean_int_slopes = (int_slope0 + int_slopes.flatten().sum())/((len(E0_values)-1)*(len(r0_values)-1) + 1)

        print "mean int ints", mean_int_ints
        print "mean int slopes", mean_int_slopes

        E0_int_all = 0
        E0_slopes = []
        E0_ints = []
        E0_poly = []

        # E0_orig = E0_values.copy()
        for l in range(len(r0_values) - 1):
            # l_d = r0_values[l+1]
            # fac = l_d*(1 - exp(-1./l_d))
            # E0_values = E0_orig/fac
            # print "l=", l
            # print "slope slopes", slope_slope0, slope_slopes[:, l]
            # raw_input()
            # print "slope ints", slope_int0,  slope_ints[:, l]
            # raw_input()
            # print "int slopes", int_slope0, int_slopes[:, l]
            # raw_input()
            # print "int ints", int_int0, int_ints[:, l]
            # raw_input()
            # print "done"

            l_slope_ints = [slope_int0] + list(slope_ints[:, l])

            self.subfigure5.plot(E0_values, l_slope_ints, "--x")

            #claim that slope_int = E0_slope*E_0 + E0_int
            E0_slope, E0_int, _, _, err = linregress(E0_values, l_slope_ints)

            #claim that E0_int is independent of l_d (= -1 / 0)
            E0_int_all += E0_int

            #E0_slope is found to be exponentially decaying with lambda from 0.24 to 0.08 on l_d 1 to 5
            E0_slopes.append(E0_slope)


            l_slope_slopes = [slope_slope0] + list(slope_slopes[:, l])
            self.subfigure4.plot(E0_values, l_slope_slopes, "--x")
            f_fit = lambda x, a, b, c: a + b*x**0.5
            (a, b, c), _ = curve_fit(f_fit, E0_values, l_slope_slopes, (0.2, 1, 0.5))
            print r0_values[l+1], a, b, c
            #claim that slope_slope = E0_ints + E0_poly[0]*E0^E0_poly[1]
            self.subfigure4.plot(E0_values, f_fit(E0_values, a, b, c), '^')


            # E0_log_slope, E0_log_int, _, _, err = linregress(np.log(E0_values), l_slope_slopes)
            # E0_values[0] -= 0.001

            E0_poly.append([b, c])
            E0_ints.append(a)


        #find <E0_int> ~ -1
        print "e0 int all:", E0_int_all/(len(r0_values)-1)

        self.subfigure6.plot(r0_values[1:], E0_slopes)
        self.subfigure6.set_xlabel("om**(ax + b*E0 + c, here show b(ld)")


        self.subfigure7.plot(r0_values[1:], E0_ints, "--o")
        self.subfigure7.set_xlabel("om**(ax + b*E0 + c(ld), here show c(ld)")

        self.subfigure8.plot(r0_values[1:], [x[0] for x in E0_poly], ":x")
        self.subfigure9.plot(r0_values[1:], [x[1] for x in E0_poly], ":o")
        self.subfigure9.set_xlabel("om**(alp(kE0**n + b) + c, here show k(ld)")
        self.subfigure8.set_xlabel("om**(alp(kE0**n + b) + c, here show n(ld)")

        self.subfigure4.set_ylabel("lslopeslopes")
        self.subfigure5.set_ylabel("lslopeints")

        #Putting it all back together we get
        #r = <R>/c_eq*(om + 1) = (om + 1)^(a*[g(0.27, 0.23) + g(0.017, 0.009)*log(E0)] + E0*g(0.24, 0.08) - 1)
        #defining G_a(F0, l_d) = [g(0.27, 0.23) + g(0.017, 0.009)*log(E0)]
        #and G_s(F0, l_d) = E0*g(0.24, 0.08)
        #and G = (a*G_a + G_s)
        #we get that <R>_om = <R>_0 (1 + om)^G
        #where in general G != 1 such that
        #we cannot have the relation v/(omega+1) = k
        #very clean


    def find_remainder(self, E0, alpha, l_d, v, mu, d_mu, omega, plot=True):

        if l_d is not None:
            F0 = E0/(l_d*(1 - exp(-1./l_d)))
        else:
            F0 = 0

        mu_eq = (mu-d_mu).mean()

        # if F0 != 0:
        #     t = mu_eq/(alpha*F0)
        #     if abs(t - 1) > 0.1:
        #         print "Thermo fail:", t, "l_d", l_d, "E0", E0

        c_over_c0 = np.exp(mu)
        c0 = exp(-2*alpha)

        #transform to "proper" time by multiplyin with c
        v_transformed = v*c_over_c0*c0

        c_eq = exp(mu_eq)*c0
        remainder = -(v_transformed - c_eq*(1 + omega))/c_eq
        if plot:
            self.subfigure.loglog(omega+1, remainder, "--" + self.markers.next(), label=r"$E_0=%.2f$" % E0)

        return remainder

class GrowthSpeed(DCVizPlotter):

    nametag = "^growthspeed_dated_\_(.*)\.npy"

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

    def plot_and_slopify(self, E0, alpha,r0, omega, mu_shifts, mu, v):
        mu0 = (mu - mu_shifts).mean()
        c_over_c0 = exp(mu0 + mu_shifts)
        c0 = exp(-2*alpha)
        c = c_over_c0*c0

        v *= c

        if E0 == 0:
            if not (omega == (c_over_c0-1)).all():

                print "OMEGA FAIL"
                print omega
                print c_over_c0
                self.Exit()

        idx_high = np.where(omega >= 0)
        idx_vhigh = np.where(omega >= 3)
        idx_low = np.where(omega < 0)
        idx_all = np.where(omega > -2)

        try:
            if int(self.argv[3]) < 0:
                idx_chosen = idx_low
            else:
                idx_chosen = idx_high
        except IndexError:
            idx_chosen = idx_high


        F0 = E0/(r0*(1-np.exp(-1./r0)))

        slope, intercept, _, _, err = linregress(omega[idx_chosen], v[idx_chosen])

        v_over_omega = v[idx_chosen]/omega[idx_chosen]
        v_over_omega = v_over_omega[np.where(omega[idx_chosen] != 0)]
        avg = v_over_omega.mean()

        slope = avg

        if str(avg) == "nan":
            print avg
            self.Exit()

        # slope_low, intercept, _, _, _ = linregress(omega[idx_low], v[idx_low])
        #
        # print E0, slope/slope_low - 1

        # slope -= np.exp(-2*alpha + alpha*F0)
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

        # E0_values = E0_values[:6]

        omega = exp(mu_shift_values) - 1

        k_values = np.zeros(shape=(len(E0_values), len(alpha_values), len(r0_values), len(s0_values)))
        mu_zero = np.zeros_like(k_values)
        errors = np.zeros_like(k_values)

        J, L, M = [int(x) for x in self.argv[:3]]

        print "alpha=%g, s0 = %g, r0 = %g" % (alpha_values[J], s0_values[M+1], r0_values[L+1])

        unloaded_with_alpha = []
        for j in range(len(alpha_values)):
            _, slope0, v0, err = self.plot_and_slopify(0, alpha_values[j], 1, omega, mu_shift_values, mu_shift_values, v_values[0, j, :, 0, 0])

            k_values[0, j, :, :] = slope0
            mu_zero[0, j, :, :] = 0
            errors[0, j, :, :] = err

            unloaded_with_alpha.append(slope0)

            if j == J:
                self.uberplot(exp(mu_shift_values) - 1, v0, 0, 0)
                # self.subfigure.loglog(abs(omega), abs(v0), "k--" + shapes[0], label="E0=0",
                #                     linewidth=1,
                #                     fillstyle='none',
                #                     markersize=7,
                #                     markeredgewidth=1.5,
                #                     color="black")

        unloaded_alpha_slope, unloaded_alpha_cut, _, _, err = linregress(alpha_values, np.log(np.array(unloaded_with_alpha)))
        print err, unloaded_alpha_slope, unloaded_alpha_cut
        # raw_input()

        k = 1
        l_skip = []
        for i, E0 in enumerate(E0_values[1:]):
            for j, alpha in enumerate(alpha_values):
                # for k, mu_shift in enumerate(mu_shift_values):
                    for l, r0 in enumerate(r0_values[1:]):
                        for m, s0 in enumerate(s0_values[1:]):

                            mu0, slope, v, error = self.plot_and_slopify(E0,
                                                                         alpha,
                                                                         r0,
                                                                         omega,
                                                                         mu_shift_values,
                                                                         mu_values[i+1, j, :, l+1, m+1],
                                                                         v_values[i+1, j, :, l+1, m+1])

                            # if slope < 0:
                            #     print "ERROR slope=%g, E0(%d)=%g, alpha(%d)=%g, r0(%d)=%g, s0(%d)=%g" % (slope, i, E0, j, alpha, l, r0, m, s0)
                            #     sys.exit(1)

                            if slope < 0 and abs(slope) < 1E-3:
                                l_skip.append(l)

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

        #K_values is zero for all l=m=0 except for for E0 = 0
        # print np.where(k_values[:, :, 1:, 1:] == 0)
        # KSA = 0
        # for i in range(4):
        #     KSA += len(np.where(k_values == 0)[i])
        # print KSA/4., (len(E0_values)-1)*len(alpha_values)*(len(r0_values) + len(s0_values) - 1)
        # sys.exit(1)
        #
        # for i, E0 in enumerate(E0_values):
        #     for j, alpha in enumerate(alpha_values):
        #         for l, r0 in enumerate(r0_values):
        #             for m, s0 in enumerate(s0_values):
        #                 if E0 == 0:
        #                     F0 = 0
        #                     if l != 0 or m != 0:
        #                         continue
        #                 else:
        #                     if l == 0 or m == 0:
        #                         continue
        #                     F0 = E0/(r0*(1 - exp(-1/r0)))
        #                     continue
        #
        #
        #                 print alpha, k_values[i, j, l, m]-exp(-2*alpha + F0*alpha)
        #
        # self.Exit()

        log_k_E0_slopes = np.zeros(shape=(len(alpha_values), len(r0_values), len(s0_values)))
        log_k_cutz = np.zeros_like(log_k_E0_slopes)
        log_k_slope_error = np.zeros_like(log_k_E0_slopes)

        ka = 0
        na = 4
        nm = na - 2
        d = len(alpha_values)/nm
        for j, alpha in enumerate(alpha_values):
            for l, r0 in enumerate(r0_values[1:]):

                if l in l_skip:
                    continue

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

                    if str(intercept) == "nan":
                        print intercept, alpha, r0, s0
                        print k_values[:, j, l+1, m+1]
                        self.Exit()

                    log_k_E0_slopes[j, l+1, m+1] = log_k_slope
                    log_k_cutz[j, l+1, m+1] = intercept
                    log_k_slope_error[j, l+1, m+1] = err

        alpha_slopes = np.zeros(shape=(len(r0_values), len(s0_values)))
        alpha_cutz = np.zeros_like(alpha_slopes)

        for l, r0 in enumerate(r0_values[1:]):
            if l in l_skip:
                continue

            for m, s0 in enumerate(s0_values[1:]):
                alpha_slope, intercept, _, _, err = linregress(alpha_values, log_k_E0_slopes[:, l+1, m+1])

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

            if l in l_skip:
                continue

            S = 0
            count = 0

            #
            # if r0 < self.l_min:
            #     continue

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

        power, log_constant, _, _, err = linregress(np.log(alpha_values), np.log(-S_tot))
        print "Power=%g, constant=%g, err=%g" % (power, exp(log_constant), err)
        print l_skip
        if len(l_skip) == 0:
            l_incr = 0
        else:
            l_incr = max([l+1 for l in l_skip])
        # self.Exit()

        self.subfigure5.loglog(alpha_values, exp(log_constant)*alpha_values**power, "g-^",
                               label=r"$%.2f\alpha^{%.2f}$" % (exp(log_constant), power))
        self.subfigure5.loglog(alpha_values, -S_tot, "r-^", label="Avg all param")

        self.subfigure5.legend(numpoints=1, handlelength=1.2, borderaxespad=0.3, )

        for l, r0 in enumerate(r0_values[1+l_incr:]):
            self.subfigure6.plot(s0_values[1:], alpha_cutz[l+1, 1:], label="r0 = %g" % r0)

        self.subfigure6.set_xlabel(r"$\sigma_0$")
        self.subfigure6.set_ylabel(r"$\mathrm{alpha shift}$")
        # self.subfigure6.legend()

        cuts = alpha_cutz[1:, 1:].mean(axis=1)

        print r0_values[1:], cuts
        self.subfigure7.plot(r0_values[1+l_incr:], cuts[l_incr:], 'b-x')
        self.subfigure7.set_xlabel(r"$\lambda_D$")
        self.subfigure7.set_ylabel(r"$\mathrm{avg alpha shift}$")

        for l, r0 in enumerate(r0_values[1+l_incr:]):
            self.subfigure8.plot(s0_values[1:], alpha_slopes[l+1, 1:], label="r0=%g" % r0)
        self.subfigure8.legend()
        self.subfigure8.set_xlabel(r"$\sigma_0$")
        self.subfigure8.set_ylabel("alpha E0 prefac")

        comb = alpha_slopes[1:, 1:].mean(axis=1)
        I = np.where(r0_values[1:] > 2.0)
        comb = comb[I]
        x = r0_values[1:][I]

        self.subfigure9.plot(x, comb, 'b-x')
        from scipy.optimize import curve_fit
        f = lambda x, a, b: a*x + b
        (a, b), _ = curve_fit(f, x, comb, (1, 1))
        self.subfigure9.plot(x, f(x, a, b), "r-x")
        print a, b

        self.subfigure9.set_xlabel(r"$\lambda_D$")
        self.subfigure9.set_ylabel(r"avg alpha E0 slope")

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
        self.subfigure4.set_xlim(0, alpha_values.max()*1.075)


        self.subfigure5.set_xlabel(r"$\alpha = E_b/kT$")
        self.subfigure5.set_ylabel(r"$\mathrm{log k shift}$")

        self.subfigure11.set_xlabel(r"$\log (\Omega + 1)$")
        self.subfigure11.set_ylabel(r"$\log(|\dot{H}|)$")
        self.subfigure11.legend(loc="upper left")



        # self.subfigure3.set_xlabel(r"$E_0$")
        # self.subfigure3.set_ylabel(r"$\gamma_\mathrm{eq}$")


class GrowthSpeed3(DCVizPlotter):

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
              "alpha_p0": "subfigure10"}

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
    shapes = ["s", "^", "v", 'o']

    transparent = True

    def plot_and_slopify(self, E0, alpha,r0, omega, mu_shifts, mu, v):
        mu0 = (mu - mu_shifts).mean()
        c_over_c0 = exp(mu0 + mu_shifts)
        c0 = exp(-2*alpha)
        c = c_over_c0*c0

        v *= c_over_c0

        if E0 == 0:
            if not (omega == (c_over_c0-1)).all():

                print "OMEGA FAIL"
                print omega
                print c_over_c0
                self.Exit()

        idx_high = np.where(omega >= 0)
        idx_vhigh = np.where(omega >= 3)
        idx_low = np.where(omega < 0)
        idx_all = np.where(omega > -2)

        try:
            if int(self.argv[3]) < 0:
                idx_chosen = idx_low
            else:
                idx_chosen = idx_high
        except IndexError:
            idx_chosen = idx_high


        F0 = E0/(r0*(1-np.exp(-1./r0)))

        slope1, intercept, _, _, err = linregress(omega[idx_chosen], v[idx_chosen])

        v_over_omega = v[idx_chosen]/omega[idx_chosen]
        v_over_omega = v_over_omega[np.where(omega[idx_chosen] != 0)]
        avg = v_over_omega.mean()

        slope = avg


        if str(avg) == "nan":
            print avg
            self.Exit()

        # slope_low, intercept, _, _, _ = linregress(omega[idx_low], v[idx_low])
        #
        # print E0, slope/slope_low - 1

        # slope -= np.exp(-2*alpha + alpha*F0)
        return mu0, slope1, v, err

    def uberplot(self, omega, v, k, E0, kg_on_plot=False):

        idx_undersat = np.where(omega <= 0)
        idx_oversat = np.where(omega >= 0)

        omega_undersat_list = list(omega[idx_undersat])
        omega_oversat_list = list(omega[idx_oversat])

        v_undersat_list = list(v[idx_undersat])
        v_oversat_list = list(v[idx_oversat])

        self.subfigure.plot(omega_undersat_list[::3] + omega_oversat_list,
                            v_undersat_list[::3] + v_oversat_list,
                            label=r"$%s=%.1f$" % (E0_tex, E0),
                            marker=self.shapes[k],
                            markeredgecolor="k",
                            linestyle="--",
                            color="r",
                            linewidth=1,
                            fillstyle='none',
                            markersize=7,
                            markeredgewidth=1.5)

        if kg_on_plot:

            pad_under = 0.4
            xpad_under = 0.1
            pad_over = 0.75
            first_over = np.where(omega[idx_oversat] >= 0.5)[0][0]
            last_over = np.where(omega[idx_oversat] >= 1.4)[0][0]
            first_under = 0
            last_under = np.where(omega[idx_undersat] >= -0.3)[0][0]

            c = "k-"

            self.subfigure.plot([omega_undersat_list[first_under]+xpad_under, omega_undersat_list[last_under]+xpad_under],
                                [v_undersat_list[first_under]-pad_under, v_undersat_list[last_under]-pad_under], c,
                                )
            self.subfigure.plot([omega_oversat_list[first_over], omega_oversat_list[last_over]],
                                [v_oversat_list[first_over]+pad_over, v_oversat_list[last_over]+pad_over], c,
                                )

            locs = [[0.44, 0.58],
                    [0.19, 0.042]]
            rots = [25, 0]
            kgs = [0, 0]

            kgs[0] = (v_oversat_list[first_over] - v_oversat_list[last_over])/(omega_oversat_list[first_over] - omega_oversat_list[last_over])
            kgs[1] = (v_undersat_list[first_under] - v_undersat_list[last_under])/(omega_undersat_list[first_under] - omega_undersat_list[last_under])

            for loc, rot, kg in zip(locs, rots, kgs):
                self.subfigure.text(loc[0], loc[1], r"$k_g \sim %.2f$" % kg, horizontalalignment="left", verticalalignment="center", transform=self.subfigure.axes.transAxes, rotation=rot)


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

        print r0_values, s0_values

        # E0_values = E0_values[:6]

        omega = exp(mu_shift_values) - 1

        k_values = np.zeros(shape=(len(E0_values), len(alpha_values), len(r0_values), len(s0_values)))
        mu_zero = np.zeros_like(k_values)
        errors = np.zeros_like(k_values)

        J, L, M = [int(x) for x in self.argv[:3]]

        print "alpha=%g, s0 = %g, r0 = %g" % (alpha_values[J], s0_values[M+1], r0_values[L+1])

        unloaded_with_alpha = []
        for j in range(len(alpha_values)):
            _, slope0, v0, err = self.plot_and_slopify(0, alpha_values[j], 1, omega, mu_shift_values, mu_shift_values, v_values[0, j, :, 0, 0])

            k_values[0, j, :, :] = slope0
            mu_zero[0, j, :, :] = 0
            errors[0, j, :, :] = err

            unloaded_with_alpha.append(slope0)

            if j == J:
                self.uberplot(exp(mu_shift_values) - 1, v0, 0, 0)

        ula = -np.log(np.array(unloaded_with_alpha))
        self.subfigure5.plot(alpha_values, ula, label="ZERO")

        unloaded_alpha_slope, unloaded_alpha_cut, _, _, err = linregress(alpha_values, np.log(np.array(unloaded_with_alpha)))
        print err, unloaded_alpha_slope, unloaded_alpha_cut

        k = 1
        l_skip = []
        for i, E0 in enumerate(E0_values[1:]):
            for j, alpha in enumerate(alpha_values):
                    for l, r0 in enumerate(r0_values[1:]):
                        for m, s0 in enumerate(s0_values[1:]):

                            mu0, slope, v, error = self.plot_and_slopify(E0,
                                                                         alpha,
                                                                         r0,
                                                                         omega,
                                                                         mu_shift_values,
                                                                         mu_values[i+1, j, :, l+1, m+1],
                                                                         v_values[i+1, j, :, l+1, m+1])

                            if slope < 0 and abs(slope) < 1E-3:
                                l_skip.append(l)

                            k_values[i+1][j][l+1][m+1] = slope
                            errors[i+1][j][l+1][m+1] = error
                            mu_zero[i+1][j][l+1][m+1] = mu0

                            if E0 in self.plot_values and j == J and l == L and m == M:
                                self.uberplot(omega, v, k, E0, E0 == 1)

                                k += 1

        log_k_E0_slopes = np.zeros(shape=(len(alpha_values), len(r0_values), len(s0_values)))
        log_k_cutz = np.zeros_like(log_k_E0_slopes)
        log_k_slope_error = np.zeros_like(log_k_E0_slopes)

        ka = 0
        na = 5
        nm = na - 2
        d = len(alpha_values)/nm
        for j, alpha in enumerate(alpha_values):
            for l, r0 in enumerate(r0_values[1:]):

                if l in l_skip:
                    continue

                # F0_values = E0_values/(r0*(1 - np.exp(-1./r0)))

                for m, s0 in enumerate(s0_values[1:]):

                    log_k_slope, intercept, _, _, err = linregress(E0_values, np.log(k_values[:, j, l+1, m+1]))

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

                    if str(intercept) == "nan":
                        print intercept, alpha, r0, s0
                        print k_values[:, j, l+1, m+1]
                        self.Exit()

                    log_k_E0_slopes[j, l+1, m+1] = log_k_slope
                    log_k_cutz[j, l+1, m+1] = intercept
                    log_k_slope_error[j, l+1, m+1] = err


        alpha_slopes = np.zeros(shape=(len(r0_values), len(s0_values)))
        alpha_cutz = np.zeros_like(alpha_slopes)

        for l, r0 in enumerate(r0_values[1:]):
            if l in l_skip:
                continue

            for m, s0 in enumerate(s0_values[1:]):
                alpha_slope, intercept, _, _, err = linregress(alpha_values, log_k_E0_slopes[:, l+1, m+1])

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

        a_cut = np.where(alpha_values > 1.0)

        S_tot = 0
        count2 = 0
        for l, r0 in enumerate(r0_values[1:]):

            if l in l_skip:
                continue

            S = 0
            count = 0

            for m, s0 in enumerate(s0_values[1:]):
                S += log_k_cutz[:, l+1, m+1]
                count += 1


                if l == L and m == M:

                    self.subfigure5.plot(alpha_values[a_cut], -log_k_cutz[:, l+1, m+1][a_cut],
                                         linestyle="--",
                                         marker="s",
                                         linewidth=1,
                                         fillstyle='none',
                                         markersize=7,
                                         markeredgewidth=1.5)


            S_tot += S
            count2 += count

        S_tot /= count2

        power, log_constant, _, _, err = linregress(np.log(alpha_values), np.log(-S_tot))
        print "Power=%g, constant=%g, err=%g" % (power, exp(log_constant), err)
        print l_skip
        if len(l_skip) == 0:
            l_incr = 0
        else:
            l_incr = max([l+1 for l in l_skip])

        self.subfigure5.plot(alpha_values, exp(log_constant)*alpha_values**power, "g-^",
                               label=r"$%.2f\alpha^{%.2f}$" % (exp(log_constant), power))

        FUNC = lambda x, a, b, c: x*(exp(b*x) - 1)
        (a, b, c), _ = curve_fit(FUNC, alpha_values, ula, (1, 1, 1))
        self.subfigure5.plot(alpha_values, FUNC(alpha_values, a, b, c), "k--x", label="cf")
        print "--a-a--a", a, b, c

        slope, intercept, _, _, _ = linregress(alpha_values[a_cut], S_tot[a_cut])
        print intercept, slope

        self.subfigure5.plot(alpha_values, -S_tot, "r-^", label="Avg all param")

        lg = self.subfigure5.legend(numpoints=1, handlelength=1.2, borderaxespad=0.3, )
        lg.get_frame().set_fill(not (self.toFile and self.transparent))

        for l, r0 in enumerate(r0_values[1+l_incr:]):
            self.subfigure6.plot(s0_values[1:], alpha_cutz[l+1, 1:], label="r0 = %g" % r0)

        self.subfigure6.set_xlabel(r"$\sigma_0$")
        self.subfigure6.set_ylabel(r"$\mathrm{alpha shift}$")
        # self.subfigure6.legend()

        cuts = alpha_cutz[1:, 1:].mean(axis=1)

        print r0_values[1:], cuts
        self.subfigure7.plot(r0_values[1+l_incr:], cuts[l_incr:], 'b-x')
        self.subfigure7.set_xlabel(r"$\lambda_D$")
        self.subfigure7.set_ylabel(r"$\mathrm{avg alpha shift}$")

        for l, r0 in enumerate(r0_values[1+l_incr:]):
            self.subfigure8.plot(s0_values[1:], alpha_slopes[l+1, 1:], label="r0=%g" % r0)
        lg = self.subfigure8.legend()
        lg.get_frame().set_fill(not (self.toFile and self.transparent))
        self.subfigure8.set_xlabel(r"$\sigma_0$")
        self.subfigure8.set_ylabel("alpha E0 prefac")

        comb = alpha_slopes[1:, 1:].mean(axis=1)
        I = np.where(r0_values[1:] > 3.0)
        comb = comb[I]
        x = r0_values[1:][I]

        self.subfigure9.plot(x, comb, 'b-x')
        f = lambda x, a, b: a*x + b
        (a, b), _ = curve_fit(f, x, comb, (1, 1))
        self.subfigure9.plot(x, f(x, a, b), "r-x")
        print a, b

        self.subfigure9.set_xlabel(r"$\lambda_D$")
        self.subfigure9.set_ylabel(r"avg alpha E0 slope")

        self.subfigure.set_xlabel(r"$\Omega = c/c_\mathrm{eq} - 1$")
        self.subfigure.set_ylabel(r"$\dot{H} = \Delta \langle h\rangle / \Delta t$")
        lg = self.subfigure.legend(loc="upper left", numpoints=1, handlelength=1.2, borderaxespad=0.3 )
        lg.get_frame().set_fill(not (self.toFile and self.transparent))
        self.subfigure.set_xlim(-1, 3)
        self.subfigure.axes.xaxis.set_ticks([-1, 0, 1, 2, 3])
        self.subfigure.axes.xaxis.set_major_formatter(FormatStrFormatter(r"$%g$"))
        self.subfigure.set_ylim(-3, 12)
        self.subfigure.axes.yaxis.set_ticks([-2, 0, 2, 4, 6, 8, 10])

        lg = self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1.2, ncol=2, columnspacing=0.3, handletextpad=0.5, borderaxespad=0.3)
        lg.get_frame().set_fill(not (self.toFile and self.transparent))
        self.subfigure2.axes.set_yscale('log')
        self.subfigure2.set_xlim(-0.05, 1.05)
        self.subfigure2.set_ylim(0.6, 23)

        self.subfigure2.set_xlabel(r"$%s \propto F_0/L$" % E0_tex)
        self.subfigure2.set_ylabel(r"$k_g = \dot{H} / \Omega$")

        self.subfigure4.set_xlabel(r"$\alpha = E_b/kT$")
        self.subfigure4.set_ylabel(r"$\log (k_g/k_0) / %s$" % E0_tex)
        self.subfigure4.set_ybound(0)
        self.subfigure4.set_xlim(0, alpha_values.max()*1.075)


        self.subfigure5.set_xlabel(r"$\alpha = E_b/kT$")
        self.subfigure5.set_ylabel(r"$\mathrm{log k shift}$")

        CV = 0
        NN = 0
        p0_mat = np.zeros(shape=(len(alpha_values)))
        log_k_averaged = np.zeros(shape=(len(E0_values), len(alpha_values)))
        for j, alpha in enumerate(alpha_values):
            N = 0
            opt_func = lambda x, p0: alpha*(x - p0)

            for l, r0 in enumerate(r0_values[1:]):

                if l in l_skip or r0 < 3:
                    continue

                for m, s0 in enumerate(s0_values[1:]):


                    pj, cv = curve_fit(opt_func, E0_values, np.log(k_values[:, j, l+1, m+1]), (0.5,))

                    # self.subfigure10.plot(E0_values, np.log(k_values[:, j, l+1, m+1]))
                    # self.subfigure10.plot(E0_values, opt_func(E0_values, pj))
                    CV += cv[0][0]

                    N += 1
                    NN += 1

                    phax = (E0_values - np.log(k_values[:, j, l+1, m+1])/alpha/alpha_slopes[l+1, m+1]).mean()

                    print phax, pj
                    p0_mat[j] += phax

                    log_k_averaged[:, j] += np.log(k_values[:, j, l+1, m+1])

            p0_mat[j] /= N
            log_k_averaged[:, j] /= N

        amin = 0

        for j, alpha in enumerate(alpha_values):
            if alpha < amin:
                continue

            self.subfigure10.plot(E0_values, log_k_averaged[:, j], label=r"$\alpha=%.2f$" % alpha)


        print "cov", CV/NN
        # self.subfigure10.plot(alpha_values, p0_mat)
        lg = self.subfigure10.legend()
        lg.get_frame().set_fill(not (self.toFile and self.transparent))
        self.subfigure10.set_xlabel("E0")
        self.subfigure10.set_ylabel("logk")

        log_k_start = log_k_averaged.min()
        log_k_end = log_k_averaged.max()

        global_match = 100000000000
        win = []
        for E0_ref in np.linspace(E0_values.min(), E0_values.max()):
            for log_k_ref in np.linspace(log_k_start, log_k_end):
                match = 0
                nm = 0

                x0 = E0_ref
                x1 = E0_values[0]
                x2 = E0_values[-1]

                y0 = log_k_ref

                for j, alpha in enumerate(alpha_values):

                    if alpha < amin:
                        continue

                    y1 = log_k_averaged[0, j]
                    y2 = log_k_averaged[-1, j]

                    # match += self.distance(x1, y1, x2, y2, x0, y0)


                    a, c, _, _, _ = linregress(E0_values, log_k_averaged[:, j])
                    b = -1

                    match += self.distance2(a, b, c, x0, y0)
                    nm+=1

                match /= nm

                if match < global_match:
                    global_match = match
                    win = [x0, y0]

        self.subfigure10.scatter(win[0], win[1], s=40)

        print win

        ##labels

        labelfigs = [self.subfigure, self.subfigure2, self.subfigure4]
        labels = ["a", "b", "c"]
        locs = [(0.4, 0.8),
                (0.075, 0.4),
                (0.075, 0.775)]

        for sfig, label, loc in zip(labelfigs, labels, locs):
            sfig.text(loc[0], loc[1], r"$\mathrm{(%s)}$" % label,
                      horizontalalignment="left",
                      verticalalignment="center",
                      transform=sfig.axes.transAxes,
                      fontsize=30)

    def distance(self, x1, y1, x2, y2, x0, y0):
        return abs((y2 - y1)*x0 - (x2 - x1)*y0 + x2*y1 - y2*x1)/np.sqrt((y2-y1)**2 + (x2-x1)**2)
    def distance2(self, a, b, c, x0, y0):
        return abs(a*x0 + b*y0 + c)/np.sqrt(a**2 + b**2)



class Quasi2D_slopes_and_stuff(DCVizPlotter):

    nametag = "linearplots\_(.*)\_?\d*\.npy"

    numpyBin = True

    isFamilyMember = True

    figMap = {"fig" : "gammaslopes", "fig2" : "E0slopes"}

    hugifyFonts = True

    transparent = True

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
        rotations = [3, 14.5, 22.75]
        plotslopes = []
        xloc = 0.65
        ylocs = [0.15, 0.475, 0.8]
        texts = []

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
                                      label=r"$%s=%1.2f$" % (E0_tex, E0s[n]),
                                      markersize=7,
                                      markeredgewidth=1.5,
                                      linewidth=1,
                                      color="black")

            kws = {}
            if n_plots == 0:
                kws["label"] = r"$\mathrm{Linear\,\,fit}$"

            self.gammaslopes.plot([0, alphas.max()], [0, slopes[n]*alphas.max()], "r--", **kws)

            texts.append(r"$%s=%1.2f$" % (E0_tex, E0s[n]))
            plotslopes.append(slopes[n])

            n_plots += 1




        self.gammaslopes.set_xbound(0)
        self.gammaslopes.set_ybound(0)
        # self.gammaslopes.legend(loc="upper left", numpoints=1, handlelength=0.8, ncol=2, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3)

        for text, yloc, rotation in zip(texts, ylocs, rotations):
            # rotation = np.arctan(slope)*180/np.pi
            print rotation, rotations, plotslopes
            self.gammaslopes.text(xloc, yloc, text,
                                  verticalalignment="center",
                                  horizontalalignment="left",
                                  transform=self.gammaslopes.axes.transAxes,
                                  rotation=rotation)

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

        self.E0slopes.set_xlabel(r"$%s \propto F_0/L$" % E0_tex)
        self.E0slopes.set_ylabel(r"$\gamma_\mathrm{eq}/\alpha$")

        self.gammaslopes.text(fig_4_x, fig_4_y, r"$\mathrm{(a)}$", verticalalignment="center", horizontalalignment="left", transform=self.gammaslopes.axes.transAxes, fontsize=fig_4_fs)
        self.E0slopes.text(fig_4_x, fig_4_y, r"$\mathrm{(b)}$", verticalalignment="center", horizontalalignment="left", transform=self.E0slopes.axes.transAxes, fontsize=fig_4_fs)

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
            print E0_value
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
                                 label="$%s = %.2f$" % (E0_tex, E0_value),
                                 markersize=7,
                                 markeredgewidth=1.5,
                                 linewidth=1)

            self.subfigure2.loglog(1./alpha_array, var_s_array, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$%s = %.2f$" % (E0_tex, E0_value),
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            self.subfigure3.loglog(1./alpha_array, (var_s_array*alpha_array)**p, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$%s = %.2f$" % (E0_tex, E0_value),
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            nplots += 1

        xmin = 0.48
        xmax = 5

        # self.subfigure.set_xlabel(r"$1/\alpha$")
        self.subfigure.set_xlim(xmin, xmax)
        # self.subfigure.axes.xaxis.set_visible(False)
        self.subfigure.axes.xaxis.set_ticklabels([])
        self.subfigure.set_ylabel(r"$\langle s_{\uparrow\downarrow} \rangle / L$")
        ax = self.subfigure.axes.twinx()
        ax.set_ylabel(r"$\langle E \rangle /E_bL$", labelpad=15)
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])
        lg = self.subfigure.legend(loc="upper left", numpoints=1, handlelength=0.5, ncol=3, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3,  bbox_to_anchor=(0, 1.55))
        lg.get_frame().set_fill(not (self.toFile and self.transparent))
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

        #Descriptions
        self.subfigure.text(0.05, 1-0.15, r"$\mathrm{surface\,\,size}$", verticalalignment="center", horizontalalignment="left", transform=ax.transAxes)
        self.subfigure2.text(0.05, 1-0.15, r"$\mathrm{surface\,\,fluctuations}$", verticalalignment="center", horizontalalignment="left", transform=ax2.transAxes)
        self.subfigure3.text(0.05, 0.15, r"$\mathrm{heat\,\,capacity}$", verticalalignment="center", horizontalalignment="left", transform=ax3.transAxes)

        x_loc = 0.725
        y_start = 0.1
        y_len = 0.75
        self.subfigure.annotate("",
                                xy=(x_loc, y_start),
                                xycoords="axes fraction",
                                xytext=(x_loc, y_start + y_len),
                                textcoords='axes fraction',
                                arrowprops=dict(facecolor='black', shrink=0.1, width=2, headwidth=10),
                                horizontalalignment='left', verticalalignment='center')
        self.subfigure.text(x_loc + 0.025, 1-(y_start + y_len - 0.2), r"$\mathrm{increasing}$", verticalalignment="center", horizontalalignment="left", transform=ax.transAxes)
        self.subfigure.text(x_loc + 0.025, 1-(y_start + y_len - 0.075), r"$\mathrm{confinement}$", verticalalignment="center", horizontalalignment="left", transform=ax.transAxes)

        self.subfigure.text(0.4, 0.725, r"$\mathrm{zero\,\,load\,\,limit}$", verticalalignment="center", horizontalalignment="left", transform=ax.transAxes, rotation=8.5)

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

    transparent = True

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

        self.mean_figure.plot(r0, self.scale(r0), "r-", label=r"$\lambda_D\left(1 - \exp\left(-1/\lambda_D\right)\right)$", linewidth=3)

        self.mean_figure.plot(r0, r0_mean, 'ks',
                              label=r"$\mathrm{KMC}$",
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5)

        lg = self.mean_figure.legend(loc="lower right", numpoints=1, handlelength=0.8, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3)
        lg.get_frame().set_fill(not (self.toFile and self.transparent))

        self.mean_figure.set_xlabel(r"$\lambda_D$")
        self.mean_figure.set_ylabel(r"$\alpha %s/\gamma_\mathrm{eq}$" % E0_tex)
        self.mean_figure.set_xbound(0)
        # self.mean_figure.set_ybound(-0.45)
        # self.mean_figure.axes.set_yticks([0, 1, 2, 3, 4])

        self.mean_figure.text(fig_4_x, fig_4_y, r"$\mathrm{(c)}$", verticalalignment="center", horizontalalignment="left", transform=self.mean_figure.axes.transAxes, fontsize=fig_4_fs)



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


class AutoCorrelation(DCVizPlotter):

    nametag = "autocorr\.arma"

    figMap = {"figure": ["subfigure"], "figure3D": [], "figure_projections" : "proj"}

    tight = False

    def plot(self, data):

        data -= data.min()


        Lm, Wm = data.shape

        L = (Lm + 1)/2
        W = (Wm + 1)/2

        ax = self.subfigure.imshow(data, extent=[-L, L, -W, W], origin="lower", aspect="auto", interpolation="none", cmap="Blues")
        self.subfigure.set_xlabel(r"$\delta x$")
        self.subfigure.set_ylabel(r"$\delta y$")
        #self.figure.colorbar(ax)
        #
        ax = Axes3D(self.figure3D)

        elements = Lm*Wm

        xpos, ypos = np.meshgrid(np.arange(Lm), np.arange(Wm))

        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(elements)
        dx = np.ones_like(zpos)
        dy = dx.copy()
        dz = data.flatten()

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, linewidth=0)

        d1 = np.diag(data)
        d2 = np.diag(np.flipud(data))

        for i in range(Lm):
            for j in range(Wm):
                print "%.2f" % data[i, j],
            print

        print
        print

        for d in d1:
            print "%.2f" % d,
        print

        for d in d2:
            print "%.2f" % d,
        print


        xl = np.linspace(-L, L, Lm)
        xw = np.linspace(-W, W, Wm)

        x1 = np.sqrt(2)*xl
        x2 = np.sqrt(2)*xw

        self.proj.plot(xl, data[L-1, :], label="X")
        self.proj.plot(xw, data[:, W-1], label="Y")
        self.proj.plot(x1, d1, label="11 diag")
        self.proj.plot(x2, d2, label="-11 diag")
        self.proj.legend()
        self.proj.set_xlabel(r"$\delta r$")
        self.proj.set_ylabel("corr")


class LatticediffSpeeds(DCVizPlotter):

    nametag = "confined_\w+_(\w+)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    ceq = exp(-3)
    h0 = 19

    figMap = {"f0": "subfigure", "f1" : "subfigure2"}

    def plot(self, data):

        supersaturations = self.get_family_member_data(data, "supersaturations")
        all_heights = self.get_family_member_data(data, "heights")
        all_times = self.get_family_member_data(data, "times")
        all_conc = self.get_family_member_data(data, "concentrations")
        lengths = self.get_family_member_data(data, "lengths")

        print all_conc.shape

        all_supersaturations = all_conc/self.ceq - 1

        print supersaturations

        M = 0
        for i, supersaturation in enumerate(supersaturations):

            s0 = all_conc[i, 0]/self.ceq - 1
            if supersaturation != s0:
                print "ERRAH", supersaturation, s0

            if len(self.argv) != 0:
                sfac = np.sign(supersaturation)
            else:
                sfac = 1

            point = sfac*s0*self.h0/(1/self.ceq - 1)

            l = lengths[i]

            T = all_times[i, :l]
            H = sfac*all_heights[i, :l]
            label = r"$\Omega=%.2f$" % s0
            m = T[-1]*1.05

            self.subfigure.scatter(0.9*m, point)
            self.subfigure.plot(T, H, label=label)
            self.subfigure.text(m, H[l/3:-l/3].mean(), label)

            self.subfigure2.plot(T, all_supersaturations[i, :l])

            if m > M:
                M = m

        self.subfigure.set_xlim(0, M*1.3)
        self.subfigure.set_xlabel("t")
        self.subfigure.set_ylabel("h")
        #lg = self.subfigure.legend(loc="upper left", numpoints=1, handlelength=0.5, ncol=3, columnspacing=0.5, handletextpad=0.5, borderaxespad=0.3)
        #lg.get_frame().set_fill(not (self.toFile and self.transparent))

    def kf(self, s, a):
        return -a*np.vectorize(self.kw)(s)

    def kw(self, s):

        if s >= 0:
            return 0.14
        else:
            return 0.53

class cconcdiffSpeeds(DCVizPlotter):

    nametag = "cconc_(\w+)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"figure" : "subfigure", "figure2" : "subfigure2"}

    ceq = exp(-3)
    h0 = 20

    def plot(self, data):

        supersaturations = self.get_family_member_data(data, "supersaturations")
        all_heights = self.get_family_member_data(data, "heights")
        all_times = self.get_family_member_data(data, "times")
        conc = self.get_family_member_data(data, "conc")
        lengths = self.get_family_member_data(data, "lengths")

        M = 0
        pad = 0.5
        for i, supersaturation in enumerate(supersaturations):

            if len(self.argv) != 0:
                sfac = np.sign(supersaturation)
            else:
                sfac = 1

            l = lengths[i]

            ss = conc[i, :l]/self.ceq - 1
            T = all_times[i, :l]/(ss + 1)
            H = sfac*all_heights[i, :l]
            label = r"$\Omega(0)=%.2f$" % supersaturation
            m = T[-1]*1.05

            point = supersaturation*self.h0/(1/self.ceq - 1)

            self.subfigure.scatter(0.9*m, point, s=20)
            self.subfigure.plot(T, H, label=label)
            self.subfigure.text(m, H[l/3:-l/3].mean(), label)

            self.subfigure2.plot(T, ss, label=label)

            if m > M:
                M = m

            self.subfigure2.text(-0.9*pad*T[-1], ss[0], label)

        self.subfigure.set_xlim(0, M*1.3)
        self.subfigure.set_xlabel("t")
        self.subfigure.set_ylabel("h")

        self.subfigure2.set_xlabel("t")
        self.subfigure2.set_ylabel(r"$\Omega$")
        self.subfigure2.set_xlim(-pad*T[-1], 1.1*T[-1])


class ignisSOS(DCVizPlotter):

    nametag = "ignisSOS\.ign"


    def plot(self, data):

        print self.loader.get_metadata()

        name_string, l, w = self.loader.get_metadata()

        names = [desc.split("@")[0] for desc in name_string.split()]
        targetname = self.argv[0]

        if targetname in names:
            self.scan_for_name(data, targetname, names)
        else:

            if targetname == "concentration":
                hi = names.index("AverageHeight")
                li = names.index("latticeDiffusion")

                c0 = exp(-3)
                L, W, H = 20, 20, 20

                conc = (data[li]/(L*W*(H - data[hi] - 2)))/c0
                cut = int(self.argv[1])

                self.subfigure.plot(conc)
                trail = np.cumsum(conc[cut:])/np.cumsum(np.ones_like(conc)[cut:])
                self.subfigure.plot(np.arange(len(conc))[cut:], trail, "r--")
                print trail[-1]
                self.subfigure.set_ylabel("concentration")

            else:

                raise RuntimeError("Invalid option %s." % targetname)

    def scan_for_name(self, data, targetname, names):

        for _data, name in zip(data, names):

            if name == targetname:
                self.subfigure.plot(_data)
                self.subfigure.set_ylabel(name)

                return

        raise RuntimeError("No match: %s in %s" % (targetname, str(names)))


class SizeTest(DCVizPlotter):

    nametag = "size_test_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    def plot(self, data):

        system_sizes = self.get_family_member_data(data, "system_sizes")
        alpha_values = self.get_family_member_data(data, "alpha_values")
        sizes = self.get_family_member_data(data, "sizes")

        for i, L in enumerate(system_sizes):

            I = np.where(sizes[i, :] != 0)
            self.subfigure.plot(alpha_values[I], sizes[i, :][I], "-"+self.markers.next(), label=r"$L=%d$" % L)

        self.subfigure.legend(loc="upper right", numpoints=1, handlelength=0.5, ncol=4, columnspacing=0.25, handletextpad=0.25, borderaxespad=0.1)

















