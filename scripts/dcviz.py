#DCVIZ

from DCViz_sup import DCVizPlotter

import numpy
import numpy as np
import os
import re
from numpy import exp, floor, ceil
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.special import lambertw
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib import rcParams

from mpl_toolkits.mplot3d import Axes3D, proj3d

E0_tex = r"P_\lambda"
fig_4_x = 0.05
fig_4_y = 0.85
fig_4_fs = 30

my_props = {

    "fmt" : {
            "markersize" : 7,
            "markeredgewidth" : 1.5,
            "linewidth" : 1,
            "fillstyle" : "none"
    }

}

def int_format_func(v, _):
    if int(v) == v:
        return r"$%d$" % v
    else:
        return ""

INTFORMATTER = FuncFormatter(int_format_func)

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

        rmslabel=r"\sigma(\mathbf{h})"
        sslabel=r"\sigma(s)"
        whlabel=r"\langle d_i \rangle"

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


        self.dh_fig.set_xlabel(r"$t\, [1/\nu]$")

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
            sfig.text(1-0.1, 0.5, r"$\mathrm{(%s)}$" % labels2[i],
                      horizontalalignment="left",
                      verticalalignment="center",
                      fontsize=30,
                      transform=sfig.axes.transAxes)

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
        ax3.set_ylabel(r"$C_V/kL^2$", labelpad=15)
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

        shifts = data[self.get_family_index_from_name("shifts.arma")]
        values = data[self.get_family_index_from_name("values.arma")]

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

    figMap = {"figure": ["subfigure"], "figure3D": [], "figure_projections" : ["proj", "proj2"]}

    # isFamilyMember = True
    # loadSequential = True
    #loadLatest = True
    #ziggyMagicNumber = 100

    tight = False

    fig_size = [10, 10]
    specific_fig_size = {"figure_projections": [15, 10]}

    def plot(self, data):

        # data = sum(_data)/len(_data)

        if not self.argv:
            self.lim = exp(-2)
        else:
            self.lim = float(self.argv[0])

        data/= data.max()

        Lm, Wm = data.shape

        L = (Lm - 1)/2
        W = (Wm - 1)/2

        ax = self.subfigure.imshow(data, extent=[-L, L, -W, W], origin="lower", aspect="auto", interpolation="none", cmap="Blues")
        self.subfigure.set_xlabel(r"$\delta x$")
        self.subfigure.set_ylabel(r"$\delta y$")

        ax = Axes3D(self.figure3D)

        xpos, ypos = np.meshgrid(np.arange(Lm), np.arange(Wm))

        ax.plot_surface(xpos, ypos, data, cstride=1, rstride=1)

        d1 = np.diag(data)
        d2 = np.diag(np.flipud(data))

        xl = np.linspace(-L, L, Lm)
        xw = np.linspace(-W, W, Wm)

        x1 = np.sqrt(2)*xl
        x2 = np.sqrt(2)*xw

        dl = data[L, :]
        dw = data[:, W]

        # f = lambda x, a, b: a*exp(-abs(x)/b)
        #
        # p0 = (1., 1.)
        # pl, _ = curve_fit(f, xl, dl, p0)
        # pw, _ = curve_fit(f, xw, dw, p0)
        # p1, _ = curve_fit(f, x1, d1, p0)
        # p2, _ = curve_fit(f, x2, d2, p0)

        self.proj2.plot(xl, dl, '-o', label="X")
        self.proj2.plot(xw, dw, '-^', label="Y")
        self.proj2.plot(x1, d1, '-d', label="11 diag")
        self.proj2.plot(x2, d2, '-s', label="-11 diag")

        print self.find_corrl(xl, dl, '-o', "X")
        print self.find_corrl(xw, dw, '-^', "Y")
        print self.find_corrl(x1, d1, '-d', "11 diag")
        print self.find_corrl(x2, d2, '-s', "-11 diag")

        self.proj.legend()
        self.proj.set_xlabel(r"$\delta r$")
        self.proj.set_ylabel("corr")
        self.proj.plot(x1, np.zeros_like(x1), 'k--')

    def find_corrl(self, x, y, c, l):
        I = np.where(y > self.lim)
        J = np.where(x[I] > 0)
        K = slice(0,4)

        X = x[I][J][K]

        # shift = y.min()
        shift = 0
        Y = np.log((y-shift)[I][J][K])

        self.proj.plot(X, Y, c, label=l)

        f = lambda x, a, b: a*x + b

        p0 = (-1., -1.)

        pl, cl = curve_fit(f, X, Y, p0)
        print pl

        _x = np.linspace(X.min(), X.max(), 10000)
        self.proj.plot(_x, f(_x, *pl))

        return -1./pl[0]

class LatticediffSpeeds(DCVizPlotter):

    nametag = "confined_\w+_(\w+)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    alpha = 1
    h0 = 10
    h0c = h0 - 2
    plot_analytical = True

    figMap = {"f0" : ["subfigure", "subfigure2", "errfigH"],
              "f2" : ["subfigure3","subfigure6"],
              "f3" : "subfigure4",
              "f4" : "subfigure5",
              "f6" : "subfigure7"}

    # plotOnly = ["errfig"]

    legend=True

    def adjust(self):

        # if self.legend or False:
        #     self.adjust_maps["f0"]["top"] = 0.89
        # else:
        self.adjust_maps["f0"]["top"] = 0.97
        self.adjust_maps["f0"]["bottom"] = 0.09
        # self.adjust_maps["f0"]["top"] = 0.86
        self.adjust_maps["f0"]["right"] = 0.95
        self.adjust_maps["f0"]["hspace"] = 0.15
        self.adjust_maps["f0"]["left"] = 0.185

        self.adjust_maps["f2"]["top"] = 0.9
        self.adjust_maps["f2"]["bottom"] = 0.13
        self.adjust_maps["f2"]["right"] = 0.95
        self.adjust_maps["f2"]["hspace"] = 0.16
        self.adjust_maps["f2"]["left"] = 0.15

        # self.adjust_maps["errfig"]["top"] = 0.93
        # self.adjust_maps["errfig"]["bottom"] = 0.25
        # self.adjust_maps["errfig"]["right"] = 0.97
        # self.adjust_maps["errfig"]["left"] = 0.14


    fig_size = [6, 6]
    specific_fig_size = {"f0": [6, 9],
                         "f2": [6, 7]}

    ceq_override = None

    def ceq_analytical(self):
        return exp(-self.alpha*3)

    def ceq(self):

        if self.ceq_override is not None:
            return self.ceq_override
        else:
            return self.ceq_analytical()

    def asympt(self, s0):
        return s0*self.h0c/(1/self.ceq() - 1)

    def ceq_from_asympt(self, a, c0):

        return (self.h0c*c0 - a)/(self.h0c - a)

    def K_over_c(self, supersaturation):

        hs = self.h0c*(1 - self.ceq()*(1 + supersaturation))

        return hs/(1 - self.ceq())

    def c(self, k):
        return k*(1/self.ceq() - 1)

    def xi(self, K_over_c):
        kappa = self.h0c/K_over_c - 1
        return kappa*exp(kappa)

    def k_convert(self, heights, times, supersaturation):

        if supersaturation == 0:
            return np.zeros_like(times)

        K_over_c = self.K_over_c(supersaturation)

        G = (self.h0c - heights)/K_over_c - 1

        xi = self.xi(K_over_c)

        expterm = G*exp(G)/xi

        core = np.log(expterm)

        c = -K_over_c*core/times

        k = c/(1/self.ceq() - 1)

        return k

    def analytical(self, times, supersaturation, k):

        K_over_c = self.K_over_c(supersaturation)
        c = self.c(k)

        core = self.xi(K_over_c)*exp(-c*times/K_over_c)

        return self.h0c - K_over_c*(1 + lambertw(core).real)

    def analytical_ss(self, times, supersaturation, k):

        if abs(supersaturation) < 1E-10:
            return np.zeros_like(times)

        K_over_c = self.K_over_c(supersaturation)

        c_over_k = 1/self.ceq() - 1

        c = self.c(k)

        core = self.xi(K_over_c)*exp(-c*times/K_over_c)

        return c_over_k/(1 + 1/lambertw(core).real)


    def find_k(self, h, t, supersaturation):

        if supersaturation == 0:
            return np.nan

        analytical = lambda time, k: self.analytical(time, supersaturation, k)

        p, _ = curve_fit(analytical, t, h, (0.1))

        return p[0]

    def find_k2(self, h, t, supersaturation):

        if supersaturation == 0:
            return np.nan

        analytical = lambda time, k: self.analytical_ss(time, supersaturation, k)

        p, _ = curve_fit(analytical, t, h, (0.1))

        return p[0]

    def find_ks(self, h, t, supersaturation):

        if supersaturation == 0:
            return np.nan

        analytical = lambda time, k, ss: self.analytical(time, ss, k)

        p, _ = curve_fit(analytical, t, h, (1, supersaturation))

        return p

    def find_ks2(self, h, t, supersaturation):

        if supersaturation == 0:
            return np.nan

        analytical = lambda time, k, ss: self.analytical_ss(time, ss, k)

        p, _ = curve_fit(analytical, t, h, (1, supersaturation))

        return p

    def reshape(self, a):
        if np.rank(a) == 1:
            return a.reshape(1, len(a))
        else:
            return a

    def plot(self, data):

        supersaturations = self.get_family_member_data(data, "supersaturations")
        all_heights = self.get_family_member_data(data, "heights")
        all_times = self.get_family_member_data(data, "times")
        all_conc = self.get_family_member_data(data, "concentrations")
        lengths = self.get_family_member_data(data, "lengths")

        all_heights = self.reshape(all_heights)
        all_times = self.reshape(all_times)
        all_conc = self.reshape(all_conc)

        if len(self.argv) > 1:
            self.h0 = float(self.argv[0])
            self.alpha = float(self.argv[1])

            self.h0c = self.h0 - 2

        nolabels = "-nolabels" in self.argv

        if "ls" in self.argv:
            idx = self.argv.index("ls")
            lsfac = float(self.argv[idx+1])
        else:
            lsfac = 2./3

        if "tmax" in self.argv:
            idx = self.argv.index("tmax") + 1
            tmax = float(self.argv[idx])
        else:
            tmax = None

        if "title" in self.argv:
            idx = self.argv.index("title") + 1
            self.subfigure3.set_title(r"$\mathrm{%s}$" % self.argv[idx], y=1.1, fontsize=30)

        if "spans" in self.argv:
            idx = self.argv.index("spans") + 1
            spans = eval(self.argv[idx]) #e.g. [[1,2],[3,4],[5,6],[7,8],[9,10]]
        else:
            spans = None

        all_supersaturations = all_conc/self.ceq() - 1

        if "onlysim" in self.argv:
            self.plot_analytical = False

        all_asympts = []
        all_ss = []
        all_k = []
        all_kss = []
        all_ceq = []

        max_l = lengths.max()

        for i in range(len(lengths)):

            l = lengths[i]

            t_end = all_times[i, l-1]

            if i == 0:
                t_min = t_end
            else:
                if t_end < t_min:
                    t_min = t_end

        if tmax is not None:
            if t_min > tmax:
                t_min = tmax

        mean_errh = np.zeros(max_l)
        mean_errss = np.zeros(max_l)
        mean_Ts = np.zeros(max_l)
        mean_cs = np.zeros(max_l)

        sfac = 1

        for i, supersaturation in enumerate(supersaturations):

            l = lengths[i]

            T = all_times[i, :l]
            I = np.where(T <= t_min)
            T = T[I]

            if "corr" in self.argv:
                l_c = len(all_conc[i, :l][I])
                self.ceq_override = all_conc[i, :l][I][(3*l_c)/4:].mean()

            all_ceq.append(100*abs(self.ceq()-self.ceq_analytical())/self.ceq_analytical())

            all_supersaturations[i, :l] = all_conc[i, :l]/self.ceq() - 1

            s0 = all_supersaturations[i, 0]

            all_ss.append(s0)


            H = sfac*all_heights[i, :l][I]
            label = r"$\Omega(0)=%.2f$" % s0

            all_asympts.append(sfac*H[(2*len(H))/3:].mean())

            if i == 0:
                _label = r"$\mathrm{KMC}$"
            else:
                _label = None

            self.subfigure.plot(T, H, "r-", label=_label)

            SS = all_supersaturations[i, :l][I]
            self.subfigure2.plot(T, SS, "r-", label=_label)

            if abs(supersaturation) > 0.1:
                self.subfigure4.plot(self.k_convert(sfac*H[1:], T[1:], s0)[1:], label=label)
                self.subfigure4.set_ylim(0, 1)

            if not self.plot_analytical:
                continue

            if i == 0:
                _labela = r"$\mathrm{Analytical}$"
            else:
                _labela = None

            ls = int(lsfac*l)
            kval1 = self.find_k(sfac*H[1:ls], T[1:ls], s0)
            kval2 = self.find_k2(SS[1:ls], T[1:ls], s0)

            kval = 0.5*(kval1 + kval2)

            if self.plot_analytical:
                analytical_h = self.analytical(T, s0, kval)
                analytical_ss = self.analytical_ss(T, s0, kval)
                self.subfigure.plot(T, analytical_h, "k--", label=_labela)
                self.subfigure2.plot(T, analytical_ss, "k--", label=_labela)

                errH = H - analytical_h
                errSS = SS - analytical_ss

                mean_errh[:l][I] +=  abs(errH)
                mean_errss[:l][I] += abs(errSS)
                mean_Ts[:l][I] += T
                mean_cs[:l][I] += 1


            if supersaturation != sorted(supersaturations)[len(supersaturations)/2]:
                all_k.append(kval)
                all_kss.append(s0)

        time_label = r"$\nu t$"

        cI = np.where(mean_cs >= 1)
        mean_errh = mean_errh[cI]/mean_cs[cI]
        mean_errss = mean_errss[cI]/mean_cs[cI]
        mean_Ts = mean_Ts[cI]/mean_cs[cI]


        self.errfigH.plot(mean_Ts[:len(mean_errh)], 1000*mean_errh/len(all_supersaturations), 'r-', label=r"$h(t)/l_0$")
        self.errfigH.plot(mean_Ts[:len(mean_errss)], 1000*mean_errss/len(all_supersaturations), color="0.75",
                          linestyle='-', label=r"$\Omega(t)$")

        y0, y1 = self.errfigH.get_ylim()
        self.errfigH.set_ybound(y0-0.25*(y1 - y0))

        if not nolabels:
            leg = self.errfigH.legend(loc="center",
                                      numpoints=1,
                                      ncol=2,
                                      handlelength=1.0,
                                      markerscale=20.0,
                                      borderpad=0.2,
                                      labelspacing=0.2,
                                      columnspacing=1.0,
                                      handletextpad=0.5,
                                      borderaxespad=0.0,
                                      frameon=False,
                                      fontsize=20,
                                      bbox_to_anchor=(0.7, 0.10))

            leg.get_frame().set_fill(not (self.toFile and self.transparent))
            for legobj in leg.legendHandles:
                legobj.set_linewidth(4.0)

        self.subfigure7.plot(all_ss, all_ceq, 'k^')
        self.subfigure7.set_xlabel(r"$\Omega(0)$")
        if not nolabels:
            self.subfigure7.set_ylabel(r"$\epsilon_\mathrm{eq} [\%]$")
        else:
            self.subfigure7.set_yticklabels([])

        self.subfigure4.legend()

        self.subfigure.set_xlim(0, t_min)
        #self.subfigure.set_ylim(self.subfigure.get_ylim()[0]*1.1, self.subfigure.get_ylim()[1]*1.1)
        #self.subfigure.set_xlabel(time_label)
        if not nolabels:
            self.subfigure.set_ylabel(r"$h(t)/l_0$")
        else:
            self.subfigure.set_yticklabels([])

        self.subfigure.xaxis.set_ticks([])

        if self.legend and not nolabels:
            lg = self.subfigure2.axes.legend(loc="upper center",
                                             numpoints=1,
                                             ncol=1,
                                             handlelength=1.0,
                                             borderpad=0.2,
                                             labelspacing=0.2,
                                             columnspacing=1.25,
                                             handletextpad=0.25,
                                             borderaxespad=0.0,
                                             frameon=False,
                                             fontsize=20,
                                             bbox_to_anchor=(0.77, 0.9))
            #bbox_to_anchor=(0.5, 1.25)
            lg.get_frame().set_fill(not (self.toFile and self.transparent))


        self.subfigure2.set_xlim(0, t_min)
        # self.subfigure2.set_xlabel(time_label)
        if not nolabels:
            self.subfigure2.set_ylabel(r"$\Omega(t)$")
        else:
            self.subfigure2.set_yticklabels([])

        self.subfigure2.xaxis.set_ticks([])

        if spans is None:
            s1ylim = max([abs(x) for x in self.subfigure.get_ylim()])
            self.subfigure.set_ylim([-s1ylim, s1ylim])

            lowSS = round(min(all_supersaturations[:, 0]), 1)
            highSS = round(max(all_supersaturations[:, 0]), 1)
            lim2 = max(abs(lowSS), abs(highSS))
            self.subfigure2.set_ylim(-lim2, lim2)
        else:
            self.subfigure.set_ylim(spans[2])
            self.subfigure2.set_ylim(spans[3])
            self.errfigH.set_ylim(spans[4])

        self.errfigH.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "" if v < 0 else r"$%g$" % v))

        self.errfigH.set_xlabel(time_label)
        if not nolabels:
            self.errfigH.set_ylabel(r"$10^3\epsilon$")
        else:
            self.errfigH.set_yticklabels([])

        if not self.plot_analytical:
            return

        ss = np.linspace(supersaturations.min(), supersaturations.max())

        self.subfigure3.plot(ss, self.asympt(ss)/self.h0c*100, "r--", label=r"$(c(0) - \langle c_\mathrm{eq}\rangle)/(1 -  \langle c_\mathrm{eq}\rangle)$")
        self.subfigure3.plot(all_ss, np.array(all_asympts)/self.h0c*100, "ks", label=r"$\langle h_\mathrm{eq} \rangle /h_l $", **my_props["fmt"])

        #self.subfigure3.yaxis.set_major_formatter(FuncFormatter(self.axisFormat))
        self.subfigure3.xaxis.set_ticks([-1, -0.5, 0, 0.5, 1])
        self.subfigure3.set_xbound(-1)
        self.subfigure3.xaxis.set_ticklabels([])
        if not nolabels:
            self.subfigure3.set_ylabel(r"$\mathrm{Eq.~(11)\,\,sides\,\,[\%]}$")
        else:
            self.subfigure3.set_yticklabels([])

        if not nolabels:
            lg = self.subfigure3.axes.legend(loc="center",
                                                 numpoints=1,
                                                 ncol=1,
                                                 handlelength=1.0,
                                                 borderpad=0.2,
                                                 labelspacing=0.2,
                                                 columnspacing=1.25,
                                                 handletextpad=0.25,
                                                 borderaxespad=0.0,
                                                 frameon=False,
                                                 fontsize=20,
                                                 bbox_to_anchor=(0.35, 0.75))
            lg.get_frame().set_fill(not (self.toFile and self.transparent))


        m = sum(all_k)/len(all_k)
        kscaled = np.array(all_k)/m

        self.subfigure6.plot(all_kss, kscaled, "ks", **my_props["fmt"])
        self.subfigure6.set_xlabel(r"$\Omega(0)$")
        if not nolabels:
            self.subfigure6.set_ylabel(r"$k/\langle k\rangle$")
        else:
            self.subfigure6.set_yticklabels([])

        if spans is None:
            self.subfigure6.set_ylim(0, max(kscaled)*1.1)
        else:
            self.subfigure3.set_ylim(spans[0])
            self.subfigure6.set_ylim(spans[1])

        self.subfigure6.set_xbound(-1)
        self.subfigure6.xaxis.set_ticks([-1, -0.5, 0, 0.5, 1])

    @staticmethod
    def axisFormat(value, index):

        if int(value) == value:
            return r'$%d$' % value
        else:
            return ''


class ignisSOS(DCVizPlotter):

    nametag = "ignisSOS(.*)\.ign"

    def plot(self, data):

        name_string, l, w = self.loader.get_metadata()

        names = [desc.split("@")[0] for desc in name_string.split()]
        yname = self.argv[0]

        if "every" in self.argv:
            i = self.argv.index("every")
            every = int(self.argv[i+1])
        else:
            every = 1

        if "-nonzero" in self.argv:
            nonzero = True
        else:
            nonzero = False

        if yname in names:
            data_y = self.scan_for_name(data, yname, names)

            try:
                xname = self.argv[1]

                data_x = self.scan_for_name(data, xname, names)

                self.subfigure.set_xlabel(xname)

            except:
                data_x = np.arange(0, len(data_y))

            data_x = data_x[::every]
            data_y = data_y[::every]

            if nonzero:
                I = np.where(data_y != 0)
                data_x = data_x[I]
                data_y = data_y[I]

            self.subfigure.plot(data_x, data_y)
            self.subfigure.set_ylabel(yname)

        else:
            raise RuntimeError("Invalid option %s." % yname)

    def scan_for_name(self, data, targetname, names):

        if len(data.shape) == 1:
            if targetname == names[0]:
                return data
        else:
            for _data, name in zip(data, names):

                if name == targetname:
                    return _data

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


class MechEq(DCVizPlotter):

    nametag = "mechEq(.*)\.arma"

    isFamilyMember = True

    loadSequential = True

    hugifyFonts = True

    def plot(self, data):
        hls, Fs = data
        self.subfigure.plot(hls, Fs)
        self.subfigure.plot([hls[0], hls[-1]], [0, 0], "k--")


class mftpc(DCVizPlotter):

    nametag = "mfptc_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    fig_size = [6, 4]

    c = r"\kappa"

    def adjust(self):

        self.adjust_maps["figure"]["top"] = 0.95
        self.adjust_maps["figure"]["bottom"] = 0.20
        self.adjust_maps["figure"]["right"] = 0.96
        self.adjust_maps["figure"]["left"] = 0.18

    def plot(self, data):

        cs = self.get_family_member_data(data, "cs")
        speeds = self.get_family_member_data(data, "growthspeeds")

        X = -np.log(cs)

        self.subfigure.plot(X, 100*speeds[0], "ks", label=r"$\mathrm{Radial}$", **my_props["fmt"])
        self.subfigure.plot(X, 100*speeds[1], "ro", label=r"$\mathrm{Pathfinding}$", **my_props["fmt"])


        ymin = self.subfigure.get_ylim()[0]
        if len(self.argv) > 0:
            cinvEq = float(self.argv[0])

            self.subfigure.plot([0, cinvEq], [0, 0], "k-.")
            self.subfigure.plot([cinvEq, cinvEq], [ymin, 0], "k-.")
            self.subfigure.text(cinvEq, ymin+1.5, r"$\log \left(1/%s'\right) \sim %.2f$" % (self.c, round(cinvEq, 2)), verticalalignment="bottom", horizontalalignment="left", fontsize=20)

            # textlabel = r"\kappa'=%.2f" % cinvEq
            # self.subfigure.text(0.75, 0.75, r"$\mathrm{%s}$" % textlabel, transform=self.subfigure.axes.transAxes, horizontalalignment="center")

        self.subfigure.set_xlabel(r"$\log \left(1/%s\right)$" % self.c)
        self.subfigure.set_ylabel(r"$\Delta \langle h\rangle/\Delta T [l_0 \nu/10^2]$")

        self.subfigure.legend(loc="upper right",
                              numpoints=1,
                              ncol=1,
                              handlelength=1.0,
                              borderpad=0.2,
                              labelspacing=0.2,
                              columnspacing=1.25,
                              handletextpad=0.25,
                              borderaxespad=0.0,
                              frameon=False,
                              fontsize=20,
                              bbox_to_anchor=(0.95, 0.9))


class cconv(DCVizPlotter):

    nametag = "cconv_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    stack = "H"
    figMap = {"f1": ["subfigure1", "subfigure2"]}

    fig_size = [8, 5]

    labelSize= 25
    ticklabelSize=21
    fontSize = 20  # Only invoked by hugifyFonts = True
    tickSize = 2

    def adjust(self):
        for f in self.figure_names:
            self.adjust_maps[f]["top"] = 0.9
            self.adjust_maps[f]["bottom"] = 0.15
            self.adjust_maps[f]["right"] = 0.98
            self.adjust_maps[f]["left"] = 0.11
            self.adjust_maps[f]["wspace"] = 0.05

    def plotsingle(self, subfigure, alphas, heights, set, do_legend=False):

        subfigure.plot([alphas.min(), alphas.max()], [2.85, 2.85], "r--", linewidth=3, label=r"$\log(1/\kappa')$")

        xvals = np.zeros((len(heights), len(set)))
        crads = np.zeros_like(xvals)
        err = np.zeros_like(xvals)

        heightFMT = ['o', '^', 's', 'd']

        for j, (a, b, i, k) in enumerate(set):
            alpha = alphas[k]

            xvals[i, j] = alpha
            crads[i, j] = (a+b)/2
            err[i, j] = -np.log(b/a)/2

        heights = sorted(heights)

        for i, height in enumerate(heights):

            I = np.where(xvals[i, :] != 0)

            subfigure.errorbar(xvals[i, :][I],
                                -np.log(crads[i, :][I]),
                                err[i, :][I],
                                fmt="k"+heightFMT[int(i)],
                                ecolor="k",
                                linestyle="none",
                                label=r"$\Delta h = %d$" % height,
                                **my_props["fmt"])

            print height, heightFMT[int(i)]

        subfigure.set_xlabel(r"$\alpha$")

        subfigure.text(0.5, 1.05, r"$\mathrm{%s}$" % ("Radial" if do_legend else "Pathfinding"),
                       horizontalalignment="center",
                       verticalalignment="center",
                       fontsize=25,
                       transform=subfigure.axes.transAxes)

        if do_legend:
            subfigure.legend(loc="lower center",
                                  numpoints=1,
                                  ncol=1,
                                  handlelength=1.0,
                                  borderpad=0.2,
                                  labelspacing=0.2,
                                  columnspacing=1.25,
                                  handletextpad=0.25,
                                  borderaxespad=0.0,
                                  frameon=False,
                                  fontsize=24,
                                  bbox_to_anchor=(0.3, 0.05)) \
                .get_frame().set_fill(not (self.toFile and self.transparent))

            subfigure.set_ylabel(r"$\log(1/\kappa)$")

        else:
            subfigure.yaxis.set_ticklabels([])

        index_list = [[1, 2], [3, 4], [1, 1], [5, 4], [3, 2], [7, 4], [2, 1]]
        def f(v, i):

            a, b = index_list[i]

            if b == 1:
                return r"$%d$" % a
            else:
                return r"\Large{$%d/%d$}" % (a, b)

        subfigure.axes.xaxis.set_major_formatter(FuncFormatter(f))

        #   subfigure.set_xlim(0.4, 2.1)
        subfigure.set_xlim(0.45, alphas.max()+0.05)
        subfigure.set_ylim(0, 3.5)
        subfigure.xaxis.set_ticks([x[0]/float(x[1]) for x in index_list])
        # subfigure.axes.tick_params(axis='x', which='major', pad=15)


    def plot(self, data):

        heights = self.get_family_member_data(data, "heights")
        alphas = self.get_family_member_data(data, "alphas")

        radial = self.get_family_member_data(data, "radial")
        self.plotsingle(self.subfigure1, alphas, heights, radial, True)

        pathfind = self.get_family_member_data(data, "pathfind")
        self.plotsingle(self.subfigure2, alphas, heights, pathfind)



class AutoCorrWoot(DCVizPlotter):

    nametag = "acf_(.*)\.npy"

    isFamilyMember = True

    figMap = {"rmsfigure" : ['RMSfig1', 'RMSfig2', 'RMSfig4'],
              "cfigure" :   ['CFig1', 'CFig2', 'CFig4']}

    styles = ["s", "^", "o", "d"]
    colors = ["k", "b", "r", "g"]

    names = {"uniform"  : "Uniform",
             "lattice"  : "Lattice",
             "radial"   : "Radial MFPT",
             "pathfind" : "Pathfind MFPT"}

    prios = {"uniform"  : 0,
             "lattice"  : 1,
             "radial"   : 2,
             "pathfind" : 3}

    styles = {"uniform"  : "ks",
              "lattice"  : "ro",
              "radial"   : "b^",
              "pathfind" : "gd"}

    hugifyFonts = True

    def adjust(self):
        for figname in self.figure_names:
            self.adjust_maps[figname]["top"] = 0.98
            self.adjust_maps[figname]["bottom"] = 0.1
            self.adjust_maps[figname]["hspace"] = 0.11

        self.adjust_maps["cfigure"]["right"] = 0.87
        self.adjust_maps["cfigure"]["left"] = 0.14

        self.adjust_maps["rmsfigure"]["right"] = 0.93
        self.adjust_maps["rmsfigure"]["left"] = 0.20

    fig_size = [5, 7]

    def plot(self, data):

        heights = self.get_family_member_data(data, "heights")

        types = []
        for name in self.familyFileNames:

            if "alphas" in name:
                _ih, type = re.findall(r"acf_h(\d+)\_(\w+)\_alphas\.npy", name)[0]

                if int(_ih) != 0:
                    continue

                types.append(type)

        types = sorted(types, key=lambda x: self.prios[x])

        print types

        #Rfs = [self.RMSfig1, self.RMSfig2, self.RMSfig4]
        Rfs = self.figures[self.figure_names.index("rmsfigure")][1:][::-1]
        #Cfs = [self.CFig1, self.CFig2, self.CFig4]
        Cfs = self.figures[self.figure_names.index("cfigure")][1:][::-1]

        for _ih, height in enumerate(heights):
            for it, type in enumerate(types):
                alphas = self.get_family_member_data(data, "h%d_%s_alphas" % (_ih, type))
                RMSes = self.get_family_member_data(data, "h%d_%s_RMSes" % (_ih, type))
                Cs = self.get_family_member_data(data, "h%d_%s_Cs" % (_ih, type))

                s = self.styles[type]

                if self.argv:
                    if _ih == 2:
                        continue
                    elif _ih == 3:
                        ih = 2
                    else:
                        ih = _ih
                else:
                    ih = _ih

                Rfs[ih].plot(alphas, RMSes, s, label=self.names[type], **my_props["fmt"])

                # if self.argv:
                #     c = Cs[:, int(self.argv[0])]
                # else:
                c = 0.5*(Cs[:, 0] + Cs[:, 1])

                Cfs[ih].plot(alphas, c, s, label=self.names[type], **my_props["fmt"])
                # self.CFig_diag.plot(alphas, Cs[:, 1], s, label=type, **my_props["fmt"])

                if type == "lattice" and height > 5:
                    Rfs[ih].text(0.75, 2.5, r"$\uparrow %.2f$" % max(RMSes))

        label = r"\Delta h"

        def format_func0(v, i):
            return ""

        def format_func(v, i):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return ""

        null_formatter = FuncFormatter(format_func0)
        int_formatter = FuncFormatter(format_func)

        for _ih, h in enumerate(heights):

            if self.argv:
                if _ih == 2:
                    continue
                elif _ih == 3:
                    ih = 2
                else:
                    ih = _ih
            else:
                ih = _ih


            rf = Rfs[ih]
            cf = Cfs[ih]

            rf.set_ylabel(r"$\sigma(\mathbf{h})/l_0$")
            rf.xaxis.set_ticks([0.5, 1, 1.5, 2])
            rf.set_xlim(0.4, 2.1)
            rf.set_ylim(0, 3)
            rf.axes.yaxis.set_major_formatter(int_formatter)

            # ax = rf.axes.twinx()
            # ax.set_ylabel(r"$%s=%.1f$" % (label, h), labelpad=15)
            # ax.yaxis.set_ticks([])
            # ax.yaxis.set_ticklabels([])

            cf.set_ylabel(r"$\xi/l_0$")
            cf.yaxis.set_ticks([0, 1, 2, 3])
            cf.xaxis.set_ticks([0.5, 1, 1.5, 2])
            cf.set_xlim(0.4, 2.1)
            cf.set_ylim(0, 3.5)
            cf.axes.yaxis.set_major_formatter(int_formatter)

            ax = cf.axes.twinx()
            ax.set_ylabel(r"$%s/l_0=%d$" % (label, h), labelpad=15)
            ax.yaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])

            if ih == 0:

                cf.set_xlabel(r"$\alpha$")
                rf.set_xlabel(r"$\alpha$")

            else:
                rf.axes.xaxis.set_major_formatter(null_formatter)
                cf.axes.xaxis.set_major_formatter(null_formatter)

        # self.CFig.text(0.05, 0.75, r"$\langle d_i \rangle = %.1f$" % float(height),
        #                horizontalalignment="left",
        #                verticalalignment="center",
        #                fontsize=24,
        #                transform=self.CFig.axes.transAxes)



        l = self.CFig1.axes.legend(loc="center",
                                   numpoints=1,
                                   ncol=2,
                                   handlelength=1.0,
                                   borderpad=0.2,
                                   labelspacing=0.2,
                                   columnspacing=0.3,
                                   handletextpad=0.25,
                                   borderaxespad=0.0,
                                   frameon=False,
                                   fontsize=12,
                                   bbox_to_anchor=(0.66, 0.15))

        l.get_frame().set_fill(not (self.toFile and self.transparent))


class ExtraNeighborPrev(DCVizPlotter):

    nametag = "^DEPR_extraneighbor_(.*)\.npy"

    isFamilyMember = True

    figMap = {"f2" : ["subfigure11",
                      "subfigure12",
                      "subfigure13"],

              "f3" : ["subfigure21",
                      "subfigure22",
                      "subfigure23"],

              "f4" : ["subfigure31",
                      "subfigure32",
                      "subfigure33"],

              "f5" : ["csubfigure11",
                      "csubfigure12",
                      "csubfigure13"],

              "f6" : ["csubfigure21",
                      "csubfigure22",
                      "csubfigure23"],

              "f7" : ["csubfigure31",
                      "csubfigure32",
                      "csubfigure33"]}

    hugifyFonts = True

    specific_fig_size = {"f2": [4.5, 8],
                         "f3": [3.6, 8],
                         "f4": [5, 8],
                         "f5": [4.5, 8],
                         "f6": [3.6, 8],
                         "f7": [5, 8]}
    def adjust(self):
        for figname in self.figure_names:

            if figname in ["f2", "f5"]:
                self.adjust_maps[figname]["right"] = 0.94
                self.adjust_maps[figname]["left"] = 0.2
            elif figname in ["f3", "f6"]:
                self.adjust_maps[figname]["right"] = 0.94
                self.adjust_maps[figname]["left"] = 0.03
            else:
                self.adjust_maps[figname]["right"] = 0.68
                self.adjust_maps[figname]["left"] = 0.03

            self.adjust_maps[figname]["top"] = 0.95
            self.adjust_maps[figname]["bottom"] = 0.10
            self.adjust_maps[figname]["hspace"] = 0.06

    plotOnly = ["f4", "f2", "f3"]

    def trans_mat(self, mat):
        return mat

        return mat.transpose()

        t = np.zeros_like(mat)

        l, w = t.shape

        for i in range(l):
            for j in range(w):
                t[i, j] = i*j/float(l*w)

        t[0, :] = 0.1
        t[:, 0] = 0.1
        t[-1, :] = 0.1
        t[:, -1] = 0.1

        return t

    tight = False

    def plot(self, data):

        omegas = self.get_family_member_data(data, "omegas")
        s0s = self.get_family_member_data(data, "s0s")
        alphas = self.get_family_member_data(data, "alphas")
        Pls = self.get_family_member_data(data, "Pls")
        ld = 5.0
        conv = ld*(1-exp(-1./ld))
        F0s = Pls/conv

        X, Y = np.meshgrid(alphas, F0s, indexing='ij')

        for io, omega in enumerate(omegas):
            cmat = self.get_family_member_data(data, "cmat_omega%d" % io)

            null_formatter = FuncFormatter(lambda v, i: "")

            for is0, s0 in enumerate(s0s):

                C = self.trans_mat(cmat[is0, :, :])

                sfig = eval("self.subfigure%d%d" % (io+1, is0+1))
                csfig = eval("self.csubfigure%d%d" % (io+1, is0+1))

                im = sfig.pcolor(C.transpose(), vmin=0, vmax=1, cmap="gist_earth_r")

                d = 0.1
                levs = [-d] + [d*i for i in range(0, int(1/d + 2))]

                #levs = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
                im2 = csfig.contourf(X, Y, C, vmin=0, vmax=1, cmap="gist_earth_r", levels=levs)

                fail = np.where(C == -1)

                # succ = np.where(C != -1)
                # Cs = C.copy()
                # Cs[fail] = 1
                # sfig.contour(X, Y, Cs)

                for y, x in zip(*fail):
                    sfig.scatter(x+0.5, y+0.5, c = 'k', marker='x')

                if is0 == 0:
                    sfig.set_title(r"$\Omega = %.1f$" % omega)
                    csfig.set_title(r"$\Omega = %.1f$" % omega)

                if io == 0:
                    label = r"$F_0/E_bA$"
                    sfig.set_ylabel(label)
                    csfig.set_ylabel(label)
                else:
                    sfig.yaxis.set_ticklabels([])
                    csfig.yaxis.set_ticklabels([])

                    if io == 2:
                        sfig.yaxis.set_label_position("right")
                        csfig.yaxis.set_label_position("right")
                        sfig.set_ylabel(r"$\sigma_0 = %.1f$" % s0, labelpad=10)
                        csfig.set_ylabel(r"$\sigma_0 = %.1f$" % s0, labelpad=10)


                # sfig.text(0.1, 0.5, r"$\sigma_0 = %.2f$" % s0, verticalalignment="bottom", horizontalalignment="left", fontsize=20)

                xi = np.arange(1, len(alphas), 2)
                yi = np.arange(1, len(Pls), 3)

                ax = sfig.axes
                xtics = [r"$%.1f$" % a for a in alphas[xi]]
                ytics = [r"$%.2f$" % p for p in Pls[yi]]

                ax.set_xticklabels(xtics, minor=False)
                ax.set_yticklabels(ytics, minor=False)

                if is0 < len(s0s) - 1:
                    sfig.axes.xaxis.set_major_formatter(null_formatter)
                    csfig.axes.xaxis.set_major_formatter(null_formatter)
                else:
                    label = r"$E_b/kT$"
                    sfig.set_xlabel(label)
                    csfig.set_xlabel(label)

                sfig.set_xlim(0, len(alphas))
                sfig.set_ylim(0, len(F0s))

        cbar_ax = self.f4.add_axes([0.85, 0.2, 0.05, 0.7])
        cbar_ax2 = self.f7.add_axes([0.85, 0.2, 0.05, 0.7])

        clabel = r'$\Theta$'
        cbar_ax.set_title(clabel)
        cbar_ax2.set_title(clabel)

        self.f4.colorbar(im, cax=cbar_ax)
        self.f7.colorbar(im, cax=cbar_ax2)


class ExtraNeighborPrev2(DCVizPlotter):

    nametag = "^DEPRextraneighbor_(.*)\.npy"

    isFamilyMember = True

    figMap = {"f2" : ["subfigure11",
                      "subfigure12"],

              "f3" : ["subfigure21",
                      "subfigure22"],

              "f4" : ["subfigure31",
                      "subfigure32"],

              "f5" : ["csubfigure11",
                      "csubfigure12"],

              "f6" : ["csubfigure21",
                      "csubfigure22"],

              "f7" : ["csubfigure31",
                      "csubfigure32"]}

    hugifyFonts = True

    specific_fig_size = {"f2": [4.5, 8],
                         "f3": [3.6, 8],
                         "f4": [5, 8],
                         "f5": [4.5, 8],
                         "f6": [3.6, 8],
                         "f7": [5, 8]}

    def adjust(self):
        for figname in self.figure_names:
            if figname in ["f2", "f5"]:
                self.adjust_maps[figname]["right"] = 0.97
                self.adjust_maps[figname]["left"] = 0.23
            elif figname in ["f3", "f6"]:
                self.adjust_maps[figname]["right"] = 0.94
                self.adjust_maps[figname]["left"] = 0.03
            else:
                self.adjust_maps[figname]["right"] = 0.68
                self.adjust_maps[figname]["left"] = 0.03

            self.adjust_maps[figname]["top"] = 0.93
            self.adjust_maps[figname]["bottom"] = 0.13
            self.adjust_maps[figname]["hspace"] = 0.10

    plotOnly = ["f2", "f3", "f4"]
    #plotOnly = "f3"

    def plot(self, data):

        omegas = self.get_family_member_data(data, "omegas")
        s0s = self.get_family_member_data(data, "s0s")
        alphas = self.get_family_member_data(data, "alphas")
        Pls = self.get_family_member_data(data, "Pls")
        cmat = self.get_family_member_data(data, "cmat")

        ld = 5.0
        conv = ld*(1-exp(-1./ld))
        F0s = Pls/conv

        X, Y = np.meshgrid(alphas, F0s, indexing='ij')

        null_formatter = FuncFormatter(lambda v, i: "")
        xes = []
        for x in range(1, len(alphas), 2):
            xes.append(x + 0.5)

        print (cmat[0, :] - cmat[1, :]).sum()

        fes = []
        for y in range(0, len(F0s), 3):
            fes.append(y+0.5)

        for io, omega in enumerate(omegas):

            for is0, s0 in enumerate(s0s):

                C = cmat[io, is0, :, :]

                sfig = eval("self.subfigure%d%d" % (is0+1, io+1))
                csfig = eval("self.csubfigure%d%d" % (is0+1, io+1))

                im = sfig.pcolor(C.transpose(), vmin=0, vmax=1, cmap="gist_earth_r")

                d = 0.1
                levs = [-d] + [d*i for i in range(0, int(1/d + 2))]

                csfig.contourf(X, Y, C, vmin=0, vmax=1, cmap="gist_earth_r", levels=levs)

                fail = np.where(C == -1)

                for y, x in zip(*fail):
                    sfig.scatter(x + 0.5, y + 0.5, c='k', marker='x')

                if io == 0:
                    sfig.set_title(r"$\sigma_0 = %.1f$" % s0)
                    csfig.set_title(r"$\sigma_0 = %.1f$" % s0)

                ax = sfig.axes
                ax.set_xticks(xes)
                ax.set_xticklabels([r"$%.1f$" % alphas[np.floor(ai)] for ai in ax.get_xticks()], minor=False)
                ax.set_yticks(fes)

                if is0 == 0:
                    label = r"$F_0/E_bA$"
                    sfig.set_ylabel(label)
                    csfig.set_ylabel(label)

                    ax.set_yticklabels([r"$%.2f$" % F0s[np.floor(fi)] for fi in ax.get_yticks()], minor=False)

                else:
                    sfig.yaxis.set_ticklabels([])
                    csfig.yaxis.set_ticklabels([])

                    if is0 == 2:
                        sfig.yaxis.set_label_position("right")
                        csfig.yaxis.set_label_position("right")
                        sfig.set_ylabel(r"$\Omega = %.1f$" % omega, labelpad=10)
                        csfig.set_ylabel(r"$\Omega = %.1f$" % omega, labelpad=10)

                if io == 0:
                    sfig.axes.xaxis.set_major_formatter(null_formatter)
                    csfig.axes.xaxis.set_major_formatter(null_formatter)
                else:
                    label = r"$E_b/kT$"
                    sfig.set_xlabel(label)
                    csfig.set_xlabel(label)

                sfig.set_xlim(0, len(alphas))
                sfig.set_ylim(0, len(F0s))
        print Pls[-2], F0s[-2]

        pos = [0.85, 0.2, 0.05, 0.7]
        cbar_ax = self.f4.add_axes(pos)
        cbar_ax2 = self.f7.add_axes(pos)

        clabel = r'$\rho_\mathrm{WV}$'
        cbar_ax.set_title(clabel, position=(0, -0.05), horizontalalignment="left", verticalalignment="center", transform=cbar_ax.transAxes)
        cbar_ax2.set_title(clabel)

        self.f4.colorbar(im, cax=cbar_ax)
        self.f7.colorbar(im, cax=cbar_ax2)

class ExtraNeighbor(DCVizPlotter):

    nametag = "^extraneighbor_(.*)\.npy"

    isFamilyMember = True

    minorfigs = ["f", "p", "t"]#, "b", "i", "P"]

    figMap = dict([["%s%d" % (n, i), "%ssubfigure%d1" % (n, i)]
                   for n in minorfigs for i in [0, 1, 2]] + [["fmjau", ["mjau0", "mjau1", "mjau2"]]])

    hugifyFonts = True

    sizes = [[4.5, 5],
             [3.6, 5],
             [5, 5]]

    specific_fig_size = dict([["%s%d" % (n, i), sizes[i]]
                              for n in minorfigs for i in [0, 1, 2]] + [["fmjau", [5, 7]]])

    #specific_fig_size["t2"] = sizes[1]

    def adjust(self):
        for figname in self.figure_names:
            if "mjau" in figname:
                self.adjust_maps[figname]["top"] = 0.98
                self.adjust_maps[figname]["bottom"] = 0.11
                self.adjust_maps[figname]["left"] = 0.17
                self.adjust_maps[figname]["right"] = 0.83
                self.adjust_maps[figname]["hspace"] = 0.17

            else:
                if "0" in figname:
                    self.adjust_maps[figname]["right"] = 0.97
                    self.adjust_maps[figname]["left"] = 0.23
                elif "1" in figname: #or figname == "t2":
                    self.adjust_maps[figname]["right"] = 0.96
                    self.adjust_maps[figname]["left"] = 0.04
                else:
                    self.adjust_maps[figname]["right"] = 0.68
                    self.adjust_maps[figname]["left"] = 0.03

                self.adjust_maps[figname]["top"] = 0.91
                self.adjust_maps[figname]["bottom"] = 0.15

    #plotOnly = "f3"

    def plot(self, data):
        s0s = self.get_family_member_data(data, "s0s")
        alphas = self.get_family_member_data(data, "alphas")
        F0s = self.get_family_member_data(data, "F0s")
        cmat = self.get_family_member_data(data, "cmat")
        pmat = self.get_family_member_data(data, "pmat")
        tcov = self.get_family_member_data(data, "typecov")
        tmat = np.zeros_like(cmat)
        avg_pmat = np.zeros_like(tmat)
        bmat = np.zeros_like(tmat)
        imat = np.zeros_like(tmat)
        Pmat = np.zeros_like(tmat)

        nbins = 10
        hist = np.zeros(shape=(nbins, 3, 4))
        # counts = np.zeros(shape=(3, 4))
        dx = 1./nbins


        K, L, M, _ = pmat.shape

        for k in range(K):
            for l, alpha in enumerate(alphas):
                for m in range(M):
                    t = np.argmax(pmat[k, l, m, :])
                    tmat[k, l, m] = t

                    p = 0.0
                    c = 0.0
                    for i in range(4):
                        p += pmat[k, l, m, i]*i
                        c += pmat[k, l, m, i]
                    avg_pmat[k, l, m] = p/c
                    imat[k, l, m] = pmat[k, l, m, 1]/c
                    bmat[k, l, m] = pmat[k, l, m, 2]/c
                    Pmat[k, l, m] = pmat[k, l, m, 3]/c
                    #
                    # asd = min(pmat[k, l, m, 1:])
                    # if asd > 2:
                    #     avg_pmat[k, l, m] = 0
                    #     print asd


        at = np.where(alphas >= 1.)[0][0]
        for t, ia, cov in tcov:
            if t == 0:
                continue

            idx = int(cov/dx)
            idx2 = 0 if ia <= at else 1

            hist[idx, idx2, 0] += 1
            hist[idx, idx2, t] += 1

            hist[idx, 2, 0] += 1
            hist[idx, 2, t] += 1

            # counts[idx2, 0] += 1
            # counts[idx2, t] += 1
            #
            # counts[2, 0] += 1
            # counts[2, t] += 1


        labels = [r"$\mathrm{I}$", r"$\mathrm{B}$", r"$\mathrm{P}$"]
        twinlabels = [r"$E_b/kT\in (%g, %g]$" % (0, alphas[at]),
                      r"$E_b/kT\in (%g, %g]$" % (alphas[at], alphas[-1]),
                      r"$E_b/kT\in (%g, %g]$" % (0, alphas[-1])]

        def f(v, i):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return r"$%.1f$" % v
        def g(v, i):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return r"$%.2f$" % v

        markers = ["s", "o", "d"]
        colors = ["r", "k", "0.5"]
        ymaxmjau = 1.05
        sfigs = [self.mjau0, self.mjau1, self.mjau2]
        for i in range(3):
            # hist[:, i, 0] /= counts[i, 0]

            sfigs[i].set_ylabel(r"$n/n_\mathrm{tot}$")
            # sfigs[i].set_yticklabels([])
            sfigs[i].set_ylim(0, ymaxmjau)
            sfigs[i].yaxis.set_major_formatter(FuncFormatter(f))
            P_tot = 0
            for j in range(1,4):
                j = 4 - j
                # hist[:, i, j] /= counts[i, j]
                P = hist[:, i, j]/hist[:, i, 0]
                P_tot += P
                sfigs[i].plot([dx*(k+0.5) for k in range(nbins)],
                              P,
                              "-" + markers[j-1],
                              label=labels[j-1],
                              color=colors[j-1],
                              **my_props["fmt"])

            sfigs[i].plot([0.25, 0.25], [0, ymaxmjau], "k:")
            sfigs[i].plot([0.75, 0.75], [0, ymaxmjau], "k:")

            ax = sfigs[i].axes.twinx()
            ax.set_ylabel(twinlabels[i], labelpad=15, fontsize=20)
            ax.yaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])

        self.mjau0.set_xticklabels([])
        self.mjau1.set_xticklabels([])
        self.mjau2.xaxis.set_major_formatter(FuncFormatter(g))

        tix = [0, 0.25, 0.5, 0.75]
        self.mjau0.set_xticks(tix)
        self.mjau1.set_xticks(tix)
        self.mjau2.set_xticks(tix)
        self.mjau2.set_xlabel(r"$\rho_\mathrm{WV}$")
        lg = self.mjau0.legend(loc="center",
                               numpoints=1,
                               ncol=1,
                               handlelength=0.9,
                               borderpad=0.2,
                               labelspacing=0.2,
                               columnspacing=0.3,
                               handletextpad=0.25,
                               borderaxespad=0.0,
                               frameon=False,
                               bbox_to_anchor=(0.14, 0.45))

        lg.get_frame().set_fill(not (self.toFile and self.transparent))

        #self.mjau.plot([dx*i for i in range(nbins)], hist[:, 0], "r-s", label=r"$\alpha \le %g$" % alphas[at])
        #self.mjau.plot([dx*i for i in range(nbins)], hist[:, 1], "k-s", label=r"$\alpha > %g$" % alphas[at])
        #self.mjau.plot([dx*i for i in range(nbins)], hist[:, 2], "r--")
        #self.mjau.plot([dx*i for i in range(nbins)], hist[:, 2], "ks", label=r"all" % alphas[at])
        #self.mjau.legend(loc="upper right")
        #
        # high_t_islands_cov = cmat[:, :at, :][np.where(tmat[:, :at, :] == 1)]
        # high_t_bands_cov = cmat[:, :at, :][np.where(tmat[:, :at, :] == 2)]
        # high_t_pits_cov = cmat[:, :at, :][np.where(tmat[:, :at, :] == 3)]
        #
        # low_t_islands_cov = cmat[:, at:, :][np.where(tmat[:, at:, :] == 1)]
        # low_t_bands_cov = cmat[:, at:, :][np.where(tmat[:, at:, :] == 2)]
        # low_t_pits_cov = cmat[:, at:, :][np.where(tmat[:, at:, :] == 3)]
        #
        # high_lower_limit = 0.5*(high_t_islands_cov.max() + high_t_bands_cov.min())
        # high_upper_limit = 0.5*(high_t_pits_cov.min() + high_t_bands_cov.max())
        #
        # low_lower_limit = 0.5*(low_t_islands_cov.max() + low_t_bands_cov.min())
        # low_upper_limit = 0.5*(low_t_pits_cov.min() + low_t_bands_cov.max())
        #
        # print "high: [%g,%g] low: [%g, %g]" % (high_lower_limit,
        #                                        high_upper_limit,
        #                                        low_lower_limit,
        #                                        low_upper_limit)
        #
        # sliz3 = cmat[:, :at, :][np.where(tmat[:, at:, :] == 3)]
        # print sliz3.min()
        #
        # t_locs = [[], [], []]
        #
        # for is0 in range(len(s0s)):
        #     for ia in range(len(alphas)):
        #         last_island = None
        #         first_pit = None
        #
        #         for iF0 in range(len(F0s)):
        #             t = tmat[is0, ia, iF0]
        #
        #             if t == 3 and first_pit is None:
        #                 first_pit = iF0
        #             elif t == 1 or t == 0:
        #                 last_island = iF0
        #
        #         if last_island is not None and first_pit is not None:
        #             t_locs[is0].append(0.5*(last_island + first_pit + 1))
        #         else:
        #             t_locs[is0].append(None)
        #
        # print t_locs

        desc = [r"$\mathrm{Pits}$",
                r"$\mathrm{Bands}$",
                r"$\mathrm{Islands}$",
                r"$\mathrm{No\,\,shape}$"]

        x_loc = len(alphas)-9
        y_locs = [8, 6, 4, 2]

        for i, y_loc in enumerate(y_locs):
            avg_pmat[-1, x_loc, y_loc] = 3 - i
            tmat[-1, x_loc, y_loc] = 3 - i

        xes = []
        for x in range(0, len(alphas), 3):
            xes.append(x + 0.5)

        fes = []
        for y in range(0, len(F0s), 3):
            fes.append(y+0.5)

        mats = [cmat, avg_pmat, tmat, bmat, imat, Pmat]

        ims = [None for _ in mats]
        for is0, s0 in enumerate(s0s):
            for i, id in enumerate(self.minorfigs):
                M = mats[i][is0, :, :]

                sfig = eval("self.%ssubfigure%d1" % (id, is0))

                sfig.set_title(r"$\sigma_0 = %.1f$" % s0)

                ims[i] = sfig.pcolor(M.transpose(), vmin=0, vmax=mats[i].max()/cmat.max(), cmap="gist_earth_r")


                #
                # fail = np.where(M == -1)
                #
                # for y, x in zip(*fail):
                #     sfig.scatter(x + 0.5, y + 0.5, c='k', marker='x')

                ax = sfig.axes
                ax.set_xticks(xes)
                ax.set_xticklabels([r"$%.1f$" % alphas[np.floor(ai)] for ai in ax.get_xticks()], minor=False)
                ax.set_yticks(fes)

                if is0 == 0:
                    label = r"$F_0/E_bA$"
                    sfig.set_ylabel(label)

                    ax.set_yticklabels([r"$%.2f$" % F0s[np.floor(fi)] for fi in ax.get_yticks()], minor=False)

                else:
                    sfig.yaxis.set_ticklabels([])

                label = r"$E_b/kT$"
                sfig.set_xlabel(label)

                sfig.set_xlim(0, len(alphas))
                sfig.set_ylim(0, len(F0s))


        # for is0, t_loc in enumerate(t_locs):
        #     if t_loc is None:
        #         continue
        #
        #     sfig = eval("self.psubfigure%d1" % (is0+1))
        #     sfig.plot(np.arange(len(alphas)) + 0.5, t_loc, 'w', linewidth=3)
        #     sfig.plot(np.arange(len(alphas)) + 0.5, t_loc, 'k--', linewidth=3)

        def stuffify(sfig):
            x0 = x_loc
            x1 = x_loc + 1
            args = ['k-']

            for i, y_loc in enumerate(y_locs):
                sfig.text(x_loc + 1.5, y_loc, desc[i])#, horizontalalignment="left", verticalalignment="center")

                y0 = y_loc
                y1 = y_loc + 1

                sfig.plot([x0, x0], [y0, y1], *args)
                sfig.plot([x1, x1], [y0, y1], *args)

                sfig.plot([x0, x1], [y0, y0], *args)
                sfig.plot([x0, x1], [y1, y1], *args)


        stuffify(self.psubfigure21)
        stuffify(self.tsubfigure21)


        pos = [0.78, 0.18, 0.05, 0.7]
        clabel = r'$\rho_\mathrm{WV}$'

        cbar_ax = self.f2.add_axes(pos)
        cbar_ax.set_title(clabel, position=(0, -0.09), horizontalalignment="left", verticalalignment="center", transform=cbar_ax.transAxes)

        ax = self.f2.colorbar(ims[0], cax=cbar_ax)

        pbar_ax = self.p2.add_axes(pos)
        ax = self.p2.colorbar(ims[1], cax=pbar_ax)
        ax.set_ticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
        ax.set_ticklabels([r"", "", r"$\mathrm{I}$", "", r"$\mathrm{B}$", "", r"$\mathrm{P}$"])


class GFPlots(DCVizPlotter):

    nametag = "GPlots.txt"

    figMap = {"gf_figure": ["G_fig","F_fig"],
              "gc_figure": "GC_fig",
              "allf_figure": ["all_F_fig_low", "all_F_fig_high"]}

    hugifyFonts = True

    specific_fig_size =  {"gf_figure": [5,7],
                          "gc_figure": [5,3],
                          "allf_figure": [5,7]}

    stack = "V"

    def adjust(self):
        self.adjust_maps["gf_figure"]["right"] = 0.86
        self.adjust_maps["gf_figure"]["left"] = 0.14
        self.adjust_maps["gf_figure"]["top"] = 0.89
        self.adjust_maps["gf_figure"]["bottom"] = 0.11
        self.adjust_maps["gf_figure"]["hspace"] = 0.13

        self.adjust_maps["gc_figure"]["right"] = 0.84
        self.adjust_maps["gc_figure"]["left"] = 0.16
        self.adjust_maps["gc_figure"]["top"] = 0.96
        self.adjust_maps["gc_figure"]["bottom"] = 0.26

        self.adjust_maps["allf_figure"]["right"] = 0.86
        self.adjust_maps["allf_figure"]["left"] = 0.14
        self.adjust_maps["allf_figure"]["top"] = 0.96
        self.adjust_maps["allf_figure"]["bottom"] = 0.11
        self.adjust_maps["allf_figure"]["hspace"] = 0.13

    def f_attraction_corr(self, di):

        f = 0
        if di >= 2:
            f += 0
        else:
            f += -(128./di**7 - 1.)/20

        return f

    def attraction_orig(self, di):

        g = 0
        g += -1./di**6

        return g

    def attraction_corr(self, di):

        g = 0
        if di >= 2:
            g += 0
        else:
            g += -(3*di + 64./di**6 - 7)/60.

        return g

    def get_repulz(self, s0, hlvec, ld, f0_over_EbA):

        xi = 1-np.exp(-1./ld)
        Z0_over_A = s0/xi

        g_rep = Z0_over_A*np.exp(-(hlvec - 1)/ld)
        g_rep[0] = 0

        H = 1 + ld*np.log(Z0_over_A/f0_over_EbA/ld)

        return H, g_rep

    def plot(self, data):

        s0 = 1.0
        ld = 5.0

        hlvec = np.linspace(1, 9, 1000)
        f0_over_EbA = 0.5

        g_f0 = (hlvec-1)*f0_over_EbA
        H, g_rep = self.get_repulz(s0, hlvec, ld, f0_over_EbA)
        f_rep = g_rep/ld

        g_attr = []
        g_attr_6 = []
        f_attr = []
        for hl in hlvec:
            g_attr.append(self.attraction_corr(hl))
            g_attr_6.append(self.attraction_orig(hl))
            f_attr.append(self.f_attraction_corr(hl))

        g_attr = np.array(g_attr)
        g_attr_6 = np.array(g_attr_6)

        """
        Attractive vs. 1/r^6
        """

        self.GC_fig.plot(hlvec, g_attr_6, "k--", label=r"$-1/d_i^6$")
        self.GC_fig.plot(hlvec, g_attr, "r-", label=r"$\tilde G_\mathrm{WV}$")
        self.GC_fig.set_xlabel(r"$d_i$")
        self.GC_fig.set_ylabel(r"$G(d_i)/E_b$", labelpad=-10)
        self.GC_fig.set_yticks([-1, -0.5, 0])
        self.GC_fig.plot([2, 2], [-1.05, 0.2], 'k:')
        self.GC_fig.set_xlim(0.9, 3)
        self.GC_fig.set_ylim(-1.05, 0.1)

        def gc_f(v, i):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return ""

        self.GC_fig.xaxis.set_major_formatter(FuncFormatter(gc_f))
        self.GC_fig.yaxis.set_major_formatter(FuncFormatter(gc_f))


        lg = self.GC_fig.legend(loc="center",
                               numpoints=1,
                               ncol=1,
                               handlelength=0.9,
                               borderpad=0.2,
                               labelspacing=0.2,
                               columnspacing=0.3,
                               handletextpad=0.25,
                               borderaxespad=0.0,
                               frameon=False,
                               bbox_to_anchor=(0.775, 0.65))

        lg.get_frame().set_fill(not (self.toFile and self.transparent))

        """
        Free energies
        """

        ymin_g = -1.5
        ymax_g = 5.5

        self.G_fig.plot(hlvec, g_f0, 'k--', label=r"$d\cdot F_0$", linewidth=1)
        self.G_fig.plot(hlvec, g_rep, 'g-.', label=r"$G_\lambda$", linewidth=2)
        self.G_fig.plot(hlvec, g_attr, 'r:', label=r"$\tilde G_\mathrm{WV}$", linewidth=2)
        self.G_fig.plot(hlvec, g_rep + g_attr + g_f0, 'k-',
                        label=r"$G_\mathrm{tot}$", linewidth=1)
        self.G_fig.plot([H, H], [ymin_g, ymax_g], "k:")

      #  self.G_fig.set_xlabel(r"$d$")
        self.G_fig.set_xticklabels([])
        self.G_fig.set_ylabel(r"$G(d)/E_bA$", labelpad=-10)

        # for Z0 in [0.5, 1, 1.5]:
        #     for f0 in [0.25, 0.5]:
        #         self.all_fig.plot(di, self.g_repulsion(di, Z0, ld) + self.attraction_corr(di) + (di-1)*f0)
        #
        # self.all_fig.set_xbound(1-0.1)
        # self.all_fig.set_ybound(-1.1)

        lg = self.G_fig.legend(loc="center",
                               numpoints=1,
                               ncol=4,
                               handlelength=0.9,
                               borderpad=0.3,
                               labelspacing=0.2,
                               columnspacing=0.3,
                               handletextpad=0.25,
                               borderaxespad=0.0,
                               frameon=True,
                               bbox_to_anchor=(0.5, 1.15))

        lg.get_frame().set_fill(not (self.toFile and self.transparent))

        self.G_fig.set_xbound(1-0.25)
        self.G_fig.set_ylim(ymin_g, ymax_g)


        """
        Forces
        """


        tot = (f_rep + f_attr - f0_over_EbA)/f0_over_EbA
        ymin = -1.1
        ymax = 1.0

        self.F_fig.plot(hlvec, tot, "r-",
                            label=r"$F_0 + F_l + F_\lambda$", linewidth=3, linestyle="-", zorder=1)
        self.F_fig.plot([H, H], [ymin, ymax], "k:", zorder=2)
        #self.F_fig.plot([di[0], di[-1]], [f0, f0], styles[1],
        #                    label=r"$F_0$", linewidth=2)
        self.F_fig.plot([0.75, hlvec[-1]], [0, 0], "k-", zorder=0)
        #self.F_fig.plot([r0], [0], 'go', markersize=10, label=r"$r_{(\mathrm{near/far})}$")
        #self.F_fig.plot([r1], [0], 'go', markersize=10)
        #self.F_fig.plot([1], [0], 'go', markersize=10)
        #self.F_fig.scatter(1, 0, s=30, c="r", edgecolors='none', zorder=2)

        self.F_fig.set_xlabel(r"$d$")
        self.F_fig.set_ylabel(r"$F_\mathrm{tot}(d)/F_0$", labelpad=-10)
        self.F_fig.set_yticks([-1, 0, 1])

        self.F_fig.set_ylim(ymin, ymax)
        #self.F_fig.set_xlim(0.9, di[-1])

        #self.F_fig.xaxis.set_ticklabels([r"$h_c$", r"$h_\lambda$"])
        #self.F_fig.axes.tick_params(axis='x', which='major', pad=10)

        self.F_fig.set_xlim(1-0.25, hlvec[-1])
      #
        def f(v, i):
            if int(v) == v:
                return "$%d$" % v
            return r"$%.1f$" % v

        self.F_fig.yaxis.set_major_formatter(FuncFormatter(f))

        x0, x1 = self.F_fig.get_xlim()
        y0, y1 = self.F_fig.get_ylim()
        self.F_fig.text((3*x1+x0)/4, 3*y0/4.,  r"$\mathrm{Attraction}$", verticalalignment="center", horizontalalignment="center", fontsize=20)
        self.F_fig.text((3*x1+x0)/4, 3*y1/4., r"$\mathrm{Repulsion}$", verticalalignment="center", horizontalalignment="center", fontsize=20)
        self.F_fig.text(H + 0.25, 0.01, r"$h_\lambda \sim %.2f$" % H, verticalalignment="bottom", horizontalalignment="left", fontsize=20)

        F0s = [0.5, 1.0]
        s0s = [0.5, 1.0, 1.5]
        styles = ["r-", "k--", "g-."]
        widths = [1, 1, 2]

        for i, f0_over_EbA in enumerate(F0s):
            all_F_fig = [self.all_F_fig_low, self.all_F_fig_high][i]
            for j, s0 in enumerate(s0s):
                H, g_rep = self.get_repulz(s0, hlvec, ld, f0_over_EbA)
                f_rep = g_rep/ld

                tot = (f_rep + f_attr - f0_over_EbA)/f0_over_EbA
                all_F_fig.plot(hlvec, tot, styles[j], label=r"$\sigma_0=%.1f$" % s0, linewidth=widths[j])

            all_F_fig.set_ylim(ymin, 2.0)
            all_F_fig.plot([1, hlvec[-1]], [0, 0], "k-", zorder=0)
            all_F_fig.set_yticks([-1, 0, 1, 2])
            all_F_fig.set_ylabel(r"$F_\mathrm{tot}(d)/F_0$", labelpad=-10)

            ax = all_F_fig.axes.twinx()
            ax.set_ylabel(r"$F_0/E_bA=%.1f$" % f0_over_EbA, labelpad=10)
            ax.yaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])

        self.all_F_fig_high.set_xlabel(r"$d$")

        self.all_F_fig_low.yaxis.set_major_formatter(FuncFormatter(f))
        self.all_F_fig_high.yaxis.set_major_formatter(FuncFormatter(f))
        self.all_F_fig_low.set_xticklabels([])


        lg = self.all_F_fig_high.legend(loc="center",
                                        numpoints=1,
                                        ncol=1,
                                        handlelength=0.9,
                                        borderpad=0.2,
                                        labelspacing=0.2,
                                        columnspacing=0.3,
                                        handletextpad=0.25,
                                        borderaxespad=0.0,
                                        frameon=False,
                                        bbox_to_anchor=(0.725, 0.7))

        lg.get_frame().set_fill(not (self.toFile and self.transparent))


class FPlots(DCVizPlotter):

    nametag = "FPlots.txt"

    hugifyFonts = True

    fig_size = [7, 5]

    def f_repulsion(self, di, s0, ld):
        return s0*np.exp(-(di-1)/ld)/ld*(di > 1)

    @staticmethod
    def f_attraction_corr_serial(di):

        if di >= 2:
            return 0

        return -(128./di**7 - 1.)/20

    f = None
    def f_attraction_corr(self, di):

        if not self.f:
            self.f = np.vectorize(self.f_attraction_corr_serial)

        return self.f(di)

    def adjust(self):
        figname = "figure"

        self.adjust_maps[figname]["right"] = 0.85
        self.adjust_maps[figname]["left"] = 0.15
        self.adjust_maps[figname]["top"] = 0.96
        self.adjust_maps[figname]["bottom"] = 0.17

    def plot(self, data):

        di = np.linspace(1, 6, 1000)
        ld = 1.0
        s0 = 2.0

        r0 = 1.29714
        r1 = 2.38629

        f_rep = self.f_repulsion(di, s0, ld)
        f_attr = self.f_attraction_corr(di)

        styles = ['g-.', 'k--', 'r-', "0.0"]

        f0 = 0.5
        tot = f_rep + f_attr - f0
        ymin = -0.1

        self.subfigure.plot(di, tot, styles[2],
                            label=r"$F_0 + F_l + F_\lambda$", linewidth=3, linestyle="-", zorder=1)
        self.subfigure.plot([r1, r1], [-2*f0 + ymin, tot.max() - ymin], "k--", zorder=2)
        #self.subfigure.plot([di[0], di[-1]], [f0, f0], styles[1],
        #                    label=r"$F_0$", linewidth=2)
        self.subfigure.plot([0, di[-1]], [0, 0], "k-", zorder=0)
        #self.subfigure.plot([r0], [0], 'go', markersize=10, label=r"$r_{(\mathrm{near/far})}$")
        #self.subfigure.plot([r1], [0], 'go', markersize=10)
        #self.subfigure.plot([1], [0], 'go', markersize=10)
        self.subfigure.scatter(1, 0, s=50, c="r", edgecolors='none', zorder=2)




        """
        lg = self.subfigure.legend(loc="center",
                                   numpoints=1,
                                   ncol=1,
                                   handlelength=0.9,
                                   borderpad=0.2,
                                   labelspacing=0.2,
                                   columnspacing=0.3,
                                   handletextpad=0.25,
                                   borderaxespad=0.0,
                                   frameon=False,
                                   bbox_to_anchor=(0.725, 0.55))

        lg.get_frame().set_fill(not (self.toFile and self.transparent))
        """

        self.subfigure.set_xlabel(r"$h_l$", labelpad=10)
        self.subfigure.set_ylabel(r"$F_\mathrm{tot}(h_l)$")

        self.subfigure.set_ylim(-2*f0  + ymin, tot.max() - ymin)
        self.subfigure.set_xlim(0.9, di[-1])

        def f(v, i):
            if v == -f0:
                return r"$-\mathrm{F_0}$"
            elif int(v) == v:
                return r"$%d$" % v

            return r"$%.1f$" % v

        self.subfigure.yaxis.set_major_formatter(FuncFormatter(f))
        self.subfigure.xaxis.set_ticks([1, r1])
        self.subfigure.xaxis.set_ticklabels([r"$h_c$", r"$h_\lambda$"])
        self.subfigure.yaxis.set_ticks([-f0, 0])
        self.subfigure.axes.tick_params(axis='x', which='major', pad=10)

        x0, x1 = self.subfigure.get_xlim()
        y0, y1 = self.subfigure.get_ylim()
        self.subfigure.text(4*x1/5, -y0/5,  r"$\mathrm{Repulsion}$", verticalalignment="center", horizontalalignment="center", fontsize=30)
        self.subfigure.text(4*x1/5, +y0/5, r"$\mathrm{Attraction}$", verticalalignment="center", horizontalalignment="center", fontsize=30)




class Extraneighbor_cluster(DCVizPlotter):

    nametag = "extran_cluster_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True
    labelSize = 35
    ticklabelSize = 25

    figMap = {"figure3D" : [],
              "figure0": ["subfigure0", "subfigure4"],
              "figure1": "subfigure1",
              "figure2": "subfigure5"}

    plotOnly = ["figure3D", "figure0"]

    tight=False

    specific_fig_size = {"figure3D": [7, 5],
                         "figure0": [7, 8]}

    def make3Dplot(self, ax, L, W, tmax, time, trace):
        max_size = 20
        min_size = 5

        proj_props = {
            "s" : max_size,
            "marker" : "x",
            "c" : "k",
            "linewidth" : 1
        }

        colors = ['r', 'g', 'b', 'k', 'm', 'y']

        for xi, yi, zi, ri, ci in trace:
            si = min_size + ri*(max_size-min_size)

            ax.scatter(time[zi], xi, yi, s=si, c=colors[int(ci)%len(colors)], edgecolors='none')


            if zi % 20 == 0:
                #ax.scatter(len(data), xi, yi, s=sproj, c='k', edgecolors='none')
                # ax.scatter(zi, x1, yi, **proj_props)
                ax.scatter(time[zi], xi, 0, **proj_props)

        #
        # def orthogonal_proj(zfront, zback):
        #     a = (zfront+zback)/(zfront-zback)
        #     b = -2*(zfront*zback)/(zfront-zback)
        #     return numpy.array([[1,0,0,0],
        #                         [0,1,0,0],
        #                         [0,0,a,b],
        #                         [0,0,0,zback]])
        #
        # proj3d.persp_transformation = orthogonal_proj

        ax.set_xlim(0, tmax)
        ax.set_zlim(0, L)
        ax.set_zlim(0, W)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])

        pad = 1
        ax.xaxis._axinfo['label']['space_factor'] = pad
        ax.yaxis._axinfo['label']['space_factor'] = pad
        ax.zaxis._axinfo['label']['space_factor'] = pad

        fs = self.labelSize
        ax.set_xlabel(r"$\nu t$", size=fs)
        ax.set_ylabel(r"$x$", size=fs)
        ax.set_zlabel(r"$y$", size=fs)

    def adjust(self):

        self.do_legend = True
        if "no_legend" in self.argv:
            self.argv.remove("no_legend")
            self.do_legend=False

        delta = 0.06
        print self.do_legend

        self.adjust_maps["figure0"]["left"] = 0.15 - delta*(not self.do_legend)
        self.adjust_maps["figure0"]["right"] = 0.85 + delta*self.do_legend
        self.adjust_maps["figure0"]["top"] = 0.93
        self.adjust_maps["figure0"]["bottom"] = 0.12
        self.adjust_maps["figure0"]["hspace"] = 0.15


    def plot(self, data):

        time = self.get_family_member_data(data, "time")
        covs = self.get_family_member_data(data, "covs")

        nbroken = self.get_family_member_data(data, "nbroken")
        ngained = self.get_family_member_data(data, "ngained")
        cumnbroken = np.cumsum(nbroken)
        cumngained = np.cumsum(ngained)

        if "force_zero" in self.argv:
            self.argv.remove("force_zero")
            idx = np.argmin(covs)
            covs[idx] = 0

        F0 = self.argv[0]
        self.subfigure0.set_title(r"$F_0/E_bA = %s$" % F0)

        L, W = [int(x) for x in self.argv[1:3]]
        A_corr = 1

        help_lines = []
        if len(self.argv) > 3:
            tmax = int(self.argv[3])

            if tmax > time[-1]:
                tmax = time[-1]

            if len(self.argv) > 4:
                for x in self.argv[4:]:
                    help_lines.append([float(x), float(x)])
        else:
            tmax = time[-1]

        def y(*args):
            _min = min([min(v[np.where(time < tmax)]) for v in args])
            _max = max([max(v[np.where(time < tmax)]) for v in args])

            return [_min, _max]

        try:
            n = self.get_family_member_data(data, "n")
            self.subfigure1.plot(time, n)
            y1 = y(n)

            for x in help_lines:
                self.subfigure1.plot(x, y1, "k--")

            self.subfigure1.set_ylim(y1)
        except:
            pass

        try:
            circs = self.get_family_member_data(data, "circs")
            size = self.get_family_member_data(data, "size")
            sphericity = np.sqrt(4*np.pi*(size+circs/2+np.pi/4))/circs
            box = np.sqrt(np.pi)/2.
            self.subfigure5.plot(time, sphericity)
            self.subfigure5.plot([time[0], time[-1]], [box, box], "k--", linewidth=3)
            y5 = y(sphericity)
            for x in help_lines:
                self.subfigure5.plot(x, y5, "k--")
            self.subfigure5.set_ylim(y5)

        except:
            pass

        try:
            trace = self.get_family_member_data(data, "trace")
            ax3d = Axes3D(self.figure3D)
            ax3d.view_init(elev=20., azim=-65)
            self.make3Dplot(ax3d, L, W, tmax, time, trace)
        except:
            pass

        self.subfigure0.plot(time, covs)
        self.subfigure4.plot(time[:-1], cumngained/A_corr, "g", linewidth=2, label=r"$n_\mathrm{Formed}$")
        self.subfigure4.plot(time[:-1], cumnbroken/A_corr, "r--", linewidth=2, label=r"$n_\mathrm{Broken}$")

        y0 = [0, np.ceil(10*max(covs[np.where(time<tmax)]))/10.]
        y4 = y(cumnbroken/A_corr, cumngained/A_corr)

        for x in help_lines:
            self.subfigure0.plot(x, y0, "k--")
            self.subfigure4.plot(x, y4, "k--")

        #self.subfigure0.set_xlabel(r"$\nu t$")
        self.subfigure1.set_xlabel(r"$\nu t$")
        self.subfigure4.set_xlabel(r"$\nu t$")
        self.subfigure5.set_xlabel(r"$\nu t$")

        if self.do_legend:
            self.subfigure0.set_ylabel(r"$\rho_\mathrm{WV}$")
            self.subfigure4.set_ylabel(r"$n/A$")
        self.subfigure1.set_ylabel("n clusters")
        self.subfigure5.set_ylabel("sphericity")

        self.subfigure0.set_ylim(y0)
        self.subfigure4.set_ylim(y4)

        self.subfigure0.set_xlim(0, tmax)
        self.subfigure1.set_xlim(0, tmax)
        self.subfigure4.set_xlim(0, tmax)
        self.subfigure5.set_xlim(0, tmax)


        self.subfigure0.set_xticks([])

        if self.do_legend:
            self.subfigure4.legend(loc="center",
                                   handlelength=1.5,
                                   markerscale=20.0,
                                   borderpad=0.2,
                                   labelspacing=0.2,
                                   columnspacing=1.0,
                                   handletextpad=0.2,
                                   borderaxespad=0.0,
                                   frameon=False,
                                   fontsize=self.labelSize,
                                   bbox_to_anchor=(0.75, 0.25))

        def format0(v, _):
            if abs(v - round(v, 1)) < 1e-3:
                return r"$%.1f$" % v
            else:
                return r""

        self.subfigure0.yaxis.set_major_formatter(FuncFormatter(format0))

        def f_x(v, _):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return ""

        self.subfigure4.yaxis.set_major_formatter(FuncFormatter(f_x))

class ExtraneighborTest(DCVizPlotter):

    nametag = "test_extraneighbor_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True

    def plot(self, data):

        Pls = self.get_family_member_data(data, "Pls")
        covs = self.get_family_member_data(data, "covs")

        self.subfigure.plot(Pls, covs, "ks")

        self.subfigure.set_xlabel("Pl")
        self.subfigure.set_ylabel("cov")


def RHS(ld, s0, F0EbA):
    xi = (1 - np.exp(-1./ld))

    return ld*np.log(s0/(ld*xi*F0EbA))

def find_integer_crossings(f):
    n0 = int(floor(f[0]))
    n1 = int(ceil(f[-1]))

    ints = range(n1, n0 + 1)
    idx = np.zeros(len(ints))

    i = 0
    for j, fi in enumerate(f):
        if fi < ints[len(ints) - i - 1]:
            idx[i] = j
            i += 1

            if i == len(ints):
                break

    return idx[::-1]

def resonance_points(s0, ld, F0, F1, n=1000):
    Fs = np.linspace(F0, F1, n)
    return [Fs[x] for x in find_integer_crossings(RHS(ld, s0, Fs))]



class ExtraneighborResonance(DCVizPlotter):

    nametag = "resonance_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True
    labelSize= 50
    ticklabelSize=40
    fontSize = 40  # Only invoked by hugifyFonts = True
    tickSize = 10

    figMap = {
        "fig": ["sfig0",
                "sfig1"],
        "rfig": ["rfig0",
                 "rfig1"]
    }

    stack = "H"

    specific_fig_size = {"fig": [16, 6],
                         "rfig": [16, 5.5]}

    def adjust(self):
        for fig in ["fig", "rfig"]:
            self.adjust_maps[fig]["left"] = 0.11
            self.adjust_maps[fig]["right"] = 0.975
            self.adjust_maps[fig]["wspace"] = 0.1

        self.adjust_maps["fig"]["top"] = 0.95
        self.adjust_maps["fig"]["bottom"] = 0.22

        self.adjust_maps["rfig"]["top"] = 0.87
        self.adjust_maps["rfig"]["bottom"] = 0.08

        rcParams['ytick.major.size'] = 8
        rcParams['ytick.major.width'] = 1

        rcParams['xtick.major.size'] = 8
        rcParams['xtick.major.width'] = 1

    def plot(self, data):

        s0s = self.get_family_member_data(data, "s0s")

        xmax1 = float(self.argv[0])
        xmaxes = [xmax1, 1000]

        ld = 5.

        for i in [0, 1]:
            subfigure = eval("self.sfig%d" % i)
            rsubfigure = eval("self.rfig%d" % i)

            F0s = self.get_family_member_data(data, "F0s%d" % i)
            cvec = self.get_family_member_data(data, "cvec%d" % i)

            start = np.where(cvec != 0)[0][0]

            F0s = F0s[start:]
            cvec = cvec[start:]

            F0EbA = np.linspace(F0s[0], F0s[-1], 1000)

            rhs = RHS(ld, s0s[i], F0EbA)

            n0 = int(floor(rhs[0]))
            n1 = int(ceil(rhs[-1]))

            idx = find_integer_crossings(rhs)

            txtshift = 0.01
            for j, n in enumerate(range(n1, n0 + 1)):
                p = F0EbA[idx[j]]

                if p > xmaxes[i]:
                    continue

                rsubfigure.plot([F0EbA[0], p], [n, n], "k--")
                rsubfigure.plot([p, p], [1, n], "k--", linewidth=3)

                textx = p+txtshift
                texty = n

                if i == 1:
                    texty += 0.1

                rsubfigure.text(textx, texty, r"$%.2f$" % p)

            rsubfigure.plot(F0EbA, rhs, "r-", linewidth=3)


            rsubfigure.set_xlim(F0s[0], F0s[-1])
            rsubfigure.set_ylim(1, 6)

            def intformatter(v, i):
                if int(v) == v:
                    return r"$%d$" % v
                else:
                    return ""

            rsubfigure.yaxis.set_major_formatter(FuncFormatter(intformatter))

            subfigure.plot(F0s, cvec, "r-", linewidth=3)

            subfigure.set_xlabel(r"$F_0/E_bA$")

            resonances = resonance_points(s0s[i], 5.0, F0s.min(), F0s.max())
            for res in resonances:
                subfigure.plot([res, res], [0, 1], "k--", linewidth=3)

            if i == 1:
                f0max = F0s.max()
            else:
                f0max = 0.5

            subfigure.set_ylim(0, 1)
            subfigure.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            subfigure.set_xlim(F0s.min(), f0max)
            subfigure.set_xticks(np.arange(round(10*F0s.min())/10.,
                                           f0max + 0.001,
                                           0.1))

            rsubfigure.set_title(r"$\sigma_0 = %.2f$" % s0s[i])

            if i == 0:
                subfigure.set_ylabel(r"$\rho_\mathrm{WV}$")
                rsubfigure.set_ylabel(r"$\lambda_D \ln \left(\frac{\sigma_0/\lambda_D\xi }{F_0/E_bA}\right)$")
            else:
                subfigure.set_yticklabels([])
                rsubfigure.set_yticklabels([])

            rsubfigure.set_xticklabels([])

        self.sfig0.set_xlim(self.sfig0.get_xlim()[0], xmax1)
        self.rfig0.set_xlim(self.sfig0.get_xlim()[0], xmax1)




class ResonancePlots(DCVizPlotter):

    nametag = "resonance\.txt"

    hugifyFonts = True

    def plot(self, data):

        F0 = 0.7
        F1 = 1.3
        s0 = 1.5
        ld = 5.

        F0EbA = np.linspace(F0, F1, 1000)
        F0EbAFull = np.linspace(0.1, 1.5, 10000)

        rhs = RHS(ld, s0, F0EbA)
        rhsFull = RHS(ld, s0, F0EbAFull)

        n0 = int(floor(rhs[0]))
        n1 = int(ceil(rhs[-1]))

        print [round(F0EbAFull[x], 2) for x in find_integer_crossings(rhsFull)]

        idx = find_integer_crossings(rhs)

        txtshift = 0.01
        for i, n in enumerate(range(n1, n0 + 1)):
            p = F0EbA[idx[i]]
            self.subfigure.plot([F0EbA[0], p], [n, n], "k--")
            self.subfigure.plot([p, p], [rhs.min(), n], "k--")
            self.subfigure.text(p+txtshift, n, r"$%.2f$" % p)

        self.subfigure.plot(F0EbA, rhs, "r-", linewidth=3)

        self.subfigure.set_xlabel(r"$F_0/E_bA$")
        self.subfigure.set_ylabel(r"$\lambda_D \log \left(\frac{\sigma_0/\lambda_D\xi }{F_0/E_bA}\right)$")

        self.subfigure.set_xlim(F0, F1)
        self.subfigure.set_ylim(rhs.min(), rhs.max())

        def intformatter(v, i):
            if int(v) == v:
                return r"$%d$" % v
            else:
                return ""

        self.subfigure.yaxis.set_major_formatter(FuncFormatter(intformatter))


class NonEqNeigz(DCVizPlotter):

    nametag = "noneq_neigz_(.*)\.npy"

    hugifyFonts = True

    isFamilyMember = True

    figMap = {"figure_high" : ["dissfig_high", "growthfig_high", "zoom_diss"],
              "figure_low": ["dissfig_low", "growthfig_low", "zoom_growth"]}

    def adjust(self):
        for figure in self.figure_names:
            self.adjust_maps[figure]["wspace"] = 0.16
            self.adjust_maps[figure]["left"] = 0.055
            self.adjust_maps[figure]["right"] = 0.965

        self.adjust_maps["figure_high"]["top"] = 0.94
        self.adjust_maps["figure_high"]["bottom"] = 0.25

        self.adjust_maps["figure_low"]["top"] = 0.87
        self.adjust_maps["figure_low"]["bottom"] = 0.12

    stack = "H"

   # fig_size = [15, 3]
    specific_fig_size = {
        "figure_high": [15, 3.25],
        "figure_low": [15, 3]
    }

    def plot(self, data):

        from matplotlib import rcParams
        print rcParams["figure.figsize"]

        logscale=3
        tscale = 10**logscale

        do_zoom = not self.argv
        growth_start = 4.0
        growth_end = 4.5
        diss_start = 2.25

        F0s = self.get_family_member_data(data, "F0s")
        omegas = self.get_family_member_data(data, "omegas")

        figs = [[self.dissfig_low,  self.growthfig_low],
                [self.dissfig_high, self.growthfig_high]]

        ylims = [0.35, 1.0]
        xlims = [5.1, 5.75]

        zylims = [0.35, 0.15]

        xzl = r"$10^%d\,\nu t$" % logscale
        args = ['k-']
        kwargs = {"linewidth": 2, "zorder": 10}

        for i in range(2):
            F0 = F0s[i]
            for j in range(2):
                omega = omegas[j]
                fig = figs[i][j]
                eqcov = self.get_family_member_data(data, "eqcov%d%d" % (i, j))
                neqcov = self.get_family_member_data(data, "neqcov%d%d" % (i, j))

                eqtime = self.get_family_member_data(data, "eqtime%d%d" % (i, j))/tscale
                neqtime = eqtime[-1] + self.get_family_member_data(data, "neqtime%d%d" % (i, j))/tscale

                fig.plot([eqtime[-1], eqtime[-1]], [0, 1], "k--")

                fig.set_ylim(0, ylims[i])
                fig.set_xlim(0, xlims[i])

                fig.yaxis.set_major_formatter(FuncFormatter(lambda v, i: r"$%.1f$" % round(v, 1) if abs(round(v, 1) - v) < 1E-3 else ""))

                #high pressure growth
                if i == 1 and j == 1:
                    fig.plot([neqtime[-1], neqtime[-1]], [neqcov[-1], 1.0], "r-")

                if do_zoom:
                    #high pressure dissolution
                    if i == 1 and j == 0:
                        k = np.where(neqtime > diss_start)[0][0]
                        self.zoom_diss.plot(neqtime[k:], neqcov[k:], 'r-')
                        self.zoom_diss.set_xlim(neqtime[k], xlims[i])
                        self.zoom_diss.set_ylim(0, zylims[i])
                        self.zoom_diss.set_xlabel(xzl)
                        self.zoom_diss.yaxis.set_label_position("right")
                        self.zoom_diss.set_ylabel(r"$F_0/E_bA = %g$" % F0, labelpad=10)

                        x0 = neqtime[k]
                        x1 = neqtime[-1]
                        args = ['k-']

                        self.dissfig_high.plot([x0, x1], [0, 0], *args, **kwargs)
                        self.dissfig_high.plot([x0, x1], [zylims[i], zylims[i]], *args, **kwargs)

                        self.dissfig_high.plot([x0, x0], [0, zylims[i]], *args, **kwargs)
                        self.dissfig_high.plot([x1, x1], [0, zylims[i]], *args, **kwargs)


                    #low pressure growth
                    elif i == 0 and j == 1:
                        ks = np.where(neqtime > growth_start)[0][0]
                        k = np.where(neqtime < growth_end)[0][-1]
                        self.zoom_growth.plot(neqtime[ks:k], neqcov[ks:k], 'r-')
                        self.zoom_growth.set_xlim(neqtime[ks], neqtime[k])
                        self.zoom_growth.set_ylim(0, zylims[i])
                        self.zoom_growth.yaxis.set_label_position("right")
                        self.zoom_growth.set_ylabel(r"$F_0/E_bA = %g$" % F0, labelpad=10)
                        self.zoom_growth.set_yticklabels([])
                        self.zoom_growth.set_title(r"$\mathrm{Zoom}$")

                        x0 = neqtime[ks-1]
                        x1 = neqtime[k-1]

                        self.growthfig_low.plot([x0, x1], [0.001, 0.001], *args, **kwargs)
                        self.growthfig_low.plot([x0, x1], [zylims[i]*0.985, zylims[i]*0.985], *args, **kwargs)

                        self.growthfig_low.plot([x0, x0], [0, zylims[i]], *args, **kwargs)
                        self.growthfig_low.plot([x1, x1], [0, zylims[i]], *args, **kwargs)


                fig.plot(eqtime, eqcov, "k-")
                fig.plot(neqtime, neqcov, "r-")

                if j == 1:
                    fig.set_yticklabels([])

                else:
                    fig.set_ylabel(r"$\rho_\mathrm{WV}$")

                if i == 0:
                    fig.set_title(r"$\Omega = %g$" % omega)
                else:
                    fig.set_xlabel(xzl)
                    # fig.set_xticklabels([])


class ExtraN_cluster_yo(DCVizPlotter):

    nametag = "cov_clusters_(.*)\.npy"

    isFamilyMember = True

    hugifyFonts = True
    labelSize = 35
    ticklabelSize = 25


    figMap = {
        "fLeft":  ["std0", "rho0", "h0"],
        "fMid":   ["std1", "rho1", "h1"],
        "fRight": ["std2", "rho2", "h2"],
        "lonelyfig1": "sfig1",
        "lonelyfig2": "sfig2",
        "hfig": "shfig",
        "testfig": "test"
    }

    specific_fig_size = {
        "fLeft":  [6.95, 8],
        "fMid":   [6, 8],
        "fRight": [6, 8],
        "lonelyfig1": [7, 6.75],
        "lonelyfig2": [7, 6.75]
       # "lonelyfig2": [6, 6.75]
    }

    def adjust(self):

        for figure in self.figure_names:

            if "lonelyfig" in figure:
                self.adjust_maps[figure]["top"] = 0.90
                self.adjust_maps[figure]["bottom"] = 0.15

                if figure == "lonelyfig1":
                    self.adjust_maps[figure]["left"] = 0.16
                    self.adjust_maps[figure]["right"] = 0.98
                else:
                    self.adjust_maps[figure]["left"] = 0.02
                    self.adjust_maps[figure]["right"] = 0.84

            else:
                self.adjust_maps[figure]["hspace"] = 0.15
                self.adjust_maps[figure]["top"] = 0.93
                self.adjust_maps[figure]["right"] = 0.965
                self.adjust_maps[figure]["bottom"] = 0.12

                if figure == "fLeft":
                    lval = 0.16
                else:
                    lval = 0.04

                self.adjust_maps[figure]["left"] = lval

    #plotOnly = "hfig"

    def plot(self, data):

        alphas = self.get_family_member_data(data, "alphas")
        F0s = self.get_family_member_data(data, "F0")
        stds = self.get_family_member_data(data, "stds")
        covs = self.get_family_member_data(data, "covs")
        heights = self.get_family_member_data(data, "heights")

        sorted_alphas = sorted(alphas)

        if self.argv:
            ymax = float(self.argv.pop(0))
            y1max = float(self.argv.pop(0))
        else:
            ymax = stds[:, :, 0].max()*101
            y1max = None

        yhmax = (heights.max()+1)#*covs.max()

        markers = ['s', '^', 'o']
        colors = ['k', 'r', '0.3']
        stdlabel = r"$\sigma(\rho_\mathrm{WV})[\%]$"
        stdnlabel = r"$\sigma(\Delta n/A)[\%]$"

        if self.argv:
            res = resonance_points(1., 5., F0s.min(), F0s.max(), 10000)
        else:
            res = []

        ypr = [0, 1]
        yps = [0, ymax]
        yph = [0, yhmax]
        print res

        for ia, a in enumerate(alphas):
            ias = sorted_alphas.index(a)
            stdfig = eval("self.std%d" % ias)
            rhofig = eval("self.rho%d" % ias)
            hfig = eval("self.h%d" % ias)

            tresh = 0.05
            I = np.where(covs[ias, :] > tresh)
            J = np.where(covs[ias, :] <= tresh)

            stdfig.plot(F0s[I], stds[ia, :, 0][I]*100, 'kd', label=r"$\mathrm{Stable}$", **my_props["fmt"])
            stdfig.plot(F0s[J], stds[ia, :, 0][J]*100, 'rx', label=r"$\mathrm{Unstable}$", **my_props["fmt"])
            rhofig.plot(F0s[I], covs[ia, :][I], 'kd', **my_props["fmt"])
            rhofig.plot(F0s[J], covs[ia, :][J], 'rx', **my_props["fmt"])
            hfig.plot(F0s[I], heights[ia, :][I], 'kd', **my_props["fmt"])
            hfig.plot(F0s[J], heights[ia, :][J], 'rx', **my_props["fmt"])

            self.test.plot(stds[ia, :, 0][I], heights[ia, :][I], markers[ias], color=colors[ias], **my_props['fmt'])

            for r in res:
                xp = [r, r]
                stdfig.plot(xp, yps, "k:")
                rhofig.plot(xp, ypr, "k:")
                hfig.plot(xp, yph, "k:")

            K = np.where(covs[ias, :] != 0)
            for i in [1, 2]:
                sfig = eval("self.sfig%d" % i)
                std_yo = stds[ias, :, 3-i]*100

                sfig.plot(F0s[K], std_yo[K], markers[ias],
                          color=colors[ias],
                          label=r"$E_b/kT=%g$" % alphas[ias],
                          **my_props["fmt"])

            hfig.set_xlabel(r"$F_0/E_b A$")

            stdfig.xaxis.set_ticklabels([])
            rhofig.xaxis.set_ticklabels([])
            stdfig.set_title(r"$E_b/kT = %g$" % a)

            stdfig.set_ylim(0, ymax)
            rhofig.set_ylim(0, 1)
            hfig.set_ylim(0, yhmax)

            stdfig.set_xlim(0.23, 1.02)
            rhofig.set_xlim(0.23, 1.02)
            hfig.set_xlim(0.23, 1.02)

            if ias == 0:
                stdfig.set_ylabel(stdlabel, labelpad=22)
                rhofig.set_ylabel(r"$\rho_\mathrm{WV}$")
                hfig.set_ylabel(r"$\langle H_c \rangle/l_0$", labelpad=16)
            else:
                stdfig.yaxis.set_ticklabels([])
                rhofig.yaxis.set_ticklabels([])
                hfig.yaxis.set_ticklabels([])

                if ias == 2:
                     stdfig.legend(loc="upper left",
                                  numpoints=1,
                                  ncol=1,
                                  handlelength=1.0,
                                  borderpad=0.2,
                                  labelspacing=0.2,
                                  columnspacing=1.25,
                                  handletextpad=0.25,
                                  borderaxespad=0.0,
                                  frameon=True,
                                  fontsize=30,
                                  bbox_to_anchor=(0.04, 0.9)) \
                    .get_frame().set_fill(not (self.toFile and self.transparent))


        if not y1max or y1max == 0:
            y1 = self.sfig1.get_ylim()[1]
            y2 = self.sfig2.get_ylim()[1]
        else:
            y1 = y1max
            y2 = 0
            self.sfig1.set_ylim(0, y1)

        self.sfig1.set_ybound(0)
        self.sfig2.set_ybound(0)

        if y1 > y2:
            self.sfig2.set_ylim(0, y1)
            ym = y1
        else:
            self.sfig1.set_ylim(0, y2)
            ym = y2

        for r in res:
            xp = [r, r]
            self.sfig1.plot(xp, [0, ym], "k:")
            self.sfig2.plot(xp, [0, ym], "k:")

            self.sfig1.set_ylim(0, ym)
            self.sfig2.set_ylim(0, ym)

        self.sfig1.legend(loc="center",
                          numpoints=1,
                          handlelength=1.0,
                          markerscale=1.0,
                          borderpad=0.3,
                          labelspacing=0.2,
                          columnspacing=1.0,
                          handletextpad=0.0,
                          borderaxespad=0.1,
                          frameon=True,
                          fontsize=self.ticklabelSize,
                          bbox_to_anchor=(0.225, 0.835))

        self.sfig1.set_xlabel(r"$F_0/E_bA$")
        self.sfig1.set_xlim(0.23, 1.02)
        self.sfig1.set_ylabel(stdnlabel)

        self.sfig2.set_xlabel(r"$F_0/E_bA$")
        self.sfig2.set_xlim(0.23, 1.02)
        self.sfig2.set_yticklabels([])

        self.sfig1.set_title(r"$\mathrm{Gained\,\,Bonds}\,\,(\Delta n_+)$")
        self.sfig2.set_title(r"$\mathrm{Broken\,\,Bonds}\,\,(\Delta n_-)$")
        self.shfig.set_title(r"$\mathrm{Contact\,\,Height}$")

        self.shfig.set_xlabel(r"$F_0/E_bA$")
        self.shfig.set_ylabel(r"$\rho_\mathrm{WV}$")
        self.shfig.set_xlim(0.23, 1.02)

        self.std0.yaxis.set_major_formatter(INTFORMATTER)


class FelixParticleH(DCVizPlotter):

    nametag = "felix_cav_phist_(.*)\.npy"
    isFamilyMember = True
    hugifyFonts = True

    def plot(self, data):

        omegas = self.get_family_member_data(data, "omegas")
        X = self.get_family_member_data(data, "X")
        H = self.get_family_member_data(data, "H")

        omega0 = 0.25

        for i, omega1 in enumerate(omegas):
            h = H[i]

            omega = (omega1 + 1)/(omega0 + 1) - 1

            self.subfigure.plot(X[:-1], h[:-1], '-s', label=r"$\Omega = %.2f$" % omega)


        self.subfigure.legend(loc="lower center",
                              numpoints=1,
                              ncol=3,
                              handlelength=1.0,
                              borderpad=0.2,
                              labelspacing=0.2,
                              columnspacing=1.0,
                              handletextpad=0.5,
                              borderaxespad=0.0,
                              frameon=False,
                              fontsize=25)

class FelixParticleHDyn(DCVizPlotter):

    nametag = "fcav_evo_(\d+)\.*"
    # transpose = True

    isFamilyMember = True
    loadSequential = True
    hugifyFonts = True
    ziggyMagicNumber = 1

    ymins = [None, None]
    ymaxs = [None, None]

    figMap = {
        "f1": "xfig",
        "f2": "yfig",
        "f3": "xyfig",
    }

    def setlims(self, fig, data, idx):
        ymin = min(data)
        ymax = max(data)

        if not self.ymaxs[idx]:
            self.ymaxs[idx] = ymax
            self.ymins[idx] = ymin
        else:
            if ymin < self.ymins[idx]:
                self.ymins[idx] = ymin
            if ymax > self.ymaxs[idx]:
                self.ymaxs[idx] = ymax

        fig.set_ylim(self.ymins[idx], self.ymaxs[idx])

    def plot(self, data):

        xdata = data.mean(axis=0)
        ydata = data.mean(axis=1)

        self.xfig.plot(xdata, "k-s")
        self.xfig.set_title("step %d" % self.getNumberForSort(self.filename))
        self.xfig.set_xlabel(r"$x$")
        self.xfig.set_ylabel(r"$P(x)$")

        self.yfig.plot(ydata, "k-s")
        self.yfig.set_title("step %d" % self.getNumberForSort(self.filename))
        self.yfig.set_xlabel(r"$y$")
        self.yfig.set_ylabel(r"$P(y)$")

        self.xyfig.pcolor(data)
        self.xyfig.set_xlabel(r"$x$")
        self.xyfig.set_ylabel(r"$y$")

        self.setlims(self.xfig, xdata, 0)
        self.setlims(self.yfig, ydata, 1)



class FelixParticleHDynCav(DCVizPlotter):

    nametag = "Fcav_evo_(\d+)\.*"
    # transpose = True

    isFamilyMember = True
    loadSequential = True
    hugifyFonts = True
    ziggyMagicNumber = 100

    ymins = [None, None]
    ymaxs = [None, None]

    figMap = {
        "f1": ["xfig", "hfig"]
    }

    def setlims(self, fig, data, idx, maxpad=0, minpad=0):
        ymin = min(data)
        ymax = max(data)

        if not self.ymaxs[idx]:
            self.ymaxs[idx] = ymax
            self.ymins[idx] = ymin
        else:
            if ymin < self.ymins[idx]:
                self.ymins[idx] = ymin
            if ymax > self.ymaxs[idx]:
                self.ymaxs[idx] = ymax

        fig.set_ylim(self.ymins[idx] - minpad, self.ymaxs[idx] + maxpad)

    def plot(self, data):

        cdata, hdata = data

        self.xfig.plot(cdata, "k-s")
        self.xfig.set_title("step %d" % self.getNumberForSort(self.filename))
        self.xfig.set_ylabel(r"$P(x)$")
        self.xfig.set_xticklabels([])

        self.hfig.plot(hdata, "k-s")
        self.hfig.set_xlabel(r"$x$")
        self.hfig.set_ylabel(r"$h(x)$")


        self.setlims(self.xfig, cdata, 0)
        self.setlims(self.hfig, hdata, 1, 1)

        # def f(i, v):
        #     if i == len(cdata) - 1:
        #         return r"$W$"
        #     else:
        #         return ""
        #
        # self.hfig.xaxis.set_major_formatter(FuncFormatter(f))

class FelixSeqC(DCVizPlotter):
    nametag = "FelixSeqC\_(.*)\.npy"
    isFamilyMember = True
    hugifyFonts = True
    labelSize = 35
    ticklabelSize = 25

    figMap = {"fig": ["subfigure",
                      "subfigurev"],
              "logfig": ["lsubfigure",
                      "lsubfigurev"],
              "xvfig": ["subfigureC",
                        "subfigurexv"]}

    specific_fig_size = {"fig": [6, 6],
                         "xvfig": [6, 7],
                         "logfig": [6, 6]}


    def adjust(self):

        for fig in self.figure_names:

            if fig == "logfig" or fig == "fig":
                self.adjust_maps[fig]['top'] = 0.95
                self.adjust_maps[fig]['bottom'] = 0.16
                self.adjust_maps[fig]['left'] = 0.22
                self.adjust_maps[fig]['right'] = 0.96
                self.adjust_maps[fig]['hspace'] = 0.3
            else:
                self.adjust_maps[fig]['top'] = 0.97
                self.adjust_maps[fig]['bottom'] = 0.16
                self.adjust_maps[fig]['left'] = 0.2
                self.adjust_maps[fig]['right'] = 0.95
                self.adjust_maps[fig]['hspace'] = 0.5

    @staticmethod
    def get_vs(t, ys):
        dt = t[1] - t[0]
        print dt

        vs = np.zeros_like(ys)

        vs[0] = (ys[1] - ys[0])/dt
        vs[-1] = (ys[-1] - ys[-2])/dt

        vs[1] = (ys[2] - ys[0])/(2*dt)
        vs[-2] = (ys[-1] - ys[-3])/(2*dt)

        for i in range(2, len(ys) - 2):
            vs[i] = (ys[i-2]-ys[i+2]+8*(ys[i+1] - ys[i-1]))/(12*dt)

        return vs

    def plot(self, data):

        if self.argv:
            cutfac = float(self.argv[0])
        else:
            cutfac = 1

        t = self.get_family_member_data(data, "t")
        ys = self.get_family_member_data(data, "ys")
        C = self.get_family_member_data(data, "C")

        imax = int(len(t)*cutfac)
        t = t[:imax]
        ys = ys[:imax]

        iterm = int(len(t)*0.25)
        ys_eq = ys[iterm:].mean()

        t /= t[-1]
        vs = self.get_vs(t, ys)

        W = np.where(t < 0.3)
        sy, i, _, _, _ = linregress(t[W], np.log(abs(1-ys[W]/ys_eq)))
        sv, i, _, _, _ = linregress(t[W], np.log(abs(vs[W])))
        print sv/sy, sv, sy

        self.subfigure.plot(t, ys, "r-", **my_props['fmt'])
        self.subfigure.set_xticklabels([])
        self.subfigure.set_ylabel(r'$y_s/L_y$')
        self.subfigure.text(0.65, ys_eq - 0.15, r"$\langle y_s\rangle/L_y$", fontsize=30)

        self.subfigurev.plot(t, vs, "r-", **my_props['fmt'])
        self.subfigurev.set_xlabel(r'$t/t_\mathrm{end}$')
        self.subfigurev.set_ylabel(r'$v_yt_\mathrm{end}/L_y$')
        self.subfigurev.set_ybound(-0.5)

        self.lsubfigure.semilogy(t, abs(1-ys/ys_eq), "r-", **my_props['fmt'])
        self.lsubfigure.set_xticklabels([])
        self.lsubfigure.set_ylabel(r'$|1-y_s/\langle y_s\rangle|$')
        self.lsubfigure.set_xbound(t.min())
        self.lsubfigure.set_ybound(1E-4)

        self.lsubfigurev.semilogy(t, abs(vs), "r-", **my_props['fmt'])
        self.lsubfigurev.set_xlabel(r'$t/t_\mathrm{end}$')
        self.lsubfigurev.set_ylabel(r'$|v_yt_\mathrm{end}/L_y|$')
        self.lsubfigurev.set_xbound(t.min())

        shift = 0.01
        self.subfigurexv.plot(ys/ys_eq, vs, "r--", **my_props['fmt'])
        self.subfigurexv.plot(ys/ys_eq, vs, "ks", **my_props['fmt'])
        self.subfigurexv.set_xlabel(r'$y_s/\langle y_s\rangle$')
        self.subfigurexv.set_ylabel(r'$v_yt_\mathrm{end}/L_y$')
        self.subfigurexv.set_xlim(min(ys/ys_eq)-shift, 1)
        self.subfigurexv.set_ybound(0)

        self.subfigureC.plot(np.arange(len(C))/float(len(C)-1), C*100, "r--", **my_props['fmt'])
        self.subfigureC.plot(np.arange(len(C))/float(len(C)-1), C*100, "ks", **my_props['fmt'])
        self.subfigureC.set_xlabel(r'$y/L_y$')
        self.subfigureC.set_ylabel(r'$10^2\langle c(y/L_y)\rangle$')
        self.subfigureC.set_xlim(-shift, 1)

        self.subfigure.set_ylim(0, 0.55)

        # mxv = self.subfigurexv.get_ylim()[1]
        # self.subfigurexv.plot([ys_eq, ys_eq], [0, mxv], "k:")
        # self.subfigurexv.set_ylim([0, mxv])
        #
        # cmin, cmax = self.subfigureC.get_ylim()
        # self.subfigureC.plot([ys_eq, ys_eq], [cmin, cmax], "k:")
        # self.subfigureC.set_ylim([cmin, cmax])




class ss_front(DCVizPlotter):
    nametag = "ss_(.*)\.npy"
    hugifyFonts = True
    labelSize = 35
    ticklabelSize = 25

    isFamilyMember = True

    figMap = {'fig': 'subfigure',
              'allfig': 'subfigure_all',
              "xfig": "subfigure_x",
              "vfig": "subfigure_v"}

    specific_fig_size = {'fig' : [7, 6],
                         'xfig': [7, 6],
                         'vfig': [7, 6],
                         'allfig': [13, 6]}

    def adjust(self):

        l = 0.165
        r = 0.98
        for f in self.figure_names:
            self.adjust_maps[f]['bottom'] = 0.16
            self.adjust_maps[f]['top'] = 0.96

            if f == "allfig":
                self.adjust_maps[f]['left'] = 0.09
                self.adjust_maps[f]['right'] = 0.99
            elif f == "xfig":
                self.adjust_maps[f]['left'] = l + l/2
                self.adjust_maps[f]['right'] = r
            else:    
                self.adjust_maps[f]['left'] = l
                self.adjust_maps[f]['right'] = r - l/2

    #plotOnly = "fig"

    def plot(self, data):

        t = self.get_family_member_data(data, "felix_t")
        x = self.get_family_member_data(data, "felix_x")
        tslice = self.get_family_member_data(data, "felix_tslice")
        xslice = self.get_family_member_data(data, "felix_fit")
        v = self.get_family_member_data(data, "felix_v")
        vslice = self.get_family_member_data(data, "felix_vfit")

        ss_edges = self.get_family_member_data(data, "edges")
        all_edges = self.get_family_member_data(data, "all_edges")

        W, L = all_edges.shape

        self.subfigure_all.pcolor(all_edges, cmap='gist_earth_r')
        self.subfigure_all.set_xlabel(r"$x/L_x$")
        self.subfigure_all.set_ylabel(r"$y/L_y$")

        self.subfigure.pcolor(ss_edges, cmap='gist_earth_r', vmin=all_edges.min(), vmax=all_edges.max())
        self.subfigure.set_xlabel(r"$(x - n\Delta T v_s)/L_x$")
        self.subfigure.set_ylabel(r"$y/L_y$")

        for subfigure in [self.subfigure, self.subfigure_all]:
            subfigure.xaxis.set_major_formatter(FuncFormatter(lambda f, _: r'$%g$' % (f/float(L))))
            subfigure.yaxis.set_major_formatter(FuncFormatter(lambda f, _: r'$%g$' % (f/float(W))))


        if False:
            Ry = 0.74*W
            Rx = 0.75*L
            dy = 0.17*W
            dx = 0

            theta = np.linspace(0, np.pi/2)

            xv = Rx*np.cos(theta) - dx
            yv = Ry*np.sin(theta) - dy

            self.subfigure.plot(xv, yv, 'k--', linewidth=4, label=r"$\mathrm{Ellipse}$")
            Ry = 0.6*W
            Rx = 0.6*L
            xv = np.linspace(0, Rx, 1000)
            yv = Ry*np.log((Rx + 1) - xv)/np.log(Rx + 1)

            self.subfigure.plot(xv, yv, 'r-', linewidth=4, label=r"$\mathrm{Exponential}$")


            leg = self.subfigure.legend(loc="center",
                                      numpoints=1,
                                      ncol=1,
                                      handlelength=1.0,
                                      markerscale=20.0,
                                      borderpad=0.2,
                                      labelspacing=0.2,
                                      columnspacing=1.0,
                                      handletextpad=0.5,
                                      borderaxespad=0.0,
                                      frameon=False,
                                      fontsize=20,
                                      bbox_to_anchor=(0.75, 0.875))

            leg.get_frame().set_fill(not (self.toFile and self.transparent))

        else:
            ly = 0.58
            a = 1

            u = np.linspace(0.01, 1, 1000)

            fac = np.sqrt(a**2-u**2)
            x1 =  a*np.log((a-fac)/u)
            xfuck = fac + x1

            xshift = 1.05
            self.subfigure.plot((xshift+xfuck/4)*ly*L, (1-u)*ly*W,'k--', linewidth=4, label=r"$\mathrm{Ellipse}$")

        self.subfigure.set_xlim(0, 0.7*L)
        self.subfigure.set_ylim(0, 0.7*W)

        #-- xvfigs

        rs = 0.02
        rsd = 0.05
        self.subfigure_x.plot(t, x/float(L), '--ks')
        self.subfigure_x.plot(tslice-rs, xslice/L+rs, "r-", linewidth=2)
        self.subfigure_x.set_xlabel(r'$t/t_\mathrm{end}$')
        self.subfigure_x.set_ylabel(r'$x_s/L_x$')
        self.subfigure_x.text(tslice.mean()-rsd, xslice.mean()/L+rsd,
                              r"$\sim v_s t$",
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontsize=30,
                              rotation=37.5)


        self.subfigure_v.plot(t, v/L, '--ks')
        self.subfigure_v.plot(tslice, vslice/L, 'r-', linewidth=2)
        self.subfigure_v.set_xlabel(r'$t/t_\mathrm{end}$')
        self.subfigure_v.set_ylabel(r'$v_xt_\mathrm{end}/L_x$')
        self.subfigure_v.text(tslice[-1]+0.05, vslice.mean()/L,
                              r"$v_s$",
                              verticalalignment='center',
                              horizontalalignment='left',
                              fontsize=30)

        for subfigure in [self.subfigure_v, self.subfigure_x]:
            subfigure.xaxis.set_major_formatter(FuncFormatter(lambda f, _: r'$%g$' % f))
            subfigure.yaxis.set_major_formatter(FuncFormatter(lambda f, _: r'$%g$' % f))

        self.subfigure_x.set_ylim(0, 1.05)


class diff_model_times(DCVizPlotter):
    nametag = "diff_model_times_(.*)\.npy"

    isFamilyMember = True
    hugifyFonts = True

    figMap = {
        "fig1":["s0", "s1", "s2", "s3"],
        "fig2": ["t0", "t1", "t2", "t3"],
        "fig": ["f0", "f1"]
    }

    stack = "H"
    specific_fig_size = {"fig": [10, 5]}

    plotOnly = ["fig"]

    ids = ["uniform", "lattice", "radial", "pathfind"]

    names = {"uniform"  : "Uniform",
             "lattice"  : "Lattice",
             "radial"   : "Radial",
             "pathfind" : "Pathfinding"}

    styles = {"uniform"  : "s",
             "lattice"  : "o",
             "radial"   : "^",
             "pathfind" : "d"}

    colors = {"uniform"  : "k",
              "lattice"  : "r",
              "radial"   : "b",
              "pathfind" : "g"}

    def adjust(self):
        left = 0.11

        self.adjust_maps["fig"]["left"] = left
        self.adjust_maps["fig"]["right"] = 1-left
        self.adjust_maps["fig"]["top"] = 0.92
        self.adjust_maps["fig"]["bottom"] = 0.15
        self.adjust_maps["fig"]["wspace"] = 0.14

    def plot(self, data):
        alphas = self.get_family_member_data(data, "alphas")
        heights = self.get_family_member_data(data, "heights")
        results = self.get_family_member_data(data, "results")

        for i, height in enumerate(heights):
            f = eval("self.s%d" % i)
            t = eval("self.t%d" % i)

            r = None
            for j in range(4):

                if not r:
                    r = [results[j, :, i, 0], results[j, :, i, 1]]

                r0 = results[j, :, i, 0]/r[0]
                r1 = results[j, :, i, 1]/r[1]

                f.semilogy(alphas, r0, label="h=%g" % height)
                t.semilogy(alphas, r1, label="h=%g" % height)

        titles = [r"$\mathrm{Per\,\,cycle}$", r"$\mathrm{Per\,\,surface\,\,event}$"]

        index_list = [[1, 2], [3, 4], [1, 1], [5, 4], [3, 2], [7, 4], [2, 1]]
        def fff(v, i):

            if not i:
                print "wtf", v

            a, b = index_list[i]

            if b == 1:
                return r"$%d$" % a
            else:
                return r"\Large{$%d/%d$}" % (a, b)

        bbox = [(0.75, 0.87), (0.75, 0.87)]

        I = 1
        for k in range(2):
            f = eval("self.f%d" % k)
            f.axes.xaxis.set_major_formatter(FuncFormatter(fff))
            f.set_title(titles[k])
            f.set_xlabel(r"$\alpha\equiv E_b/kT$")
            f.set_xticks(alphas)
            f.set_xlim(0.45, 2.05)
            f.set_ylim(0.1, 15000)
            for j in range(4):
                id = self.ids[j]

                if k == 0 and j < 2:
                    label = r"$\mathrm{%s}$" % self.names[id]
                elif k == 1 and j >= 2:
                    label = r"$\mathrm{%s}$" % self.names[id]
                else:
                    label = None

                f.semilogy(alphas,
                           results[j, :, I, 1-k]/results[0, 0, I, 1-k],
                           label = label,
                           marker=self.styles[id],
                           markeredgecolor=self.colors[id],
                           linestyle="--",
                           color=self.colors[id],
                           linewidth=1,
                           fillstyle='none',
                           markersize=7,
                           markeredgewidth=1.5)


            l = f.axes.legend(loc="center",
                               numpoints=1,
                               ncol=1,
                               handlelength=1.0,
                               borderpad=0.2,
                               labelspacing=0.2,
                               columnspacing=0.3,
                               handletextpad=0.25,
                               borderaxespad=0.0,
                               frameon=False,
                               fontsize=20,
                               bbox_to_anchor=bbox[k])

            l.get_frame().set_fill(not (self.toFile and self.transparent))

        self.f0.set_ylabel(r"$t_\mathrm{CPU}/t_\mathrm{CPU}^\ast$")
        self.f1.set_yticklabels([])
