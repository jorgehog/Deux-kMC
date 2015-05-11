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

    figMap = {"converge_figure": "subfigure", "value_figure": "subfigure2", "slope_figure": "subfigure3"}

    def trans(self, v):
        return (v/v.min())**10

    shapes = ["s", "^", "v"]

    def plot(self, data):

        clip = 1.1
        start = 1
        transform = False

        rmslabel=r"$\sigma(h)$"
        sslabel=r"$\sigma(s)$"
        whlabel=r"$\langle \delta h_l(t) \rangle$"


        I, J, K = [int(x) for x in self.argv]

        E0_values = data[self.get_family_index_from_name("steadystate_E0.npy")]
        alpha_values = data[self.get_family_index_from_name("steadystate_alpha.npy")]
        mu_shift_values = data[self.get_family_index_from_name("steadystate_mu_shift.npy")]

        count = 0
        slopes = []
        for i, E0 in enumerate(E0_values):
            rms_values = []
            print i, E0

            for j, alpha in enumerate(alpha_values):
                for k, mu_shift in enumerate(mu_shift_values):

                    if i == I and j == J and k == K:
                        print "E0 = %g, alpha= %g, mus = %g" % (E0, alpha, mu_shift)

                    time = data[self.get_family_index_from_name("steadystate_Time_%d.npy" % count)]
                    L = int(len(time)/clip)
                    time = time[start:L]

                    rms = data[self.get_family_index_from_name("steadystate_HeightRMS_%d.npy" % count)][start:L]

                    ss = data[self.get_family_index_from_name("steadystate_SurfaceSize_%d.npy" % count)][start:L]

                    wh = data[self.get_family_index_from_name("steadystate_PressureWall_%d.npy" % count)][start:L]

                    if i == I and j == J and k == K:

                        if transform:
                            self.subfigure.loglog(time, rms, label=rmslabel)
                            self.subfigure.loglog(time, (ss - ss[0])/(ss.max() - ss[0]), label=sslabel)
                            self.subfigure.loglog(time, (wh - wh[0])/(wh.max() - wh[0]), label=whlabel)
                        else:
                            self.subfigure.loglog(time, self.trans(rms), "-", label=rmslabel,
                              linewidth=3,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="red")

                            self.subfigure.loglog(time, self.trans(wh), "-.", label=whlabel,
                              linewidth=3,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="black")

                            self.subfigure.loglog(time, self.trans(ss), "--", label=sslabel,
                              linewidth=3,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="green")



                    L2 = len(rms)/4

                    rms_value = rms[L2:].mean()
                    ss_value = ss[L2:].mean()
                    wh_value = wh[L2:].mean()

                    if j == J:

                        rms_values.append(rms_value)

                    count += 1

            self.subfigure2.plot(exp(mu_shift_values), rms_values, "k--" + self.shapes[i], label="E0=%.2f" % E0,
                                 linewidth=1,
                                 fillstyle='none',
                                 markersize=7,
                                 markeredgewidth=1.5,
                                 color="black")

            slope, intercept, _, _, _ = linregress(exp(mu_shift_values), rms_values)
            print i, E0, intercept, slope
            slopes.append(slope)

        self.subfigure3.plot(E0_values, slopes)

        self.subfigure.set_xlabel(r"$t$")
        self.subfigure.axes.get_yaxis().set_ticklabels([])
        self.subfigure.legend(loc="upper left")
        # self.subfigure.set_xlim(0, self.subfigure.get_xlim()[1])
        # self.subfigure.set_ylim(0, 1.1)

        self.subfigure2.set_xlabel(r"$c/c_0$")
        self.subfigure2.set_ylabel(rmslabel)
        self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1.2)


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
        self.subfigure.set_xlabel(r"$\Omega$")
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
              "neighborstuff": "subfigure10"}

    # figMap = {"omega_vs_v": ["subfigure",
    #           "subfigure2",
    #           "subfigure4",
    #           "subfigure5",
    #           "subfigure6",
    #           "subfigure7",
    #           "subfigure8",
    #           "subfigure9"]}

    # figMap = {"asd": "subfigure"}

    plot_values = [0.2, 0.4]
    shapes = ["s", "^", "v"]

    def plot_and_slopify(self, E0, omega, mu_shifts, mu, v):
        mu0 = (mu - mu_shifts).mean()
        c_over_c0 = exp(mu0 + mu_shifts)
        v *= c_over_c0

        if E0 == 0:
            if not (omega == (c_over_c0-1)).all():

                print "OMEGA FAIL"
                print omega
                print c_over_c0
                self.Exit()

        idx_high = np.where(omega >= 0)
        idx_low = np.where(omega < 0)

        idx_chosen = idx_high

        slope, intercept, _, _, _ = linregress(omega[idx_chosen], v[idx_chosen])
        # slope_low, intercept, _, _, _ = linregress(omega[idx_low], v[idx_low])
        #
        # print E0, slope/slope_low - 1

        return mu0, slope, v

    def uberplot(self, omega, v, k, E0):

        idx_undersat = np.where(omega <= 0)
        idx_oversat = np.where(omega >= 0)

        omega_undersat = abs(omega[idx_undersat])
        omega_oversat = omega[idx_oversat]

        v_undersat = abs(v[idx_undersat])
        v_oversat = v[idx_oversat]

        v = v_undersat
        omega = omega_undersat

        v_log_under = np.sign(v)*np.log(abs(v))
        o_log_under = np.sign(omega)*np.log(abs(omega))
        self.subfigure.plot(o_log_under, v_log_under, "k--" + self.shapes[k], label="E0=%.2f" % E0,
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


        self.subfigure.plot(o_log_over - xshift, v_log_over - yshift, "k-." + self.shapes[k],
                              linewidth=1,
                              fillstyle='none',
                              markersize=7,
                              markeredgewidth=1.5,
                              color="black")
        #


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

        slopes = np.zeros(shape=(len(E0_values), len(alpha_values), len(r0_values), len(s0_values)))
        mu_zero = np.zeros_like(slopes)

        J, L, M = [int(x) for x in self.argv]

        print "alpha=%g, s0 = %g, r0 = %g" % (alpha_values[J], s0_values[M+1], r0_values[L+1])

        for j in range(len(alpha_values)):
            _, slope0, v0 = self.plot_and_slopify(0, omega, mu_shift_values, mu_shift_values, v_values[0, j, :, 0, 0])

            slopes[0, j, :, :] = slope0
            mu_zero[0, j, :, :] = 0

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

                            mu0, slope, v = self.plot_and_slopify(E0,
                                                                  omega,
                                                                  mu_shift_values,
                                                                  mu_values[i+1, j, :, l+1, m+1],
                                                                  v_values[i+1, j, :, l+1, m+1])

                            if slope < 0:
                                print "ERROR slope=%g, E0(%d)=%g, alpha(%d)=%g, r0(%d)=%g, s0(%d)=%g" % (slope, i, E0, j, alpha, l, r0, m, s0)

                            slopes[i+1][j][l+1][m+1] = slope
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

        ka = 0
        na = 4
        nm = na - 2
        d = len(alpha_values)/nm
        for j, alpha in enumerate(alpha_values):
            for l, r0 in enumerate(r0_values[1:]):
                for m, s0 in enumerate(s0_values[1:]):

                    if l == L and m == M:
                        if j == 0 or j == len(alpha_values) - 1 or j % d == 0:

                            self.subfigure2.plot(E0_values, slopes[:, j, l+1, m+1], "k--" + self.shapes[ka], label=r"$\alpha=%g$" % round(alpha, 1),
                                                 linewidth=1,
                                                 fillstyle='none',
                                                 markersize=7,
                                                 markeredgewidth=1.5)

                            ka += 1

                    log_k_slope, intercept, _, _, _ = linregress(E0_values, np.log(slopes[:, j, l+1, m+1]))

                    log_k_E0_slopes[j, l+1, m+1] = log_k_slope
                    log_k_cutz[j, l+1, m+1] = intercept

        alpha_slopes = np.zeros(shape=(len(r0_values), len(s0_values)))
        alpha_cutz = np.zeros_like(alpha_slopes)

        for l, r0 in enumerate(r0_values[1:]):
                for m, s0 in enumerate(s0_values[1:]):
                    alpha_slope, intercept, _, _, _ = linregress(alpha_values, log_k_E0_slopes[:, l+1, m+1])

                    alpha_slopes[l+1, m+1] = alpha_slope
                    alpha_cutz[l+1, m+1] = intercept



        self.subfigure4.plot(alpha_values, log_k_E0_slopes[:, L+1, M+1], "ks--",
                                 linewidth=1,
                                 fillstyle='none',
                                 markersize=7,
                                 markeredgewidth=1.5)

        S_tot = 0
        count2 = 0
        for l, r0 in enumerate(r0_values[1:]):
            S = 0
            count = 0

            if r0 == 0.5 or r0 == 1:
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

        self.subfigure5.legend(numpoints=1, handlelength=1.2)

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
        self.subfigure8.set_ylabel(r"$\mathrm{slopes}$")

        comb = alpha_slopes[1:, 1:].mean(axis=1)

        self.subfigure9.plot(r0_values[1:], comb)
        self.subfigure9.set_xlabel(r"$\lambda_D$")
        self.subfigure9.set_ylabel(r"$\mathrm{avg slope}$")

        self.subfigure.set_xlabel(r"$|\Omega + 1|$")
        self.subfigure.set_ylabel(r"$|\dot{H}|$")
        self.subfigure.legend(loc="upper left", numpoints=1, handlelength=1.2)
        # self.subfigure.axes.set_xscale('log')
        # self.subfigure.axes.set_yscale('log')

        self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1.2)
        self.subfigure2.axes.set_yscale('log')

        self.subfigure2.set_xlabel(r"$E_0$")
        self.subfigure2.set_ylabel(r"$k_g = \partial \dot{H} / \partial \Omega$")

        self.subfigure4.set_xlabel(r"$\alpha$")
        self.subfigure4.set_ylabel(r"$\mathrm{slope}$")
        self.subfigure4.set_ybound(0)

        self.subfigure5.set_xlabel(r"$\alpha$")
        self.subfigure5.set_ylabel(r"$\mathrm{log k shift}$")



        # self.subfigure3.set_xlabel(r"$E_0$")
        # self.subfigure3.set_ylabel(r"$\gamma_\mathrm{eq}$")


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

    plot_values = [0.01, 0.1, 0.19]

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
            self.gammaslopes.plot([0, alphas.max()], [0, slopes[n]*alphas.max()], "r--")

            n_plots += 1

        self.gammaslopes.set_xbound(0)
        self.gammaslopes.set_ybound(0)
        self.gammaslopes.legend(loc="upper left")

        self.gammaslopes.set_xlabel(r"$\alpha$")
        self.gammaslopes.set_ylabel(r"$\gamma_\mathrm{eq}$")

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

        self.E0slopes.set_xlabel(r"$E_0$")
        self.E0slopes.set_ylabel(r"$K(E_0)$")

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
            self.subfigure.loglog(1./analytical[:, 0], analytical[:, 1], 'r-',
                                linewidth=3,
                                label="$E_0 = 0.00$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure2.loglog(1./analytical[:, 0], analytical[:, 2], 'r-',
                                linewidth=3,
                                label="$E_0 = 0.00$",
                                fillstyle='none',
                                markersize=5)
            self.subfigure3.loglog(1./analytical[:, 0], analytical[:, 0]**p*analytical[:, 2], 'r-',
                                linewidth=3,
                                label="$E_0 = 0.00$",
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

            self.subfigure3.loglog(1./alpha_array, var_s_array*alpha_array**p, 'k%s' % shapes[nplots],
                     fillstyle='none',
                     label="$E_0 = %.2f$" % E0_value,
                     markersize=7,
                     markeredgewidth=1.5,
                     linewidth=1)

            nplots += 1

        self.subfigure.set_xlabel(r"$1/\alpha$")
        self.subfigure.set_ylabel(r"$\langle s \rangle / L$")
        ax = self.subfigure.axes.twinx()
        ax.set_ylabel(r"$\langle E \rangle /E_bL$")
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])
        self.subfigure.set_xbound(0)
        self.subfigure.legend(loc="lower right", numpoints=1, handlelength=1)


        self.subfigure2.set_xlabel(r"$1/\alpha$")
        self.subfigure2.set_ylabel(r"$\sigma (s) / L$")
        ax2 = self.subfigure2.axes.twinx()
        ax2.set_ylabel(r"$\sigma (E) / E_bL$")
        ax2.yaxis.set_ticks([])
        ax2.yaxis.set_ticklabels([])
        self.subfigure2.set_xbound(0)
        # self.subfigure2.legend(loc="upper left", numpoints=1, handlelength=1)

        self.subfigure3.set_xlabel(r"$1/\alpha$")
        self.subfigure3.set_ylabel(r"$\alpha^2\sigma(s)^2/L^2$")
        ax3 = self.subfigure3.axes.twinx()
        ax3.set_ylabel(r"$C_V/k$")
        ax3.yaxis.set_ticks([])
        ax3.yaxis.set_ticklabels([])
        self.subfigure3.set_xbound(0)
        # self.subfigure3.legend(numpoints=1, handlelength=1, loc="upper right")


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