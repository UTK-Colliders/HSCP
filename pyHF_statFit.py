# First run `/usr/bin/python3.12 -m venv pyHF`, followed by `source pyHF/bin/activate` to access the pyhf library
# For your first time, run `pip install pyhf` and `pip install matplotlib`
# Run `deactivate` when you are done with the venv

import pyhf
import argparse
import numpy as np
from matplotlib import pyplot as plt
from pyhf.contrib.viz import brazil
import mplhep as hep
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase


class HandlerBand(HandlerBase):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent,
        width, height, fontsize, trans
    ):
        thin_top, thick, thin_bottom = orig_handle

        y_top = ydescent + 0.80 * height
        y_mid = ydescent + 0.50 * height
        y_bottom = ydescent + 0.20 * height

        # Shaded band between the two thin lines
        band = Rectangle(
            (xdescent, y_bottom),
            width,
            y_top - y_bottom,
            facecolor=thin_top.get_color(),
            alpha=thin_top.get_alpha(),
            edgecolor="none",
            transform=trans
        )

        # Top thin line
        top_line = Line2D(
            [xdescent, xdescent + width],
            [y_top, y_top],
            color=thin_top.get_color(),
            linewidth=thin_top.get_linewidth(),
            linestyle=thin_top.get_linestyle(),
            transform=trans
        )

        # Center thick line
        middle_line = Line2D(
            [xdescent, xdescent + width],
            [y_mid, y_mid],
            color=thick.get_color(),
            linewidth=thick.get_linewidth(),
            linestyle=thick.get_linestyle(),
            transform=trans
        )

        # Bottom thin line
        bottom_line = Line2D(
            [xdescent, xdescent + width],
            [y_bottom, y_bottom],
            color=thin_bottom.get_color(),
            linewidth=thin_bottom.get_linewidth(),
            linestyle=thin_bottom.get_linestyle(),
            transform=trans
        )

        return [band, top_line, middle_line, bottom_line]

class StatFit:
    def __init__(self, n_signal, n_background, background_uncertainty, gluino_mass=1000, neutralino_mass=100, tau=0.1, quark_decay="uds", do_debug_print=False):
        self.do_debug_print = do_debug_print
        self.gluino_mass = gluino_mass
        self.neutralino_mass = neutralino_mass
        self.tau = tau
        self.quark_decay = quark_decay

        self.model = pyhf.simplemodels.uncorrelated_background(
            signal=[n_signal], bkg=[n_background], bkg_uncertainty=[background_uncertainty]
        )

        if self.do_debug_print:
            self.debug_print()

        #pp -> gogo 13.6TeV theory xsecs in fb with error
        self.xsecDict = {1000: [484.2, 46.48],
                         1200: [128.8, 13.83],
                         1400: [38.83, 4.780],
                         1600: [12.82, 1.891],
                         1800: [4.524, 0.8342],
                         2000: [1.684, 0.4095],
                         2200: [0.6535, 0.2158],
                         2400: [0.2627, 0.1192],
                         2600: [0.1089, 0.06709],
                         2800: [0.04623, 0.03916]}

        # Look up the theory cross section (and its uncertainty) for this mass point
        if self.gluino_mass not in self.xsecDict:
            raise ValueError(
                f"No theory cross section entry for gluino mass {self.gluino_mass}. "
                f"Available masses: {sorted(self.xsecDict.keys())}"
            )
        self.theory_xsec, self.theory_xsec_unc = self.xsecDict[self.gluino_mass]

    def background_only_parameters(self, pars):
        bkg_pars = pars.copy()
        bkg_pars[self.model.config.poi_index] = 0
        return bkg_pars

    def compute_log_likelihood(self, observations, sigPlusBkg_pars, bkg_pars):
        log_likelihood_sigPlusBkg = self.model.logpdf(pars=sigPlusBkg_pars, data=observations)
        log_likelihood_bkgOnly = self.model.logpdf(pars=bkg_pars, data=observations)
        if self.do_debug_print:
            print(f"Log likelihood for signal + background model: {log_likelihood_sigPlusBkg}")
            print(f"Log likelihood for background only model: {log_likelihood_bkgOnly}")

        return log_likelihood_sigPlusBkg, log_likelihood_bkgOnly

    def perform_fit(self, observations):
        fit = pyhf.infer.mle.fit(data=observations, pdf=self.model)
        if self.do_debug_print:
            print(f"Model parameter best fits: {fit}")
            print("-----END DEBUG INFO-----\n")

        return fit

    def perform_SM_only_hypothesis_test(self, observations):
        CLs_obs, CLs_exp = pyhf.infer.hypotest(
            1.0,  # test the nominal signal hypothesis (mu=1) against background-only
            observations, self.model, test_stat="qtilde", return_expected_set=True,
        )

        print(f"      Observed CLs: {CLs_obs:.4f}")
        for expected_value, n_sigma in zip(CLs_exp, np.arange(-2, 3)):
            print(f"Expected CLs({n_sigma:2d} σ): {expected_value:.4f}")

    def mu_to_xsec(self, mu_value):
        return mu_value * self.theory_xsec

    def set_limits(self, observations):
        poi_values = np.linspace(0.1, 5, 50)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
            observations, self.model, poi_values, level=0.05, return_results=True
        )

        print(f"Upper limit (obs): μ = {obs_limit:.4f}")
        print(f"Upper limit (exp): μ = {exp_limits[2]:.4f}")

        x_sec_scan = [self.mu_to_xsec(mu) for mu in poi_values]
        obs_xsec_limit = self.mu_to_xsec(obs_limit)
        exp_xsec_limits = [self.mu_to_xsec(mu) for mu in exp_limits]
        print(f"Upper limit (obs): σ = {obs_xsec_limit:.4f} fb "
                f"(theory σ = {self.theory_xsec:.4f} +- {self.theory_xsec_unc:.4f} fb)")
        print(f"Upper limit (exp): σ = {exp_xsec_limits[2]:.4f} fb")

        return x_sec_scan, results, obs_xsec_limit, exp_xsec_limits

    def plot_limits(self, x_sec_scan, results, output, gluino_mass):
        fig, ax = plt.subplots()
        fig.set_size_inches(10.5, 7)

        # Define the title and labels
        if self.quark_decay == "uds":
            decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + (u\bar{{u}},d\bar{{d}},s\bar{{s}})$"
        else:
            decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + t\bar{{t}}$"
        plt.text(.55, .99, fr'$m_{{\tilde{{g}}}}={self.gluino_mass}$ [GeV]', ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)
        plt.text(.55, .94, fr'$\tau_{{\tilde{{g}}}}={self.tau}$ [ns]', ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)
        plt.text(.55, .89, decayString, ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)
        plt.text(.55, .84, fr'$m_{{\tilde{{\chi}}^0_1}}={self.neutralino_mass}$ [GeV]', ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)

        # Fake results for testing warning
        plt.text(0.01, 0.94, "Fake results for testing", transform=ax.transAxes, fontname="TeX Gyre Heros", size=20, fontweight='bold')

        artists = brazil.plot_results(x_sec_scan, results, ax=ax)
        first_legend = ax.get_legend()
        ax.set_xlabel(r"$\sigma$ [fb]")
        #ax.legend(loc="best")

        # CMS specific
        hep.cms.label("Work in progress", loc=0, ax=ax, fontsize=18, fontproperties="TeX Gyre Heros:italic", com=13.6, lumi=283.8)

        # Plot a vertical band for the theory xsec
        theory_xsec = self.xsecDict[gluino_mass][0]
        theory_xsec_error = self.xsecDict[gluino_mass][1]
        plt.axvline(x=theory_xsec, color='indianred')
        plt.axvspan(theory_xsec - theory_xsec_error, theory_xsec + theory_xsec_error, color='lightcoral', alpha=0.3)

        exp_thin_a = Line2D([0], [0], color="lightcoral", linewidth=2, alpha=0.3)
        exp_thick = Line2D([0], [0], color="indianred", linewidth=2)
        exp_thin_b = Line2D([0], [0], color="lightcoral", linewidth=2, alpha=0.3)

        legend_handles = [
            (exp_thin_a, exp_thick, exp_thin_b),
        ]

        legend_labels = [
            r"$\sigma_{theory} \pm$ uncertainty",
        ]

        handles = first_legend.legend_handles if first_legend is not None else []
        labels = [text.get_text() for text in first_legend.get_texts()] if first_legend is not None else []
        ax.legend(
            handles + legend_handles,
            labels + legend_labels,
            handler_map={tuple: HandlerBand()},
            loc="upper right"
        )

        plt.savefig(output)

    def debug_print(self):
        print("-----DEBUG INFO-----")
        print(f"  channels: {self.model.config.channels}")
        print(f"     nbins: {self.model.config.channel_nbins}")
        print(f"   samples: {self.model.config.samples}")
        print(f" modifiers: {self.model.config.modifiers}")
        print(f"parameters: {self.model.config.parameters}")
        print(f"  nauxdata: {self.model.config.nauxdata}")
        print(f"   auxdata: {self.model.config.auxdata}\n")

        print(f"Suggested initial parameters:          {self.model.config.suggested_init()}")
        print(f"Suggested parameter bounds:            {self.model.config.suggested_bounds()}")
        print(f"Should parameters be fixed during fit: {self.model.config.suggested_fixed()}\n")


def main():
    parser = argparse.ArgumentParser(description="Perform a single-channel statistical fit with uncorrelated background")
    parser.add_argument("--n_signal", type=float, default=10.0, help="Number of signal events expected in the signal region")
    parser.add_argument("--n_background", type=float, default=1.0, help="Number of background events expected in the signal region")
    parser.add_argument("--background_uncertainty", type=float, default=1.0, help="Uncertainty on number of background events expected in the signal region. For example, if the number of background events is 50 +- 5, you would pass 5 into this argument.")
    parser.add_argument("--n_observed", type=int, default=1, help="Number of events observed in data in the signal region")
    parser.add_argument("--gluino_mass", type=int, default=1000, help="Gluino mass point in GeV, used to look up the theory cross section for mu->xsec conversion and plotting")
    parser.add_argument("--neutralino_mass", type=int, default=100)
    parser.add_argument("--tau", type=float, default=0.1)
    parser.add_argument("--quark_decay", default="uds")
    parser.add_argument("--limitPlotOutput", default="pdfs/output.pdf")
    parser.add_argument("--do_debug_print", type=bool, default=False)
    args = parser.parse_args()

    stat_fit = StatFit(args.n_signal, args.n_background, args.background_uncertainty,
                        gluino_mass=args.gluino_mass, do_debug_print=args.do_debug_print)
    observations = [args.n_observed] + stat_fit.model.config.auxdata

    # For a single channel, the parameters are [mu, gamma_b] (see https://pyhf.github.io/pyhf-tutorial/helloworld/)
    # Mu is the signal strength, gamma_b is a multiplicative factor calculated by the background uncertainty given an uncorrelated background
    # Here by default the model automatically calculates what mu and gamma_b should be
    init_pars = stat_fit.model.config.suggested_init()
    stat_fit.model.expected_actualdata(init_pars)

    # Check the background only model
    bkg_pars = stat_fit.background_only_parameters(init_pars)
    stat_fit.model.expected_actualdata(bkg_pars)

    # Compute the loglikelihood
    stat_fit.compute_log_likelihood(observations, init_pars, bkg_pars)

    # Perform the model parameter fit
    stat_fit.perform_fit(observations)

    # Perform the SM only hypothesis test
    stat_fit.perform_SM_only_hypothesis_test(observations)

    # Set upper limits (mu, and xsec if gluino_mass was given) and plot them
    poi_values, results, obs_limit, exp_limits = stat_fit.set_limits(observations)
    stat_fit.plot_limits(poi_values, results, args.limitPlotOutput, args.gluino_mass)

if __name__ == "__main__":
    main()