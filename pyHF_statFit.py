# First run `/usr/bin/python3.12 -m venv pyHF`, followed by `source pyHF/bin/activate` to access the pyhf library
# For your first time, run `pip install pyhf` and `pip install matplotlib`
# Run `deactivate` when you are done with the venv

import pyhf
import argparse
import numpy as np
from matplotlib import pyplot as plt
from pyhf.contrib.viz import brazil

class StatFit:
    def __init__(self, n_signal, n_background, background_uncertainty, do_debug_print=False):
        self.do_debug_print = do_debug_print

        self.model = pyhf.simplemodels.uncorrelated_background(
            signal=[n_signal], bkg=[n_background], bkg_uncertainty=[background_uncertainty]
        )

        if self.do_debug_print:
            self.debug_print()

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

    def set_limits(self, observations):
        poi_values = np.linspace(0.1, 5, 50)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
            observations, self.model, poi_values, level=0.05, return_results=True
        )

        print(f"Upper limit (obs): μ = {obs_limit:.4f}")
        print(f"Upper limit (exp): μ = {exp_limits[2]:.4f}")

        return poi_values, results, obs_limit, exp_limits

    def plot_limits(self, poi_values, results, output):
        fig, ax = plt.subplots()
        fig.set_size_inches(10.5, 7)
        ax.set_title("Hypothesis Tests")

        artists = brazil.plot_results(poi_values, results, ax=ax)
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
    parser.add_argument("--n_signal", type=float, default=1.0, help="Number of signal events expected in the signal region")
    parser.add_argument("--n_background", type=float, default=1.0, help="Number of background events expected in the signal region")
    parser.add_argument("--background_uncertainty", type=float, default=1.0, help="Uncertainty on number of background events expected in the signal region. For example, if the number of background events is 50 +- 5, you would pass 5 into this argument.")
    parser.add_argument("--n_observed", type=int, default=1, help="Number of events observed in data in the signal region")
    parser.add_argument("--limitPlotOutput", default="output.pdf")
    parser.add_argument("--do_debug_print", type=bool, default=False)
    args = parser.parse_args()

    stat_fit = StatFit(args.n_signal, args.n_background, args.background_uncertainty, args.do_debug_print)
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

    # Set upper limits and plot them
    poi_values, results, obs_limit, exp_limits = stat_fit.set_limits(observations)
    stat_fit.plot_limits(poi_values, results, args.limitPlotOutput)

if __name__ == "__main__":
    main()