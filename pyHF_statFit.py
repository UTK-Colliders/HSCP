# First run `/usr/bin/python3.12 -m venv pyHF`, followed by `source pyHF/bin/activate` to access the pyhf library
# For your first time, run `pip install pyhf`
# Run `deactivate` when you are done with the venv

import pyhf
import argparse


def build_model(n_signal, n_background, background_uncertainty):
    model = pyhf.simplemodels.uncorrelated_background(
        signal=[n_signal], bkg=[n_background], bkg_uncertainty=[background_uncertainty]
    )

    if do_debug_print:
        debug_print(model)

    return model

def background_only_parameters(model, pars):
    bkg_pars = pars.copy()
    bkg_pars[model.config.poi_index] = 0
    return bkg_pars

def compute_log_likelihood(model, n_observed, sigPlusBkg_pars, bkg_pars):
    observations = [n_observed] + model.config.auxdata
    if do_debug_print:
        print(f"Log likelihood for signal + background model: {model.logpdf(pars=sigPlusBkg_pars, data=observations)}")
        print(f"Log likelihood for background only model: {model.logpdf(pars=bkg_pars, data=observations)}")

    return model.logpdf(pars=sigPlusBkg_pars, data=observations), model.logpdf(pars=bkg_pars, data=observations)

def debug_print(model):
    print(f"  channels: {model.config.channels}")
    print(f"     nbins: {model.config.channel_nbins}")
    print(f"   samples: {model.config.samples}")
    print(f" modifiers: {model.config.modifiers}")
    print(f"parameters: {model.config.parameters}")
    print(f"  nauxdata: {model.config.nauxdata}")
    print(f"   auxdata: {model.config.auxdata}\n")

    print(f"Suggested initial parameters:          {model.config.suggested_init()}")
    print(f"Suggested parameter bounds:            {model.config.suggested_bounds()}")
    print(f"Should parameters be fixed during fit: {model.config.suggested_fixed()}")

def main():
    parser = argparse.ArgumentParser(description="Perform a single-channel statistical fit with uncorrelated background")
    parser.add_argument("--n_signal", type=float, default=1.0, help="Number of signal events expected in the signal region")
    parser.add_argument("--n_background", type=float, default=1.0, help="Number of background events expected in the signal region")
    parser.add_argument("--background_uncertainty", type=float, default=1.0, help="Uncertainty on number of background events expected in the signal region. For example, if the number of background events is 50 +- 5, you would pass 5 into this argument.")
    parser.add_argument("--n_observed", type=int, default=1, help="Number of events observed in data in the signal region")
    parser.add_argument("--do_debug_print", type=bool, default=False)
    args = parser.parse_args()
    global do_debug_print
    do_debug_print = args.do_debug_print

    model = build_model(args.n_signal, args.n_background, args.background_uncertainty)

    # For a single channel, the parameters are [mu, gamma_b] (see https://pyhf.github.io/pyhf-tutorial/helloworld/)
    # Mu is the signal strength, gamma_b is a multiplicative factor calculated by the background uncertainty given an uncorrelated background
    # Here by default the model automatically calculates what mu and gamma_b should be
    init_pars = model.config.suggested_init()
    model.expected_actualdata(init_pars)

    # Check the background only model
    bkg_pars = background_only_parameters(model, init_pars)
    model.expected_actualdata(bkg_pars)

    # Compute the loglikelihood
    compute_log_likelihood(model, args.n_observed, init_pars, bkg_pars)

if __name__ == "__main__":
    main()