# First run `/usr/bin/python3.12 -m venv pyHF`, followed by `source pyHF/bin/activate` to access the pyhf library
# For your first time, run `pip install pyhf`
# Run `deactivate` when you are done with the venv

import pyhf
import argparse


def build_model(n_signal, n_background, background_uncertainty, do_debug_print):
    model = pyhf.simplemodels.uncorrelated_background(
        signal=[n_signal], bkg=[n_background], bkg_uncertainty=[background_uncertainty]
    )

    if do_debug_print:
        debug_print(model)

    return model

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
    parser.add_argument("--n_signal", type=float, default=1.0)
    parser.add_argument("--n_background", type=float, default=1.0)
    parser.add_argument("--background_uncertainty", type=float, default=1.0)
    parser.add_argument("--do_debug_print", type=bool, default=True)
    args = parser.parse_args()

    model = build_model(args.n_signal, args.n_background, args.background_uncertainty, args.do_debug_print)

    # For a single channel, the parameters are [mu, gamma_b] (see https://pyhf.github.io/pyhf-tutorial/helloworld/)
    # Mu is the signal strength, gamma_b is a multiplicative factor calculated by the background uncertainty given an uncorrelated background
    # Here by default the model automatically calculates what mu and gamma_b should be
    init_pars = model.config.suggested_init()
    model.expected_actualdata(init_pars)

if __name__ == "__main__":
    main()