# First run `source pyHF/bin/activate` to access the pyhf library
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
    print(f"   auxdata: {model.config.auxdata}")

def main():
    parser = argparse.ArgumentParser(description="Perform a single-channel statistical fit with uncorrelated background")
    parser.add_argument("--n_signal", default=1.0)
    parser.add_argument("--n_background", default=1.0)
    parser.add_argument("--background_uncertainty", default=1.0)
    parser.add_argument("--do_debug_print", default=False)
    args = parser.parse_args()

    model = build_model(args.n_signal, args.n_background, args.background_uncertainty, args.do_debug_print)

if __name__ == "__main__":
    main()