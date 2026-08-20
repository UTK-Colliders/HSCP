# Run inside the same pyHF venv as pyhf_stat_fit.py:
# source pyHF/bin/activate

import argparse
import csv
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from pyHF_statFit import StatFit

def read_signal_points(csv_path):
    signal_points = []
    with open(csv_path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            signal_points.append({
                "gluinoMass": float(row["gluinoMass"]),
                "tau": float(row["tau"]),
                "neutralinoMass": float(row["neutralinoMass"]),
                "quarkDecay": row["quarkDecay"],
                "n_signal": float(row["n_signal"]),
                "n_background": float(row["n_background"]),
                "background_uncertainty": float(row["background_uncertainty"]),
                "n_observed": int(row["n_observed"]),
            })

    return signal_points

def compute_upper_limits(signal_points, do_debug_print=False):
    obs_limits = {}
    exp_limits_median = {}

    os.makedirs("pdfs", exist_ok=True)

    for point in signal_points:
        gluinoMass = point["gluinoMass"]
        tau = point["tau"]
        neutralinoMass = point["neutralinoMass"]
        quarkDecay = point["quarkDecay"]
        print(f"Running hypothesis test for gluinoMass={gluinoMass:.0f} GeV, tau={tau} ns")

        stat_fit = StatFit(point["n_signal"], point["n_background"], point["background_uncertainty"], do_debug_print)
        observations = [point["n_observed"]] + stat_fit.model.config.auxdata

        poi_values, results, obs_limit, exp_limits = stat_fit.set_limits(observations)

        brazil_output = f"pdfs/gluinoMass{gluinoMass:.0f}_tau{tau:g}_neutralinoMass{neutralinoMass:.0f}_quarkDecay{quarkDecay}_limit.pdf"
        stat_fit.plot_limits(poi_values, results, brazil_output)

        obs_limits[(gluinoMass, tau)] = obs_limit
        exp_limits_median[(gluinoMass, tau)] = exp_limits[2]

    return obs_limits, exp_limits_median

def build_grid(obs_limits, exp_limits_median):
    gluinoMasses = sorted(set(gluinoMass for gluinoMass, tau in obs_limits))
    taus = sorted(set(tau for gluinoMass, tau in obs_limits))

    Z_obs = np.zeros((len(taus), len(gluinoMasses)))
    Z_exp = np.zeros((len(taus), len(gluinoMasses)))

    for i, tau in enumerate(taus):
        for j, gluinoMass in enumerate(gluinoMasses):
            Z_obs[i, j] = obs_limits[(gluinoMass, tau)]
            Z_exp[i, j] = exp_limits_median[(gluinoMass, tau)]

    return gluinoMasses, taus, Z_obs, Z_exp

def plot_exclusion(gluinoMasses, taus, Z_obs, Z_exp, output):
    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_title(r"TEST ($\mu$ Upper Limit)")
    ax.set_ylabel(r"$m_{\tilde{g}}$ (GeV)")
    ax.set_xlabel(r"$\tau$ (ns)")
    ax.set_xscale("log")

    mesh = ax.pcolormesh(taus, gluinoMasses, Z_obs, shading="nearest", norm=LogNorm(), cmap="RdYlBu_r")
    fig.colorbar(mesh, ax=ax, label=r"95% CL upper limit on $\mu$")

    # mu = 1 is where the observed/expected limit equals the nominal signal prediction,
    # i.e. the exclusion boundary
    ax.contour(taus, gluinoMasses, Z_obs, levels=[1.0], colors="black", linewidths=2)
    ax.contour(taus, gluinoMasses, Z_exp, levels=[1.0], colors="red", linestyles="dashed", linewidths=2)

    # Only tick the actual gluino mass / tau values that came from the csv, not
    # whatever matplotlib would pick on its own for a log axis
    ax.set_xticks(taus)
    ax.set_yticks(gluinoMasses)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(plt.NullFormatter())

    plt.savefig(output)

def main():
    parser = argparse.ArgumentParser(description="Scan gluino mass/lifetime signal points and build a 2D exclusion plot")
    parser.add_argument("--csv_input", default="sampleCounts.csv", help="CSV file with signal points to scan")
    parser.add_argument("--exclusionPlotOutput", default="pdfs/exclusion.pdf")
    parser.add_argument("--do_debug_print", type=bool, default=False)
    args = parser.parse_args()

    signal_points = read_signal_points(args.csv_input)
    obs_limits, exp_limits_median = compute_upper_limits(signal_points, args.do_debug_print)
    gluinoMasses, taus, Z_obs, Z_exp = build_grid(obs_limits, exp_limits_median)
    plot_exclusion(gluinoMasses, taus, Z_obs, Z_exp, args.exclusionPlotOutput)

if __name__ == "__main__":
    main()