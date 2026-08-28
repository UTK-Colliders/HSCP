# Run inside the same pyHF venv as pyhf_stat_fit.py:
# source pyHF/bin/activate

import argparse
import csv
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from pyHF_statFit import StatFit
import mplhep as hep
from scipy.interpolate import griddata

TAUS = [0.1, 0.3, 1, 3, 10, 30]
GLUINO_MASSES = [1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]

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

def compute_upper_limits(signal_points, plot=True, do_debug_print=False):
    obs_limits = {}
    exp_limits_median = {}
    exp_limits_full = {}

    os.makedirs("pdfs", exist_ok=True)

    for point in signal_points:
        gluinoMass = point["gluinoMass"]
        tau = point["tau"]
        neutralinoMass = point["neutralinoMass"]
        quarkDecay = point["quarkDecay"]
        print(f"Running hypothesis test for gluinoMass={gluinoMass:.0f} GeV, tau={tau} ns")

        stat_fit = StatFit(point["n_signal"], point["n_background"], point["background_uncertainty"], gluinoMass, neutralinoMass, tau, quarkDecay, do_debug_print)
        observations = [point["n_observed"]] + stat_fit.model.config.auxdata

        poi_values, results, obs_limit, exp_limits = stat_fit.set_limits(observations)

        if plot:
            brazil_output = f"pdfs/gluinoMass{gluinoMass:.0f}_tau{tau:g}_neutralinoMass{neutralinoMass:.0f}_quarkDecay{quarkDecay}_limit.pdf"
            stat_fit.plot_limits(poi_values, results, brazil_output, gluinoMass)

        obs_limits[(tau, gluinoMass, neutralinoMass, quarkDecay)] = obs_limit
        exp_limits_median[(tau, gluinoMass, neutralinoMass, quarkDecay)] = exp_limits[2]
        exp_limits_full[(tau, gluinoMass, neutralinoMass, quarkDecay)] = exp_limits

    return obs_limits, exp_limits_median, exp_limits_full

def filter_limits(limits_dict, predicate, fixed_tau=None, tau_tol=1e-6):
    """
    Reduce a {(tau, gluinoMass, neutralinoMass, quarkDecay): value} dict down to a
    {(tau, gluinoMass): value} dict, keeping only entries where
    predicate(gluinoMass, neutralinoMass, quarkDecay) is True.

    If fixed_tau is provided, reduce instead to a {gluinoMass: value} dict at that
    lifetime.

    predicate receives floats for gluinoMass/neutralinoMass and a string for
    quarkDecay, so conditions like a fixed neutralino mass, or a fixed mass
    splitting (gluinoMass - neutralinoMass), can both be expressed here.

    value can be a scalar (e.g. obs_limits) or an array (e.g. exp_limits_full,
    the full [-2,-1,0,+1,+2 sigma] band) -- filter_limits just carries it through.
    """
    filtered = {}
    for (tau, gluinoMass, neutralinoMass, quarkDecay), value in limits_dict.items():
        if predicate(gluinoMass, neutralinoMass, quarkDecay) and (
            fixed_tau is None or np.isclose(tau, fixed_tau, atol=tau_tol)
        ):
            key = (tau, gluinoMass) if fixed_tau is None else gluinoMass
            if key in filtered:
                if fixed_tau is None:
                    raise ValueError(
                        f"Duplicate (tau={tau}, gluinoMass={gluinoMass}) after filtering -- "
                        f"the predicate isn't narrowing down to a single neutralinoMass/quarkDecay "
                        f"per point, so the 2D grid would be ambiguous."
                    )
                raise ValueError(
                    f"Duplicate (tau~={fixed_tau}, gluinoMass={gluinoMass}) after filtering -- "
                    f"the predicate isn't narrowing down to a single neutralinoMass/quarkDecay "
                    f"per point, so the 1D mass scan would be ambiguous."
                )
            filtered[key] = value
    return filtered

def plot_xsec_limit_vs_mass(mass_to_results, xsecDict, tau, quarkDecay, massSplitting, output="xsec_limit_vs_mass.pdf"):
    """
    mass_to_results: dict mapping gluino_mass -> (obs_limit_mu, exp_limits_mu array-like of 5 values [-2,-1,0,+1,+2 sigma])
    """
    masses = sorted(mass_to_results.keys())

    obs_xsec = []
    exp_xsec = {sigma: [] for sigma in range(-2, 3)}
    theory_xsec = []
    theory_xsec_unc = []

    for m in masses:
        obs_mu, exp_mu = mass_to_results[m]
        xsec, xsec_unc = xsecDict[m]
        theory_xsec.append(xsec)
        theory_xsec_unc.append(xsec_unc)
        obs_xsec.append(obs_mu)
        for sigma_index, sigma in enumerate(range(-2, 3)):
            exp_xsec[sigma].append(exp_mu[sigma_index])

    masses = np.array(masses)
    theory_xsec = np.array(theory_xsec)
    theory_xsec_unc = np.array(theory_xsec_unc)

    if quarkDecay == "uds":
        decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + (u\bar{{u}},d\bar{{d}},s\bar{{s}})$"
    else:
        decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + t\bar{{t}}$"
    if massSplitting == "small":
        neutralinoMassString = fr'$m_{{\tilde{{\chi}}^0_1}}=m_{{\tilde{{g}}}}-100$ [GeV]'
    else:
        neutralinoMassString = fr'$m_{{\tilde{{\chi}}^0_1}}=100$ [GeV]'

    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)

    plt.text(.55, .99, fr"$\tau_{{\tilde{{g}}}}={tau}$ [ns]", ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)    
    plt.text(.55, .94, decayString, ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)
    plt.text(.55, .89, neutralinoMassString, ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)

    # Fake results for testing warning
    plt.text(0.01, 0.94, "Fake results for testing", transform=ax.transAxes, fontname="TeX Gyre Heros", size=20, fontweight='bold')

    # Expected band (2 sigma, then 1 sigma on top)
    ax.fill_between(masses, exp_xsec[-2], exp_xsec[2], color="#ffcc00", label="Expected ± 2σ")
    ax.fill_between(masses, exp_xsec[-1], exp_xsec[1], color="#228b22", label="Expected ± 1σ")
    ax.plot(masses, exp_xsec[0], color="black", linestyle="--", label="Expected limit")

    # Observed limit
    ax.plot(masses, obs_xsec, color="black", marker="o", label="Observed limit")

    # Theory curve with its own uncertainty band
    ax.plot(masses, theory_xsec, color="red", linewidth=2, label="Theory")
    ax.fill_between(masses, theory_xsec - theory_xsec_unc, theory_xsec + theory_xsec_unc,
                    color="red", alpha=0.2)

    ax.set_yscale("log")
    ax.set_xlabel(r"$m_{\tilde{g}}$ [GeV]")
    ax.set_ylabel("95% CL upper limit on σ [fb]")
    ax.set_xticks(masses)
    # CMS specific
    hep.cms.label("Work in progress", loc=0, ax=ax, fontsize=18, fontproperties="TeX Gyre Heros:italic", com=13.6, lumi=283.8)
    ax.legend(loc="upper right")

    plt.savefig(output)

def extract_sigma_dict(exp_full_dict, sigma_index):
    """
    exp_full_dict maps keys -> the full 5-element [-2,-1,0,+1,+2 sigma]
    expected-limit array (as stored in exp_limits_full). Pulls out one sigma
    level as a plain scalar dict, so it can be run through build_single_grid()
    exactly like obs_limits/exp_limits_median.
    """
    return {key: value[sigma_index] for key, value in exp_full_dict.items()}

def build_single_grid(limits_dict):
    """Build a (len(GLUINO_MASSES) x len(TAUS)) grid from a
    {(tau, gluinoMass): value} dict -- see filter_limits(). Missing
    combinations are left as 0, matching the "missing scan point" convention
    used everywhere else in this script."""
    Z = np.zeros((len(GLUINO_MASSES), len(TAUS)))
    for i, gluinoMass in enumerate(GLUINO_MASSES):
        for j, tau in enumerate(TAUS):
            try:
                Z[i, j] = limits_dict[(tau, gluinoMass)]
            except KeyError:
                continue
    return Z

def build_grid(obs_limits, exp_limits_median):
    """obs_limits/exp_limits_median here are already filtered down to
    {(tau, gluinoMass): value} — see filter_limits()."""
    Z_obs = build_single_grid(obs_limits)
    Z_exp = build_single_grid(exp_limits_median)
    return TAUS, GLUINO_MASSES, Z_obs, Z_exp

def interpolate_upperlimits(*Zs, n_points=500):
    """
    Upsample any number of (gluinoMass x tau) grids -- all sharing the
    GLUINO_MASSES/TAUS axes -- onto a common fine n_points x n_points grid,
    for smooth pcolormesh/contour rendering.

    Both the interpolation coordinate for tau and the interpolated quantity
    itself are handled in log space:
      - tau is interpolated/queried in log10(tau), since the plot displays tau
        on a log axis and the limits vary smoothly in log(lifetime).
      - the values themselves are interpolated as log10(value), then
        exponentiated back, since they span multiple orders of magnitude and
        linear interpolation across that range behaves poorly.

    Returns (taus_fine_grid, gluinoMasses_fine_grid, [interpolated_Z, ...]),
    with the interpolated list in the same order as the Zs passed in.
    """
    taus = np.array(TAUS)
    gluinoMasses = np.array(GLUINO_MASSES)

    log_taus = np.log10(taus)
    gluinoMasses_fine = np.linspace(gluinoMasses.min(), gluinoMasses.max(), n_points)
    log_taus_fine = np.linspace(log_taus.min(), log_taus.max(), n_points)

    gluinoMasses_grid, log_taus_grid = np.meshgrid(gluinoMasses, log_taus, indexing="ij")
    gluinoMasses_fine_grid, log_taus_fine_grid = np.meshgrid(
        gluinoMasses_fine, log_taus_fine, indexing="ij"
    )
    taus_fine_grid = 10 ** log_taus_fine_grid

    # Zero entries represent missing scan points, not limits.
    def interpolate(values):
        valid = values > 0
        points = (gluinoMasses_grid[valid], log_taus_grid[valid])
        log_values = np.log10(values[valid])
        interpolated_log = griddata(
            points,
            log_values,
            (gluinoMasses_fine_grid, log_taus_fine_grid),
            method="linear",
        )
        return 10 ** interpolated_log

    interpolated = [interpolate(Z) for Z in Zs]
    return taus_fine_grid, gluinoMasses_fine_grid, interpolated

def plot_massVstau_exclusion(taus, gluinoMasses, massSplitting, quarkDecay,
                              Z_obs, Z_obs_theory_up, Z_obs_theory_down,
                              Z_exp, Z_exp_m1sigma, Z_exp_p1sigma, output):
    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_xlabel(r"$\tau$ [ns]")
    ax.set_ylabel(r"$m_{\tilde{g}}$ [GeV]")
 
    if quarkDecay == "uds":
        decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + (u\bar{{u}},d\bar{{d}},s\bar{{s}})$"
    else:
        decayString = fr"$\tilde{{g}} \rightarrow \tilde{{\chi}}^0_1 + t\bar{{t}}$"
    if massSplitting == "small":
        neutralinoMassString = fr'$m_{{\tilde{{\chi}}^0_1}}=m_{{\tilde{{g}}}}-100$ [GeV]'
    else:
        neutralinoMassString = fr'$m_{{\tilde{{\chi}}^0_1}}=100$ [GeV]'
 
    plt.text(.01, .91, decayString, ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)
    plt.text(.01, .86, neutralinoMassString, ha='left', va='top', transform=ax.transAxes, fontname="TeX Gyre Heros", size=10)

    # Fake results for testing warning
    plt.text(0.01, 0.94, "Fake results for testing", transform=ax.transAxes, fontname="TeX Gyre Heros", size=20, fontweight='bold')

    zmin, zmax = 1e-2, 1e3
    norm = LogNorm(vmin=zmin, vmax=zmax)
    mesh = ax.pcolormesh(taus, gluinoMasses, Z_obs, shading="nearest", norm=norm, cmap="viridis")
 
    cbar_ticks = np.logspace(np.log10(zmin), np.log10(zmax), num=6)
    cbar = fig.colorbar(mesh, ax=ax, label=r"95% CL upper limit on $\sigma$ [fb]", ticks=cbar_ticks)
    cbar.ax.yaxis.set_major_formatter(ScalarFormatter())
    cbar.ax.yaxis.set_minor_formatter(plt.NullFormatter())

    # Nominal contours (mu=1, i.e. obs_xsec == theory_xsec)
    ax.contour(taus, gluinoMasses, Z_obs, levels=[1.0], colors="black", linewidths=2)
    ax.contour(taus, gluinoMasses, Z_exp, levels=[1.0], colors="red", linestyles="dashed", linewidths=2)

    # +-1 sigma_theory band: since mu = obs_xsec/theory_xsec, the crossing
    # against a theory prediction shifted by its own +-1 sigma uncertainty is
    # equivalent to contouring obs_mu/(1+-relunc) at the same level=1.0 --
    # algebraically the same as finding where the observed cross section
    # crosses the upper/lower edge of the theory uncertainty band.
    ax.contour(taus, gluinoMasses, Z_obs_theory_up, levels=[1.0], colors="black", linewidths=1)
    ax.contour(taus, gluinoMasses, Z_obs_theory_down, levels=[1.0], colors="black", linewidths=1)

    # +-1 sigma_experiment band: the -1/+1 sigma quantiles of the expected
    # limit distribution (pulled straight from pyhf's exp_limits array), not
    # an approximation -- these are exact, unlike the theory band above.
    ax.contour(taus, gluinoMasses, Z_exp_m1sigma, levels=[1.0], colors="red", linestyles="dashed", linewidths=1)
    ax.contour(taus, gluinoMasses, Z_exp_p1sigma, levels=[1.0], colors="red", linestyles="dashed", linewidths=1)

    # Composite thin/thick/thin legend handles, matching the reference plot's
    # style, using HandlerTuple(ndivide=None) to overlay all three proxy lines
    # in the same legend icon rather than splitting the icon into segments.
    obs_thin_a = Line2D([0], [0], color="black", linewidth=1)
    obs_thick = Line2D([0], [0], color="black", linewidth=2)
    obs_thin_b = Line2D([0], [0], color="black", linewidth=1)
    exp_thin_a = Line2D([0], [0], color="red", linewidth=1, linestyle="dashed")
    exp_thick = Line2D([0], [0], color="red", linewidth=2, linestyle="dashed")
    exp_thin_b = Line2D([0], [0], color="red", linewidth=1, linestyle="dashed")

    legend_handles = [
        (obs_thin_a, obs_thick, obs_thin_b),
        (exp_thin_a, exp_thick, exp_thin_b),
    ]
    legend_labels = [
        r"Observed $\pm 1\sigma_{theory}$",
        r"Expected $\pm 1\sigma_{experiment}$",
    ]
    ax.legend(legend_handles, legend_labels, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="lower right")
 
    ax.set_xlim(0.1, 30)
    ax.set_xscale("log")
    ax.set_yticks(GLUINO_MASSES)

    # CMS specific
    hep.cms.label("Work in progress", loc=0, ax=ax, fontsize=18, fontproperties="TeX Gyre Heros:italic", com=13.6, lumi=283.8)
 
    plt.savefig(output)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Scan gluino mass/lifetime signal points and build a 2D exclusion plot")
    parser.add_argument("--csv_input", default="sampleCounts.csv", help="CSV file with signal points to scan")
    parser.add_argument("--createIndividualPlots", type=bool, default=True)
    parser.add_argument("--exclusionPlotOutputDir", default="pdfs", help="Directory to write the three exclusion plots into")
    parser.add_argument("--do_debug_print", type=bool, default=False)
    args = parser.parse_args()

    signal_points = read_signal_points(args.csv_input)
    obs_limits, exp_limits_median, exp_limits_full = compute_upper_limits(
        signal_points, args.createIndividualPlots, args.do_debug_print
    )

    os.makedirs(args.exclusionPlotOutputDir, exist_ok=True)

    # Create exclusion as a function of mass plot. One for each lifetime and decay
    mass_split = 100.0
    fixed_neutralino_mass = 100.0
    all_taus = TAUS
    mass_scan_scenarios = [
        (
            "neutralino100_uds",
            "large",
            "uds",
            lambda gm, nm, qd: np.isclose(nm, fixed_neutralino_mass) and qd == "uds",
        ),
        (
            "neutralino100_ttbar",
            "large",
            "ttbar",
            lambda gm, nm, qd: np.isclose(nm, fixed_neutralino_mass) and qd == "ttbar",
        ),
        (
            "massSplit100_uds",
            "small",
            "uds",
            lambda gm, nm, qd: np.isclose(gm - nm, mass_split) and qd == "uds",
        ),
    ]

    # Grab the shared theory cross section dictionary from any valid StatFit point.
    first_point = signal_points[0]
    xsec_dict = StatFit(
        first_point["n_signal"],
        first_point["n_background"],
        first_point["background_uncertainty"],
        first_point["gluinoMass"],
        first_point["neutralinoMass"],
        first_point["tau"],
        first_point["quarkDecay"],
        False,
    ).xsecDict

    # Relative theory cross-section uncertainty per gluino mass, in the same
    # order as GLUINO_MASSES, used below to build the observed +-1 sigma_theory
    # contour bands.
    rel_theory_unc = np.array([xsec_dict[m][1] / xsec_dict[m][0] for m in GLUINO_MASSES])

    for tau in all_taus:
        for scenario_label, massSplitting, quarkDecay, predicate in mass_scan_scenarios:
            filtered_obs = filter_limits(obs_limits, predicate, fixed_tau=tau)
            filtered_exp_full = filter_limits(exp_limits_full, predicate, fixed_tau=tau)

            masses = sorted(set(filtered_obs.keys()) & set(filtered_exp_full.keys()))
            if not masses:
                continue

            mass_to_results = {
                mass: (filtered_obs[mass], filtered_exp_full[mass]) for mass in masses
            }

            output_filename = f"exclusion_tau{tau:g}_{scenario_label}.pdf"
            output_path = os.path.join(args.exclusionPlotOutputDir, output_filename)
            plot_xsec_limit_vs_mass(mass_to_results, xsec_dict, tau, quarkDecay, massSplitting, output_path)
            print(f"Wrote {output_path}")

    # Create 2D exclusion as a function of mass and lifetime. One for each decay
    # Each entry: (output filename, plot title, predicate over (gluinoMass, neutralinoMass, quarkDecay))
    scenarios = [
        (
            "exclusion_neutralino100_uds.pdf",
            "large",
            "uds",
            lambda gm, nm, qd: np.isclose(nm, fixed_neutralino_mass) and qd == "uds",
        ),
        (
            "exclusion_neutralino100_ttbar.pdf",
            "large",
            "ttbar",
            lambda gm, nm, qd: np.isclose(nm, fixed_neutralino_mass) and qd == "ttbar",
        ),
        (
            "exclusion_massSplit100_uds.pdf",
            "small",
            "uds",
            lambda gm, nm, qd: np.isclose(gm - nm, mass_split) and qd == "uds",
        ),
    ]

    for filename, massSplitting, quarkDecay, predicate in scenarios:
        filtered_obs = filter_limits(obs_limits, predicate)
        filtered_exp_full = filter_limits(exp_limits_full, predicate)
        filtered_exp_median = {key: value[2] for key, value in filtered_exp_full.items()}
        filtered_exp_m1 = {key: value[1] for key, value in filtered_exp_full.items()}
        filtered_exp_p1 = {key: value[3] for key, value in filtered_exp_full.items()}

        Z_obs = build_single_grid(filtered_obs)
        Z_exp_median = build_single_grid(filtered_exp_median)
        Z_exp_m1 = build_single_grid(filtered_exp_m1)
        Z_exp_p1 = build_single_grid(filtered_exp_p1)

        # Rescale mu by (1 +- relative theory uncertainty) per mass row, then
        # contour at the same level=1.0 -- equivalent to crossing the
        # observed cross section against the shifted theory band edges.
        Z_obs_theory_up = np.where(Z_obs > 0, Z_obs / (1 + rel_theory_unc[:, None]), 0.0)
        Z_obs_theory_down = np.where(Z_obs > 0, Z_obs / (1 - rel_theory_unc[:, None]), 0.0)

        taus_fine, gluinoMasses_fine, interpolated = interpolate_upperlimits(
            Z_obs, Z_obs_theory_up, Z_obs_theory_down, Z_exp_median, Z_exp_m1, Z_exp_p1
        )
        Z_obs_i, Z_obs_up_i, Z_obs_down_i, Z_exp_i, Z_exp_m1_i, Z_exp_p1_i = interpolated

        output_path = os.path.join(args.exclusionPlotOutputDir, filename)
        plot_massVstau_exclusion(
            taus_fine, gluinoMasses_fine, massSplitting, quarkDecay,
            Z_obs_i, Z_obs_up_i, Z_obs_down_i, Z_exp_i, Z_exp_m1_i, Z_exp_p1_i,
            output_path,
        )
        print(f"Wrote {output_path}")

if __name__ == "__main__":
    main()