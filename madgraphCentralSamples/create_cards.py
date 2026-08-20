#!/usr/bin/env python3


from __future__ import annotations

import re
import shutil
from pathlib import Path


GLUINO_MASSES = [1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
TAUS = ["100ps", "300ps", "1ns", "3ns", "10ns", "30ns"]
WIDTHS = ["6.58211957E-15", "2.19403986E-15", "6.58211957E-16", "2.19403986E-16", "6.58211957E-17", "2.19403986E-17"]
DECAYS = ["udsQuarkDecay", "ttbarQuarkDecay"]
NEVENTS = {1000: 200000,
           1200: 100000,
           1400: 100000,
           1600: 20000,
           1800: 20000,
           2000: 20000,
           2200: 20000,
           2400: 20000,
           2600: 20000,
           2800: 20000}

SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_PARAM_CARD = SCRIPT_DIR / "LLgluinoTemplate_param_card.dat"
TEMPLATE_RUN_CARD = SCRIPT_DIR / "LLgluinoTemplate_run_card.dat"
TEMPLATE_PROC_CARD = SCRIPT_DIR / "LLgluinoTemplate_proc_card.dat"
OUTPUT_DIR = Path("genproductions_scripts/bin/MadGraph5_aMCatNLO/InputCards")
GLUINO_TEMPLATE_MASS = 1000.0
NEUTRALINO_TEMPLATE_MASS = 100.0

def replace_gluino_mass(line: str, gluino_mass: int) -> str:
    return re.sub(
        r"^(\s*1000021\s+)([0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?\d+)?)(\s+# ~g\b.*)$",
        rf"\g<1>{gluino_mass:.2f}\g<3>",
        line,
    )


def replace_rhadron_masses(line: str, gluino_mass: int) -> str:
    match = re.match(r"^(\s*-?\d+\s+)([0-9]+)(\.[0-9]+(?:[Ee][+-]?\d+)?)(\s+# ~g_.*)$", line)
    if not match:
        return line
    prefix, integer_part, fractional_part, suffix = match.groups()
    shifted_integer = int(integer_part) + (gluino_mass - int(GLUINO_TEMPLATE_MASS))
    line_ending = "\n" if line.endswith("\n") else ""
    return f"{prefix}{shifted_integer}{fractional_part}{suffix}{line_ending}"


def replace_neutralino_mass(line: str, neutralino_mass: int) -> str:
    match = re.match(r"^(\s*-?\d+\s+)([0-9]+)(\.[0-9]+(?:[Ee][+-]?\d+)?)(\s+# ~chi_[1]0\b.*)$", line)
    if not match:
        return line
    prefix, integer_part, fractional_part, suffix = match.groups()
    shifted_integer = int(integer_part) + (neutralino_mass - int(NEUTRALINO_TEMPLATE_MASS))
    line_ending = "\n" if line.endswith("\n") else ""
    return f"{prefix}{shifted_integer}{fractional_part}{suffix}{line_ending}"


def replace_gluino_width(line: str, width: str) -> str:
    return re.sub(
        r"^DECAY\s+1000021\s+[0-9.E+-]+\s+# gluino decays$",
        f"DECAY   1000021     {width}   # gluino decays",
        line,
    )


def replace_uds_branching_ratios(line: str) -> str:
    replacements = {
        "# BR(~g -> ~chi_10 d db)": "3.33333333E-01",
        "# BR(~g -> ~chi_10 u ub)": "3.33333333E-01",
        "# BR(~g -> ~chi_10 s sb)": "3.33333333E-01",
        "# BR(~g -> ~chi_10 t tb)": "0.00000000E+00",
    }
    for marker, value in replacements.items():
        if marker in line:
            return re.sub(r"^\s*[0-9.E+-]+", f"     {value}", line, count=1)
    return line


def update_param_card(param_card_path: Path, gluino_mass: int, neutralino_mass: int, width: str, decay: str) -> None:
    updated_lines = []
    for line in param_card_path.read_text().splitlines(keepends=True):
        line = replace_gluino_mass(line, gluino_mass)
        line = replace_rhadron_masses(line, gluino_mass)
        line = replace_neutralino_mass(line, neutralino_mass)
        line = replace_gluino_width(line, width)
        if decay == "udsQuarkDecay":
            line = replace_uds_branching_ratios(line)
        updated_lines.append(line)
    param_card_path.write_text("".join(updated_lines))


def main() -> None:
    if not TEMPLATE_PARAM_CARD.exists():
        raise FileNotFoundError(f"Missing template card: {TEMPLATE_PARAM_CARD}")
    if not TEMPLATE_RUN_CARD.exists():
        raise FileNotFoundError(f"Missing template card: {TEMPLATE_RUN_CARD}")
    if not TEMPLATE_PROC_CARD.exists():
        raise FileNotFoundError(f"Missing template card: {TEMPLATE_PROC_CARD}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    COUNTER=0

    for gluino_mass in GLUINO_MASSES:
        neutralino_masses = [100, gluino_mass - 100]
        for neutralino_mass in neutralino_masses:
            for tau in TAUS:
                for decay in DECAYS:
                    name = f"LLgluino_M-{gluino_mass}_tau-{tau}_{decay}_chi10_M-{neutralino_mass}"
                    width_index = TAUS.index(tau)
                    WIDTH = WIDTHS[width_index]
                    if skipSample(gluino_mass, neutralino_mass, tau, decay):
                        continue

                    sample_dir = OUTPUT_DIR / name
                    sample_dir.mkdir(parents=True, exist_ok=True)
                    COUNTER+=1

                    param_card = sample_dir / f"{name}_param_card.dat"
                    run_card = sample_dir / f"{name}_run_card.dat"
                    proc_card = sample_dir / f"{name}_proc_card.dat"

                    shutil.copyfile(TEMPLATE_PARAM_CARD, param_card)
                    shutil.copyfile(TEMPLATE_RUN_CARD, run_card)
                    shutil.copyfile(TEMPLATE_PROC_CARD, proc_card)

                    update_param_card(param_card, gluino_mass, neutralino_mass, WIDTH, decay)

                    proc_text = proc_card.read_text()
                    proc_text = re.sub(
                        r"^output\s+LLgluino_M-1000_tau-1ns_ttbarQuarkDecay_chi10_M-100\s+-nojpeg\s*$",
                        f"output {name} -nojpeg",
                        proc_text,
                        flags=re.MULTILINE,
                    )
                    proc_card.write_text(proc_text)

                    run_text = run_card.read_text()
                    run_text = re.sub(
                        r"^LLgluino_M-1000_tau-1ns_ttbarQuarkDecay_chi10_M-100\s+= run_tag ! name of the run\s*$",
                        f"{name} = run_tag ! name of the run",
                        run_text,
                        flags=re.MULTILINE,
                    )

                    if (gluino_mass == 1000 and (tau == "10ns" or tau == "30ns")):
                        nevents = 100000
                    else:
                        nevents = NEVENTS[gluino_mass]
                    run_text = re.sub(
                        r"^200000 = nevents ! Number of unweighted events requested\s*$",
                        f"{nevents} = nevents ! Number of unweighted events requested",
                        run_text,
                        flags=re.MULTILINE,
                    )

                    run_card.write_text(run_text)

                    print(f"Created {sample_dir}")
    print(f"Total samples created: {COUNTER}")

def skipSample(gluino_mass: int, neutralino_mass: int, tau: str, decay: str) -> bool:
    if ((decay == "ttbarQuarkDecay" and gluino_mass-neutralino_mass == 100) or
        (tau == "100ps" and gluino_mass >= 1400) or
        (tau == "300ps" and gluino_mass >= 2200) or
        (tau == "1ns" and gluino_mass >= 2600) or
        (tau == "3ns" and gluino_mass >= 2800) or
        (tau == "300ps" and gluino_mass == 2000 and gluino_mass-neutralino_mass == 100) or
        (tau == "1ns" and gluino_mass == 2400 and gluino_mass-neutralino_mass == 100) or
        (tau == "3ns" and gluino_mass == 2600 and gluino_mass-neutralino_mass == 100) or
        ((tau == "10ns" or tau == "30ns") and gluino_mass == 2800 and gluino_mass-neutralino_mass == 100)):
        return True
    return False


if __name__ == "__main__":
    main()
