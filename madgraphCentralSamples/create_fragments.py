from pathlib import Path
import shutil


template = "LLgluino_fragmet_forRequests.py"
directory = "fragments"

GLUINO_MASSES = [1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
TAUS = ["100ps", "300ps", "1ns", "3ns", "10ns", "30ns"]
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
QCUT =    {1000: 140,
           1200: 143,
           1400: 147,
           1600: 150,
           1800: 156,
           2000: 156,
           2200: 160,
           2400: 162,
           2600: 162,
           2800: 168}

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

def main():
    base_directory = Path(__file__).resolve().parent
    template_path = base_directory / template
    output_directory = base_directory / directory
    output_directory.mkdir(exist_ok=True)

    for gluino_mass in GLUINO_MASSES:
        neutralino_masses = [100, gluino_mass - 100]
        for neutralino_mass in neutralino_masses:
            for tau in TAUS:
                for decay in DECAYS:
                    if skipSample(gluino_mass, neutralino_mass, tau, decay):
                        continue
                    new_filename = f"HSCP-Gluino_Par-M-{gluino_mass}-tau-{tau}-chi10M-{neutralino_mass}-{decay}_TuneCP5_13p6TeV_madgraphMLM-pythia8_fragment.py"
                    new_file_path = output_directory / new_filename
                    shutil.copy2(template_path, new_file_path)

                    content = new_file_path.read_text()
                    content = content.replace("gluino_mass = 1000", f"gluino_mass = {gluino_mass}", 1)
                    content = content.replace('tau = "100ps"', f'tau = "{tau}"', 1)
                    content = content.replace("neutralino_mass = 100", f"neutralino_mass = {neutralino_mass}", 1)
                    content = content.replace('decay = "uds"', f'decay = "{decay.removesuffix("QuarkDecay")}"', 1)
                    content = content.replace("JetMatching:qCut = 140.", f"JetMatching:qCut = {QCUT[gluino_mass]}.")
                    new_file_path.write_text(content)

if __name__ == "__main__":
    main()
