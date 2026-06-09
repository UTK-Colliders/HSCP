"""
Blind an HSCP data ntuple: drop any event in which ANY track simultaneously satisfies ALL THREE
signal-region cuts.

    IsoTrack_pt       > 200 GeV   -- hard HSCP candidate pT threshold
    DeDx_FiPixelNoL1  > 0.9       -- fine-pixel dE/dx, no L1 correction ("fpix")
    DeDx_GiStrip      > 0.45      -- Gi estimator on strip detector ("gstips")

Selection (keep an event only if NO track passes all three):
    Sum$(IsoTrack_pt > 200 && DeDx_FiPixelNoL1 > 0.9 && DeDx_GiStrip > 0.45) == 0
        - 3 tracks, 1 passes all three -> Sum$()=1 -> fails ==0 -> event DROPPED
        - 3 tracks, 0 pass  all three -> Sum$()=0 -> passes ==0 -> event KEPT
"""

import os
import time
import argparse
import ROOT
INPUT_FILE  = ("/eos/uscms/store/user/avendras/HSCP/data/2024C/JetMET1_Run2024C_Collisions2024_378981_386951_Golden/JetMET1/JetMET1_Run2024C-MINIv6NANOv15-v1_MINIAOD_Cert_Collisions2024_378981_386951_Golden/260604_180557/0000/2024C_JETMet1_unblinded.root")
TREE_PATH   = "HSCPMiniAODAnalyzer/Events"
OUTPUT_FILE = "/eos/uscms/store/user/avendras/HSCP/data/2024C/JetMET1_blinded.root"

SIGNAL_REGION_VETO = ("Sum$(IsoTrack_pt > 200  && DeDx_FiPixelNoL1 > 0.9 && DeDx_GiStrip > 0.45) == 0")
MAX_TREE_SIZE = 1_000_000_000_000  # bytes (1 TB)

def _open_with_retry(path, mode, tries=20, wait=30):
    last = ""
    for attempt in range(1, tries + 1):
        f = ROOT.TFile.Open(path, mode)
        if f and not f.IsZombie():
            if attempt > 1:
                print(f"  opened {path} on attempt {attempt}")
            return f
        if f:
            f.Close()
        last = f"attempt {attempt}/{tries} failed"
        if attempt < tries:
            print(f"  {last} opening {path} -- retrying in {wait}s ...")
            time.sleep(wait)
    raise RuntimeError(f"could not open {path} after {tries} attempts ({last})")


def run_blinder(input_file=INPUT_FILE, output_file=OUTPUT_FILE, tree_path=TREE_PATH):
    ROOT.TTree.SetMaxTreeSize(MAX_TREE_SIZE)
    if not output_file.startswith("root://"):
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    src = _open_with_retry(input_file, "READ")
    in_tree = src.Get(tree_path)
    if not in_tree:
        raise RuntimeError(f"tree '{tree_path}' not found in {input_file}")
    n_before = in_tree.GetEntries()

    print(f"Input:  {n_before} events")
    print(f"File:   {input_file}")
    print(f"\nkill cuts (remove event if ANY track passes ALL):")
    print(f"  IsoTrack_pt       > 200 GeV")
    print(f"  DeDx_FiPixelNoL1  > 0.9")
    print(f"  DeDx_GiStrip      > 0.45")
    print(f"\nexcluding events if any track passes:")
    print(f"  {SIGNAL_REGION_VETO}\n")

    out = _open_with_retry(output_file, "RECREATE")
    out.cd()
    print("Copying selected events (PyROOT CopyTree) ...")
    out_tree = in_tree.CopyTree(SIGNAL_REGION_VETO)
    n_after = out_tree.GetEntries()
    out_tree.Write()
    out.Close()
    src.Close()

    n_blinded = n_before - n_after
    frac_blinded = n_blinded / n_before if n_before > 0 else 0.0
    print(f"\nOutput: {n_after} events  ({n_blinded} blinded, {frac_blinded:.4%} of input)")
    print(f"File:   {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Blind an HSCP data ntuple (signal-region veto).")
    parser.add_argument("--input",  default=INPUT_FILE,  help="input unblinded .root (local path or root:// URL)")
    parser.add_argument("--output", default=OUTPUT_FILE, help="output blinded .root (local path or root:// URL)")
    parser.add_argument("--tree",   default=TREE_PATH,   help="path to the Events tree inside the input file")
    args = parser.parse_args()
    run_blinder(args.input, args.output, args.tree)