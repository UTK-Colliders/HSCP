from __future__ import annotations
import operator
import os
import awkward as ak
import numpy as np
import uproot
import uproot.source.xrootd
import ROOT
import cmsstyle as CMS
from dataclasses import dataclass
from typing import Tuple


ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

INPUT_FILES = os.environ["PLOT_INPUTS"].split() if os.environ.get("PLOT_INPUTS") else [
    "/eos/uscms/store/user/avendras/HSCP/data/2024C/JetMET0_blinded.root",
    "/eos/uscms/store/user/avendras/HSCP/data/2024C/JetMET1_blinded.root",
]
TREE_PATH  = "Events"
OUTPUT_DIR = os.environ.get(
    "PLOT_OUTDIR",
    "/uscms/home/avendras/nobackup/HSCP/scripts/custom_Framework_scripts/Data/plots_combined2024_noTrg/",
)
CHUNK_SIZE = "50 MB"    # uproot.iterate step size (decompressed). Had to add a step size.

@dataclass
class Cut:
    """One selection cut: keep tracks where `op(branch, threshold)` holds (|branch| if use_abs)."""
    branch:    str
    threshold: float
    op:        str  = ">"
    use_abs:   bool = False


@dataclass
class Plot1DSpec:
    branch:      str
    xlabel:      str
    xrange:      Tuple[float, float]
    nbins:       int  = 50
    jagged:      bool = True          # per-track branch (vs. per-event like MET)
    log_y:       bool = False
    skip_presel: bool = False         # bypass the preselection mask for this branch


@dataclass
class Plot2DSpec:
    branch_x:  str
    branch_y:  str
    xlabel:    str
    ylabel:    str
    xrange:    Tuple[float, float]
    yrange:    Tuple[float, float]
    nbins_x:   int  = 20
    nbins_y:   int  = 20
    log_z:     bool = True
    abcd:      bool = False           # use ABCD_CUTS instead of the global preselection


# The FPix=FPIX_BOUNDARY / GStrips=GSTRIPS_BOUNDARY lines split the plane into A/B/C/D. The
# high-high corner (D) is the signal region and matches the blinding veto exactly.s
FPIX_BOUNDARY     = 0.9
GSTRIPS_BOUNDARY  = 0.45
ABCD_AXES         = ("DeDx_FiPixelNoL1", "DeDx_GiStrip") 
ABCD_PT_THRESHOLD = 200.0

# Global preselection — &ed together and applied to every per-track plot. Per-event plots (jagged=False, e.g. MET) bypass it automatically.
PRESELECTION_CUTS = [
    Cut("IsoTrack_pt",                  55.0,  ">"),                  # p_T > 55 GeV
    Cut("IsoTrack_eta",                 2.4,    "<",  use_abs=True),  # |eta| < 2.4
    Cut("DeDx_PixelNoL1NOM",            2,      ">="),                # valid pixel hits L2-L4
    Cut("IsoTrack_fractionOfValidHits", 0.8,    ">"),                 # fraction of valid hits
    Cut("IsoTrack_numberOfValidHits",   10,     ">="),                # # of dE/dx measurements
    Cut("IsoTrack_isHighPurityTrack",   True,   "=="),                # high-purity track
    Cut("IsoTrack_normChi2",            5,      "<"),                 # chi2/dof < 5
    Cut("IsoTrack_dz",                  0.1,    "<",  use_abs=True),  # |dz| < 0.1 cm
    Cut("IsoTrack_dxy",                 0.02,   "<",  use_abs=True),  # |dxy| < 0.02 cm
    Cut("IsoTrack_pfMiniRelIsoAll",     0.02,   "<"),                 # I_PF^rel < 0.02
    Cut("IsoTrack_IsoSumPt_dr03",       15.0,   "<"),                 # I_trk < 15 GeV (tracker iso)
    Cut("IsoTrack_pfEnergyOverP",       0.3,    "<"),                 # PF E/p < 0.3
    Cut("IsoTrack_ptErrOverPt2",        0.0008, "<"),                 # sigma(pT)/pT^2
    Cut("DeDx_FiPixelNoL1",             0.3,    ">"),                 # F_i^Pixels > 0.3

  #TO ADD TO  PRESELECTION_CUTS
    #hscp_candidate veto
    #add triggers
    #pseudoMET

]

# ABCD selection = preselection with IsoTrack_pt tightened to ABCD_PT_THRESHOLD, so region D
# coincides with the blinded corner.
# FiPix>0.3 is used so  A/C span 0.3 < FPix < 0.9;
ABCD_CUTS = [c for c in PRESELECTION_CUTS if c.branch != "IsoTrack_pt"]
ABCD_CUTS.append(Cut("IsoTrack_pt", ABCD_PT_THRESHOLD, ">"))

# FPix preselection edge — the left boundary of sidebands A and C.
FPIX_PRESEL_CUT = next((c.threshold for c in ABCD_CUTS if c.branch == "DeDx_FiPixelNoL1"), None)

PLOT1D_SPECS = [
    # dE/dx
    Plot1DSpec("DeDx_FiPixelNoL1", "FPix",    (0, 1),  nbins=50),
    Plot1DSpec("DeDx_GiStrip",     "GStrips", (0, 1),  nbins=50),
    Plot1DSpec("DeDx_Ih",          "I_{h}",   (0, 15), nbins=50),

    # Isolated tracks
    Plot1DSpec("IsoTrack_IsoSumPt_dr03",   "#Sigma p_{T}^{iso} (#DeltaR<0.3) [GeV]", (0, 100),       nbins=50),
    Plot1DSpec("IsoTrack_PseudoTrack_eta", "pseudo-track #eta",                      (-5, 5),        nbins=50),
    Plot1DSpec("IsoTrack_PseudoTrack_p",   "pseudo-track p [GeV]",                   (0, 2000),      nbins=50),
    Plot1DSpec("IsoTrack_PseudoTrack_phi", "pseudo-track #phi [rad]",                (-4, 4),        nbins=50),
    Plot1DSpec("IsoTrack_PseudoTrack_pt",  "pseudo-track p_{T} [GeV]",               (0, 2000),      nbins=50),
    Plot1DSpec("IsoTrack_charge",          "track charge",                            (-2, 2),       nbins=50),
    Plot1DSpec("IsoTrack_dxy",             "d_{xy} [cm]",                             (-0.5, 0.5),   nbins=50),
    Plot1DSpec("IsoTrack_dxyError",        "#sigma(d_{xy}) [cm]",                     (0, 0.1),      nbins=50),
    Plot1DSpec("IsoTrack_dz",              "d_{z} [cm]",                              (-30, 30),     nbins=50),
    Plot1DSpec("IsoTrack_eta",             "#eta",                                    (-5, 5),       nbins=50),
    Plot1DSpec("IsoTrack_phi",             "#phi [rad]",                              (-4, 4),       nbins=50),
    Plot1DSpec("IsoTrack_pt",              "p_{T} [GeV]",                             (0, 2000),     nbins=50),
    Plot1DSpec("IsoTrack_dzError",         "#sigma(d_{z}) [cm]",                      (0, 1),        nbins=50),
    Plot1DSpec("IsoTrack_ptErrOverPt",     "#sigma(p_{T})/p_{T}",                     (0, 1),        nbins=50),
    Plot1DSpec("IsoTrack_ptErrOverPt2",    "#sigma(p_{T})/p_{T}^{2} [GeV^{-1}]",      (0, 0.01),     nbins=50),
    Plot1DSpec("IsoTrack_ptError",         "#sigma(p_{T}) [GeV]",                     (0, 100),      nbins=50),
    Plot1DSpec("IsoTrack_px",              "p_{x} [GeV]",                             (-2000, 2000), nbins=50),
    Plot1DSpec("IsoTrack_py",              "p_{y} [GeV]",                             (-2000, 2000), nbins=50),
    Plot1DSpec("IsoTrack_pz",              "p_{z} [GeV]",                             (-2000, 2000), nbins=50),

    # MET (per-event scalars)
    Plot1DSpec("RecoPFMET",        "Reco PF E_{T}^{miss} [GeV]",    (0, 2000), jagged=False),
    Plot1DSpec("RecoPFMET_phi",    "Reco PF E_{T}^{miss} #phi",     (-4, 4),   jagged=False),
    Plot1DSpec("RecoPuppiMET",     "Reco PUPPI E_{T}^{miss} [GeV]", (0, 2000), jagged=False),
    Plot1DSpec("RecoPuppiMET_phi", "Reco PUPPI E_{T}^{miss} #phi",  (-4, 4),   jagged=False),
    Plot1DSpec("HLTPFMET",         "HLT PF E_{T}^{miss} [GeV]",     (0, 2000), jagged=False),
    Plot1DSpec("HLTPFMET_phi",     "HLT PF E_{T}^{miss} #phi",      (-4, 4),   jagged=False),
]

PLOT2D_SPECS = [
    Plot2DSpec("DeDx_FiPixelNoL1", "DeDx_GiStrip", "FPix", "GStrips",                        (0, 1),  (0, 1),    nbins_x=40, nbins_y=40, abcd=True),
    Plot2DSpec("DeDx_FiPixelNoL1", "IsoTrack_pt",  "FPix", "IsoTrack p_{T} [GeV]",           (0, 1),  (0, 2000), nbins_x=50, nbins_y=50),
    Plot2DSpec("DeDx_FiPixelNoL1", "DeDx_Ih",      "FPix", "I_{h}",                          (0, 1),  (0, 15),   nbins_x=50, nbins_y=50),
    Plot2DSpec("DeDx_GiStrip",     "IsoTrack_pt",  "GStrips", "IsoTrack p_{T} [GeV]",        (0, 1),  (0, 2000), nbins_x=50, nbins_y=50),
    Plot2DSpec("DeDx_GiStrip",     "DeDx_Ih",      "GStrips", "I_{h}",                       (0, 1),  (0, 15),   nbins_x=50, nbins_y=50),
    Plot2DSpec("DeDx_Ih",          "IsoTrack_pt",  "I_{h} [MeV/cm]", "IsoTrack p_{T} [GeV]", (0, 15), (0, 2000), nbins_x=50, nbins_y=50),
]

_OPERATORS = {">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "==": operator.eq}

def _flatten(values) -> np.ndarray:
    """Flatten a jagged awkward array to a 1D float64 numpy array."""
    return ak.to_numpy(ak.flatten(values)).astype(np.float64)


def build_track_mask(data, cuts=PRESELECTION_CUTS):
    """
    Per-track boolean mask = AND of every cut in `cuts`.
    All IsoTrack_*/DeDx_* branches share one per-track layout (length nIsoTrack per event), so the
    same mask applies to every per-track plot. Pass ABCD_CUTS for the ABCD selection.
    """
    mask = None
    for cut in cuts:
        values = abs(data[cut.branch]) if cut.use_abs else data[cut.branch]
        passes = _OPERATORS[cut.op](values, cut.threshold)
        mask = passes if mask is None else (mask & passes)
    return mask


def branches_to_read() -> list:
    """Every branch the run touches: preselection + ABCD cuts + ABCD axes + all 1D + all 2D x/y."""
    branches = {cut.branch for cut in PRESELECTION_CUTS}
    branches.update(cut.branch for cut in ABCD_CUTS)
    branches.update(ABCD_AXES)
    branches.add("IsoTrack_pt")
    for spec in PLOT1D_SPECS:
        branches.add(spec.branch)
    for spec in PLOT2D_SPECS:
        branches.update((spec.branch_x, spec.branch_y))
    return sorted(branches)


def make_hist_1d(spec: Plot1DSpec) -> ROOT.TH1F:
    hist = ROOT.TH1F(f"h_{spec.branch}", "", spec.nbins, spec.xrange[0], spec.xrange[1])
    hist.SetDirectory(0)
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetFillStyle(0)
    hist.SetMarkerStyle(0)
    hist.SetMarkerColor(0)
    return hist

def make_hist_2d(spec: Plot2DSpec) -> ROOT.TH2F:
    hist = ROOT.TH2F(
        f"h_{spec.branch_x}_vs_{spec.branch_y}", "",
        spec.nbins_x, spec.xrange[0], spec.xrange[1],
        spec.nbins_y, spec.yrange[0], spec.yrange[1],
    )
    hist.SetDirectory(0)
    return hist


def fill_hist_1d(hist: ROOT.TH1F, spec: Plot1DSpec, data, preselection_mask):
    """Fill one chunk into a 1D hist. Per-track plots apply preselection_mask; scalars don't."""
    raw = data[spec.branch]
    if spec.jagged and not spec.skip_presel and preselection_mask is not None:
        values = _flatten(raw[preselection_mask])
    elif spec.jagged:
        values = _flatten(raw)
    else:
        values = ak.to_numpy(raw).astype(np.float64)

    values = values[np.isfinite(values)]
    if values.size:
        hist.FillN(len(values), values, np.ones(len(values), dtype=np.float64))


def fill_hist_2d(hist: ROOT.TH2F, spec: Plot2DSpec, data, track_mask) -> int:
    """Fill one chunk into a 2D hist; returns the number of tracks filled."""
    x = _flatten(data[spec.branch_x][track_mask])
    y = _flatten(data[spec.branch_y][track_mask])
    if x.size:
        hist.FillN(len(x), np.ascontiguousarray(x), np.ascontiguousarray(y),
                   np.ones(len(x), dtype=np.float64))
    return len(x)


def set_cms_style():
    CMS.SetExtraText("Preliminary")
    CMS.SetLumi(None, run="Run 3")
    CMS.SetEnergy(13.6)


def _style_axes(obj):
    obj.GetXaxis().SetTitleSize(0.045)
    obj.GetYaxis().SetTitleSize(0.045)
    obj.GetXaxis().SetLabelSize(0.040)
    obj.GetYaxis().SetLabelSize(0.040)


def _dashed_line(x1, y1, x2, y2, color, style=ROOT.kDashed) -> ROOT.TLine:
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineStyle(style)
    line.SetLineWidth(2)
    line.SetLineColor(color)
    line.Draw("SAME")
    return line


def draw_1d(spec: Plot1DSpec, hist: ROOT.TH1F, outpath: str):
    set_cms_style()
    xmin, xmax = spec.xrange
    ymin = 0.1 if spec.log_y else 0.0
    ymax = max(hist.GetMaximum(), 1.0) * (10.0 if spec.log_y else 1.35)

    canv = CMS.cmsCanvas(
        spec.branch, xmin, xmax, ymin, ymax,
        spec.xlabel, "Events", square=False, iPos=0,
    )
    if spec.log_y:
        canv.SetLogy(True)
        hist.SetMinimum(0.1)

    frame = CMS.GetCmsCanvasHist(canv)
    _style_axes(frame)
    frame.GetYaxis().SetTitleOffset(1.2)
    canv.Modified()
    canv.Update()

    CMS.cmsDraw(hist, "HIST", lwidth=2, lcolor=ROOT.kBlue + 1, fstyle=0)
    canv.SaveAs(outpath)
    canv.Close()


def _abcd_label(x, y, text, size=0.042):
    """
    TLatex centered at (x, y), drawn as a black 'shadow' + white text so it stays legible
    over any COLZ palette color. Returns both TLatex objects (keep alive until SaveAs).
    """
    objs = []
    for dx, dy, color in ((0.005, -0.005, ROOT.kBlack), (0.0, 0.0, ROOT.kWhite)):
        latex = ROOT.TLatex(x + dx, y + dy, text)
        latex.SetTextAlign(22)          # horizontally + vertically centered
        latex.SetTextFont(62)           # bold
        latex.SetTextSize(size)
        latex.SetTextColor(color)
        latex.Draw()
        objs.append(latex)
    return objs


def _draw_abcd_regions(spec: Plot2DSpec, abcd: "AbcdCount"):
    """
    Mark the four ABCD quadrants (letter + candidate count) on the FPix-vs-GStrips plot, with
    the B*C/A prediction for the blinded D region in a header. Returns kept-alive TObjects.
    Quadrant centers sit midway between each boundary and the axis edge, so they follow
    FPIX_BOUNDARY / GSTRIPS_BOUNDARY automatically.
    """
    gstr_lo = (spec.yrange[0] + GSTRIPS_BOUNDARY) / 2.0
    # A/C live in the sideband 0.3 < FPix < FPIX_BOUNDARY; center the label in that band so it sits
    # to the RIGHT of the FiPix>0.3 preselection line rather than on top of it.
    fpix_ac_left = FPIX_PRESEL_CUT if FPIX_PRESEL_CUT is not None else spec.xrange[0]
    fpix_lo = (fpix_ac_left + FPIX_BOUNDARY) / 2.0
    gstr_hi = (GSTRIPS_BOUNDARY + spec.yrange[1]) / 2.0
    fpix_hi = (FPIX_BOUNDARY + spec.xrange[1]) / 2.0

    drawn = []
    # ABCD quadrant labels temporarily disabled (uncomment to restore the A/B/C/D letters + counts):
    # drawn += _abcd_label(fpix_lo, gstr_lo, "#splitline{#scale[1.6]{A}}{%s}" % f"{abcd.A:,}")
    # drawn += _abcd_label(fpix_hi, gstr_lo, "#splitline{#scale[1.6]{B}}{%s}" % f"{abcd.B:,}")
    # drawn += _abcd_label(fpix_lo, gstr_hi, "#splitline{#scale[1.6]{C}}{%s}" % f"{abcd.C:,}")
    # drawn += _abcd_label(fpix_hi, gstr_hi, "#splitline{#scale[1.6]{D}}{blind}")

    pred, unc = abcd.predicted_D()
    if pred is not None:
        header_text = "D_{pred} = B#timesC/A = %.1f #pm %.1f" % (pred, unc)
    else:
        header_text = "D_{pred}: need A, B, C > 0"
    header = ROOT.TLatex(0.14, 0.945, header_text)
    header.SetNDC(True)
    header.SetTextFont(42)
    header.SetTextSize(0.038)
    header.SetTextColor(ROOT.kBlack)
    header.Draw()
    drawn.append(header)
    return drawn


def draw_2d(spec: Plot2DSpec, hist: ROOT.TH2F, abcd: "AbcdCount" = None):
    set_cms_style()
    ROOT.gStyle.SetPalette(ROOT.kViridis)

    canv = ROOT.TCanvas(f"c_{spec.branch_x}_vs_{spec.branch_y}", "", 800, 700)
    canv.SetRightMargin(0.15)
    canv.SetLeftMargin(0.12)
    if spec.log_z:
        canv.SetLogz(True)

    hist.GetXaxis().SetTitle(spec.xlabel)
    hist.GetYaxis().SetTitle(spec.ylabel)
    hist.GetZaxis().SetTitle("Events")
    _style_axes(hist)
    hist.Draw("COLZ")

    overlays = [] 
    if (spec.branch_x, spec.branch_y) == ABCD_AXES:
        overlays.append(_dashed_line(FPIX_BOUNDARY, spec.yrange[0], FPIX_BOUNDARY, spec.yrange[1], ROOT.kRed + 1))
        overlays.append(_dashed_line(spec.xrange[0], GSTRIPS_BOUNDARY, spec.xrange[1], GSTRIPS_BOUNDARY, ROOT.kOrange + 1))
        if FPIX_PRESEL_CUT is not None:   # left edge of A/C (FiPix>0.3 preselection), dotted gray
            overlays.append(_dashed_line(FPIX_PRESEL_CUT, spec.yrange[0], FPIX_PRESEL_CUT, spec.yrange[1], ROOT.kGray + 3, style=ROOT.kDotted))
        if spec.abcd and abcd is not None:
            overlays += _draw_abcd_regions(spec, abcd)

    outpath = os.path.join(OUTPUT_DIR, f"{spec.branch_x}_vs_{spec.branch_y}.png")
    canv.SaveAs(outpath)
    canv.Close()
    print(f"    Saved: {outpath}")

# ABCD background estimate
@dataclass
class AbcdCount:
    """ ABCD resuls accumulated across all chunks loaded into the Fpix vs Gstrip hist.

    A/B/C/D count candidates in the four quadrants. The control-region
    FPix-GStrips correlation is accumulated as running sums (n, Sx, Sy, Sxx, Syy, Sxy) instead of
    keeping every candidate's coordinates, so memory stays flat over arbitrarily large inputs.
    """
    A: int = 0
    B: int = 0
    C: int = 0
    D: int = 0
    # Running sums over control-region (A+B+C) candidates, for the Pearson correlation.
    ctrl_n:   int   = 0
    ctrl_sx:  float = 0.0
    ctrl_sy:  float = 0.0
    ctrl_sxx: float = 0.0
    ctrl_syy: float = 0.0
    ctrl_sxy: float = 0.0

    @property
    def total(self) -> int:
        return self.A + self.B + self.C + self.D

    def predicted_D(self):
        """(prediction, stat. uncertainty) for D = B*C/A, or (None, None) if a region is empty."""
        if self.A > 0 and self.B > 0 and self.C > 0:
            pred = self.B * self.C / self.A
            return pred, pred * np.sqrt(1.0 / self.A + 1.0 / self.B + 1.0 / self.C)
        return None, None

    def control_correlation(self):
        """ Pearson FPix-GStrips correlation over the control region (A+B+C), computed from the
        running sums."""
        n = self.ctrl_n
        if n < 2:
            return None
        cov  = self.ctrl_sxy - self.ctrl_sx * self.ctrl_sy / n
        varx = self.ctrl_sxx - self.ctrl_sx * self.ctrl_sx / n
        vary = self.ctrl_syy - self.ctrl_sy * self.ctrl_sy / n
        if varx <= 0.0 or vary <= 0.0:
            return None
        return cov / np.sqrt(varx * vary)


def accumulate_abcd(data, abcd_mask, abcd: AbcdCount):
    """Add chunk's candidates to the tally."""
    fpix = _flatten(data["DeDx_FiPixelNoL1"][abcd_mask])
    gstr = _flatten(data["DeDx_GiStrip"][abcd_mask])
    if fpix.size == 0:
        return

    high_fpix = fpix > FPIX_BOUNDARY
    high_gstr = gstr > GSTRIPS_BOUNDARY
    abcd.A += int(np.sum(~high_fpix & ~high_gstr))
    abcd.B += int(np.sum( high_fpix & ~high_gstr))
    abcd.C += int(np.sum(~high_fpix &  high_gstr))
    abcd.D += int(np.sum( high_fpix &  high_gstr))

    # Accumulate running sums over the control region (everything except the blinded D quadrant)
    # so the FPix-GStrips correlation needs no per-track buffers that grow with the input size.
    control = ~(high_fpix & high_gstr)
    cf = fpix[control]
    cg = gstr[control]
    abcd.ctrl_n   += int(cf.size)
    abcd.ctrl_sx  += float(cf.sum())
    abcd.ctrl_sy  += float(cg.sum())
    abcd.ctrl_sxx += float(np.dot(cf, cf))
    abcd.ctrl_syy += float(np.dot(cg, cg))
    abcd.ctrl_sxy += float(np.dot(cf, cg))


def report_abcd(abcd: AbcdCount):
    """Print the A/B/C/D yields, the B*C/A prediction (stat. only), and the FPix-GStrips
    correlation in the control region (everything except the blinded D quadrant)."""
    A, B, C, D = abcd.A, abcd.B, abcd.C, abcd.D

    print("\nABCD background estimate (per-track candidates)")
    print(f"  selection : PRESELECTION_CUTS with IsoTrack_pt tightened to > {ABCD_PT_THRESHOLD:g} (FiPix>0.3 kept)")
    print(f"  boundaries: FPix > {FPIX_BOUNDARY}, GStrips > {GSTRIPS_BOUNDARY}")
    print(f"  ABCD-preselected tracks (candidates): {abcd.total:,}")
    print(f"    A (FPix<{FPIX_BOUNDARY}, GStrips<{GSTRIPS_BOUNDARY}) = {A:,}")
    print(f"    B (FPix>{FPIX_BOUNDARY}, GStrips<{GSTRIPS_BOUNDARY}) = {B:,}")
    print(f"    C (FPix<{FPIX_BOUNDARY}, GStrips>{GSTRIPS_BOUNDARY}) = {C:,}")
    print(f"    D (FPix>{FPIX_BOUNDARY}, GStrips>{GSTRIPS_BOUNDARY}) = {D:,}   <- signal region (blinded -> expect 0)")

    pred, unc = abcd.predicted_D()
    if pred is not None:
        print(f"    predicted D = B*C/A = {pred:.2f} +/- {unc:.2f} (stat)")
    else:
        print("    predicted D: undefined (need A, B, C > 0)")

    corr = abcd.control_correlation()
    if corr is not None:
        print(f"    FPix-GStrips Pearson corr in control (A+B+C) = {corr:+.3f}  (ABCD assumes ~0)")


def make_plots():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Inputs ({len(INPUT_FILES)}):")
    for path in INPUT_FILES:
        print(f"  {path}")
    print(f"Output: {OUTPUT_DIR}\n")


    hists_1d        = {spec.branch: make_hist_1d(spec) for spec in PLOT1D_SPECS}
    hists_2d        = {(spec.branch_x, spec.branch_y): make_hist_2d(spec) for spec in PLOT2D_SPECS}
    track_counts_2d = {(spec.branch_x, spec.branch_y): 0 for spec in PLOT2D_SPECS}
    abcd            = AbcdCount()

    print(f"Filling histograms (single pass over {len(INPUT_FILES)} file(s)) ...")
    branches = branches_to_read()
    n_chunks = 0
    for path in INPUT_FILES:
        # For root:// URLs use uproot's  XRootD
        handler = uproot.source.xrootd.XRootDSource if path.startswith("root://") else None
        with uproot.open(path, handler=handler) as f:
            tree = f[TREE_PATH]
            for data in tree.iterate(branches, library="ak", step_size=CHUNK_SIZE):
                n_chunks += 1
                preselection_mask = build_track_mask(data) if PRESELECTION_CUTS else None
                abcd_mask         = build_track_mask(data, ABCD_CUTS)

                for spec in PLOT1D_SPECS:
                    fill_hist_1d(hists_1d[spec.branch], spec, data, preselection_mask)

                for spec in PLOT2D_SPECS:
                    key  = (spec.branch_x, spec.branch_y)
                    mask = abcd_mask if spec.abcd else preselection_mask
                    track_counts_2d[key] += fill_hist_2d(hists_2d[key], spec, data, mask)

                accumulate_abcd(data, abcd_mask, abcd)
    print(f"  processed {n_chunks} chunk(s)\n")
    print(f"2-D histograms ({len(PLOT2D_SPECS)}) ...")


    for spec in PLOT2D_SPECS:
        key = (spec.branch_x, spec.branch_y)
        print(f"  {spec.branch_x} vs {spec.branch_y}: {track_counts_2d[key]:,} tracks after preselection")
        draw_2d(spec, hists_2d[key], abcd)

    print(f"\n1-D histograms ({len(PLOT1D_SPECS)} branches) ...")
    for spec in PLOT1D_SPECS:
        hist = hists_1d[spec.branch]
        print(f"  {spec.branch}: {int(hist.GetEntries()):,} events")
        draw_1d(spec, hist, os.path.join(OUTPUT_DIR, f"{spec.branch}.png"))

    report_abcd(abcd)
    print("\nDone.")


if __name__ == "__main__":
    make_plots()
