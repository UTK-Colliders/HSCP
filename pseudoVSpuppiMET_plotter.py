import ROOT


def _get_events(filepath):
    root_file = ROOT.TFile.Open(filepath)
    if not root_file or root_file.IsZombie():
        raise OSError(f"Unable to open ROOT file: {filepath}")

    analyzer_dir = root_file.Get("HSCPMiniAODAnalyzer")
    if not analyzer_dir:
        root_file.Close()
        raise KeyError("Directory 'HSCPMiniAODAnalyzer' not found")

    events = analyzer_dir.Get("Events")
    if not events:
        root_file.Close()
        raise KeyError("Tree 'Events' not found in 'HSCPMiniAODAnalyzer'")

    return root_file, events


def _configure_style():
    ROOT.gStyle.SetStatX(0.90)
    ROOT.gStyle.SetStatY(0.90)
    ROOT.gStyle.SetStatW(0.20)
    ROOT.gStyle.SetStatH(0.15)

def fullPreselection_2d(filepath, tau, neutralinoMass):
    root_file, events = _get_events(filepath)
    canvas = ROOT.TCanvas("c1", "c1", 800, 600)
    _configure_style()
    events.Draw(
        "PseudoMET_viaCaloJets:RecoPuppiMET >> h2(50,0,2500,50,0,2500)",
        "HLT_FilterOR==1 && IsoTrack_pt>55 && abs(IsoTrack_eta)<2.4 && DeDx_PixelNoL1NOM>=2 && IsoTrack_fractionOfValidHits>0.8 && IsoTrack_numberOfValidHits>=10 && IsoTrack_isHighPurityTrack && IsoTrack_normChi2<5 && abs(IsoTrack_dz)<0.1 && abs(IsoTrack_dxy)<0.02 && IsoTrack_pfMiniRelIsoAll<0.02 && IsoTrack_IsoSumPt_dr03<15 && IsoTrack_pfEnergyOverP<0.3 && IsoTrack_ptErrOverPt2<0.0008 && DeDx_FiPixelNoL1>0.3",
        "COLZ",
    )
    h2 = ROOT.gDirectory.Get("h2")
    h2.GetXaxis().SetTitle("PUppi MET")
    h2.GetYaxis().SetTitle("Pseudo MET")
    ROOT.gPad.Update()
    canvas.SaveAs(f"/eos/user/c/cthompso/analysis/CMSSW_15_0_19/src/plots/pseudoMET_vs_puppiMET/gluino1800_neutralino{neutralinoMass}_{tau}_uds_pseudoVSpuppi.pdf")
    root_file.Close()

def triggerPreselection_2d(filepath, tau, neutralinoMass):
    root_file, events = _get_events(filepath)
    canvas = ROOT.TCanvas("c1", "c1", 800, 600)
    _configure_style()
    events.Draw("PseudoMET_viaCaloJets:RecoPuppiMET >> h2(50,0,2500,50,0,2500)", "HLT_FilterOR==1", "COLZ")
    h2 = ROOT.gDirectory.Get("h2")
    h2.GetXaxis().SetTitle("PUppi MET")
    h2.GetYaxis().SetTitle("Pseudo MET")
    ROOT.gPad.Update()
    canvas.SaveAs(f"/eos/user/c/cthompso/analysis/CMSSW_15_0_19/src/plots/pseudoMET_vs_puppiMET/gluino1800_neutralino{neutralinoMass}_{tau}_uds_pseudoVSpuppi_triggerLevel.pdf")
    root_file.Close()

def fullPreselection_1d_pseudo(filepath, tau, neutralinoMass):
    root_file, events = _get_events(filepath)
    canvas = ROOT.TCanvas("c1", "c1", 800, 600)
    _configure_style()
    events.Draw(
        "PseudoMET_viaCaloJets >> h1(50,0,2500)",
        "HLT_FilterOR==1 && IsoTrack_pt>55 && abs(IsoTrack_eta)<2.4 && DeDx_PixelNoL1NOM>=2 && IsoTrack_fractionOfValidHits>0.8 && IsoTrack_numberOfValidHits>=10 && IsoTrack_isHighPurityTrack && IsoTrack_normChi2<5 && abs(IsoTrack_dz)<0.1 && abs(IsoTrack_dxy)<0.02 && IsoTrack_pfMiniRelIsoAll<0.02 && IsoTrack_IsoSumPt_dr03<15 && IsoTrack_pfEnergyOverP<0.3 && IsoTrack_ptErrOverPt2<0.0008 && DeDx_FiPixelNoL1>0.3",
    )
    h1 = ROOT.gDirectory.Get("h1")
    h1.GetXaxis().SetTitle("Pseudo MET")
    ROOT.gPad.Update()
    canvas.SaveAs(f"/eos/user/c/cthompso/analysis/CMSSW_15_0_19/src/plots/pseudoMET_vs_puppiMET/gluino1800_neutralino{neutralinoMass}_{tau}_uds_pseudoMET.pdf")
    root_file.Close()

def fullPreselection_1d_puppi(filepath, tau, neutralinoMass):
    root_file, events = _get_events(filepath)
    canvas = ROOT.TCanvas("c1", "c1", 800, 600)
    _configure_style()
    events.Draw(
        "RecoPuppiMET >> h1(50,0,2500)",
        "HLT_FilterOR==1 && IsoTrack_pt>55 && abs(IsoTrack_eta)<2.4 && DeDx_PixelNoL1NOM>=2 && IsoTrack_fractionOfValidHits>0.8 && IsoTrack_numberOfValidHits>=10 && IsoTrack_isHighPurityTrack && IsoTrack_normChi2<5 && abs(IsoTrack_dz)<0.1 && abs(IsoTrack_dxy)<0.02 && IsoTrack_pfMiniRelIsoAll<0.02 && IsoTrack_IsoSumPt_dr03<15 && IsoTrack_pfEnergyOverP<0.3 && IsoTrack_ptErrOverPt2<0.0008 && DeDx_FiPixelNoL1>0.3",
    )
    h1 = ROOT.gDirectory.Get("h1")
    h1.GetXaxis().SetTitle("PUppi MET")
    ROOT.gPad.Update()
    canvas.SaveAs(f"/eos/user/c/cthompso/analysis/CMSSW_15_0_19/src/plots/pseudoMET_vs_puppiMET/gluino1800_neutralino{neutralinoMass}_{tau}_uds_puppiMET.pdf")
    root_file.Close()


taus = ['300ps', '1ns', '3ns', '10ns', '30ns']
neutralinoMasses = ['100', '1700']
for tau in taus:
    for neutralinoMass in neutralinoMasses:
        filepath = f"/eos/user/c/cthompso/analysis/CMSSW_15_0_19/src/ntuples/gluino1800_chi10{neutralinoMass}_{tau}_lightDecay_130X_mcRun3_2023_realistic_postBPix_v2_noPUweight.root"
        fullPreselection_2d(filepath, tau, neutralinoMass)
        triggerPreselection_2d(filepath, tau, neutralinoMass)
        fullPreselection_1d_pseudo(filepath, tau, neutralinoMass)
        fullPreselection_1d_puppi(filepath, tau, neutralinoMass)