#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

class MET_efficiency {
public:
    explicit MET_efficiency(std::string inputFile, std::string outputDir);
    ~MET_efficiency() = default;

    struct HistStore {
        std::unordered_map<std::string, TH1F*> h1;
        std::unordered_map<std::string, TH2F*> h2;
        void FillHisto1F(const std::string& name, double value, double weight) {
            if (h1.find(name) != h1.end()) {
                h1[name]->Fill(value, weight);
            } else {
                std::cerr << "Histogram " << name << " not found!" << std::endl;
            }
        }
        void FillHisto2F(const std::string& name, double xvalue, double yvalue, double weight) {
            if (h2.find(name) != h2.end()) {
                h2[name]->Fill(xvalue, yvalue, weight);
            } else {
                std::cerr << "Histogram " << name << " not found!" << std::endl;
            }
        }
    };
    std::vector<HistStore> hists;

private:
    void bindBranches();
    void createHistograms();
    void calculateMETs();
    void fillHistograms();
    void writeHistograms();

    std::string inputFile;
    std::string outputDir;
    std::unique_ptr<TFile> inputFileHandle;
    std::unique_ptr<TFile> outputFileHandle;
    std::unique_ptr<TTreeReader> reader;
    TTree* tree = nullptr;

    TTreeReaderValue<bool>* HLT_PFMET120_PFMHT120_IDTight = nullptr;
    TTreeReaderValue<bool>* HLT_PFHT500_PFMET100_PFMHT100_IDTight = nullptr;
    TTreeReaderValue<bool>* HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = nullptr;
    TTreeReaderValue<bool>* HLT_MET105_IsoTrk50 = nullptr;
    TTreeReaderValue<double>* weightPU = nullptr;

    TTreeReaderArray<double>* L1MET = nullptr;
    TTreeReaderArray<float>* HLTCaloMET = nullptr;
    TTreeReaderArray<float>* HLTCaloMHT = nullptr;
    TTreeReaderArray<float>* HLTPFMHT = nullptr;
    TTreeReaderArray<float>* HLTPFMET = nullptr;
    TTreeReaderArray<double>* RecoPFMET = nullptr;
    TTreeReaderArray<double>* PseudoMET_viaCaloJets = nullptr;
    TTreeReaderArray<bool>* Flag_allMETFilters = nullptr;
    TTreeReaderArray<double>* RecoPUppiMET = nullptr;
    TTreeReaderArray<double>* RecoPUppiMET_phi = nullptr;
    TTreeReaderArray<bool>* HSCP_hasTrack = nullptr;
    TTreeReaderArray<double>* Pt = nullptr;
    TTreeReaderArray<double>* Jet_px = nullptr;
    TTreeReaderArray<double>* Jet_py = nullptr;
    TTreeReaderArray<double>* Jet_pt = nullptr;
    TTreeReaderArray<bool>* Jet_passJetID = nullptr;
    TTreeReaderArray<double>* muon_pt = nullptr;
    TTreeReaderArray<double>* muon_phi = nullptr;
    TTreeReaderArray<double>* muon_eta = nullptr;
    TTreeReaderArray<bool>* muon_isPFMuon = nullptr;
    TTreeReaderArray<bool>* muon_isTight = nullptr;
    TTreeReaderArray<float>* muon_pfMiniRelIsoAll = nullptr;

    TTreeReaderArray<double>* Pt_pseudo = nullptr;
    TTreeReaderArray<double>* Eta = nullptr;
    TTreeReaderArray<double>* Px = nullptr;
    TTreeReaderArray<double>* Py = nullptr;
    TTreeReaderArray<int>* NbPixelHit_noL1 = nullptr;
    TTreeReaderArray<double>* FracOfValidHit = nullptr;
    TTreeReaderArray<int>* NOM_noL1 = nullptr;
    TTreeReaderArray<bool>* isHighPurityTrack = nullptr;
    TTreeReaderArray<double>* normChi2 = nullptr;
    TTreeReaderArray<double>* dz = nullptr;
    TTreeReaderArray<double>* dxy = nullptr;
    TTreeReaderArray<float>* miniRelIsoAll = nullptr;
    TTreeReaderArray<double>* EoP = nullptr;
    TTreeReaderArray<float>* IsoSumPt_dr03 = nullptr;
    TTreeReaderArray<double>* ptOverptErrptErr = nullptr;
    TTreeReaderArray<double>* ptOverptErr = nullptr;
    TTreeReaderArray<float>* Ih_Strip = nullptr;
    TTreeReaderArray<float>* Fpix = nullptr;

    std::vector<std::string> selLabels = {"TriggerEffCalib", "TriggerEffCalib_PseudoMETrescaled", "TriggerEffCalib__Signal", "TriggerEffCalib_PseudoMETrescaled__Signal"};

    double PFHT = 0., PFMHT = 0., PFMHTNoMu = 0.;
    double PFMHTx = 0., PFMHTy = 0.;
    double PFMHTNoMux = 0., PFMHTNoMuy = 0.;
    double PUppiMET_NoMu = 0., PUppiMET_x = 0., PUppiMET_y = 0.;
};

namespace {

template <typename T>
void bindValue(TTreeReader& reader, TTreeReaderValue<T>*& handle, const char* branchName) {
    handle = new TTreeReaderValue<T>(reader, branchName);
}

template <typename T>
void bindArray(TTreeReader& reader, TTreeReaderArray<T>*& handle, const char* branchName) {
    handle = new TTreeReaderArray<T>(reader, branchName);
}

double axisMax1D(const std::string& name) {
    if (name.find("eta") != std::string::npos) return 3.0;
    if (name.find("phi") != std::string::npos) return 3.2;
    if (name.find("dxy") != std::string::npos) return 0.1;
    if (name.find("dz") != std::string::npos) return 0.5;
    if (name.find("IsoSumPt") != std::string::npos) return 50.0;
    if (name.find("PT") != std::string::npos || name.find("pt") != std::string::npos) return 500.0;
    return 1000.0;
}

} // namespace

MET_efficiency::MET_efficiency(std::string inputFile, std::string outputDir)
    : inputFile(std::move(inputFile)), outputDir(std::move(outputDir)) {

    inputFileHandle.reset(TFile::Open(this->inputFile.c_str(), "READ"));
    if (!inputFileHandle || inputFileHandle->IsZombie()) {
        std::cerr << "cannot open input file " << this->inputFile << std::endl;
        return;
    }

    tree = dynamic_cast<TTree*>(inputFileHandle->Get("HSCPMiniAODAnalyzer/Events"));
    if (!tree) {
        inputFileHandle->GetObject("Events", tree);
    }
    if (!tree) {
        std::cerr << "cannot find Events tree in " << this->inputFile << std::endl;
        return;
    }

    reader = std::make_unique<TTreeReader>(tree);
    bindBranches();
    createHistograms();

    if (!this->outputDir.empty()) {
        gSystem->mkdir(this->outputDir.c_str(), true);
        outputFileHandle.reset(TFile::Open((this->outputDir + "/MET_efficiency.root").c_str(), "RECREATE"));
    }

    while (reader->Next()) {
        calculateMETs();
        fillHistograms();
    }

    writeHistograms();
}

void MET_efficiency::bindBranches() {
    bindValue(*reader, HLT_PFMET120_PFMHT120_IDTight, "HLT_PFMET120_PFMHT120_IDTight");
    bindValue(*reader, HLT_PFHT500_PFMET100_PFMHT100_IDTight, "HLT_PFHT500_PFMET100_PFMHT100_IDTight");
    bindValue(*reader, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60");
    bindValue(*reader, HLT_MET105_IsoTrk50, "HLT_MET105_IsoTrk50");
    bindValue(*reader, weightPU, "weightPU");

    bindArray(*reader, L1MET, "L1MET");
    bindArray(*reader, HLTCaloMET, "HLTCaloMET");
    bindArray(*reader, HLTCaloMHT, "HLTCaloMHT");
    bindArray(*reader, HLTPFMHT, "HLTPFMHT");
    bindArray(*reader, HLTPFMET, "HLTPFMET");
    bindArray(*reader, RecoPFMET, "RecoPFMET");
    bindArray(*reader, PseudoMET_viaCaloJets, "PseudoMET_viaCaloJets");
    bindArray(*reader, RecoPUppiMET, "RecoPuppiMET");
    bindArray(*reader, RecoPUppiMET_phi, "RecoPuppiMET_phi");
    bindArray(*reader, Jet_px, "puppijet_px");
    bindArray(*reader, Jet_py, "puppijet_py");
    bindArray(*reader, Jet_pt, "puppijet_pt");
    bindArray(*reader, Jet_passJetID, "puppijet_passJetID");
    bindArray(*reader, muon_pt, "muon_pt");
    bindArray(*reader, muon_eta, "muon_eta");
    bindArray(*reader, muon_phi, "muon_phi");
    bindArray(*reader, muon_isPFMuon, "muon_isPFMuon");
    bindArray(*reader, muon_isTight, "muon_isTight");
    bindArray(*reader, muon_pfMiniRelIsoAll, "muon_pfMiniRelIsoAll");
    bindArray(*reader, HSCP_hasTrack, "HSCP_hasTrack");
    bindArray(*reader, Pt, "IsoTrack_pt");
    bindArray(*reader, Pt_pseudo, "IsoTrack_PseudoTrack_pt");
    bindArray(*reader, Eta, "IsoTrack_eta");
    bindArray(*reader, Px, "IsoTrack_px");
    bindArray(*reader, Py, "IsoTrack_py");
    bindArray(*reader, NbPixelHit_noL1, "IsoTrack_numberOfValidPixelHits");
    bindArray(*reader, FracOfValidHit, "IsoTrack_fractionOfValidHits");
    bindArray(*reader, NOM_noL1, "IsoTrack_numberOfTrackerLayers");
    bindArray(*reader, isHighPurityTrack, "IsoTrack_isHighPurityTrack");
    bindArray(*reader, normChi2, "IsoTrack_normChi2");
    bindArray(*reader, dz, "IsoTrack_dz");
    bindArray(*reader, dxy, "IsoTrack_dxy");
    bindArray(*reader, miniRelIsoAll, "IsoTrack_pfMiniRelIsoAll");
    bindArray(*reader, EoP, "IsoTrack_pfEnergyOverP");
    bindArray(*reader, IsoSumPt_dr03, "IsoTrack_IsoSumPt_dr03");
    bindArray(*reader, ptOverptErrptErr, "IsoTrack_ptErrOverPt2");
    bindArray(*reader, ptOverptErr, "IsoTrack_ptErrOverPt");
    bindArray(*reader, Ih_Strip, "DeDx_IhStrip");
    bindArray(*reader, Fpix, "DeDx_FiPixelNoL1");
}

void MET_efficiency::createHistograms() {
    hists.resize(selLabels.size());
}

void MET_efficiency::calculateMETs() {
    const auto& Jet_pt_ = *Jet_pt;
    const auto& Jet_px_ = *Jet_px;
    const auto& Jet_py_ = *Jet_py;
    const auto& Jet_passJetID_ = *Jet_passJetID;
    const auto& RecoPUppiMET_ = *RecoPUppiMET;
    const auto& RecoPUppiMET_phi_ = *RecoPUppiMET_phi;
    const auto& muon_pt_ = *muon_pt;
    const auto& muon_phi_ = *muon_phi;
    const auto& muon_isPFMuon_ = *muon_isPFMuon;

    PFHT = 0., PFMHT = 0., PFMHTNoMu = 0.;
    PFMHTx = 0., PFMHTy = 0.;
    PFMHTNoMux = 0., PFMHTNoMuy = 0.;
    PUppiMET_NoMu = 0., PUppiMET_x = 0., PUppiMET_y = 0.;
    for ( unsigned int j = 0; j < Jet_pt_.GetSize(); j++ ) {
        if (Jet_passJetID_[j]) {
            PFHT += std::abs(Jet_pt_[j]);
            PFMHTx += Jet_px_[j];
            PFMHTy += Jet_py_[j];
        }
    }
    PFMHTNoMux = PFMHTx;
    PFMHTNoMuy = PFMHTy;

    PUppiMET_x = RecoPUppiMET_[0] * std::cos(RecoPUppiMET_phi_[0]);
    PUppiMET_y = RecoPUppiMET_[0] * std::sin(RecoPUppiMET_phi_[0]);
    for (unsigned m = 0; m < muon_pt_.GetSize(); m++) {
        if (muon_isPFMuon_[m]) {
            PFMHTNoMux -= muon_pt_[m] * std::cos(muon_phi_[m]);
            PFMHTNoMuy -= muon_pt_[m] * std::sin(muon_phi_[m]);

            PUppiMET_x -= muon_pt_[m] * std::cos(muon_phi_[m]);
            PUppiMET_y -= muon_pt_[m] * std::sin(muon_phi_[m]);
        }
    }
    PFMHT = std::sqrt(PFMHTx*PFMHTx + PFMHTy*PFMHTy);
    PFMHTNoMu = std::sqrt(PFMHTNoMux*PFMHTNoMux + PFMHTNoMuy*PFMHTNoMuy);
    PUppiMET_NoMu = std::sqrt(PUppiMET_x*PUppiMET_x + PUppiMET_y*PUppiMET_y);
}

void MET_efficiency::fillHistograms() {
    const bool HLT_PFMET120_PFMHT120_IDTight_ = **HLT_PFMET120_PFMHT120_IDTight;
    const bool HLT_PFHT500_PFMET100_PFMHT100_IDTight_ = **HLT_PFHT500_PFMET100_PFMHT100_IDTight;
    const bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ = **HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
    const bool HLT_MET105_IsoTrk50_ = **HLT_MET105_IsoTrk50;
    const double weightPU_ = **weightPU;

    const auto& RecoPUppiMET_ = *RecoPUppiMET;
    const auto& RecoPUppiMET_phi_ = *RecoPUppiMET_phi;
    const auto& PseudoMET_viaCaloJets_ = *PseudoMET_viaCaloJets;
    const auto& RecoPFMET_ = *RecoPFMET;
    const auto& HSCP_hasTrack_ = *HSCP_hasTrack;
    const auto& Pt_ = *Pt;
    const auto& Pt_pseudo_ = *Pt_pseudo;
    const auto& Eta_ = *Eta;
    const auto& NbPixelHit_noL1_ = *NbPixelHit_noL1;
    const auto& FracOfValidHit_ = *FracOfValidHit;
    const auto& NOM_noL1_ = *NOM_noL1;
    const auto& isHighPurityTrack_ = *isHighPurityTrack;
    const auto& normChi2_ = *normChi2;
    const auto& dz_ = *dz;
    const auto& dxy_ = *dxy;
    const auto& miniRelIsoAll_ = *miniRelIsoAll;
    const auto& EoP_ = *EoP;
    const auto& IsoSumPt_dr03_ = *IsoSumPt_dr03;
    const auto& ptOverptErrptErr_ = *ptOverptErrptErr;
    const auto& ptOverptErr_ = *ptOverptErr;
    const auto& Ih_Strip_ = *Ih_Strip;
    const auto& Fpix_ = *Fpix;
    const auto& muon_pt_ = *muon_pt;
    const auto& muon_eta_ = *muon_eta;
    const auto& muon_phi_ = *muon_phi;
    const auto& muon_isTight_ = *muon_isTight;
    const auto& muon_pfMiniRelIsoAll_ = *muon_pfMiniRelIsoAll;

    const auto& Flag_allMETFilters_ = *Flag_allMETFilters;
    const auto& L1MET_ = *L1MET;
    const auto& HLTCaloMET_ = *HLTCaloMET;
    const auto& HLTCaloMHT_ = *HLTCaloMHT;
    const auto& HLTPFMHT_ = *HLTPFMHT;
    const auto& HLTPFMET_ = *HLTPFMET;

    for(unsigned int s=0;s<selLabels.size();s++) {
            // Trigger Eff vs PUppiMET and vs PseudoMET
        if (selLabels[s] == "TriggerEffCalib" || selLabels[s] == "TriggerEffCalib_PseudoMETrescaled") {
            bool BasicSel = false;
            unsigned int hasPassedMuon = 0;
            unsigned int hadPassedHSCP = 0;

            if (muon_pt_.GetSize()==1 && Flag_allMETFilters_[0]) BasicSel = true;

            for (unsigned int im = 0; im < muon_pt_.GetSize(); im++) {
                if (muon_isTight_[im] && muon_pt_[im]>30 && muon_pfMiniRelIsoAll_[im]<0.15) hasPassedMuon++;
            }

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack_.GetSize(); j++){    

                if (!HSCP_hasTrack_[j]) continue;

                if ( (Pt_[i_track] > 50.0) && (Pt_pseudo_[i_track] > 50.0) && (std::abs(Eta_[i_track]) < 2.4) && (NbPixelHit_noL1_[i_track] >= 2) && (FracOfValidHit_[i_track] > 0.8) && 
                (NOM_noL1_[i_track] >= 10) && (isHighPurityTrack_[i_track] == true) && (normChi2_[i_track] < 5.0) && (std::abs(dz_[i_track]) < 0.1) && (std::abs(dxy_[i_track]) < 0.02) && 
                (miniRelIsoAll_[i_track] < 0.02) && (EoP_[i_track] < 0.3) && (IsoSumPt_dr03_[i_track] < 15) && (ptOverptErrptErr_[i_track] < 0.0008) && 
                (Fpix_[i_track] < 0.9)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (BasicSel && hasPassedMuon==1 && hadPassedHSCP>0) {

                for (unsigned int im = 0; im < muon_pt_.GetSize(); im++) {
                    if (muon_isTight_[im] && muon_pt_[im]>30 && muon_pfMiniRelIsoAll_[im]<0.15) {
                        hists[s].FillHisto1F(selLabels[s] + "_muon_pt", muon_pt_[im], weightPU_);
                        hists[s].FillHisto1F(selLabels[s] + "_muon_eta", muon_eta_[im], weightPU_);
                        hists[s].FillHisto1F(selLabels[s] + "_muon_phi", muon_phi_[im], weightPU_);
                    }
                }

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = weightPU_;
                if (selLabels[s] == "TriggerEffCalib_PseudoMETrescaled") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) hists[s].FillHisto1F(selLabels[s] + "_PFtrackPT", Pt_[i_candPartialHSCP[iHSCP]], weightOnCalib);
                hists[s].FillHisto2F(selLabels[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto2F(selLabels[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_RecoPFMET", RecoPFMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_L1MET", L1MET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTCaloMET", HLTCaloMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTCaloMHT", HLTCaloMHT_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTPFMHT", HLTPFMHT_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTPFMET", HLTPFMET_[0], weightOnCalib);


                if (HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_MET105_IsoTrk50_) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt_[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                    || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60) || HLT_MET105_IsoTrk50_) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                    || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                     || HLT_MET105_IsoTrk50_) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)
                     || HLT_MET105_IsoTrk50_) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100) || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)
                     || HLT_MET105_IsoTrk50_) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }

            }
        }

            // Trigger Eff vs PUppiMET and vs PseudoMET for Signal MC
        if (selLabels[s] == "TriggerEffCalib__Signal" || selLabels[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") {
            unsigned int hadPassedHSCP = 0;

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack_.GetSize(); j++){    

                if (!HSCP_hasTrack_[j]) continue;

                if ( (Flag_allMETFilters_[0] == true) && (Pt_[i_track] > 50.0) && (Pt_pseudo_[i_track] > 50.0) && (std::abs(Eta_[i_track]) < 2.4) && (NbPixelHit_noL1_[i_track] >= 2) && 
                (FracOfValidHit_[i_track] > 0.8) && (NOM_noL1_[i_track] >= 10) && (isHighPurityTrack_[i_track] == true) && (normChi2_[i_track] < 5.0) && (std::abs(dz_[i_track]) < 0.1) &&
                (std::abs(dxy_[i_track]) < 0.02) && (miniRelIsoAll_[i_track] < 0.02) && (EoP_[i_track] < 0.3) && (IsoSumPt_dr03_[i_track] < 15) && (ptOverptErrptErr_[i_track] < 0.0008) && 
                (Fpix_[i_track] > 0.3) && (ptOverptErr_[i_track] < 1) && (Ih_Strip_[i_track] > 2.9784)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (hadPassedHSCP>0) {

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = weightPU_;
                if (selLabels[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) hists[s].FillHisto1F(selLabels[s] + "_PFtrackPT", Pt_[i_candPartialHSCP[iHSCP]], weightOnCalib);
                hists[s].FillHisto2F(selLabels[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto2F(selLabels[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_RecoPFMET", RecoPFMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_L1MET", L1MET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTCaloMET", HLTCaloMET_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTCaloMHT", HLTCaloMHT_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTPFMHT", HLTPFMHT_[0], weightOnCalib);
                hists[s].FillHisto1F(selLabels[s] + "_HLTPFMET", HLTPFMET_[0], weightOnCalib);

                if (HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if (HLT_MET105_IsoTrk50_) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt_[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                    || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60) || (HLT_MET105_IsoTrk50_)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                    || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100)
                     || (HLT_MET105_IsoTrk50_)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFMET120_PFMHT120_IDTight_ && PFMHT>120) || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)
                     || (HLT_MET105_IsoTrk50_)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
                if ((HLT_PFHT500_PFMET100_PFMHT100_IDTight_ && PFHT>500 && PFMHT>100) || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ && PFMHTNoMu>120 && PFHT>60)
                     || (HLT_MET105_IsoTrk50_)) {
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET_[0], weightOnCalib);
                    hists[s].FillHisto1F(selLabels[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET_[0], PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                    hists[s].FillHisto2F(selLabels[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoMET_viaCaloJets_[0]*isRescaled, weightOnCalib);
                }
            }
        }
    }
}

void MET_efficiency::writeHistograms() {
    if (!outputFileHandle || outputFileHandle->IsZombie()) {
        return;
    }
    outputFileHandle->cd();
    for (auto& store : hists) {
        for (auto& item : store.h1) {
            if (item.second) {
                item.second->Write();
            }
        }
        for (auto& item : store.h2) {
            if (item.second) {
                item.second->Write();
            }
        }
    }
    outputFileHandle->Write();
    outputFileHandle->Close();
}