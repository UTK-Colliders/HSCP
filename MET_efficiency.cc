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
    explicit MET_efficiency(std::string inputFile, std::string outputDir,
                            bool HLT_PFMET120_PFMHT120_IDTight = true,
                            bool HLT_PFHT500_PFMET100_PFMHT100_IDTight = true,
                            bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = true,
                            bool HLT_MET105_IsoTrk50 = true);
    ~MET_efficiency() = default;

    struct HistStore {
        std::unordered_map<std::string, TH1F*> h1;
        std::unordered_map<std::string, TH2F*> h2;
    };
    std::vector<HistStore> hists;

private:
    void bindBranches();
    void createHistograms();
    void calculateMETs();
    void fillHistograms();
    void writeHistograms();

    std::string inputFile_;
    std::string outputDir_;
    std::unique_ptr<TFile> inputFileHandle_;
    std::unique_ptr<TFile> outputFileHandle_;
    std::unique_ptr<TTreeReader> reader_;
    TTree* tree_ = nullptr;

    TTreeReaderValue<bool>* HLT_PFMET120_PFMHT120_IDTight_ = nullptr;
    TTreeReaderValue<bool>* HLT_PFHT500_PFMET100_PFMHT100_IDTight_ = nullptr;
    TTreeReaderValue<bool>* HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ = nullptr;
    TTreeReaderValue<bool>* HLT_MET105_IsoTrk50_ = nullptr;
    TTreeReaderValue<bool>* Flag_allMETFilters_ = nullptr;
    TTreeReaderValue<double>* weightPU_ = nullptr;

    TTreeReaderValue<double>* L1MET_ = nullptr;
    TTreeReaderValue<float>* HLTCaloMET_ = nullptr;
    TTreeReaderValue<float>* HLTCaloMHT_ = nullptr;
    TTreeReaderValue<float>* HLTPFMHT_ = nullptr;
    TTreeReaderValue<float>* HLTPFMET_ = nullptr;
    TTreeReaderValue<double>* RecoPFMET_ = nullptr;
    TTreeReaderValue<double>* PseudoMET_viaCaloJets_ = nullptr;
    TTreeReaderValue<double>* RecoPuppiMET_ = nullptr;
    TTreeReaderValue<double>* RecoPuppiMET_phi_ = nullptr;

    TTreeReaderArray<double>* muon_pt_ = nullptr;
    TTreeReaderArray<double>* muon_eta_ = nullptr;
    TTreeReaderArray<double>* muon_phi_ = nullptr;
    TTreeReaderArray<bool>* muon_isTight_ = nullptr;
    TTreeReaderArray<bool>* muon_isPFcand_ = nullptr;
    TTreeReaderArray<float>* muon_pfMiniRelIsoAll_ = nullptr;

    TTreeReaderArray<bool>* HSCP_hasTrack_ = nullptr;
    TTreeReaderArray<double>* Pt_ = nullptr;
    TTreeReaderArray<double>* Pt_pseudo_ = nullptr;
    TTreeReaderArray<double>* Eta_ = nullptr;
    TTreeReaderArray<double>* Px_ = nullptr;
    TTreeReaderArray<double>* Py_ = nullptr;
    TTreeReaderArray<int>* NbPixelHit_noL1_ = nullptr;
    TTreeReaderArray<double>* FracOfValidHit_ = nullptr;
    TTreeReaderArray<int>* NOM_noL1_ = nullptr;
    TTreeReaderArray<bool>* isHighPurityTrack_ = nullptr;
    TTreeReaderArray<double>* normChi2_ = nullptr;
    TTreeReaderArray<double>* dz_ = nullptr;
    TTreeReaderArray<double>* dxy_ = nullptr;
    TTreeReaderArray<float>* miniRelIsoAll_ = nullptr;
    TTreeReaderArray<double>* EoP_ = nullptr;
    TTreeReaderArray<float>* IsoSumPt_dr03_ = nullptr;
    TTreeReaderArray<double>* ptOverptErrptErr_ = nullptr;
    TTreeReaderArray<double>* ptOverptErr_ = nullptr;
    TTreeReaderArray<float>* Ih_Strip_ = nullptr;

    std::vector<double> Fpix_;
    std::vector<std::string> selections_;
    std::vector<std::string> selLabels_;

    double PFHT_ = 0.0;
    double PFMHT_ = 0.0;
    double PFMHTNoMu_ = 0.0;
    double PFMHT_x_ = 0.0;
    double PFMHT_y_ = 0.0;
    double PFMHTNoMu_x_ = 0.0;
    double PFMHTNoMu_y_ = 0.0;
    double PUppiMET_NoMu_ = 0.0;
    double PUppiMET_x_ = 0.0;
    double PUppiMET_y_ = 0.0;
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

template <typename T>
const T& value(TTreeReaderValue<T>* handle) {
    return **handle;
}

template <typename T>
const TTreeReaderArray<T>& Array(const TTreeReaderArray<T>* handle) {
    return *handle;
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

MET_efficiency::MET_efficiency(std::string inputFile, std::string outputDir,
                               bool HLT_PFMET120_PFMHT120_IDTight,
                               bool HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                               bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                               bool HLT_MET105_IsoTrk50)
    : inputFile_(std::move(inputFile)), outputDir_(std::move(outputDir)) {
    (void)HLT_PFMET120_PFMHT120_IDTight;
    (void)HLT_PFHT500_PFMET100_PFMHT100_IDTight;
    (void)HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
    (void)HLT_MET105_IsoTrk50;

    inputFileHandle_.reset(TFile::Open(inputFile_.c_str(), "READ"));
    if (!inputFileHandle_ || inputFileHandle_->IsZombie()) {
        std::cerr << "cannot open input file " << inputFile_ << std::endl;
        return;
    }

    tree_ = dynamic_cast<TTree*>(inputFileHandle_->Get("HSCPMiniAODAnalyzer/Events"));
    if (!tree_) {
        inputFileHandle_->GetObject("Events", tree_);
    }
    if (!tree_) {
        std::cerr << "cannot find Events tree in " << inputFile_ << std::endl;
        return;
    }

    reader_ = std::make_unique<TTreeReader>(tree_);
    bindBranches();
    createHistograms();

    if (!outputDir_.empty()) {
        gSystem->mkdir(outputDir_.c_str(), true);
        outputFileHandle_.reset(TFile::Open((outputDir_ + "/MET_efficiency.root").c_str(), "RECREATE"));
    }

    while (reader_->Next()) {
        calculateMETs();
        fillHistograms();
    }

    writeHistograms();
}

void MET_efficiency::bindBranches() {
    bindValue(*reader_, HLT_PFMET120_PFMHT120_IDTight_, "HLT_PFMET120_PFMHT120_IDTight");
    bindValue(*reader_, HLT_PFHT500_PFMET100_PFMHT100_IDTight_, "HLT_PFHT500_PFMET100_PFMHT100_IDTight");
    bindValue(*reader_, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60");
    bindValue(*reader_, HLT_MET105_IsoTrk50_, "HLT_MET105_IsoTrk50");
    bindValue(*reader_, Flag_allMETFilters_, "Flag_allMETFilters");
    bindValue(*reader_, weightPU_, "weightPU");

    bindValue(*reader_, L1MET_, "L1MET");
    bindValue(*reader_, HLTCaloMET_, "HLTCaloMET");
    bindValue(*reader_, HLTCaloMHT_, "HLTCaloMHT");
    bindValue(*reader_, HLTPFMHT_, "HLTPFMHT");
    bindValue(*reader_, HLTPFMET_, "HLTPFMET");
    bindValue(*reader_, RecoPFMET_, "RecoPFMET");
    bindValue(*reader_, PseudoMET_viaCaloJets_, "PseudoMET_viaCaloJets");
    bindValue(*reader_, RecoPuppiMET_, "RecoPuppiMET");
    bindValue(*reader_, RecoPuppiMET_phi_, "RecoPuppiMET_phi");

    bindArray(*reader_, muon_pt_, "muon_pt");
    bindArray(*reader_, muon_eta_, "muon_eta");
    bindArray(*reader_, muon_phi_, "muon_phi");
    bindArray(*reader_, muon_isTight_, "muon_isTight");
    bindArray(*reader_, muon_isPFcand_, "Muon_isPFcand");
    bindArray(*reader_, muon_pfMiniRelIsoAll_, "muon_pfMiniRelIsoAll");

    bindArray(*reader_, HSCP_hasTrack_, "HSCP_hasTrack");
    bindArray(*reader_, Pt_, "IsoTrack_pt");
    bindArray(*reader_, Pt_pseudo_, "IsoTrack_PseudoTrack_pt");
    bindArray(*reader_, Eta_, "IsoTrack_eta");
    bindArray(*reader_, Px_, "IsoTrack_px");
    bindArray(*reader_, Py_, "IsoTrack_py");
    bindArray(*reader_, NbPixelHit_noL1_, "IsoTrack_numberOfValidPixelHits");
    bindArray(*reader_, FracOfValidHit_, "IsoTrack_fractionOfValidHits");
    bindArray(*reader_, NOM_noL1_, "IsoTrack_numberOfTrackerLayers");
    bindArray(*reader_, isHighPurityTrack_, "IsoTrack_isHighPurityTrack");
    bindArray(*reader_, normChi2_, "IsoTrack_normChi2");
    bindArray(*reader_, dz_, "IsoTrack_dz");
    bindArray(*reader_, dxy_, "IsoTrack_dxy");
    bindArray(*reader_, miniRelIsoAll_, "IsoTrack_pfMiniRelIsoAll");
    bindArray(*reader_, EoP_, "IsoTrack_pfEnergyOverP");
    bindArray(*reader_, IsoSumPt_dr03_, "IsoTrack_IsoSumPt_dr03");
    bindArray(*reader_, ptOverptErrptErr_, "IsoTrack_ptErrOverPt2");
    bindArray(*reader_, ptOverptErr_, "IsoTrack_ptErrOverPt");
    bindArray(*reader_, Ih_Strip_, "DeDx_IhStrip");
}

void MET_efficiency::createHistograms() {
    selections_ = {"TriggerEffCalib", "TriggerEffCalib_PseudoMETrescaled", "TriggerEffCalib__Signal", "TriggerEffCalib_PseudoMETrescaled__Signal"};
    selLabels_ = selections_;
    hists.resize(selLabels_.size());
}

void MET_efficiency::calculateMETs() {
    PFHT_ = 0.0;
    PFMHT_ = 0.0;
    PFMHTNoMu_ = 0.0;
    PFMHT_x_ = 0.0;
    PFMHT_y_ = 0.0;
    PFMHTNoMu_x_ = 0.0;
    PFMHTNoMu_y_ = 0.0;
    PUppiMET_x_ = 0.0;
    PUppiMET_y_ = 0.0;
    PUppiMET_NoMu_ = 0.0;
    Fpix_.clear();

    const auto& pt = Array(Pt_);
    const auto& px = Array(Px_);
    const auto& py = Array(Py_);
    const auto& muonPt = Array(muon_pt_);
    const auto& muonPhi = Array(muon_phi_);
    const auto& muonIsPF = Array(muon_isPFcand_);

    for (unsigned int i = 0; i < pt.GetSize(); ++i) {
        PFHT_ += std::abs(pt[i]);
        PFMHT_x_ -= px[i];
        PFMHT_y_ -= py[i];
    }

    PFMHTNoMu_x_ = PFMHT_x_;
    PFMHTNoMu_y_ = PFMHT_y_;

    PUppiMET_x_ = value(RecoPuppiMET_) * std::cos(value(RecoPuppiMET_phi_));
    PUppiMET_y_ = value(RecoPuppiMET_) * std::sin(value(RecoPuppiMET_phi_));

    for (unsigned int m = 0; m < muonPt.GetSize(); ++m) {
        if (m < muonIsPF.GetSize() && muonIsPF[m]) {
            const double muPt = muonPt[m];
            const double muPhi = muonPhi[m];
            PFMHTNoMu_x_ -= muPt * std::cos(muPhi);
            PFMHTNoMu_y_ -= muPt * std::sin(muPhi);
            PUppiMET_x_ -= muPt * std::cos(muPhi);
            PUppiMET_y_ -= muPt * std::sin(muPhi);
        }
    }

    PFMHT_ = std::sqrt(PFMHT_x_ * PFMHT_x_ + PFMHT_y_ * PFMHT_y_);
    PFMHTNoMu_ = std::sqrt(PFMHTNoMu_x_ * PFMHTNoMu_x_ + PFMHTNoMu_y_ * PFMHTNoMu_y_);
    PUppiMET_NoMu_ = std::sqrt(PUppiMET_x_ * PUppiMET_x_ + PUppiMET_y_ * PUppiMET_y_);
}

void MET_efficiency::fillHistograms() {
    for (std::size_t s = 0; s < selLabels_.size(); ++s) {
        const std::string& label = selLabels_[s];
        const double weight = value(weightPU_);
        const double scale = (label.find("PseudoMETrescaled") != std::string::npos) ? 163.451 / 146.864 : 1.0;
        const double pseudoCaloMET = value(PseudoMET_viaCaloJets_) * scale;
        const double recoPFMET = value(RecoPFMET_);
        const double l1MET = value(L1MET_);
        const double hltCaloMET = value(HLTCaloMET_);
        const double hltCaloMHT = value(HLTCaloMHT_);
        const double hltPFMHT = value(HLTPFMHT_);
        const double hltPFMET = value(HLTPFMET_);
        const double recoPuppiMET = value(RecoPuppiMET_);
        const bool hltPFMET120 = value(HLT_PFMET120_PFMHT120_IDTight_);
        const bool hltPFHT500 = value(HLT_PFHT500_PFMET100_PFMHT100_IDTight_);
        const bool hltNoMu120 = value(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_);
        const bool hltIsoTrk50 = value(HLT_MET105_IsoTrk50_);

        auto& store = hists[s];
        const std::vector<std::pair<std::string, double>> hist1 = {
            {label + "_PUppiMET", recoPuppiMET},
            {label + "_PUppiMETNoMu", PUppiMET_NoMu_},
            {label + "_PseudoCaloMET", pseudoCaloMET},
            {label + "_RecoPFMET", recoPFMET},
            {label + "_L1MET", l1MET},
            {label + "_HLTCaloMET", hltCaloMET},
            {label + "_HLTCaloMHT", hltCaloMHT},
            {label + "_HLTPFMHT", hltPFMHT},
            {label + "_HLTPFMET", hltPFMET}
        };
        for (const auto& [name, value] : hist1) {
            if (!store.h1[name]) {
                store.h1[name] = new TH1F(name.c_str(), name.c_str(), 100, 0.0, axisMax1D(name));
                store.h1[name]->Sumw2();
            }
            store.h1[name]->Fill(value, weight);
        }

        if (hltPFMET120 && PFMHT_ > 120.0) {
            const std::string name = label + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET";
            if (!store.h1[name]) {
                store.h1[name] = new TH1F(name.c_str(), name.c_str(), 100, 0.0, axisMax1D(name));
                store.h1[name]->Sumw2();
            }
            store.h1[name]->Fill(pseudoCaloMET, weight);
        }

        if (hltPFHT500 && PFHT_ > 500.0 && PFMHT_ > 100.0) {
            const std::string name = label + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET";
            if (!store.h1[name]) {
                store.h1[name] = new TH1F(name.c_str(), name.c_str(), 100, 0.0, axisMax1D(name));
                store.h1[name]->Sumw2();
            }
            store.h1[name]->Fill(pseudoCaloMET, weight);
        }

        if (hltNoMu120 && PFMHTNoMu_ > 120.0 && PFHT_ > 60.0) {
            const std::string name = label + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET";
            if (!store.h1[name]) {
                store.h1[name] = new TH1F(name.c_str(), name.c_str(), 100, 0.0, axisMax1D(name));
                store.h1[name]->Sumw2();
            }
            store.h1[name]->Fill(pseudoCaloMET, weight);
        }

        if (hltIsoTrk50) {
            const std::string name = label + "_if___HLT_MET105_IsoTrk50___PFtrackPT";
            if (!store.h1[name]) {
                store.h1[name] = new TH1F(name.c_str(), name.c_str(), 100, 0.0, axisMax1D(name));
                store.h1[name]->Sumw2();
            }
            const auto& pt = Array(Pt_);
            const auto& hasTrack = Array(HSCP_hasTrack_);
            for (unsigned int i = 0; i < pt.GetSize(); ++i) {
                if (i < hasTrack.GetSize() && hasTrack[i]) {
                    store.h1[name]->Fill(pt[i], weight);
                }
            }
        }
    }
}

void MET_efficiency::writeHistograms() {
    if (!outputFileHandle_ || outputFileHandle_->IsZombie()) {
        return;
    }
    outputFileHandle_->cd();
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
    outputFileHandle_->Write();
    outputFileHandle_->Close();
}