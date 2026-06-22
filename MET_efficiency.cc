#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
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
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

class MET_efficiency {
public:
    explicit MET_efficiency(std::string inputFile, std::string outputDir);
    ~MET_efficiency() override = default;
private:
    void createHistograms();
    void calculateMETs();
    void fillHistograms();

    bool HLT_PFMET120_PFMHT120_IDTight_;
    bool HLT_PFHT500_PFMET100_PFMHT100_IDTight_;
    bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_;
    bool HLT_MET105_IsoTrk50_;
    double PFHT = 0., PFMHT = 0., PFMHTNoMu = 0.;
    double PFMHT_x = 0., PFMHT_y = 0.;
    double PFMHTNoMu_x = 0., PFMHTNoMu_y = 0.;
    double PUppiMET_NoMu = 0., PUppiMET_x = 0., PUppiMET_y = 0.;
};

void MET_efficiency::MET_efficiency(std::string inputFile, std::string outputDir, bool HLT_PFMET120_PFMHT120_IDTight, bool HLT_PFHT500_PFMET100_PFMHT100_IDTight, bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, bool HLT_MET105_IsoTrk50) {
    TFile inputFile_(inputFile.c_str(), "READ");
    if (inputFile_.IsZombie()) {
        std::cerr << "cannot open input file " << inputFile << std::endl;
        return;
    }

    TTree* tree = nullptr;
    TDirectory* stageDir = inputFile_.GetDirectory("stage");
    if (stageDir) {
        stageDir->GetObject("ttree", tree);
    }
    if (!tree) {
        inputFile_.GetObject("ttree", tree);
    }
    if (!tree) {
        std::cerr << "cannot find ttree in " << inputFile << std::endl;
        return;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<int> runNumber(reader, "runNumber");

}

void MET_efficiency::createHistograms(){
}

void MET_efficiency::calculateMETs(){
    for ( unsigned int j = 0; j < Jet_pt.GetSize(); j++ ) {
        if (Jet_passJetID[j]) {
            PFHT += std::abs(Jet_pt[j]);
            
            PFMHT_x += Jet_px[j];
            PFMHT_y += Jet_py[j];
        }
    }
    PFMHTNoMu_x = PFMHT_x;
    PFMHTNoMu_y = PFMHT_y;

    PUppiMET_x = RecoPUppiMET[0] * cos(RecoPUppiMET_phi[0]);
    PUppiMET_y = RecoPUppiMET[0] * sin(RecoPUppiMET_phi[0]);
    for (unsigned m = 0; m <muon_pt.GetSize(); m++) {
        if (muon_isPFMuon[m]) {
            PFMHTNoMu_x -= muon_pt[m] * cos(muon_phi[m]);
            PFMHTNoMu_y -= muon_pt[m] * sin(muon_phi[m]);

            PUppiMET_x -= muon_pt[m] * cos(muon_phi[m]);
            PUppiMET_y -= muon_pt[m] * sin(muon_phi[m]);
        }
    }
    PFMHT = sqrt(PFMHT_x*PFMHT_x + PFMHT_y*PFMHT_y);
    PFMHTNoMu = sqrt(PFMHTNoMu_x*PFMHTNoMu_x + PFMHTNoMu_y*PFMHTNoMu_y);
    PUppiMET_NoMu = sqrt(PUppiMET_x*PUppiMET_x + PUppiMET_y*PUppiMET_y);
}
    
void MET_efficiency::fillHistograms(){
    for(unsigned int s=0;s<selections_.size();s++) {
        // Trigger Eff vs PUppiMET and vs PseudoMET
        if (selLabels_[s] == "TriggerEffCalib" || selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled") {
            bool BasicSel = false;
            unsigned int hasPassedMuon = 0;
            unsigned int hadPassedHSCP = 0;

            if (*HLT_IsoMu27 && muon_pt.GetSize()==1 && Flag_allMETFilters[0]) BasicSel = true;

            for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) hasPassedMuon++;
            }

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && 
                (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && 
                (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && 
                (Fpix[i_track] < 0.9)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (BasicSel && hasPassedMuon==1 && hadPassedHSCP>0) {

                for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                    if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_pt", muon_pt[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_eta", muon_eta[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_phi", muon_phi[im], *weightPU);
                    }
                }

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = *weightPU;
                if (selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);


                if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

            }
        }

            // Trigger Eff vs PUppiMET and vs PseudoMET for Signal MC
        if (selLabels_[s] == "TriggerEffCalib__Signal" || selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") {
            unsigned int hadPassedHSCP = 0;

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Flag_allMETFilters[0] == true) && (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && 
                (FracOfValidHit[i_track] > 0.8) && (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) &&
                (std::abs(dxy[i_track]) < 0.02) && (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && 
                (Fpix[i_track] > 0.3) && (ptOverptErr[i_track] < 1) && (Ih_Strip[i_track] > 2.9784)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (hadPassedHSCP>0) {

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = *weightPU;
                if (selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);

                if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                        || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
            }  
        }
    }
}