# Create MET efficiency plots from ntuples

import argparse
import os
import sys
import ROOT


def load_and_run(inputFile, outputDir, 
                 HLT_PFMET120_PFMHT120_IDTight=True, 
                 HLT_PFHT500_PFMET100_PFMHT100_IDTight=True, 
                 HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=True, 
                 HLT_MET105_IsoTrk50=True):
    cmssw_base = os.environ['CMSSW_BASE']
    source_path = os.path.join(cmssw_base, 'src', 'HSCP', 'MET_efficiency.cc')
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load('libTree')
    ROOT.gROOT.ProcessLine(f'.L {source_path}+')
    ROOT.MET_efficiency(inputFile, outputDir, 
                        HLT_PFMET120_PFMHT120_IDTight, 
                        HLT_PFHT500_PFMET100_PFMHT100_IDTight, 
                        HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, 
                        HLT_MET105_IsoTrk50)


def main():
    parser = argparse.ArgumentParser(description='Create MET efficiency plots from ntuples')
    parser.add_argument('--inputFile', required=True, help='ntuple file')
    parser.add_argument('--outputDir', required=True, help='Output directory for plots')
    parser.add_argument('--HLT_PFMET120_PFMHT120_IDTight', default=True, help='Trigger option')
    parser.add_argument('--HLT_PFHT500_PFMET100_PFMHT100_IDTight', default=True, help='Trigger option')
    parser.add_argument('--HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60', default=True, help='Trigger option')
    parser.add_argument('--HLT_MET105_IsoTrk50', default=True, help='Trigger option')
    args = parser.parse_args()

    print('\n')
    print('Output directory       : {}'.format(args.outputDir))
    print('Input File       : {}'.format(args.inputFile))
    print('HLT_PFMET120_PFMHT120_IDTight       : {}'.format(args.HLT_PFMET120_PFMHT120_IDTight))
    print('HLT_PFHT500_PFMET100_PFMHT100_IDTight       : {}'.format(args.HLT_PFHT500_PFMET100_PFMHT100_IDTight))
    print('HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60       : {}'.format(args.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60))
    print('HLT_MET105_IsoTrk50       : {}'.format(args.HLT_MET105_IsoTrk50))
    print('\n')

    load_and_run(args.inputFile, args.outputDir, 
                        args.HLT_PFMET120_PFMHT120_IDTight, 
                        args.HLT_PFHT500_PFMET100_PFMHT100_IDTight, 
                        args.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, 
                        args.HLT_MET105_IsoTrk50)
    return 0


if __name__ == '__main__':
  sys.exit(main())