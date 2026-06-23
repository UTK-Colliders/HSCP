# Create MET efficiency plots from ntuples

import argparse
import os
import sys
import ROOT


def load_and_run(inputFile, outputDir):
    cmssw_base = os.environ['CMSSW_BASE']
    source_path = os.path.join(cmssw_base, 'src', 'HSCP', 'MET_efficiency.cc')
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load('libTree')
    ROOT.gROOT.ProcessLine(f'.L {source_path}+')
    ROOT.MET_efficiency(inputFile, outputDir)


def main():
    parser = argparse.ArgumentParser(description='Create MET efficiency plots from ntuples')
    parser.add_argument('--inputFile', required=True, help='ntuple file')
    parser.add_argument('--outputDir', required=True, help='Output directory for plots')
    args = parser.parse_args()

    print('\n')
    print('Output directory       : {}'.format(args.outputDir))
    print('Input File       : {}'.format(args.inputFile))
    print('\n')

    load_and_run(args.inputFile, args.outputDir)
    return 0


if __name__ == '__main__':
  sys.exit(main())