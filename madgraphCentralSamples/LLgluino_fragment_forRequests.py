import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *

gluino_mass = 1000
tau = "100ps"
neutralino_mass = 100
decay = "uds"

dirhadrongenfilter = cms.EDFilter("MCParticlePairFilter",
    MaxEta = cms.untracked.vdouble(100.0, 100.0),
    MinEta = cms.untracked.vdouble(-100, -100),
    MinP = cms.untracked.vdouble(0.0, 0.0),
    MinPt = cms.untracked.vdouble(0.0, 0.0),
    ParticleCharge = cms.untracked.int32(0),
    ParticleID1 = cms.untracked.vint32( 
        1000993, 1009213, 1009313, 1009323, 1009113, 
        1009223, 1009333, 1091114, 1092114, 1092214, 
        1092224, 1093114, 1093214, 1093224, 1093314, 
        1093324, 1093334
    ),
    ParticleID2 = cms.untracked.vint32(
        1000993, 1009213, 1009313, 1009323, 1009113, 
        1009223, 1009333, 1091114, 1092114, 1092214, 
        1092224, 1093114, 1093214, 1093224, 1093314, 
        1093324, 1093334
    ),
    Status = cms.untracked.vint32(1, 1)
)

generator = cms.EDFilter("Pythia8ConcurrentHadronizerFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CP5Settings', 
            'processParameters'
        ),
        processParameters = cms.vstring(
            'RHadrons:allow = on', 
            'RHadrons:allowDecay = off', 
            'RHadrons:setMasses = on', 
            'RHadrons:probGluinoball = 0.1',

            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 5.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 140.', #this is the actual merging scale
            'JetMatching:clFact = 1', # determines jet-to parton matching
            'JetMatching:nQmatch = 5', #5 for 5-flavour scheme (matching of b-quarks)
            'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity
            'JetMatching:doShowerKt = off'
        ),
        pythia8CP5Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:ecmPow=0.03344', 
            'MultipartonInteractions:bProfile=2', 
            'MultipartonInteractions:pT0Ref=1.41', 
            'MultipartonInteractions:coreRadius=0.7634', 
            'MultipartonInteractions:coreFraction=0.63', 
            'ColourReconnection:range=5.176', 
            'SigmaTotal:zeroAXB=off', 
            'SpaceShower:alphaSorder=2', 
            'SpaceShower:alphaSvalue=0.118', 
            'SigmaProcess:alphaSvalue=0.118', 
            'SigmaProcess:alphaSorder=2', 
            'MultipartonInteractions:alphaSvalue=0.118', 
            'MultipartonInteractions:alphaSorder=2', 
            'TimeShower:alphaSorder=2', 
            'TimeShower:alphaSvalue=0.118', 
            'SigmaTotal:mode = 0', 
            'SigmaTotal:sigmaEl = 21.89', 
            'SigmaTotal:sigmaTot = 100.309', 
            'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118'
        ),
        pythia8CommonSettings = cms.vstring(
            'Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            '1000021:mayDecay = off'
        )
    ),
    SLHAFileForPythia8 = cms.string(f'Configuration/Generator/data/SUSY/LLGluino/LLgluino_M-{gluino_mass}_tau-{tau}_{decay}QuarkDecay_chi10_M-{neutralino_mass}.slha'),
    comEnergy = cms.double(13600.0),
    crossSection = cms.untracked.double(-1),
    hscpFlavor = cms.untracked.string('gluino'),
    massPoint = cms.untracked.int32(gluino_mass),
    maxEventsToPrint = cms.untracked.int32(0),
    particleFile = cms.untracked.string(f'Configuration/Generator/data/particles_gluino_{gluino_mass}_GeV.txt'),
    pdtFile = cms.FileInPath(f'Configuration/Generator/data/hscppythiapdtgluino{gluino_mass}.tbl'),
    processFile = cms.untracked.string('SimG4Core/CustomPhysics/data/RhadronProcessList.txt'),
    useregge = cms.bool(False),
)


ProductionFilterSequence = cms.Sequence(generator * dirhadrongenfilter)
from SimG4Core.CustomPhysics.Exotica_HSCP_SIM_cfi import customise 
