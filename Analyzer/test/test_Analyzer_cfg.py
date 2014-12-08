import FWCore.ParameterSet.Config as cms

#-- Options -------------------------------------------
# can also be given from command line
# eg. cmsRun Analyzer.py useMiniAOD=True
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('standard')

options.register('useMiniAOD', False,  VarParsing.multiplicity.singleton, VarParsing.varType.bool,   "MiniAOD: True, AOD: False (Default)")
options.register('isMC',       True,   VarParsing.multiplicity.singleton, VarParsing.varType.bool,   "MC: True (Default), Data: False")
#options.register('globalTag',  'auto', VarParsing.multiplicity.singleton, VarParsing.varType.string, "GlobalTag to use (auto (Default): selected automatically, empty: given from cfg)")
options.register('globalTag',  '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "GlobalTag to use (auto: selected automatically, empty: given from cfg (Default))")
options.register('nevt',       -1,     VarParsing.multiplicity.singleton, VarParsing.varType.int,    "Spefify number of events: - Default = -1 (all events)")
options.parseArguments()

# From MiniAOD config file
if options.useMiniAOD:
    process = cms.Process("ANA")
else:
    process = cms.Process("PAT")

# Import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if options.isMC:
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')
    if not options.useMiniAOD:
        process.load('Configuration.StandardSequences.PATMC_cff')
else:
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    if not options.useMiniAOD:
        process.load('Configuration.StandardSequences.PAT_cff')

# Input file
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
# Normal AOD
if not options.useMiniAOD:
    process.source.fileNames= [
        "file:/data/store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003944B6-CEC7-E311-A9BD-002590A2CD68.root"
    ]
# MiniAOD
else:
    process.source.fileNames= [
        #"file:/data/jkarancs/CMSSW/SusyAnalysis/CMSSW_7_0_6_patch1/src/SusyAnalysis/Analyzer/test/miniAOD-prod_PAT.root" # generated from AOD (above)
        # signal sample
        "file:/data/store/mc/Spring14miniaod/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/26155C45-520E-E411-9910-20CF3027A624.root",
        # background samples
        #"file:/data/store/mc/Spring14miniaod/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/063013AD-9907-E411-8135-0026189438BD.root",
        #"file:/data/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/00245653-EC23-E411-B9BF-02163E006C73.root",
        #"file:/data/store/mc/Spring14miniaod/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/00BC2643-2709-E411-9578-90B11C08C6BF.root",
        #"file:/data/store/mc/Spring14miniaod/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/005B02AC-DE18-E411-A3C1-0025905A48EC.root",
    ]

# Set GlobalTag
if options.globalTag == 'auto':
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
elif options.globalTag:
    process.GlobalTag.globaltag = options.globalTag
else:
    process.GlobalTag.globaltag = "PLS170_V6AN1::All" if options.isMC else "FT_R_70_V1::All" # 2012D ReReco

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nevt) )

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

#-- Event Filters -------------------------------------
# Trigger Selection
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    HLTPaths = cms.vstring(),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = False
)
# SingleEle
process.hltFilter.HLTPaths.append('HLT_Ele27_WP80_v*')
# SingleMu
process.hltFilter.HLTPaths.append('HLT_IsoMu24_v*')
# MuHad
#process.hltFilter.HLTPaths.append('HLT_Mu40_PFHT350_v*')
#process.hltFilter.HLTPaths.append('HLT_Mu60_PFHT350_v*')
#process.hltFilter.HLTPaths.append('HLT_Mu40_PFNoPUHT350_v*')
#process.hltFilter.HLTPaths.append('HLT_Mu60_PFNoPUHT350_v*')
#process.hltFilter.HLTPaths.append('HLT_PFHT350_Mu15_PFMET45_v*')
#process.hltFilter.HLTPaths.append('HLT_PFHT350_Mu15_PFMET50_v*')
#process.hltFilter.HLTPaths.append('HLT_PFHT400_Mu5_PFMET45_v*')
#process.hltFilter.HLTPaths.append('HLT_PFHT400_Mu5_PFMET50_v*')
#process.hltFilter.HLTPaths.append('HLT_PFNoPUHT350_Mu15_PFMET45_v*')
#process.hltFilter.HLTPaths.append('HLT_PFNoPUHT350_Mu15_PFMET50_v*')
#process.hltFilter.HLTPaths.append('HLT_PFNoPUHT400_Mu5_PFMET45_v*')
#process.hltFilter.HLTPaths.append('HLT_PFNoPUHT400_Mu5_PFMET50_v*')

#-- Analyzer ------------------------------------------
process.Analyzer = cms.EDAnalyzer("Analyzer",
    vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons     = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus      = cms.InputTag("slimmedTaus"),
    photons   = cms.InputTag("slimmedPhotons"),
    jets      = cms.InputTag("slimmedJets"),
    #fatjets   = cms.InputTag("slimmedJetsAK8"),
    mets      = cms.InputTag("slimmedMETs"),
)

#-- Path ----------------------------------------------
if not options.useMiniAOD:
    process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
    process.Flag_goodVertices = cms.Path(process.goodVertices)
    process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
    process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
    process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
    process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
    process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
    process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
    process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
    process.Flag_METFilters = cms.Path(process.metFilters)
    process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilter)
    process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
    process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
    process.endjob_step = cms.EndPath(process.endOfProcess)
    if options.isMC:
        from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 
        process = miniAOD_customizeAllMC(process)
    else:  
        from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
        process = miniAOD_customizeAllData(process)

process.p = cms.Path(
    #process.hltFilter*
    #process.metCorrectionSequence*
    process.Analyzer
)

#if not options.useMiniAOD: process.p.insert(0, process.patDefaultSequence)


