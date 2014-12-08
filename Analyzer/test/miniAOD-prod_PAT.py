# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD-prod -s PAT --eventcontent MINIAODSIM --runUnscheduled --mc --filein /store/relval/CMSSW_7_0_0/RelValTTbar_13/GEN-SIM-RECO/PU50ns_POSTLS170_V4-v2/00000/265B9219-FF98-E311-BF4A-02163E00EA95.root --conditions PLS170_V6AN1::All --no_exec
import FWCore.ParameterSet.Config as cms

#-- Options -------------------------------------------
# can also be given from command line
# eg. cmsRun Analyzer.py useMiniAOD=True
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('standard')

options.register('useMiniAOD', False,  VarParsing.multiplicity.singleton, VarParsing.varType.bool,   "MiniAOD: True, AOD: False (Default)")
options.register('isMC',       True,   VarParsing.multiplicity.singleton, VarParsing.varType.bool,   "MC: True (Default), Data: False")
options.register('globalTag',  'auto', VarParsing.multiplicity.singleton, VarParsing.varType.string, "GlobalTag to use (auto (Default): selected automatically, empty: given from cfg)")
options.register('nevt',       -1,     VarParsing.multiplicity.singleton, VarParsing.varType.int,    "Spefify number of events: - Default = -1 (all events)")

options.parseArguments() #Read above from command line

# From MiniAOD config file
if options.useMiniAOD:
    process = cms.Process("ANA")
else:
    process = cms.Process("PAT")

# import of standard configurations
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

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PLS170_V6AN1::All', '')

#process = cms.Process('PAT')
#
## import of standard configurations
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/data/store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003944B6-CEC7-E311-A9BD-002590A2CD68.root')
)

#process.options = cms.untracked.PSet(
#    allowUnscheduled = cms.untracked.bool(True)
#)

# # Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     version = cms.untracked.string('$Revision: 1.19 $'),
#     annotation = cms.untracked.string('miniAOD-prod nevts:1'),
#     name = cms.untracked.string('Applications')
# )
# 
# # Output definition
# 
# process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
#     compressionLevel = cms.untracked.int32(4),
#     compressionAlgorithm = cms.untracked.string('LZMA'),
#     eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#     outputCommands = process.MINIAODSIMEventContent.outputCommands,
#     fileName = cms.untracked.string('miniAOD-prod_PAT.root'),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string(''),
#         dataTier = cms.untracked.string('')
#     ),
#     dropMetaData = cms.untracked.string('ALL'),
#     fastCloning = cms.untracked.bool(False),
#     overrideInputFileSplitLevels = cms.untracked.bool(True)
# )

# Additional output definition

# Path and EndPath definitions
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
#process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

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

process.p = cms.Path(
    process.Analyzer
)
# End of customisation functions
