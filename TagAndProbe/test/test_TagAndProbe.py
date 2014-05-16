#  ----------------------------------------------------
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 53X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyPatLayer1DefV12

# Import Default Pat Settings
from PhysicsTools.PatAlgos.patTemplate_cfg import *
# Use pf Isolation
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPFBasedIsolation#PAT_configuration
# Needed also for Rho corrected pfIsolation for Electrons
#from PhysicsTools.PatAlgos.tools.pfTools import *
#usePFIso( process )

# remove MC matching from the default sequence
#removeMCMatching(process, ['Muons'])
#runOnData(process)

#-- Options -------------------------------------------
# can also be given from command line
# eg. cmsRun testTagAndProbe.py isMC=True isEle=False
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('standard')

savePatTuple = False

options.register('isMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "MC: True, Data: False")
options.register('isEle', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Ele: True, Muon: False")

#   other options
options.register('mcVersion', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Currently not needed and supported")
options.register('doValidation', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Include the validation histograms from SusyDQM (needs extra tags)")
options.register('doExtensiveMatching', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Matching to simtracks (needs extra tags)")
options.register('doSusyTopProjection', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Apply Susy selection in PF2PAT to obtain lepton cleaned jets (needs validation)")
options.register('hltName', 'HLT', VarParsing.multiplicity.singleton, VarParsing.varType.string, "HLT menu to use for trigger matching")

# GlobalTag
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "GlobalTag to use (auto: selected automatically, empty: given from cfg)")

options.parseArguments() #Read above from command line

if options.globalTag == 'auto':
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
elif options.globalTag:
    process.GlobalTag.globaltag = options.globalTag
else:
    process.GlobalTag.globaltag = "START53_V7G::All" if options.isMC else "FT53_V21A_AN6::All"

#   Jet Type and corrections (specify these here instead)
options.register('jetTypes', 'AK5PF', VarParsing.multiplicity.list, VarParsing.varType.string, "Additional jet types that will be produced (AK5Calo and AK5PF, cross cleaned in PF2PAT, are included anyway)")
options.register('jetCorrections', 'L1FastJet', VarParsing.multiplicity.list, VarParsing.varType.string, "Level of jet corrections to use: Note the factors are read from DB via GlobalTag")
options.jetCorrections.append('L2Relative')
options.jetCorrections.append('L3Absolute')
if not options.isMC:
    options.jetCorrections.append('L2L3Residual')

# Input files
process.source.fileNames = [
    'file:/data/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/GEN-SIM-RECO/HLT8E33_PU_S10_START53_V7I-v1/20000/F8B40CFD-267A-E211-B13B-0025904B144E.root'
    #'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/GEN-SIM-RECO/HLT8E33_PU_S10_START53_V7I-v1/20000/F8B40CFD-267A-E211-B13B-0025904B144E.root'
]

# Number of Events
process.maxEvents.input = -1

# Message Logger
process.options.wantSummary = False # Get Summary in the end of job
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# High Level Trigger Selections
# - to filter Events (with OR)
# - match triggers to corresponding Pat objects (used by TagAndProbe)
# Plugin to examine PatTriggers:
#process.patTriggers = cms.EDProducer("PatTriggers",
#    src = cms.InputTag( 'patTrigger' ),
#)
if options.isEle:
    # 2012 SingleElectron triggers
    hltFilterSelection = cms.vstring(
        'HLT_Ele22_CaloIdL_CaloIsoVL_v*',
        'HLT_Ele27_CaloIdL_CaloIsoVL_v*',
        'HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
        'HLT_Ele30_CaloIdVT_TrkIdT_v*',
        'HLT_Ele80_CaloIdVT_GsfTrkIdT_v*',
        'HLT_Ele90_CaloIdVT_GsfTrkIdT_v*',
        #'HLT_Ele27_WP80_v*',
    )
    patObject_src = cms.InputTag('cleanPatElectrons')
    triggerObjectCuts = cms.string(
        'type( "TriggerElectron" ) && ('
        '  path("HLT_Ele22_CaloIdL_CaloIsoVL_v*"                                          ,1,0) || '
        '  path("HLT_Ele27_CaloIdL_CaloIsoVL_v*"                                          ,1,0) || '
        '  path("HLT_Ele30_CaloIdVT_TrkIdT_v*"                                            ,1,0) || '
        '  path("HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                         ,1,0) || '
        '  path("HLT_Ele80_CaloIdVT_GsfTrkIdT_v*"                                         ,1,0) || '
        '  path("HLT_Ele90_CaloIdVT_GsfTrkIdT_v*"                                         ,1,0) '
        # WP80 trigger conflicts with others, so it should only be selected alone
        #'  || path("HLT_Ele27_WP80_v*"                                                       ,1,0) '
        ')'
    )
else:
    # 2012 SingleMu (and MuHad) Triggers
    hltFilterSelection = cms.vstring(
        # SingleMu
        'HLT_Mu*',
        'HLT_IsoMu*',
        # MuHad
        'HLT_Mu40_PFHT350_v*',
        'HLT_Mu60_PFHT350_v*',
        'HLT_Mu40_PFNoPUHT350_v*',
        'HLT_Mu60_PFNoPUHT350_v*',
        'HLT_PFHT350_Mu15_PFMET45_v*',
        'HLT_PFHT350_Mu15_PFMET50_v*',
        'HLT_PFHT400_Mu5_PFMET45_v*',
        'HLT_PFHT400_Mu5_PFMET50_v*',
        'HLT_PFNoPUHT350_Mu15_PFMET45_v*',
        'HLT_PFNoPUHT350_Mu15_PFMET50_v*',
        'HLT_PFNoPUHT400_Mu5_PFMET45_v*',
        'HLT_PFNoPUHT400_Mu5_PFMET50_v*',
    )
    patObject_src = cms.InputTag('cleanPatMuons')
    triggerObjectCuts  = cms.string(
        'type("TriggerMuon") && ('
        # Mu*
        '  (hasFilterLabel("hltL3fL1sMu16L1f0L2f16QL3Filtered24Q") && path("HLT_Mu24_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16L1f0L2f16QL3Filtered30Q") && path("HLT_Mu30_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16L1f0L2f16QL3Filtered40Q") && path("HLT_Mu40_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q") && path("HLT_Mu24_eta2p1_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered30Q") && path("HLT_Mu30_eta2p1_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q") && path("HLT_Mu40_eta2p1_v*",0,0)) || '
        '  (hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered50Q") && path("HLT_Mu50_eta2p1_v*",0,0)) || '
        # IsoMu*
        '  (hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15") && path("HLT_IsoMu24_v*",0,0)) || '
        '  (hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15") && path("HLT_IsoMu30_v*",0,0)) || '
        '  ((hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20L3crIsoFiltered10") || hasFilterLabel(" hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20L3crIsoRhoFiltered0p15")) && path("HLT_IsoMu20_eta2p1_v*",0,0)) || '
        '  ((hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10") || hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15")) && path("HLT_IsoMu24_eta2p1_v*",0,0)) || '
        '  ((hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoFiltered10") || hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f30QL3crIsoRhoFiltered0p15")) && path("HLT_IsoMu30_eta2p1_v*",0,0)) || '
        '  ((hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoFiltered10") || hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f34QL3crIsoRhoFiltered0p15")) && path("HLT_IsoMu34_eta2p1_v*",0,0)) || '
        '  ((hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoFiltered10") || hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15")) && path("HLT_IsoMu40_eta2p1_v*",0,0)) '
        # Hadronic Cross Triggers (Mu + HT,MET)
        '  || '
        '  (hasFilterLabel("hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered40") && path("HLT_Mu40_PFHT350_v*",0,0)) || '
        '  (hasFilterLabel("hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered40") && path("HLT_Mu40_PFNoPUHT350_v*",0,0)) || '
        '  (hasFilterLabel("hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered60") && path("HLT_Mu60_PFHT350_v*",0,0)) || '
        '  (hasFilterLabel("hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered60") && path("HLT_Mu60_PFNoPUHT350_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered15") && path("HLT_PFHT350_Mu15_PFMET45_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered15") && path("HLT_PFNoPUHT350_Mu15_PFMET45_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered15") && path("HLT_PFHT350_Mu15_PFMET50_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered15") && path("HLT_PFNoPUHT350_Mu15_PFMET50_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered5") && path("HLT_PFHT400_Mu5_PFMET45_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered5") && path("HLT_PFNoPUHT400_Mu5_PFMET45_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered5") && path("HLT_PFHT400_Mu5_PFMET50_v*",0,0)) || '
        '  (hasFilterLabel("hltL1HTT150singleMuL3PreFiltered5") && path("HLT_PFNoPUHT400_Mu5_PFMET50_v*",0,0))'
        ')'
    )

# Trigger Matcher
#process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi') # examples
process.triggerMatcher = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
#process.triggerMatcher = cms.EDProducer("PATTriggerMatcherModified", # modified plugin
    src     = patObject_src,
    matched = cms.InputTag("patTrigger"),
    matchedCuts = triggerObjectCuts,
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

#-- Event Filters -------------------------------------
# Trigger Selection
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    HLTPaths = hltFilterSelection,
    TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
    throw = False
)

# Vertex filter (also needed for trackingFailureFilter)
process.goodVertices = cms.EDFilter("VertexSelector",
    filter = cms.bool(True),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

# Beam background (monster) removal
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

# MET Filters (Same as that previously used by RA4), copied from:
# process.load('RecoMET.METFilters.metFilters_cff') # load doesn't work for VertexSelector
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
process.load('RecoMET.METFilters.trackingPOGFilters_cff')
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.metFilters = cms.Sequence(
    process.HBHENoiseFilter *
    process.CSCTightHaloFilter *
    process.hcalLaserEventFilter *
    process.EcalDeadCellTriggerPrimitiveFilter *
    process.trackingFailureFilter *
    process.eeBadScFilter *
    process.ecalLaserCorrFilter *
    process.trkPOGFilters
)
process.filterSequence = cms.Sequence(
    process.hltFilter*
    process.goodVertices *
    process.scrapingVeto *
    process.metFilters
)

#-- SusyPat Specifics ---------------------------------
# Apply settings to SusyPat
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
addDefaultSUSYPAT(process,options.isMC,options.hltName,options.jetCorrections,options.mcVersion,options.jetTypes,options.doValidation,options.doExtensiveMatching,options.doSusyTopProjection)

# Embed matched trigger objects inside pat objects (used by TagAndProbe)
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, hltProcess = options.hltName, sequence =  "susyPatDefaultSequence") # This is optional and can be omitted.
switchOnTriggerMatchEmbedding( process, [ 'triggerMatcher' ], hltProcess = options.hltName, sequence =  "susyPatDefaultSequence")

# Save PatTuple to a separate file
if savePatTuple:
    process.out.fileName = 'susyPatTuple.root'
    SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
    process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )
else:
    process.outpath = cms.EndPath()

#-- MET Corrections -----------------------------------
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
if options.isMC:
  process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3"
  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc
else:
  process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3Residual"
  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
process.patPFMETs = process.patMETs.clone(
    metSource = cms.InputTag('pfMet'),
    addMuonCorrections = cms.bool(False),
    #genMETSource = cms.InputTag('genMetTrue'),
    #addGenMET = cms.bool(True)
)
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1') ,
    cms.InputTag('pfMEtSysShiftCorr')
)
process.patPFMETsTypeIcorrected = process.patPFMETs.clone(
    metSource = cms.InputTag('pfType1CorrectedMet'),
)

process.metCorrectionSequence = cms.Sequence(
    process.pfMEtSysShiftCorrSequence *
    process.producePFMETCorrections *
    process.patPFMETsTypeIcorrected
)

#-- TagAndProbe ---------------------------------------
if options.isEle:
    from SusyAnalysis.TagAndProbe.TagAndProbe_Electrons import *
else:
    from SusyAnalysis.TagAndProbe.TagAndProbe_Muons import *
initTP(process,options.isMC,options.hltName)  
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")

lepSuffix  = "_Ele" if options.isEle else "_Mu"
dataSuffix = "_MC" if options.isMC else "_Data" 
process.TFileService = cms.Service("TFileService", fileName = cms.string("TNP_Ntuple"+dataSuffix+lepSuffix+".root"))

#-- Path ----------------------------------------------
process.p = cms.Path(
    process.filterSequence*
    process.susyPatDefaultSequence *
    process.metCorrectionSequence *
    process.TagAndProbe
)
