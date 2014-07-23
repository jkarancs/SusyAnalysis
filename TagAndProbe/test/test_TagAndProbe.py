#  ----------------------------------------------------
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 53X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyPatLayer1DefV12

is7X = True
useSusyPat = False if is7X else True
savePatTuple = True

# Standard Pat configuration, geometry and detector conditions
# Can also do:
# from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.ParameterSet.Config as cms
process = cms.Process("PAT")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 53X
        #'file:/data/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/GEN-SIM-RECO/HLT8E33_PU_S10_START53_V7I-v1/20000/F8B40CFD-267A-E211-B13B-0025904B144E.root'
        # 70X
        'file:/data/store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003944B6-CEC7-E311-A9BD-002590A2CD68.root'
    )
)

# Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patTuple.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
    outputCommands = cms.untracked.vstring('keep *')
)

# Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
)

if useSusyPat:
    from PhysicsTools.Configuration.SUSY_pattuple_cff import getSUSY_pattuple_outputCommands
    SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
    process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )
    process.out.fileName = 'susyPatTuple.root'

if savePatTuple: process.outpath = cms.EndPath(process.out) if savePatTuple else cms.EndPath()

# This mode allows to run processes only when another process depends on it
#process.options.allowUnscheduled = cms.untracked.bool( True ) # problematic for matcher

#-- Further Options -----------------------------------
# can also be given from command line
# eg. cmsRun testTagAndProbe.py isMC=True runEle=True runMuon=False
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('standard')

options.register('isMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "MC: True, Data: False")
options.register('runEle',  True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Run on events with Electrons")
options.register('runMuon', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Run on events with Muons")

#   other options for SusyPat
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
    if is7X: process.GlobalTag.globaltag = "POSTLS170_V6::All" if options.isMC else "FT_R_70_V1::All" # 2012D ReReco
    else: process.GlobalTag.globaltag = "START53_V7G::All" if options.isMC else "FT53_V21A_AN6::All"

#   Jet Type and corrections (specify these here instead)
options.register('jetTypes', 'AK5PF', VarParsing.multiplicity.list, VarParsing.varType.string, "Additional jet types that will be produced (AK5Calo and AK5PF, cross cleaned in PF2PAT, are included anyway)")
options.register('jetMetCorrections', 'L1FastJet', VarParsing.multiplicity.list, VarParsing.varType.string, "Level of jet corrections to use: Note the factors are read from DB via GlobalTag")
options.jetMetCorrections.append('L2Relative')
options.jetMetCorrections.append('L3Absolute')
if not options.isMC:
    options.jetMetCorrections.append('L2L3Residual')

# Match triggers to corresponding Pat objects (used by TagAndProbe)
#process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi') # some examples
# Plugin to examine PatTriggers:
#process.patTriggers = cms.EDProducer("PatTriggers",
#    src = cms.InputTag( 'patTrigger' ),
#)
#process.triggerMatcher = cms.EDProducer("PATTriggerMatcherModified", # modified plugin
process.eleTriggerMatcher = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    src = cms.InputTag('cleanPatElectrons'),
    matched = cms.InputTag("patTrigger"),
    matchedCuts = cms.string(
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
    ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

process.muonTriggerMatcher = process.eleTriggerMatcher.clone(
    src = cms.InputTag('cleanPatMuons'),
    matchedCuts = cms.string(
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
    ),
)

#-- Event Filters -------------------------------------
# Trigger Selection
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    HLTPaths = cms.vstring(),
    TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
    throw = False
)
if options.runEle:
    # SingleEle
    process.hltFilter.HLTPaths.append('HLT_Ele27_WP80_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele22_CaloIdL_CaloIsoVL_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele27_CaloIdL_CaloIsoVL_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele30_CaloIdVT_TrkIdT_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele80_CaloIdVT_GsfTrkIdT_v*')
    process.hltFilter.HLTPaths.append('HLT_Ele90_CaloIdVT_GsfTrkIdT_v*')
if options.runMuon:
    # SingleMu
    process.hltFilter.HLTPaths.append('HLT_Mu*')
    process.hltFilter.HLTPaths.append('HLT_IsoMu*')
    # MuHad
    process.hltFilter.HLTPaths.append('HLT_Mu40_PFHT350_v*')
    process.hltFilter.HLTPaths.append('HLT_Mu60_PFHT350_v*')
    process.hltFilter.HLTPaths.append('HLT_Mu40_PFNoPUHT350_v*')
    process.hltFilter.HLTPaths.append('HLT_Mu60_PFNoPUHT350_v*')
    process.hltFilter.HLTPaths.append('HLT_PFHT350_Mu15_PFMET45_v*')
    process.hltFilter.HLTPaths.append('HLT_PFHT350_Mu15_PFMET50_v*')
    process.hltFilter.HLTPaths.append('HLT_PFHT400_Mu5_PFMET45_v*')
    process.hltFilter.HLTPaths.append('HLT_PFHT400_Mu5_PFMET50_v*')
    process.hltFilter.HLTPaths.append('HLT_PFNoPUHT350_Mu15_PFMET45_v*')
    process.hltFilter.HLTPaths.append('HLT_PFNoPUHT350_Mu15_PFMET50_v*')
    process.hltFilter.HLTPaths.append('HLT_PFNoPUHT400_Mu5_PFMET45_v*')
    process.hltFilter.HLTPaths.append('HLT_PFNoPUHT400_Mu5_PFMET50_v*')

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
if not is7X: # CMSSW 53X
    if useSusyPat:
        # Create default susyPattuple
        from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT
        addDefaultSUSYPAT(process,options.isMC,options.hltName,options.jetMetCorrections,options.mcVersion,options.jetTypes,options.doValidation,options.doExtensiveMatching,options.doSusyTopProjection)
        
        # Embed matched trigger objects inside pat objects (used by TagAndProbe)
        from PhysicsTools.PatAlgos.tools.trigTools import switchOnTriggerMatchEmbedding
        switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'eleTriggerMatcher', 'muonTriggerMatcher' ], hltProcess = options.hltName, sequence =  "susyPatDefaultSequence")
    else: # or use patDefaultSequence
        # Use pf Isolation
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPFBasedIsolation#PAT_configuration
        # Needed also for Rho corrected pfIsolation for Electrons
        from PhysicsTools.PatAlgos.tools.pfTools import usePFIso
        usePFIso( process )
        
        # Add Jet Collections
        from PhysicsTools.Configuration.SUSY_pattuple_cff import addSUSYJetCollection
        for jetName in options.jetTypes:
            addSUSYJetCollection(process,options.jetMetCorrections,jetName,options.mcVersion)
        
else: # SusyPat is not available for 7X yet
    
    # pf Isolation is in CMSSW7X by default, but the producers were renamed
    # eg. elPFIsoValueCharged03PFIdPFIso -> elPFIsoValueCharged03PFId (no PFIso at the end)
    
    #-- Remove MC dependencies ----------------------------
    if not options.isMC:
        from PhysicsTools.PatAlgos.tools.coreTools import *
        runOnData(process)
        # Above also adds 'L2L3Residual' corrections to Jet Correction factors
    
    #-- Jet corrections -----------------------------------
    #process.patJetCorrFactors.levels = options.jetMetCorrections
    
    #-- Trigger Match embedding ---------------------------
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger, switchOnTriggerMatchEmbedding
    switchOnTrigger(process, hltProcess=options.hltName)
    switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'eleTriggerMatcher', 'muonTriggerMatcher' ], hltProcess = options.hltName)
    process.patDefaultSequence += process.patTrigger
    process.patDefaultSequence += process.eleTriggerMatcher
    process.patDefaultSequence += process.muonTriggerMatcher
    process.patDefaultSequence += process.cleanPatElectronsTriggerMatch
    process.patDefaultSequence += process.cleanPatMuonsTriggerMatch
    
    #Additional electron ids as defined for VBTF
    process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
    process.simpleEleId80relIso.src = "gedGsfElectrons"
    process.patElectrons.electronIDSources.simpleEleId80relIso = cms.InputTag("simpleEleId80relIso")
    process.patDefaultSequence.replace(process.patElectrons,process.simpleEleId80relIso+process.patElectrons)
    
    ##-- Add jet collections ------------------------------
    #from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    ## lines below taken from (and updated to 7X):
    ## PhysicsTools/Configuration/python/SUSY_pattuple_cff.py and
    ## PhysicsTools/PatAlgos/test/patTuple_addJets_cfg.py
    #for jetName in options.jetTypes:
    #    algo = jetName[0:3]
    #    if 'IC' in algo: collection = algo.replace('IC','iterativeCone')
    #    elif 'SC' in algo: collection = algo.replace('SC','sisCone')
    #    elif 'AK' in algo: collection = algo.replace('AK','ak')
    #    elif 'KT' in algo: collection = algo.replace('KT','kt')
    #    else: raise ValueError, "Unknown jet algo: %s" % (jetName)
    #    jetType = jetName[3:len(jetName)]
    #    if jetType == 'Calo':
    #        jetCollection = '%(collection)sCaloJets' % locals()
    #    elif jetType == 'PF':
    #        jetCollection = '%(collection)sPFJets' % locals()
    #    elif jetType == 'JPT':
    #        if 'IC' in algo: collectionJPT = algo.replace('IC','Icone')
    #        elif 'SC' in algo: collectionJPT = algo.replace('SC','Siscone')
    #        elif 'AK' in algo: collectionJPT = algo.replace('AK','AntiKt')
    #        else: raise ValueError, "Unknown jet algo: %s" % (jetName)
    #        jetCollection = 'JetPlusTrackZSPCorJet%(collectionJPT)s' % locals()
    #    elif jetType == 'Track':
    #        jetCollection = '%(collection)sTrackJets' % locals()
    #    else: raise ValueError, "Unknown jet jetType: %s" % (jetName)
    #    
    #    addJetCollection(
    #        process,
    #        labelName = jetName,
    #        jetSource = cms.InputTag(jetCollection),
    #        jetCorrections = (jetName, cms.vstring(options.jetMetCorrections), 'Type-1'),
    #        genJetCollection = cms.InputTag('%(collection)sGenJets' % locals()),
    #        btagDiscriminators = [
    #            'jetBProbabilityBJetTags',
    #            'jetProbabilityBJetTags',
    #            'trackCountingHighPurBJetTags',
    #            'trackCountingHighEffBJetTags',
    #            'simpleSecondaryVertexHighEffBJetTags',
    #            'simpleSecondaryVertexHighPurBJetTags',
    #            'combinedSecondaryVertexBJetTags'
    #        ],
    #    )
    #    
    #    process.patCandidateSummary.candidates.append(cms.InputTag('patJets'+jetName))
    #    process.selectedPatCandidateSummary.candidates.append(cms.InputTag('selectedPatJets'+jetName))
    #    process.cleanPatCandidateSummary.candidates.append(cms.InputTag('cleanPatJets'+jetName))
    #    
    #    print "Jet collection added: "+jetCollection
    #
    #    #-- PF2PAT -------------------------------------------
    #    #from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
    #    ## for charged hadron subtraction (PFnoPU) use: jetName + "chs"
    #    #usePF2PAT(
    #    #    process,
    #    #    runPF2PAT=True,
    #    #    jetAlgo=algo,
    #    #    runOnMC=options.isMC,
    #    #    postfix="PF",
    #    #    jetCorrections=(jetName, cms.vstring(options.jetMetCorrections)),
    #    #    typeIMetCorrections=True
    #    #)

#-- MET Corrections -----------------------------------
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
#if options.isMC:
#  process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3"
#  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc
#else:
#  process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3Residual"
#  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
#process.patPFMETs = process.patMETs.clone(
#    metSource = cms.InputTag('pfMet'),
#    addMuonCorrections = cms.bool(False),
#    #genMETSource = cms.InputTag('genMetTrue'),
#    #addGenMET = cms.bool(True)
#)
#process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
#process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
#    cms.InputTag('pfJetMETcorr', 'type1') ,
#    cms.InputTag('pfMEtSysShiftCorr')
#)
#process.patPFMETsTypeIcorrected = process.patPFMETs.clone(
#    metSource = cms.InputTag('pfType1CorrectedMet'),
#)
#
#process.metCorrectionSequence = cms.Sequence(
#    process.pfMEtSysShiftCorrSequence *
#    process.producePFMETCorrections *
#    process.patPFMETsTypeIcorrected
#)


#-- TagAndProbe ---------------------------------------
if options.runEle:
    from SusyAnalysis.TagAndProbe.TagAndProbe_Electrons import *
    initTP_Electrons(process, options.isMC, options.hltName, is7X)
if options.runMuon:
    from SusyAnalysis.TagAndProbe.TagAndProbe_Muons import *
    initTP_Muons(process, options.isMC, options.hltName, is7X)

dataSuffix = "_MC" if options.isMC else "_Data" 
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService = cms.Service("TFileService", fileName = cms.string("TNP_Ntuple"+dataSuffix+".root"))

#-- Path ----------------------------------------------
process.p = cms.Path(
    process.filterSequence*
    process.patDefaultSequence
    #process.metCorrectionSequence *
)
if useSusyPat: process.p.replace(process.patDefaultSequence,process.susyPatDefaultSequence)

if options.runEle: process.p += process.TagAndProbe_Electrons
if options.runMuon: process.p += process.TagAndProbe_Muons
