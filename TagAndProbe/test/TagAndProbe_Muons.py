##  *****************************************************************************
##  Code taken from PhysicsTools/TagAndProbe
##  and modified to run with SUSYPat recipie for Lepton Efficiency studies
##  by Janos Karancsi (ATOMKI Debrecen/University of Debrecen)
##  and Viktor Veszpremi (WIGNER Budapest)
##  *****************************************************************************

import FWCore.ParameterSet.Config as cms

def initTP(process, mcInfo=False, hltName="HLT"):
    
    ################################################################################################################################
    #   ___ ____ ____    ____ _  _ ___     ___  ____ ____ ___  ____ 
    #    |  |__| | __    |__| |\ | |  \    |__] |__/ |  | |__] |___ 
    #    |  |  | |__]    |  | | \| |__/    |    |  \ |__| |__] |___ 
    ################################################################################################################################
    
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
    process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
    process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
    process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
    
    ################################################################################################################################
    #   ____ _  _ ___ ____ 
    #   |    |  |  |  [__  
    #   |___ |__|  |  ___] 
    
    TAG_CUTS       = "pt > 10 && abs(eta) < 2.4 && isGlobalMuon() && isPFMuon() && numberOfMatchedStations() > 1  && abs(track.d0) < 2 && abs(track.dz) < 24"
    
    STA_PROBE_CUTS = "pt > 5 && abs(eta) < 2.4"
    TRK_PROBE_CUTS = "pt > 5 && abs(eta) < 2.4 && numberOfValidHits > 5"
    GLB_PROBE_CUTS = "abs(eta) < 2.4 && isGood('GlobalMuonPromptTight') && isPFMuon() && track().hitPattern().trackerLayersWithMeasurement() > 5 && numberOfMatchedStations() >= 2 && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
    
    JET_CUTS       =  "pt() >= 40.0 && abs(eta()) <= 2.4 && neutralHadronEnergyFraction() < 0.99 && neutralEmEnergyFraction() < 0.99 && chargedHadronEnergyFraction() > 0  && chargedMultiplicity() > 0  && chargedEmEnergyFraction() < 0.99 && (chargedMultiplicity() + neutralMultiplicity() + muonMultiplicity()) > 1"
    
    # Trigger Matches and Passes
    # SingleMu (2012)
    PASS_Mu24                  ="!triggerObjectMatchesByPath('HLT_Mu24_v*',1,0).empty()"
    PASS_Mu30                  ="!triggerObjectMatchesByPath('HLT_Mu30_v*',1,0).empty()"
    PASS_Mu40                  ="!triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty()"
    PASS_Mu24_eta2p1           ="!triggerObjectMatchesByPath('HLT_Mu24_eta2p1_v*',1,0).empty()"
    PASS_Mu30_eta2p1           ="!triggerObjectMatchesByPath('HLT_Mu30_eta2p1_v*',1,0).empty()"
    PASS_Mu40_eta2p1           ="!triggerObjectMatchesByPath('HLT_Mu40_eta2p1_v*',1,0).empty()"
    PASS_Mu50_eta2p1           ="!triggerObjectMatchesByPath('HLT_Mu50_eta2p1_v*',1,0).empty()"
    PASS_IsoMu24               ="!triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty()"
    PASS_IsoMu30               ="!triggerObjectMatchesByPath('HLT_IsoMu30_v*',1,0).empty()"
    PASS_IsoMu20_eta2p1        ="!triggerObjectMatchesByPath('HLT_IsoMu20_eta2p1_v*',1,0).empty()"
    PASS_IsoMu24_eta2p1        ="!triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty()"
    PASS_IsoMu30_eta2p1        ="!triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty()"
    PASS_IsoMu34_eta2p1        ="!triggerObjectMatchesByPath('HLT_IsoMu34_eta2p1_v*',1,0).empty()"
    PASS_IsoMu40_eta2p1        ="!triggerObjectMatchesByPath('HLT_IsoMu40_eta2p1_v*',1,0).empty()"
    
    # MuHad (2012)
    PASS_Mu40_PFHT350            ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_Mu40_PFHT350_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered40').empty()"
    PASS_Mu40_PFNoPUHT350        ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_Mu40_PFNoPUHT350_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered40').empty()"
    PASS_Mu60_PFHT350            ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_Mu60_PFHT350_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered60').empty()"
    PASS_Mu60_PFNoPUHT350        ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_Mu60_PFNoPUHT350_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1Mu0HTT100ORL1Mu4HTT125L2QualL3MuFiltered60').empty()"
    PASS_PFHT350_Mu15_PFMET45    ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFHT350_Mu15_PFMET45_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered15').empty()"
    PASS_PFNoPUHT350_Mu15_PFMET45="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFNoPUHT350_Mu15_PFMET45_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered15').empty()"
    PASS_PFHT350_Mu15_PFMET50    ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFHT350_Mu15_PFMET50_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered15').empty()"
    PASS_PFNoPUHT350_Mu15_PFMET50="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFNoPUHT350_Mu15_PFMET50_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered15').empty()"
    PASS_PFHT400_Mu5_PFMET45     ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFHT400_Mu5_PFMET45_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered5').empty()"
    PASS_PFNoPUHT400_Mu5_PFMET45 ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFNoPUHT400_Mu5_PFMET45_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered5').empty()"
    PASS_PFHT400_Mu5_PFMET50     ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFHT400_Mu5_PFMET50_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered5').empty()"
    PASS_PFNoPUHT400_Mu5_PFMET50 ="!triggerObjectMatchesByType('TriggerMuon').empty() && !triggerObjectMatchesByPath('HLT_PFNoPUHT400_Mu5_PFMET50_v*',0,0).empty() && !triggerObjectMatchesByFilter('hltL1HTT150singleMuL3PreFiltered5').empty()"
    
    PASS_SINGLEELE  = "!triggerObjectMatchesByType('TriggerMuon').empty() && ("
    PASS_SINGLEELE += PASS_Mu5 + "||"
    PASS_SINGLEELE += PASS_Mu8 + "||"
    PASS_SINGLEELE += PASS_Mu12 + "||"
    PASS_SINGLEELE += PASS_Mu17 + "||"
    PASS_SINGLEELE += PASS_Mu24 + "||"
    PASS_SINGLEELE += PASS_Mu30 + "||"
    PASS_SINGLEELE += PASS_Mu40 + "||"
    PASS_SINGLEELE += PASS_Mu15_eta2p1 + "||"
    PASS_SINGLEELE += PASS_Mu24_eta2p1 + "||"
    PASS_SINGLEELE += PASS_Mu30_eta2p1 + "||"
    PASS_SINGLEELE += PASS_Mu40_eta2p1 + "||"
    PASS_SINGLEELE += PASS_Mu50_eta2p1 + "||"
    PASS_SINGLEELE += PASS_IsoMu24 + "||"
    PASS_SINGLEELE += PASS_IsoMu30 + "||"
    PASS_SINGLEELE += PASS_IsoMu20_eta2p1 + "||"
    PASS_SINGLEELE += PASS_IsoMu24_eta2p1 + "||"
    PASS_SINGLEELE += PASS_IsoMu30_eta2p1 + "||"
    PASS_SINGLEELE += PASS_IsoMu34_eta2p1 + "||"
    PASS_SINGLEELE += PASS_IsoMu40_eta2p1 + ")"
    
    PASS_ANY  = "(" + PASS_SINGLEELE
    PASS_ANY += "||"+ "(" + PASS_Mu40_PFHT350 + ")"
    PASS_ANY += "||"+ "(" + PASS_Mu40_PFNoPUHT350 + ")"
    PASS_ANY += "||"+ "(" + PASS_Mu60_PFHT350 + ")"
    PASS_ANY += "||"+ "(" + PASS_Mu60_PFNoPUHT350 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFHT350_Mu15_PFMET45 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFNoPUHT350_Mu15_PFMET45 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFHT350_Mu15_PFMET50 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFNoPUHT350_Mu15_PFMET50 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFHT400_Mu5_PFMET45 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFNoPUHT400_Mu5_PFMET45 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFHT400_Mu5_PFMET50 + ")"
    PASS_ANY += "||"+ "(" + PASS_PFNoPUHT400_Mu5_PFMET50  + ")" + ")"
    
    AllTriggerFlags = cms.PSet(
        passing_HLT_Mu5                      = cms.string(PASS_Mu5),
        passing_HLT_Mu8                      = cms.string(PASS_Mu8),
        passing_HLT_Mu12                     = cms.string(PASS_Mu12),
        passing_HLT_Mu17                     = cms.string(PASS_Mu17),
        passing_HLT_Mu24                     = cms.string(PASS_Mu24),
        passing_HLT_Mu30                     = cms.string(PASS_Mu30),
        passing_HLT_Mu40                     = cms.string(PASS_Mu40),
        passing_HLT_Mu15_eta2p1              = cms.string(PASS_Mu15_eta2p1),
        passing_HLT_Mu24_eta2p1              = cms.string(PASS_Mu24_eta2p1),
        passing_HLT_Mu30_eta2p1              = cms.string(PASS_Mu30_eta2p1),
        passing_HLT_Mu40_eta2p1              = cms.string(PASS_Mu40_eta2p1),
        passing_HLT_Mu50_eta2p1              = cms.string(PASS_Mu50_eta2p1),
        passing_HLT_IsoMu24                  = cms.string(PASS_IsoMu24),
        passing_HLT_IsoMu30                  = cms.string(PASS_IsoMu30),
        passing_HLT_IsoMu20_eta2p1           = cms.string(PASS_IsoMu20_eta2p1),
        passing_HLT_IsoMu24_eta2p1           = cms.string(PASS_IsoMu24_eta2p1),
        passing_HLT_IsoMu30_eta2p1           = cms.string(PASS_IsoMu30_eta2p1),
        passing_HLT_IsoMu34_eta2p1           = cms.string(PASS_IsoMu34_eta2p1),
        passing_HLT_IsoMu40_eta2p1           = cms.string(PASS_IsoMu40_eta2p1),
        passing_HLT_Mu40_PFHT350             = cms.string(PASS_Mu40_PFHT350),
        passing_HLT_Mu60_PFHT350             = cms.string(PASS_Mu60_PFHT350),
        passing_HLT_Mu40_PFNoPUHT350         = cms.string(PASS_Mu40_PFNoPUHT350),
        passing_HLT_Mu60_PFNoPUHT350         = cms.string(PASS_Mu60_PFNoPUHT350),
        passing_HLT_PFHT350_Mu15_PFMET45     = cms.string(PASS_PFHT350_Mu15_PFMET45),
        passing_HLT_PFHT350_Mu15_PFMET50     = cms.string(PASS_PFHT350_Mu15_PFMET50),
        passing_HLT_PFHT400_Mu5_PFMET45      = cms.string(PASS_PFHT400_Mu5_PFMET45),
        passing_HLT_PFHT400_Mu5_PFMET50      = cms.string(PASS_PFHT400_Mu5_PFMET50),
        passing_HLT_PFNoPUHT350_Mu15_PFMET45 = cms.string(PASS_PFNoPUHT350_Mu15_PFMET45),
        passing_HLT_PFNoPUHT350_Mu15_PFMET50 = cms.string(PASS_PFNoPUHT350_Mu15_PFMET50),
        passing_HLT_PFNoPUHT400_Mu5_PFMET45  = cms.string(PASS_PFNoPUHT400_Mu5_PFMET45),
        passing_HLT_PFNoPUHT400_Mu5_PFMET50  = cms.string(PASS_PFNoPUHT400_Mu5_PFMET50),
    ) 
    
    ################################################################################################################################
    #   ___  ____ ____ ___  ____ ____ 
    #   |__] |__/ |  | |__] |___ [__  
    #   |    |  \ |__| |__] |___ ___] 
                         
    # probe1: standalone muons
    process.staTracks = cms.EDProducer("TrackViewCandidateProducer", 
        src  = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
        particleType = cms.string("mu+"),
        cut = cms.string(""),
    )
    process.staProbes = cms.EDFilter("RecoChargedCandidateRefSelector",
        src = cms.InputTag("staTracks"),
        cut = cms.string(STA_PROBE_CUTS),
    )
    
    # probe2: tracker muons
    process.tkTracks = cms.EDProducer("TrackViewCandidateProducer",
        src = cms.InputTag("generalTracks"),
        particleType = cms.string('mu+'),
        cut = cms.string(TRK_PROBE_CUTS),
    )
    process.tkProbes = cms.EDFilter("RecoChargedCandidateRefSelector",
        src = cms.InputTag("tkTracks"),
        cut = cms.string("")
    )
    
    # probe3: global muons
    process.glbProbes = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("cleanPatMuonsTriggerMatch"),
        cut = cms.string(GLB_PROBE_CUTS),
    )
    
    ################################################################################################################################
    #   ___ ____ ____ ____ 
    #    |  |__| | __ [__  
    #    |  |  | |__] ___] 
    
    process.tagMuonsSingleEle = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("cleanPatMuonsTriggerMatch"),
        cut = cms.string(TAG_CUTS + " && "+ PASS_SINGLEELE), 
    )
    
    process.tagMuonsAll = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("cleanPatMuonsTriggerMatch"),
        cut = cms.string(TAG_CUTS + " && " + PASS_ANY),
    )
    
    process.allTagsAndProbes = cms.Sequence(
        process.tagMuonsSingleEle +
        process.tagMuonsAll +
        process.staTracks * process.staProbes +
        process.tkTracks * process.tkProbes +
        process.glbProbes
    )
    
    ################################################################################################################################
    #   ___ ____ ____    _  _    ___  ____ ____ ___  ____    ___  ____ _ ____ ____ 
    #    |  |__| | __    |\ |    |__] |__/ |  | |__] |___    |__] |__| | |__/ [__  
    #    |  |  | |__]    | \|    |    |  \ |__| |__] |___    |    |  | | |  \ ___] 
    
    # T&P Pair1: globalmuons & standalone muons 
    process.tpGlbSta = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("tagMuonsSingleEle@+ staProbes@-"), # charge conjugate states are implied
        cut   = cms.string("40 < mass < 140"),
    )
    
    # T&P Pair2: global muons & tracker muons
    process.tpGlbTk = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("tagMuonsSingleEle@+ tkProbes@-"), # charge conjugate states are implied
        cut   = cms.string("40 < mass < 140"),
    )
    
    # T&P Pair3: global muons & global muons
    process.tpGlbGlb = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("tagMuonsAll@+ glbProbes@-"), # charge conjugate states are implied
        cut   = cms.string("40 < mass < 140"),
    )

    process.allTPPairs = cms.Sequence(
        process.tpGlbSta +
        process.tpGlbTk +
        process.tpGlbGlb
    )
    
    ################################################################################################################################
    #   _  _ ____ ___ ____ _  _    ____ _  _ ___     ___  ____ ____ ____ 
    #   |\/| |__|  |  |    |__|    |__| |\ | |  \    |__] |__| [__  [__  
    #   |  | |  |  |  |___ |  |    |  | | \| |__/    |    |  | ___] ___] 
    
    # passing1: standalone muons passing as tracker muons
    process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
        src     = cms.InputTag("staTracks"), # all standalone muons
        matched = cms.InputTag("tkTracks"),  # to all tk tracks
        algorithm = cms.string("byDirectComparison"), # using parameters at PCA
        srcTrack = cms.string("tracker"),  # 'staTracks' is a 'RecoChargedCandidate', so it thinks
        srcState = cms.string("atVertex"), # it has a 'tracker' track, not a standalone one
        matchedTrack = cms.string("tracker"),
        matchedState = cms.string("atVertex"),
        maxDeltaR        = cms.double(1.),   # large range in DR
        maxDeltaEta      = cms.double(0.2),  # small in eta, which is more precise
        maxDeltaLocalPos = cms.double(100),
        maxDeltaPtRel    = cms.double(3),
        sortBy           = cms.string("deltaR"),
    )
    process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
        src   = cms.InputTag("staProbes"),
        match = cms.InputTag("staToTkMatch"),
    )
    
    # passing2: tracker muons passing as global muons
    process.tkToGlbMatch = cms.EDProducer("MatcherUsingTracksMatchInfo",
        src     = cms.InputTag("tkTracks"),
        matched = cms.InputTag("glbProbes"),
        algorithm = cms.string("byDirectComparison"), 
        srcTrack = cms.string("tracker"),
        srcState = cms.string("atVertex"),
        matchedTrack = cms.string("tracker"),
        matchedState = cms.string("atVertex"),
        maxDeltaR        = cms.double(0.01),
        maxDeltaLocalPos = cms.double(0.01),
        maxDeltaPtRel    = cms.double(0.01),
        sortBy           = cms.string("deltaR"),
    )
    process.tkPassingGlb = cms.EDProducer("MatchedCandidateSelector",
        src   = cms.InputTag("tkProbes"),
        match = cms.InputTag("tkToGlbMatch"),
    )
    
    process.allPassingProbes = cms.Sequence(
        process.tkToGlbMatch * process.tkPassingGlb +
        process.staToTkMatch * process.staPassingTk 
    )
    
    # passing 3: global muons with ID passing HLT paths
    # matches and passes can be found in the top CUTS section
    
    # Additional variable to check if event passed the specific cross trigger
    process.glbEventPassingHLTMu40PFHT350 = cms.EDProducer("HLTResultProducer",
        probes = cms.InputTag("cleanPatMuonsTriggerMatch"),
        TriggerResultsTag = cms.InputTag("TriggerResults","",hltName),
        HLTPaths = cms.vstring("HLT_Mu40_PFHT350_v.*"),
        andOr = cms.bool(True)
    )
    process.glbEventPassingHLTMu60PFHT350            = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_Mu60_PFHT350_v.*"))
    process.glbEventPassingHLTMu40PFNoPUHT350        = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_Mu40_PFNoPUHT350_v.*"))
    process.glbEventPassingHLTMu60PFNoPUHT350        = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_Mu60_PFNoPUHT350_v.*"))
    process.glbEventPassingHLTPFHT350Mu15PFMET45     = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFHT350_Mu15_PFMET45_v.*"))
    process.glbEventPassingHLTPFHT350Mu15PFMET50     = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFHT350_Mu15_PFMET50_v.*"))
    process.glbEventPassingHLTPFHT400Mu5PFMET45      = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFHT400_Mu5_PFMET45_v.*"))
    process.glbEventPassingHLTPFHT400Mu5PFMET50      = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFHT400_Mu5_PFMET50_v.*"))
    process.glbEventPassingHLTPFNoPUHT350Mu15PFMET45 = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFNoPUHT350_Mu15_PFMET45_v.*"))
    process.glbEventPassingHLTPFNoPUHT350Mu15PFMET50 = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFNoPUHT350_Mu15_PFMET50_v.*"))
    process.glbEventPassingHLTPFNoPUHT400Mu5PFMET45  = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFNoPUHT400_Mu5_PFMET45_v.*"))
    process.glbEventPassingHLTPFNoPUHT400Mu5PFMET50  = process.glbEventPassingHLTMu40PFHT350.clone(HLTPaths = cms.vstring("HLT_PFNoPUHT400_Mu5_PFMET50_v.*"))
    
    process.allHLTResults = cms.Sequence(
        process.glbEventPassingHLTMu40PFHT350 +
        process.glbEventPassingHLTMu60PFHT350 +
        process.glbEventPassingHLTMu40PFNoPUHT350 +
        process.glbEventPassingHLTMu60PFNoPUHT350 +
        process.glbEventPassingHLTPFHT350Mu15PFMET45 +
        process.glbEventPassingHLTPFHT350Mu15PFMET50 +
        process.glbEventPassingHLTPFHT400Mu5PFMET45 +
        process.glbEventPassingHLTPFHT400Mu5PFMET50 +
        process.glbEventPassingHLTPFNoPUHT350Mu15PFMET45 +
        process.glbEventPassingHLTPFNoPUHT350Mu15PFMET50 +
        process.glbEventPassingHLTPFNoPUHT400Mu5PFMET45 +
        process.glbEventPassingHLTPFNoPUHT400Mu5PFMET50
    )
    
    ##    __  __  ____   __  __       _       _               
    ##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
    ##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
    ##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
    ##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
    
    if mcInfo:
        # MC muon matches
        process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
            pdgId = cms.vint32(13),
            src = cms.InputTag("cleanPatMuonsTriggerMatch"),
            distMin = cms.double(0.3),
            matched = cms.InputTag("genParticles")
        )
        # MC tracker muon matches & standalone muon matches
        process.tkMcMatch  = process.muMcMatch.clone(src = "tkTracks")
        process.staMcMatch = process.muMcMatch.clone(src = "staTracks",distMin = 0.6)

        process.allMcMatches = cms.Sequence(
            process.staMcMatch +
            process.tkMcMatch +
            process.muMcMatch 
        )

    ################################################################################################################################
    ##    _____      _                        _  __     __             
    ##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
    ##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
    ##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
    ##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
    ##   
    ################################################################################################################################

    ################################################################################################################################
    #   Vertices
    
    process.staNvertices = cms.EDProducer("VertexMultiplicityCounter", 
        probes = cms.InputTag("staTracks"),
        objects = cms.InputTag("offlinePrimaryVertices"),
        objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    )
    process.trkNvertices = process.staNvertices.clone(probes = "tkTracks")
    process.glbNvertices = process.staNvertices.clone(probes = "cleanPatMuonsTriggerMatch")
    
    process.allVertices = cms.Sequence(
        process.staNvertices +
        process.trkNvertices +
        process.glbNvertices
    )
    
    ################################################################################################################################
    #   PU Reweighting
    
    process.staPileupWeight = cms.EDProducer("PileupWeightComputer",
        probes = cms.InputTag("staTracks"),
        isMC = cms.bool(mcInfo),
        dataPileupFile = cms.string("PileupHistogram_2012Data_FlatPU_50bins.root"),
        mcPileupFile   = cms.string("PileupHistogram_2012Data_FlatPU_50bins.root"),
        dataPileupHistoName = cms.string("pileup"),
        mcPileupHistoName = cms.string("mcpileup"),
        #mcPileupHistoName = cms.string("pileup"),
        mcLumiScale = cms.double(221.95),
        dataPileupInputFile = cms.string("run_ls_instlumi_pileup_2012.txt"),
     )
    process.trkPileupWeight = process.staPileupWeight.clone(probes = "tkTracks")
    process.glbPileupWeight = process.staPileupWeight.clone(probes = "cleanPatMuonsTriggerMatch")
    
    process.allPileupWeights = cms.Sequence(
        process.staPileupWeight +
        process.trkPileupWeight +
        process.glbPileupWeight
    )
    
    ################################################################################################################################
    #   _ _  _ ___  ____ ____ ___    ___  ____ ____ ____ _  _ ____ ___ ____ ____ 
    #   | |\/| |__] |__| |     |     |__] |__| |__/ |__| |\/| |___  |  |___ |__/ 
    #   | |  | |    |  | |___  |     |    |  | |  \ |  | |  | |___  |  |___ |  \ 
    
    process.staImpactParameter = cms.EDProducer("ChargedCandidateImpactParameter",
        probes = cms.InputTag("staTracks"),
    )
    process.trkImpactParameter = process.staImpactParameter.clone(probes = "tkTracks")
    process.glbImpactParameter = cms.EDProducer("PatMuonImpactParameter",
        probes = cms.InputTag("cleanPatMuonsTriggerMatch"),
    )
    
    process.allImpactParameters = cms.Sequence(
        process.staImpactParameter +
        process.trkImpactParameter +
        process.glbImpactParameter
    )
    
    ################################################################################################################################
    #    _ ____ ___ ____ 
    #    | |___  |  [__  
    #   _| |___  |  ___] 
    
    process.selectedJets = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("cleanPatJetsAK5PF"),
        cut = cms.string( JET_CUTS ),
    )

    # DeltaR
    process.stadRToNearestJet = cms.EDProducer("minCutDeltaRNearestPatJetComputer",
        probes = cms.InputTag("staTracks"),
        objects = cms.InputTag("selectedJets"),
        minDeltaR = cms.double(0.1),
        objectSelection = cms.InputTag(""),
    )
    process.trkdRToNearestJet = process.stadRToNearestJet.clone(probes = "tkTracks")
    process.glbdRToNearestJet = process.stadRToNearestJet.clone(probes = "cleanPatMuonsTriggerMatch")

    # Njet
    process.staJetMultiplicity = cms.EDProducer("PatJetMultiplicityCounter",
        probes = cms.InputTag("staTracks"),
        objects = cms.InputTag("selectedJets"),
        objectSelection = cms.InputTag(""),
    )
    process.trkJetMultiplicity = process.staJetMultiplicity.clone(probes = "tkTracks")
    process.glbJetMultiplicity = process.staJetMultiplicity.clone(probes = "cleanPatMuonsTriggerMatch")

    # HT
    process.staHT = cms.EDProducer("PatJetHTComputer",
        probes = cms.InputTag("staTracks"),
        objects = cms.InputTag("selectedJets"),
        objectSelection = cms.InputTag(""),
    )
    process.trkHT = process.staHT.clone(probes = "tkTracks")
    process.glbHT = process.staHT.clone(probes = "cleanPatMuonsTriggerMatch")
    
    #   MET
    process.staMet = cms.EDProducer("PatMetAssociator",
        probes = cms.InputTag("staTracks"),
        metTag = cms.InputTag("patMETsPF"),
    )
    process.trkMet = process.staMet.clone(probes = "tkTracks")
    process.glbMet = process.staMet.clone(probes = "cleanPatMuonsTriggerMatch")
    
    #   ST
    process.staST = cms.EDProducer("PatMetSTComputer",
        probes = cms.InputTag("staTracks"),
        metTag = cms.InputTag("patMETsPF"),
    )
    process.trkST = process.staST.clone(probes = "tkTracks")
    process.glbST = process.staST.clone(probes = "cleanPatMuonsTriggerMatch")
    
    process.allJets = cms.Sequence(
        process.selectedJets*
        ( process.stadRToNearestJet +
          process.trkdRToNearestJet +
          process.glbdRToNearestJet +
          process.staJetMultiplicity +
          process.trkJetMultiplicity +
          process.glbJetMultiplicity +
          process.staHT +
          process.trkHT +
          process.glbHT)
    )
    
    process.allMet = cms.Sequence(
        process.staMet +
        process.trkMet +
        process.glbMet
    )
    
    process.allST = cms.Sequence(
        process.staST +
        process.trkST +
        process.glbST
    )
    
    ################################################################################################################################
    #   _ ____ ____ _    ____ ___ _ ____ _  _ ____ 
    #   | [__  |  | |    |__|  |  | |  | |\ | [__  
    #   | ___] |__| |___ |  |  |  | |__| | \| ___] 
    
    # Now Using Partcle Flow Relative Isolations
    
    # from RecoMuon.MuonIsolationProducers.trackExtractorBlocks_cff import MIsoTrackExtractorBlock
    # process.trkIsoDepositTk = cms.EDProducer("CandIsoDepositProducer",
    #     src = cms.InputTag("tkTracks"),
    #     MultipleDepositsFlag = cms.bool(False),
    #     trackType = cms.string('best'),
    #     ExtractorPSet = cms.PSet(
    #         MIsoTrackExtractorBlock
    #     )
    # )
    # from RecoMuon.MuonIsolationProducers.caloExtractorByAssociatorBlocks_cff import MIsoCaloExtractorByAssociatorTowersBlock
    # process.trkIsoDepositCalByAssociatorTowers = cms.EDProducer("CandIsoDepositProducer",
    #     src = cms.InputTag("tkTracks"),
    #     MultipleDepositsFlag = cms.bool(True),
    #     trackType = cms.string('best'),
    #     ExtractorPSet = cms.PSet(
    #         MIsoCaloExtractorByAssociatorTowersBlock
    #     )
    # )
    # 
    # process.TrackIsolationForTrk = cms.EDProducer("IsolationProducerForTracks",
    #     highPtTracks = cms.InputTag("tkTracks"),
    #     tracks = cms.InputTag("tkTracks"),
    #     isoDeps = cms.InputTag("trkIsoDepositTk"),
    #     coneSize = cms.double(0.3),
    #     trackPtMin = cms.double(3.0)
    # )
    # 
    # process.TrackIsolationForTrk04 = process.TrackIsolationForTrk.clone(coneSize = cms.double(0.4))
    # 
    # process.EcalIsolationForTrk = cms.EDProducer("IsolationProducerForTracks",
    #     highPtTracks = cms.InputTag("tkTracks"),
    #     tracks = cms.InputTag("tkTracks"),
    #     isoDeps = cms.InputTag("trkIsoDepositCalByAssociatorTowers","ecal"),
    #     coneSize = cms.double(0.3),
    #     trackPtMin = cms.double(3.0)
    # )
    # 
    # process.EcalIsolationForTrk04                 = process.EcalIsolationForTrk.clone(coneSize = cms.double(0.4))
    # 
    # process.HcalIsolationForTrk = cms.EDProducer("IsolationProducerForTracks",
    #     highPtTracks = cms.InputTag("tkTracks"),
    #     tracks = cms.InputTag("tkTracks"),
    #     isoDeps = cms.InputTag("trkIsoDepositCalByAssociatorTowers","hcal"),
    #     coneSize = cms.double(0.3),
    #     trackPtMin = cms.double(3.0)
    # )
    # 
    # process.HcalIsolationForTrk04               = process.HcalIsolationForTrk.clone(coneSize = cms.double(0.4) )
    # 
    # process.RelIsolationForTrk = cms.EDProducer("RelIsolationProducerForTracks",
    #     highPtTracks = cms.InputTag("tkTracks"),
    #     tracks = cms.InputTag("tkTracks"),
    #     trkisoDeps = cms.InputTag("trkIsoDepositTk"),
    #     ecalisoDeps = cms.InputTag("trkIsoDepositCalByAssociatorTowers","ecal"),
    #     hcalisoDeps = cms.InputTag("trkIsoDepositCalByAssociatorTowers","hcal"),
    #     coneSize = cms.double(0.3),
    #     trackPtMin = cms.double(3.0)
    # )
    # 
    # process.RelIsolationForTrk04                = process.RelIsolationForTrk.clone(coneSize=cms.double(0.4))
    # 
    # process.staIsoDepositTk                     = process.trkIsoDepositTk.clone( src = cms.InputTag("staTracks") )
    # process.staIsoDepositCalByAssociatorTowers  = process.trkIsoDepositCalByAssociatorTowers.clone( src = cms.InputTag("staTracks") )
    # process.TrackIsolationForStA                = process.TrackIsolationForTrk.clone( highPtTacks = cms.InputTag("staTracks") )
    # process.EcalIsolationForStA                 = process.EcalIsolationForTrk.clone( highPtTacks = cms.InputTag("staTracks") )
    # process.HcalIsolationForStA                 = process.HcalIsolationForTrk.clone( highPtTacks = cms.InputTag("staTracks") )
    # process.RelIsolationForStA                  = process.RelIsolationForTrk.clone( highPtTacks = cms.InputTag("staTracks") )
    # 
    # process.TrackIsolationForStA.tracks = cms.InputTag("staTracks")
    # process.EcalIsolationForStA.tracks = cms.InputTag("staTracks")
    # process.HcalIsolationForStA.tracks = cms.InputTag("staTracks")
    # process.RelIsolationForStA.tracks = cms.InputTag("staTracks")
    # 
    # process.TrackIsolationForStA.isoDeps = cms.InputTag("staIsoDepositTk")
    # process.EcalIsolationForStA.isoDeps = cms.InputTag("staIsoDepositCalByAssociatorTowers","ecal")
    # process.HcalIsolationForStA.isoDeps = cms.InputTag("staIsoDepositCalByAssociatorTowers","hcal")
    # process.RelIsolationForStA.trkisoDeps = cms.InputTag("staIsoDepositTk")
    # 
    # process.RelIsolationForStA.ecalisoDeps = cms.InputTag("staIsoDepositCalByAssociatorTowers","ecal")
    # process.RelIsolationForStA.hcalisoDeps = cms.InputTag("staIsoDepositCalByAssociatorTowers","hcal")
    # 
    # process.TrackIsolationForStA04              = process.TrackIsolationForStA.clone( coneSize=cms.double(0.4) )
    # process.EcalIsolationForStA04               = process.EcalIsolationForStA.clone( coneSize = cms.double(0.4) )
    # process.HcalIsolationForStA04               = process.HcalIsolationForStA.clone( coneSize = cms.double(0.4) )
    # process.RelIsolationForStA04                = process.RelIsolationForStA.clone( coneSize = cms.double(0.4) )

    #process.allIsolations = cms.Sequence(
    #    (process.trkIsoDepositTk +
    #     process.trkIsoDepositCalByAssociatorTowers) *
    #    (process.TrackIsolationForTrk +
    #     process.EcalIsolationForTrk +
    #     process.HcalIsolationForTrk +
    #     process.RelIsolationForTrk) *
    #    (process.staIsoDepositTk +
    #     process.staIsoDepositCalByAssociatorTowers) *
    #    (process.TrackIsolationForStA +
    #     process.EcalIsolationForStA +
    #     process.HcalIsolationForStA +
    #     process.RelIsolationForStA) *
    #    (process.TrackIsolationForTrk04 +
    #     process.EcalIsolationForTrk04 +
    #     process.HcalIsolationForTrk04 +
    #     process.RelIsolationForTrk04) *
    #    (process.TrackIsolationForStA04 +
    #     process.EcalIsolationForStA04 +
    #     process.HcalIsolationForStA04 +
    #     process.RelIsolationForStA04) 
    #) 
    
    ################################################################################################################################
    #   Delta Reco-PF Muon Pt
    
    process.glbDeltaPfRecoPt = cms.EDProducer("MuonDeltaPfRecoPt",
        probes = cms.InputTag("cleanPatMuonsTriggerMatch"),
    )
    
    ################################################################################################################################
    #   ___  ____ ____ ____ _  _ ____ ___ ____ ____ ____ 
    #   |__] |__| |__/ |__| |\/| |___  |  |___ |__/ [__  
    #   |    |  | |  \ |  | |  | |___  |  |___ |  \ ___] 

    if mcInfo:
        mcTruthCommonStuff = cms.PSet(
            isMC = cms.bool(True),
            addRunLumiInfo = cms.bool(False),
            makeMCUnbiasTree = cms.bool(False),
            checkMotherInUnbiasEff = cms.bool(True),
            motherPdgId = cms.vint32(22,23), # Gamma or Z
            tagMatches = cms.InputTag("muMcMatch"),
            #mcVariables = cms.PSet(
            #    mass  = cms.string("mass"),
            #    mt  = cms.string("mt"),
            #    pt  = cms.string("pt"),
            #    et  = cms.string("et"),
            #    phi  = cms.string("phi"),
            #    eta = cms.string("eta"),
            #    e  = cms.string("energy"),
            #    p  = cms.string("p"),
            #    px  = cms.string("px"),
            #    py  = cms.string("py"),
            #    pz  = cms.string("pz"),
            #    theta  = cms.string("theta"),
            #    vx     = cms.string("vx"),
            #    vy     = cms.string("vy"),
            #    vz     = cms.string("vz"),
            #    charge = cms.string("charge"),
            #    rapidity  = cms.string("rapidity"),
            #),
            #mcFlags     =  cms.PSet(
            #    flag = cms.string("pt>0")
            #),
        )
    else:
        mcTruthCommonStuff = cms.PSet(
            isMC = cms.bool(False),
            addRunLumiInfo = cms.bool(True),
        )
    
    commonStuff = cms.PSet(
        mcTruthCommonStuff,
        addEventVariablesInfo = cms.bool(False),
        ignoreExceptions = cms.bool(False),
        arbitration = cms.string("OneProbe"),
        #pairFlags = cms.PSet(
        #    mass60to120 = cms.string("60 < mass < 120")
        #),
        #pairVariables =  cms.PSet(
        #    #mass  = cms.string("mass"), # created automatically
        #    mt  = cms.string("mt"), 
        #    pt  = cms.string("pt"),
        #    et  = cms.string("et"),
        #    phi  = cms.string("phi"),
        #    eta = cms.string("eta"),
        #    abs_eta = cms.string("abs(eta)"),
        #    #e  = cms.string("energy"),
        #    #p  = cms.string("p"),
        #    #px  = cms.string("px"),
        #    #py  = cms.string("py"),
        #    #pz  = cms.string("pz"),
        #    #theta  = cms.string("theta"),    
        #    #vx     = cms.string("vx"),
        #    #vy     = cms.string("vy"),
        #    #vz     = cms.string("vz"),
        #    #rapidity  = cms.string("rapidity"),
        #),
    )
    
    staParameters = cms.PSet(
        variables = cms.PSet(
            pt       = cms.string("pt"),
            phi      = cms.string("phi"),
            eta      = cms.string("eta"),
            abs_eta  = cms.string("abs(eta)"),
            nvtx     = cms.InputTag("staNvertices"),
            pileup   = cms.InputTag("staPileupWeight","pileup"),
            instlumi = cms.InputTag("staPileupWeight","instlumi"),
            weight   = cms.InputTag("staPileupWeight","weight"),
            d0_v     = cms.InputTag("staImpactParameter","d0v"),
            d0_b     = cms.InputTag("staImpactParameter","d0b"),
            dz_v     = cms.InputTag("staImpactParameter","dzv"),
            dz_b     = cms.InputTag("staImpactParameter","dzb"),
            drjet    = cms.InputTag("stadRToNearestJet"),
            njet     = cms.InputTag("staJetMultiplicity",),
            ht       = cms.InputTag("staHT"),
            met      = cms.InputTag("staMet"),
            st       = cms.InputTag("staST"),
            #d0 = cms.string("track.d0"),
            #dz = cms.string("track.dz"),
            #validhits = cms.string("track.numberOfValidHits"),
            #trkiso  = cms.InputTag("TrackIsolationForStA"),
            #ecaliso = cms.InputTag("EcalIsolationForStA"),
            #hcaliso = cms.InputTag("HcalIsolationForStA"),
            #reliso = cms.InputTag("RelIsolationForStA"),
            #trkiso04  = cms.InputTag("TrackIsolationForStA04"),
            #ecaliso04 = cms.InputTag("EcalIsolationForStA04"),
            #hcaliso04 = cms.InputTag("HcalIsolationForStA04"),
            #reliso04 = cms.InputTag("RelIsolationForStA04"),
        ),
        tagVariables = cms.PSet(),
        flags = cms.PSet(
            passing_trk = cms.InputTag("staPassingTk")
        ),
        tagFlags = cms.PSet(),
    )
    
    trkParameters = cms.PSet(
        variables = cms.PSet(
            pt         = cms.string("pt"),
            phi        = cms.string("phi"),
            eta        = cms.string("eta"),
            abs_eta    = cms.string("abs(eta)"),
            nvtx       = cms.InputTag("trkNvertices"),
            pileup     = cms.InputTag("trkPileupWeight","pileup"),
            instlumi   = cms.InputTag("trkPileupWeight","instlumi"),
            weight     = cms.InputTag("trkPileupWeight","weight"),
            d0_v       = cms.InputTag("trkImpactParameter","d0v"),
            d0_b       = cms.InputTag("trkImpactParameter","d0b"),
            dz_v       = cms.InputTag("trkImpactParameter","dzv"),
            dz_b       = cms.InputTag("trkImpactParameter","dzb"),
            drjet      = cms.InputTag("trkdRToNearestJet"),
            njet       = cms.InputTag("trkJetMultiplicity"),
            ht         = cms.InputTag("trkHT"),
            met        = cms.InputTag("trkMet"),
            st         = cms.InputTag("trkST"),
            absdeltapt = cms.InputTag("tkToGlbMatch", "absdeltaPtRecoPF"), # -9999: under 10 GeV, 999999: no global match
            deltapt    = cms.InputTag("tkToGlbMatch", "deltaPtRecoPF"), # -9999: under 10 GeV, 999999: no global match
            reliso     = cms.InputTag("tkToGlbMatch", "pfreliso"),
            #d0 = cms.string("track.d0"),
            #dz = cms.string("track.dz"),
            #ptErrorByPt2 = cms.InputTag("tkToGlbMatch", "ptErrByPt2"),
            #pixlayer = cms.string("track.hitPattern.pixelLayersWithMeasurement"),
            #trkiso  = cms.InputTag("TrackIsolationForTrk"),
            #ecaliso = cms.InputTag("EcalIsolationForTrk"),
            #hcaliso = cms.InputTag("HcalIsolationForTrk"),
            #reliso = cms.InputTag("RelIsolationForTrk"),
            #trkiso04  = cms.InputTag("TrackIsolationForTrk04"),
            #ecaliso04 = cms.InputTag("EcalIsolationForTrk04"),
            #hcaliso04 = cms.InputTag("HcalIsolationForTrk04"),
            #reliso04 = cms.InputTag("RelIsolationForTrk04"),
        ),
        tagVariables = cms.PSet(),
        flags = cms.PSet(
            passing_glb = cms.InputTag("tkPassingGlb")
        ),
        tagFlags = cms.PSet(),
    )
    
    glbParameters = cms.PSet(
        variables = cms.PSet(
            pt         = cms.string("pt"),
            phi        = cms.string("phi"),
            eta        = cms.string("eta"),
            abs_eta    = cms.string("abs(eta)"),
            nvtx       = cms.InputTag("glbNvertices"),
            pileup     = cms.InputTag("glbPileupWeight","pileup"),
            instlumi   = cms.InputTag("glbPileupWeight","instlumi"),
            weight     = cms.InputTag("glbPileupWeight","weight"),
            d0_v       = cms.InputTag("glbImpactParameter","d0v"),
            d0_b       = cms.InputTag("glbImpactParameter","d0b"),
            dz_v       = cms.InputTag("glbImpactParameter","dzv"),
            dz_b       = cms.InputTag("glbImpactParameter","dzb"),
            drjet      = cms.InputTag("glbdRToNearestJet"),
            njet       = cms.InputTag("glbJetMultiplicity"),
            ht         = cms.InputTag("glbHT"),
            met        = cms.InputTag("glbMet"),
            st         = cms.InputTag("glbST"),
            absdeltapt = cms.InputTag("glbDeltaPfRecoPt","absdeltapt"), # -9999: under 10 GeV
            reliso     = cms.string("(pfIsolationR03.sumChargedHadronPt + max(pfIsolationR03.sumNeutralHadronEt + pfIsolationR03.sumPhotonEt - pfIsolationR03.sumPUPt/2,0.0))/pt"),
            #d0 = cms.string("track.d0"),
            #dz = cms.string("track.dz"),
            #pterror = cms.string("globalTrack.ptError"),  
            #validhits = cms.string("globalTrack.hitPattern.numberOfValidTrackerHits"),
            #validpixhits = cms.string("innerTrack.hitPattern.numberOfValidPixelHits"),
            #trklayer = cms.string("track.hitPattern.trackerLayersWithMeasurement"),
            #pixlayer = cms.string("track.hitPattern.pixelLayersWithMeasurement"),
            #nmatched = cms.string("numberOfMatchedStations"),
            #trkiso  = cms.string("trackIso"),
            #ecaliso = cms.string("ecalIso"),
            #hcaliso = cms.string("hcalIso"),
            #reliso = cms.string("( trackIso + ecalIso + hcalIso )/pt"),
            event_passing_HLT_Mu40_PFHT350             = cms.InputTag("glbEventPassingHLTMu40PFHT350"),
            event_passing_HLT_Mu60_PFHT350             = cms.InputTag("glbEventPassingHLTMu60PFHT350"),
            event_passing_HLT_Mu40_PFNoPUHT350         = cms.InputTag("glbEventPassingHLTMu40PFNoPUHT350"),
            event_passing_HLT_Mu60_PFNoPUHT350         = cms.InputTag("glbEventPassingHLTMu60PFNoPUHT350"),
            event_passing_HLT_PFHT350_Mu15_PFMET45     = cms.InputTag("glbEventPassingHLTPFHT350Mu15PFMET45"),
            event_passing_HLT_PFHT350_Mu15_PFMET50     = cms.InputTag("glbEventPassingHLTPFHT350Mu15PFMET50"),
            event_passing_HLT_PFHT400_Mu5_PFMET45      = cms.InputTag("glbEventPassingHLTPFHT400Mu5PFMET45"),
            event_passing_HLT_PFHT400_Mu5_PFMET50      = cms.InputTag("glbEventPassingHLTPFHT400Mu5PFMET50"),
            event_passing_HLT_PFNoPUHT350_Mu15_PFMET45 = cms.InputTag("glbEventPassingHLTPFNoPUHT350Mu15PFMET45"),
            event_passing_HLT_PFNoPUHT350_Mu15_PFMET50 = cms.InputTag("glbEventPassingHLTPFNoPUHT350Mu15PFMET50"),
            event_passing_HLT_PFNoPUHT400_Mu5_PFMET45  = cms.InputTag("glbEventPassingHLTPFNoPUHT400Mu5PFMET45"),
            event_passing_HLT_PFNoPUHT400_Mu5_PFMET50  = cms.InputTag("glbEventPassingHLTPFNoPUHT400Mu5PFMET50"),
        ),
        tagVariables = cms.PSet(),
        flags    = cms.PSet( AllTriggerFlags ),
        tagFlags = cms.PSet( AllTriggerFlags ),
    )
    
    ################################################################################################################################
    #   ____ _ ___    ___ ____ ____ ____    ___  ____ ____ ___  _  _ ____ ____ ____ 
    #   |___ |  |      |  |__/ |___ |___    |__] |__/ |  | |  \ |  | |    |___ |__/ 
    #   |    |  |      |  |  \ |___ |___    |    |  \ |__| |__/ |__| |___ |___ |  \ 
    
    # tracker muon efficiency from standalone muons
    process.fitTkFromSta = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, staParameters,
        tagProbePairs = cms.InputTag("tpGlbSta"),
        probeMatches  = cms.InputTag("staMcMatch"),
        allProbes     = cms.InputTag("staProbes"),
    )
    
    # global muon efficieny from tracker muons
    process.fitGlbFromTk = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, trkParameters,
        tagProbePairs = cms.InputTag("tpGlbTk"),
        probeMatches  = cms.InputTag("tkMcMatch"),
        allProbes     = cms.InputTag("tkProbes"),
    )
    
    # HLT efficiency from global muons
    process.fitHltFromGlb = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, glbParameters,
        tagProbePairs = cms.InputTag("tpGlbGlb"),
        probeMatches  = cms.InputTag("muMcMatch"),
        allProbes     = cms.InputTag("glbProbes"),
    )
    
    process.allTPTrees = cms.Sequence(
        process.fitGlbFromTk +
        process.fitTkFromSta +
        process.fitHltFromGlb 
    )
    
    ##############################################################################################################
    ##    ____       _   _     
    ##   |  _ \ __ _| |_| |__  
    ##   | |_) / _` | __| '_ \ 
    ##   |  __/ (_| | |_| | | |
    ##   |_|   \__,_|\__|_| |_|
    ##
    
    if mcInfo:
        process.TagAndProbe = cms.Sequence(
            process.allTagsAndProbes *
            ( process.allTPPairs +
              process.allPassingProbes +
              process.allHLTResults +
              process.allMcMatches +
              process.allVertices +
              process.allPileupWeights +
              process.allImpactParameters +
              process.allJets + 
              process.allMet +
              process.allST +
              #process.allIsolations +
              process.glbDeltaPfRecoPt )
            * process.allTPTrees
        )
    else:
        process.TagAndProbe = cms.Sequence(
            process.allTagsAndProbes *
            ( process.allTPPairs +
              process.allPassingProbes +
              process.allHLTResults +
              process.allVertices +
              process.allPileupWeights +
              process.allImpactParameters +
              process.allJets + 
              process.allMet +
              process.allST +
              #process.allIsolations +
              process.glbDeltaPfRecoPt )
            * process.allTPTrees
        )
