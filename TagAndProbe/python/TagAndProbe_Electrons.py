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
    
    ################################################################################################################################
    #   ____ _  _ ___ ____ 
    #   |    |  |  |  [__  
    #   |___ |__|  |  ___] 
    
    TAG_CUTS = "p4().Pt() > 10 && abs(eta) < 2.5 && electronID('simpleEleId80cIso') == 7"
    
    SC_PROBE_CUTS = "et>5 && abs(eta)<2.5"
    GSF_PROBE_CUTS = "(ecalEnergy*sin(superClusterPosition.theta))>5 && abs(superCluster.eta) <= 2.5 && ecalDrivenSeed==1"
    PAT_PROBE_CUTS = ("(isEB||isEE) && (abs(eta())<= 2.5) "
                      "&& (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1)"
                      "&& ( (isEB"
                      "      && (sigmaIetaIeta<0.01)"
                      "      && ( abs(deltaPhiSuperClusterTrackAtVtx)<0.06 )"
                      "      && ( abs(deltaEtaSuperClusterTrackAtVtx)<0.004 )"
                      "      && (hadronicOverEm<0.12)"
                      "      )"
                      "     || (isEE"
                      "         && (sigmaIetaIeta<0.03)"
                      "         && ( abs(deltaPhiSuperClusterTrackAtVtx)<0.03 )"
                      "         && ( abs(deltaEtaSuperClusterTrackAtVtx)<0.007 )"
                      "         && (hadronicOverEm<0.1) "
                      "         )"
                      "    )"
                      "&& passConversionVeto")
    
    JET_CUTS =  ("pt() >= 40.0 && abs(eta()) <= 2.4 "
                 " && neutralHadronEnergyFraction() < 0.99 "
                 " && neutralEmEnergyFraction() < 0.99 "
                 " && chargedEmEnergyFraction() < 0.99 "
                 " && chargedHadronEnergyFraction() > 0 "
                 " && chargedMultiplicity() > 0 "
                 " && (chargedMultiplicity() + neutralMultiplicity() + muonMultiplicity()) > 1")
    
    # 2012 SingleElectron Triggers
    PASS_ANY =  '!triggerObjectMatchesByType("TriggerElectron").empty() && ('
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele22_CaloIdL_CaloIsoVL_v*"                 ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele27_CaloIdL_CaloIsoVL_v*"                 ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele27_WP80_v*"                              ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele30_CaloIdVT_TrkIdT_v*"                   ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele80_CaloIdVT_GsfTrkIdT_v*"                ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_Ele90_CaloIdVT_GsfTrkIdT_v*"                ,0,1).empty() || '
    # 2012 EleHad Triggers
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"     ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"     ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*" ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*" ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"    ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"    ,0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",0,1).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*"                              ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*"                              ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*"                          ,1,0).empty() || '
    PASS_ANY += '!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*"                          ,1,0).empty() '
    PASS_ANY += ')'
    
    AllTriggerFlags = cms.PSet(
        ############################################################
        #   DO NOT USE LONG VARIABLE NAMES - THERE'S A BUG ! ! !   #
        ############################################################
        # 2012 SingleElectron Triggers
        passing_HLT_Ele22      = cms.string('!triggerObjectMatchesByPath("HLT_Ele22_CaloIdL_CaloIsoVL_v*"                       ,1,0).empty()'),
        passing_HLT_Ele27      = cms.string('!triggerObjectMatchesByPath("HLT_Ele27_CaloIdL_CaloIsoVL_v*"                       ,1,0).empty()'),
        passing_HLT_Ele27_WP80 = cms.string('!triggerObjectMatchesByPath("HLT_Ele27_WP80_v*"                                    ,1,0).empty()'),
        passing_HLT_Ele30      = cms.string('!triggerObjectMatchesByPath("HLT_Ele30_CaloIdVT_TrkIdT_v*"                         ,1,0).empty()'),
        passing_HLT_Ele32      = cms.string('!triggerObjectMatchesByPath("HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"      ,1,0).empty()'),
        passing_HLT_Ele80      = cms.string('!triggerObjectMatchesByPath("HLT_Ele80_CaloIdVT_GsfTrkIdT_v*"                      ,0,1).empty()'),
        passing_HLT_Ele90      = cms.string('!triggerObjectMatchesByPath("HLT_Ele90_CaloIdVT_GsfTrkIdT_v*"                      ,0,1).empty()'),
        # 2012 EleHad Triggers
        passing_HLT_CleanPFHT350_Ele5_PFMET45                = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"     ,0,1).empty()'),
        passing_HLT_CleanPFHT350_Ele5_PFMET50                = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"     ,0,1).empty()'),
        passing_HLT_CleanPFNoPUHT350_Ele5_PFMET45            = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*" ,0,1).empty()'),
        passing_HLT_CleanPFNoPUHT350_Ele5_PFMET50            = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*" ,0,1).empty()'),
        passing_HLT_CleanPFHT300_Ele15_PFMET45               = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"    ,0,1).empty()'),
        passing_HLT_CleanPFHT300_Ele15_PFMET50               = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"    ,0,1).empty()'),
        passing_HLT_CleanPFNoPUHT300_Ele15_PFMET45           = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",0,1).empty()'),
        passing_HLT_CleanPFNoPUHT300_Ele15_PFMET50           = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",0,1).empty()'),
        passing_HLT_CleanPFHT300_Ele40                       = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*"                              ,0,1).empty()'),
        passing_HLT_CleanPFHT300_Ele60                       = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*"                              ,0,1).empty()'),
        passing_HLT_CleanPFNoPUHT300_Ele40                   = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*"                          ,0,1).empty()'),
        passing_HLT_CleanPFNoPUHT300_Ele60                   = cms.string('!triggerObjectMatchesByPath("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*"                          ,0,1).empty()'),
    )

    ################################################################################################################################
    #   ___ ____ ____ ____ 
    #    |  |__| | __ [__  
    #    |  |  | |__] ___] 
    
    process.tagPATElectrons = cms.EDFilter("PATElectronRefSelector",
        src = cms.InputTag("cleanPatElectronsTriggerMatch"),
        cut = cms.string( TAG_CUTS + " && " + PASS_ANY )
    )

    ################################################################################################################################
    #   ___  ____ ____ ___  ____ ____ 
    #   |__] |__/ |  | |__] |___ [__  
    #   |    |  \ |__| |__] |___ ___] 
    
    ##   ____                         ____ _           _            
    ##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
    ##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
    ##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
    ##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
    ##  
    
    # probe1: superclusters
    #  SuperClusters  ################
    process.superClusters = cms.EDProducer("SuperClusterMerger",
       src = cms.VInputTag(cms.InputTag( "correctedHybridSuperClusters" ,""),
                           cms.InputTag( "correctedMulti5x5SuperClustersWithPreshower" ,"") )  
    )
    
    process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
       src = cms.InputTag("superClusters"),
       particleType = cms.int32(11),
    )
    
    #   Get the above SC's Candidates and place a cut on their Et and eta
    process.goodSuperClusters = cms.EDFilter("CandViewSelector",
          src = cms.InputTag("superClusterCands"),
          cut = cms.string( SC_PROBE_CUTS ),
          filter = cms.bool(True)
    )
    
    #### remove real jets (with high hadronic energy fraction) from SC collection
    ##### this improves the purity of the probe sample without affecting efficiency    
    process.JetsToRemoveFromSuperCluster = cms.EDFilter("CaloJetSelector",   
        src = cms.InputTag("ak5CaloJets"),
        cut = cms.string('pt>5 && energyFractionHadronic > 0.15')
    )
    process.goodSuperClustersClean = cms.EDProducer("CandViewCleaner",
        srcObject = cms.InputTag("goodSuperClusters"),
        module_label = cms.string(''),
        srcObjectsToRemove = cms.VInputTag(cms.InputTag("JetsToRemoveFromSuperCluster")),
        deltaRMin = cms.double(0.1)
    )

    ##    ____      __ _____ _           _                   
    ##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
    ##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
    ##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
    ##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
    ##

    # probe2: GsfElectrons
    #  GsfElectron ################ 
    process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
        src = cms.InputTag("gsfElectrons"),
        cut = cms.string( GSF_PROBE_CUTS )
    )
    
    # probe3: Pat Electrons
    process.goodPATElectrons = cms.EDFilter("PATElectronRefSelector",
        src = cms.InputTag("cleanPatElectronsTriggerMatch"),
        cut = cms.string( PAT_PROBE_CUTS ),
    )
    
    process.allTagsAndProbes = cms.Sequence(
        process.tagPATElectrons +
        (process.superClusters *
         process.superClusterCands *
         (process.goodSuperClusters + process.JetsToRemoveFromSuperCluster) *
         process.goodSuperClustersClean) +
        process.goodElectrons +
        process.goodPATElectrons
    )

    ################################################################################################################################
    #   ___ ____ ____    _  _    ___  ____ ____ ___  ____    ___  ____ _ ____ ____ 
    #    |  |__| | __    |\ |    |__] |__/ |  | |__] |___    |__] |__| | |__/ [__  
    #    |  |  | |__]    | \|    |    |  \ |__| |__] |___    |    |  | | |  \ ___] 
    #
    ##    _____ ___   ____    ____       _          
    ##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
    ##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
    ##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
    ##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
    ##                                              
    ##   
    #  Tag & probe selection ######
    process.tagPATSC = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("tagPATElectrons@+ goodSuperClustersClean@-"),
        checkCharge = cms.bool(False),                           
        cut   = cms.string("40 < mass < 140"),
    )
    process.tagPATGsf = process.tagPATSC.clone()
    process.tagPATGsf.decay = cms.string("tagPATElectrons@+ goodElectrons@-")
    process.tagPATGoodPATElectron = process.tagPATSC.clone()
    process.tagPATGoodPATElectron.decay = cms.string("tagPATElectrons@+ goodPATElectrons@-")

    process.allTPPairs = cms.Sequence(
        process.tagPATSC +
        process.tagPATGsf +
        process.tagPATGoodPATElectron
    )

    ################################################################################################################################
    #   _  _ ____ ___ ____ _  _    ____ _  _ ___     ___  ____ ____ ____ 
    #   |\/| |__|  |  |    |__|    |__| |\ | |  \    |__] |__| [__  [__  
    #   |  | |  |  |  |___ |  |    |  | | \| |__/    |    |  | ___] ___] 
    
    ##    ____   ____       __     ____      __ 
    ##   / ___| / ___|      \ \   / ___|___ / _|
    ##   \___ \| |      _____\ \ | |  _/ __| |_ 
    ##    ___) | |___  |_____/ / | |_| \__ \  _|
    ##   |____/ \____|      /_/   \____|___/_|  
    ##
    # passing1: Superclusters passing as Gsf Electrons

    process.GsfMatchedSuperClusterCands = cms.EDProducer("ElectronMatchedCandidateProducer",
       src     = cms.InputTag("goodSuperClustersClean"),
       ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
       deltaR =  cms.untracked.double(0.3)
    )

    ##   ____      __       __    ___                 ___    _ 
    ##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
    ## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
    ## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
    ##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
    ##                                           |/            
    # passing2: Gsf electrons passing isolation, ID cuts (Pat electrons)
    process.GSFPassingGoodPat = cms.EDProducer("MatchGsfElectronsToPAT",
        electrons   = cms.InputTag("goodElectrons"),
        pat = cms.InputTag("goodPATElectrons"),
        patCut = cms.string("pt>0"),
        matchByReference = cms.bool(False)
    )
    
    # Same can be done for Gsf Matched Superclusters
    process.GSFtoPATMatchedSuperClusterCandsClean = cms.EDProducer("ElectronMatchedCandidateProducer",
       src     = cms.InputTag("goodSuperClustersClean"),
       ReferenceElectronCollection = cms.untracked.InputTag("GSFPassingGoodPat"),
       deltaR =  cms.untracked.double(0.3)
    )
    
    process.allPassingProbes = cms.Sequence(
        process.GsfMatchedSuperClusterCands +
        process.GSFPassingGoodPat * process.GSFtoPATMatchedSuperClusterCandsClean
    )
    
    ##    ___    _       __    _   _ _   _____ 
    ##   |_ _|__| |      \ \  | | | | | |_   _|
    ##    | |/ _` |  _____\ \ | |_| | |   | |  
    ##    | | (_| | |_____/ / |  _  | |___| |  
    ##   |___\__,_|      /_/  |_| |_|_____|_|
    ##
    
    # passing 3: ID electrons passing HLT paths
    # corresponding matches and passes can be found in the top CUTS sections
    
    # Additional variable to check if event passed the specific cross trigger
    process.patEventPassingHLTCleanPFHT350Ele5PFMET45 = cms.EDProducer("HLTResultProducer",
        probes = cms.InputTag("cleanPatElectronsTriggerMatch"),
        TriggerResultsTag = cms.InputTag("TriggerResults","",hltName),
        HLTPaths = cms.vstring("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v.*"),
        andOr = cms.bool(True)
    )
    process.patEventPassingHLTCleanPFHT350Ele5PFMET50      = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v.*"))
    process.patEventPassingHLTCleanPFHT300Ele15PFMET45     = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v.*"))
    process.patEventPassingHLTCleanPFHT300Ele15PFMET50     = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v.*"))
    process.patEventPassingHLTCleanPFHT300Ele40            = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v.*"))
    process.patEventPassingHLTCleanPFHT300Ele60            = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT350Ele5PFMET45  = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT350Ele5PFMET50  = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT300Ele15PFMET45 = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT300Ele15PFMET50 = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT300Ele40        = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v.*"))
    process.patEventPassingHLTCleanPFNoPUHT300Ele60        = process.patEventPassingHLTCleanPFHT350Ele5PFMET45.clone(HLTPaths = cms.vstring("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v.*"))
    
    process.allHLTResults = cms.Sequence(
        process.patEventPassingHLTCleanPFHT350Ele5PFMET45 +
        process.patEventPassingHLTCleanPFHT350Ele5PFMET50 +
        process.patEventPassingHLTCleanPFHT300Ele15PFMET45 +
        process.patEventPassingHLTCleanPFHT300Ele15PFMET50 +
        process.patEventPassingHLTCleanPFHT300Ele40 +
        process.patEventPassingHLTCleanPFHT300Ele60 +
        process.patEventPassingHLTCleanPFNoPUHT350Ele5PFMET45 +
        process.patEventPassingHLTCleanPFNoPUHT350Ele5PFMET50 +
        process.patEventPassingHLTCleanPFNoPUHT300Ele15PFMET45 +
        process.patEventPassingHLTCleanPFNoPUHT300Ele15PFMET50 +
        process.patEventPassingHLTCleanPFNoPUHT300Ele40 +
        process.patEventPassingHLTCleanPFNoPUHT300Ele60
    )
    
    ##    __  __  ____   __  __       _       _               
    ##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
    ##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
    ##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
    ##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
    ##                                                        
    process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId = cms.vint32(11),
        src = cms.InputTag("goodSuperClustersClean"),
        distMin = cms.double(0.3),
        matched = cms.InputTag("genParticles")
    )
    process.McMatchGsf = process.McMatchSC.clone()
    process.McMatchGsf.src = cms.InputTag("goodElectrons")
    process.McMatchGSF = process.McMatchSC.clone()
    process.McMatchGSF.src = cms.InputTag("gsfElectrons")
    process.McMatchPATElectron = process.McMatchSC.clone()
    process.McMatchPATElectron.src = cms.InputTag("goodPATElectrons")
    process.McMatchTagPATElectron = process.McMatchSC.clone()
    process.McMatchTagPATElectron.src = cms.InputTag("tagPATElectrons")
    
    process.allMcMatches = cms.Sequence(
       process.McMatchSC +
       process.McMatchGSF +
       process.McMatchGsf +
       process.McMatchPATElectron +
       process.McMatchTagPATElectron
    )
    
    ################################################################################################################################
    ##    _____      _                        _  __     __             
    ##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
    ##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
    ##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
    ##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
    ##

    ################################################################################################################################
    #   Vertices
    
    process.scNvertices = cms.EDProducer("VertexMultiplicityCounter", 
        probes = cms.InputTag("goodSuperClusters"),
        objects = cms.InputTag("offlinePrimaryVertices"),
        objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    )
    process.gsfNvertices = process.scNvertices.clone(probes = "gsfElectrons")
    process.patNvertices = process.scNvertices.clone(probes = "cleanPatElectronsTriggerMatch")
    
    process.allVertices = cms.Sequence(
        process.scNvertices +
        process.gsfNvertices +
        process.patNvertices
    )
    
    ################################################################################################################################
    #   PU Reweighting

    process.scPileup = cms.EDProducer("PileupWeightComputer",
        probes = cms.InputTag("goodSuperClusters"),
        isMC = cms.bool(mcInfo),
        dataPileupFile = cms.string("PileupHistogram_2012Data_FlatPU_50bins.root"),
        mcPileupFile   = cms.string("PileupHistogram_2012Data_FlatPU_50bins.root"),
        dataPileupHistoName = cms.string("pileup"),
        mcPileupHistoName = cms.string("mcpileup"),
        mcLumiScale = cms.double(221.95),
        dataPileupInputFile = cms.string("run_ls_instlumi_pileup_2012.txt"),
    )
    process.gsfPileup = process.scPileup.clone(probes = cms.InputTag("gsfElectrons"))
    process.patPileup = process.scPileup.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))

    process.allPileup = cms.Sequence(
        process.scPileup +
        process.gsfPileup +
        process.patPileup
    )
    
    ################################################################################################################################
    #   _ _  _ ___  ____ ____ ___    ___  ____ ____ ____ _  _ ____ ___ ____ ____ 
    #   | |\/| |__] |__| |     |     |__] |__| |__/ |__| |\/| |___  |  |___ |__/ 
    #   | |  | |    |  | |___  |     |    |  | |  \ |  | |  | |___  |  |___ |  \ 
    
    process.gsfImpactParameter = cms.EDProducer("GsfElectronImpactParameter",
        probes = cms.InputTag("gsfElectrons"),
    )
     
    process.patImpactParameter = cms.EDProducer("PatElectronImpactParameter",
        probes = cms.InputTag("cleanPatElectronsTriggerMatch"),
    )

    process.allImpactParameters = cms.Sequence(
        process.gsfImpactParameter +
        process.patImpactParameter
    )
    
    ################################################################################################################################
    #    _ ____ ___ ____ 
    #    | |___  |  [__  
    #   _| |___  |  ___] 
    
    process.selectedJets = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("cleanPatJetsAK5PF"),
        cut = cms.string( JET_CUTS ), # <= anpassen
    )

    # DeltaR
    process.scDRToNearestJet = cms.EDProducer("minCutDeltaRNearestPatJetComputer",
        probes = cms.InputTag("goodSuperClusters"),
        objects = cms.InputTag("selectedJets"),
        minDeltaR = cms.double(0.3),
        objectSelection = cms.InputTag(""),
    )
    process.gsfDRToNearestJet = process.scDRToNearestJet.clone(probes = cms.InputTag("gsfElectrons"))
    process.patDRToNearestJet = process.scDRToNearestJet.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))

    # Njet
    process.scJetMultiplicity = cms.EDProducer("PatJetMultiplicityCounter",
        probes = cms.InputTag("goodSuperClusters"),
        objects = cms.InputTag("selectedJets"),
        minDeltaR = cms.double(0.3),
        objectSelection = cms.InputTag(""),
    )
    process.gsfJetMultiplicity = process.scJetMultiplicity.clone(probes = cms.InputTag("gsfElectrons"))
    process.patJetMultiplicity = process.scJetMultiplicity.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))

    # HT
    process.scHT = cms.EDProducer("PatJetHTComputer",
        probes = cms.InputTag("goodSuperClusters"),
        objects = cms.InputTag("selectedJets"),
        objectSelection = cms.InputTag(""),
    )
    process.gsfHT = process.scHT.clone(probes = cms.InputTag("gsfElectrons"))
    process.patHT = process.scHT.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))
    
    # MET
    process.scMet = cms.EDProducer("PatMetAssociator",
        probes = cms.InputTag("goodSuperClusters"),
        metTag = cms.InputTag("patMETsPF"),
    )
    process.gsfMet = process.scMet.clone(probes = cms.InputTag("gsfElectrons"))
    process.patMet = process.scMet.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))
    
    # ST
    process.gsfST = cms.EDProducer("PatMetSTComputer",
        probes = cms.InputTag("gsfElectrons"),
        metTag = cms.InputTag("patMETsPF"),
    )
    process.patST = process.gsfST.clone(probes = cms.InputTag("cleanPatElectronsTriggerMatch"))
    
    process.allJets = cms.Sequence(
        process.selectedJets*
        (process.scDRToNearestJet +
         process.gsfDRToNearestJet +
         process.patDRToNearestJet +
         process.scJetMultiplicity +
         process.gsfJetMultiplicity +
         process.patJetMultiplicity +
         process.scHT +
         process.gsfHT +
         process.patHT)
    )
    
    process.allMet = cms.Sequence(
        process.scMet +
        process.gsfMet +
        process.patMet
    )
    
    process.allST = cms.Sequence(
        process.gsfST +
        process.patST
    )
    
    ################################################################################################################################
    #   _ ____ ____ _    ____ ___ _ ____ _  _ ____ 
    #   | [__  |  | |    |__|  |  | |  | |\ | [__  
    #   | ___] |__| |___ |  |  |  | |__| | \| ___] 

    #compute rho for 2011 effective area Egamma isolation corrections
    from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
    process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
    process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
    
    #   Rho Corrected Relative Isolation
    process.gsfRelIso = cms.EDProducer("GsfElectronRelIsoProducer",
        isMC            = cms.bool(mcInfo),
        ElectronProbes  = cms.InputTag("gsfElectrons"),
        rhoIsoInputTag  = cms.InputTag("kt6PFJetsForIsolation", "rho"),
        isoValInputTags = cms.VInputTag(
            cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
        ),
    )
    process.patRelIso = cms.EDProducer("PatElectronRelIsoProducer",
        isMC            = cms.bool(mcInfo),
        ElectronProbes  = cms.InputTag("cleanPatElectronsTriggerMatch"),
        rhoIsoInputTag  = cms.InputTag("kt6PFJetsForIsolation", "rho"),
        isoValInputTags = cms.VInputTag(
            cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
        ),
    )    
    
    process.allRelIso = cms.Sequence(
        process.kt6PFJetsForIsolation *
        ( process.gsfRelIso +
          process.patRelIso )
    )
    
    ################################################################################################################################
    #   Delta Reco-PF Electron Pt
    
    process.gsfDeltaPfRecoPt = cms.EDProducer("GsfElectronDeltaPfRecoPt",
        probes = cms.InputTag("gsfElectrons"),
    )
    process.patDeltaPfRecoPt = cms.EDProducer("ElectronDeltaPfRecoPt",
        probes = cms.InputTag("cleanPatElectronsTriggerMatch"),
    )
    
    process.allDeltaPfRecoPt = cms.Sequence(
        process.gsfDeltaPfRecoPt +
        process.patDeltaPfRecoPt
    )
    
    ################################################################################################################################
    #   Conversion Rejection
    
    process.gsfConvRejVars = cms.EDProducer("ElectronConversionRejectionVars",
        probes = cms.InputTag("gsfElectrons")
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
            tagMatches = cms.InputTag("McMatchTagPATElectron"),
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
        arbitration = cms.string("Random2"),
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
    
    scParameters = cms.PSet(
        variables = cms.PSet(
            pt      = cms.string("pt"),
            et      = cms.string("et"),
            phi     = cms.string("phi"),
            eta     = cms.string("eta"),
            abs_eta = cms.string("abs(eta)"),
            nvtx        = cms.InputTag("scNvertices"),
            pileup      = cms.InputTag("scPileup","pileup"),
            instlumi    = cms.InputTag("scPileup","instlumi"),
            weight      = cms.InputTag("scPileup","weight"),
            drjet       = cms.InputTag("scDRToNearestJet"),
            njet        = cms.InputTag("scJetMultiplicity"),
            ht          = cms.InputTag("scHT"),
            met         = cms.InputTag("scMet"),
            #e       = cms.string("energy"),
            #p       = cms.string("p"),
            #px      = cms.string("px"),
            #py      = cms.string("py"),
            #pz      = cms.string("pz"),
            #theta   = cms.string("theta"),
        ),
        tagVariables = cms.PSet(),
        flags = cms.PSet(
            passing_gsf = cms.InputTag("GsfMatchedSuperClusterCands"),
            passing_pat = cms.InputTag("GSFtoPATMatchedSuperClusterCandsClean"),
        ),
        tagFlags = cms.PSet(),
    )
    
    gsfParameters = cms.PSet(
        variables = cms.PSet(
            pt           = cms.string("pt"),
            et           = cms.string("et"),
            phi          = cms.string("phi"),
            eta          = cms.string("eta"),
            abs_eta      = cms.string("abs(eta)"),
            sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
            sc_phi   = cms.string("superCluster.phi"),
            sc_eta   = cms.string("superCluster.eta"),
            track_pt       = cms.string("gsfTrack.pt"),
            track_phi      = cms.string("gsfTrack.phi"),
            track_eta      = cms.string("gsfTrack.eta"),
            nvtx        = cms.InputTag("gsfNvertices"),
            pileup      = cms.InputTag("gsfPileup","pileup"),
            instlumi    = cms.InputTag("gsfPileup","instlumi"),
            weight      = cms.InputTag("gsfPileup","weight"),
            d0_v        = cms.InputTag("gsfImpactParameter","d0v"),
            d0_b        = cms.InputTag("gsfImpactParameter","d0b"),
            dz_v        = cms.InputTag("gsfImpactParameter","dzv"),
            dz_b        = cms.InputTag("gsfImpactParameter","dzb"),
            drjet       = cms.InputTag("gsfDRToNearestJet"),
            njet        = cms.InputTag("gsfJetMultiplicity"),
            ht          = cms.InputTag("gsfHT"),
            met         = cms.InputTag("gsfMet"),
            st          = cms.InputTag("gsfST"),
            absdeltapt  = cms.InputTag("gsfDeltaPfRecoPt","absdeltapt"), # -9999: under 10 GeV
            reliso      = cms.InputTag("gsfRelIso","reliso"),
            passConvRej = cms.InputTag("gsfConvRejVars","passConvRej"),
            #conv_dist   = cms.InputTag("gsfConvRejVars","dist"),
            #conv_dcot   = cms.InputTag("gsfConvRejVars","dcot"),
            #conv_radius = cms.InputTag("gsfConvRejVars","convradius"),
            #theta        = cms.string("theta"),    
            #rapidity     = cms.string("rapidity"),
            #e            = cms.string("energy"),
            #p            = cms.string("p"),
            #px           = cms.string("px"),
            #py           = cms.string("py"),
            #pz           = cms.string("pz"),
            #charge       = cms.string("charge"),
            #missingHits  = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
            #convDist     = cms.string("convDist"),
            #convDcot     = cms.string("convDcot"),
            #convRadius   = cms.string("convRadius"),        
            #hasL1BPixHit = cms.string("gsfTrack.hitPattern.hasValidHitInFirstPixelBarrel"),
            ## super cluster quantities
            #sc_e     = cms.string("superCluster.energy"),
            #sc_x     = cms.string("superCluster.x"),
            #sc_y     = cms.string("superCluster.y"),
            #sc_z     = cms.string("superCluster.z"),
            #sc_theta = cms.string("superClusterPosition.theta"),
            #sc_nhit  = cms.string("superCluster.size"),
            ## track quantities
            #track_p        = cms.string("gsfTrack.p"),
            #track_px       = cms.string("gsfTrack.px"),
            #track_py       = cms.string("gsfTrack.py"),
            #track_pz       = cms.string("gsfTrack.pz"),
            #track_theta    = cms.string("gsfTrack.theta"),   
            #track_vx       = cms.string("gsfTrack.vx"),
            #track_vy       = cms.string("gsfTrack.vy"),
            #track_vz       = cms.string("gsfTrack.vz"),    
            #track_dxy      = cms.string("gsfTrack.dxy"),
            #track_d0       = cms.string("gsfTrack.d0"),
            #track_dsz      = cms.string("gsfTrack.dsz"),
            #track_charge   = cms.string("gsfTrack.charge"),
            #track_qoverp   = cms.string("gsfTrack.qoverp"),
            #track_normChi2 = cms.string("gsfTrack.normalizedChi2"),
            ## isolation 
            #trackiso = cms.string("dr03TkSumPt"),
            #ecaliso  = cms.string("dr03EcalRecHitSumEt"),
            #hcaliso  = cms.string("dr03HcalTowerSumEt"),
            ## classification, location, etc.    
            #classification = cms.string("classification"),
            #numberOfBrems  = cms.string("numberOfBrems"),     
            #bremFraction   = cms.string("fbrem"),
            #mva            = cms.string("mva"),        
            #deltaEta       = cms.string("deltaEtaSuperClusterTrackAtVtx"),
            #deltaPhi       = cms.string("deltaPhiSuperClusterTrackAtVtx"),
            #deltaPhiOut    = cms.string("deltaPhiSeedClusterTrackAtCalo"),
            #deltaEtaOut    = cms.string("deltaEtaSeedClusterTrackAtCalo"),
            #isEB           = cms.string("isEB"),
            #isEE           = cms.string("isEE"),
            #isGap          = cms.string("isGap"),
            ## Hcal energy over Ecal Energy
            #HoverE         = cms.string("hcalOverEcal"),    
            #EoverP         = cms.string("eSuperClusterOverP"),
            #eSeedClusterOverP = cms.string("eSeedClusterOverP"),    
            ## Cluster shape information
            #sigmaEtaEta  = cms.string("sigmaEtaEta"),
            #sigmaIetaIeta = cms.string("sigmaIetaIeta"),
            #e1x5               = cms.string("e1x5"),
            #e2x5Max            = cms.string("e2x5Max"),
            #e5x5               = cms.string("e5x5"),
            ## is ECAL driven ? is Track driven ?
            #ecalDrivenSeed     = cms.string("ecalDrivenSeed"),
            #trackerDrivenSeed  = cms.string("trackerDrivenSeed"),
        ),
        tagVariables = cms.PSet(),
        flags = cms.PSet(
            passing_pat = cms.InputTag("GSFPassingGoodPat"),
        ),
        tagFlags = cms.PSet(),
    )
    
    patParameters = cms.PSet(
        variables = cms.PSet(
            pt        = cms.string("pt"),
            et        = cms.string("et"),
            phi       = cms.string("phi"),
            eta       = cms.string("eta"),
            abs_eta   = cms.string("abs(eta)"),
            sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
            sc_phi    = cms.string("superCluster.phi"),
            sc_eta    = cms.string("superCluster.eta"), 
            track_pt  = cms.string("gsfTrack.pt"),
            track_phi = cms.string("gsfTrack.phi"),
            track_eta = cms.string("gsfTrack.eta"),
            nvtx       = cms.InputTag("patNvertices"),
            pileup     = cms.InputTag("patPileup","pileup"),
            instlumi   = cms.InputTag("patPileup","instlumi"),
            weight     = cms.InputTag("patPileup","weight"),
            d0_v       = cms.InputTag("patImpactParameter","d0v"),
            d0_b       = cms.InputTag("patImpactParameter","d0b"),
            dz_v       = cms.InputTag("patImpactParameter","dzv"),
            dz_b       = cms.InputTag("patImpactParameter","dzb"),
            drjet      = cms.InputTag("patDRToNearestJet"),
            njet       = cms.InputTag("patJetMultiplicity"),
            ht         = cms.InputTag("patHT"),
            met        = cms.InputTag("patMet"),
            st         = cms.InputTag("patST"),
            absdeltapt = cms.InputTag("patDeltaPfRecoPt","absdeltapt"), # -9999: under 10 GeV
            reliso     = cms.InputTag("patRelIso","reliso"),
            #charge    = cms.string("charge"),
            #isEB      = cms.string("isEB"),
            #isEE      = cms.string("isEE"),
            #isGap     = cms.string("isGap"),
            #trackiso  = cms.string("dr03TkSumPt"),
            #ecaliso   = cms.string("dr03EcalRecHitSumEt"),
            #hcaliso   = cms.string("dr03HcalTowerSumEt"),
            event_passing_HLT_CleanPFHT350_Ele5_PFMET45      = cms.InputTag("patEventPassingHLTCleanPFHT350Ele5PFMET45"),
            event_passing_HLT_CleanPFHT350_Ele5_PFMET50      = cms.InputTag("patEventPassingHLTCleanPFHT350Ele5PFMET50"),
            event_passing_HLT_CleanPFHT300_Ele15_PFMET45     = cms.InputTag("patEventPassingHLTCleanPFHT300Ele15PFMET45"),
            event_passing_HLT_CleanPFHT300_Ele15_PFMET50     = cms.InputTag("patEventPassingHLTCleanPFHT300Ele15PFMET50"),
            event_passing_HLT_CleanPFHT300_Ele40             = cms.InputTag("patEventPassingHLTCleanPFHT300Ele40"),
            event_passing_HLT_CleanPFHT300_Ele60             = cms.InputTag("patEventPassingHLTCleanPFHT300Ele60"),
            event_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET45  = cms.InputTag("patEventPassingHLTCleanPFNoPUHT350Ele5PFMET45"),
            event_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET50  = cms.InputTag("patEventPassingHLTCleanPFNoPUHT350Ele5PFMET50"),
            event_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET45 = cms.InputTag("patEventPassingHLTCleanPFNoPUHT300Ele15PFMET45"),
            event_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET50 = cms.InputTag("patEventPassingHLTCleanPFNoPUHT300Ele15PFMET50"),
            event_passing_HLT_CleanPFNoPUHT300_Ele40         = cms.InputTag("patEventPassingHLTCleanPFNoPUHT300Ele40"),
            event_passing_HLT_CleanPFNoPUHT300_Ele60         = cms.InputTag("patEventPassingHLTCleanPFNoPUHT300Ele60"),            
        ),
        tagVariables = cms.PSet(),
        flags    = cms.PSet( AllTriggerFlags ),
        tagFlags = cms.PSet( AllTriggerFlags ),
    )
    
    ################################################################################################################################
    #   ____ _ ___    ___ ____ ____ ____    ___  ____ ____ ___  _  _ ____ ____ ____ 
    #   |___ |  |      |  |__/ |___ |___    |__] |__/ |  | |  \ |  | |    |___ |__/ 
    #   |    |  |      |  |  \ |___ |___    |    |  \ |__| |__/ |__| |___ |___ |  \     
    
    ##    ____   ____       __     ____      __ 
    ##   / ___| / ___|      \ \   / ___|___ / _|
    ##   \___ \| |      _____\ \ | |  _/ __| |_ 
    ##    ___) | |___  |_____/ / | |_| \__ \  _|
    ##   |____/ \____|      /_/   \____|___/_|  
    ##
    
    process.SuperClusterToGsfElectronPATTag = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, scParameters,
        tagProbePairs = cms.InputTag("tagPATSC"),
        allProbes     = cms.InputTag("goodSuperClustersClean"),    
        probeMatches  = cms.InputTag("McMatchSC"),
    )
    
    ##   ____      __       __    ___                 ___    _ 
    ##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
    ## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
    ## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
    ##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
    ##                                           |/
    
    process.GsfElectronToIdPATTag = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, gsfParameters,
        tagProbePairs = cms.InputTag("tagPATGsf"),
        allProbes     = cms.InputTag("goodElectrons"),
        probeMatches  = cms.InputTag("McMatchGsf"),
    )
    
    ##    ___    _       __    _   _ _   _____ 
    ##   |_ _|__| |      \ \  | | | | | |_   _|
    ##    | |/ _` |  _____\ \ | |_| | |   | |  
    ##    | | (_| | |_____/ / |  _  | |___| |  
    ##   |___\__,_|      /_/  |_| |_|_____|_|
    ##
    
    process.goodPATEleToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
        commonStuff, patParameters,
        tagProbePairs = cms.InputTag("tagPATGoodPATElectron"),
        allProbes     = cms.InputTag("goodPATElectrons"),
        probeMatches  = cms.InputTag("McMatchPATElectron"),
    )
    
    process.tree_sequence = cms.Sequence(
        process.SuperClusterToGsfElectronPATTag +
        process.GsfElectronToIdPATTag +
        process.goodPATEleToHLT #+
        #process.goodGsfEleProbeTree
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
              process.allPileup +
              process.allImpactParameters +
              process.allJets + 
              process.allMet +
              process.allST +
              process.allDeltaPfRecoPt +
              process.allRelIso +
              process.gsfConvRejVars ) *
            process.tree_sequence
        )
    else:
        process.TagAndProbe = cms.Sequence(
            process.allTagsAndProbes *
            ( process.allTPPairs +
              process.allPassingProbes +
              process.allHLTResults +
              process.allVertices +
              process.allPileup +
              process.allImpactParameters +
              process.allJets + 
              process.allMet +
              process.allST +
              process.allDeltaPfRecoPt +
              process.allRelIso +
              process.gsfConvRejVars ) *
            process.tree_sequence
        )
