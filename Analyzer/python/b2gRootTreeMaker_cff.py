import FWCore.ParameterSet.Config as cms
import copy

leptonssize = cms.untracked.int32(5)
jetssize = cms.untracked.int32(20)

mulabel = cms.string("muons")
elelabel = cms.string("electrons")
jetlabel = cms.string("jets")
jetAK8label = cms.string("jetsAK8")
metlabel = cms.string("met")

DMTreesDumper = cms.EDAnalyzer(
    'DMAnalysisTreeMaker',
    lhes = cms.InputTag('source'),
    muLabel = mulabel,
    eleLabel = elelabel,
    jetsLabel = jetlabel,
    boostedTopsLabel = jetlabel,
    metLabel = metlabel,
    physicsObjects = cms.VPSet(
        cms.PSet(
            label = elelabel,
            prefix = cms.string("el"),
            maxInstances = leptonssize,
            variablesF = cms.VInputTag(
                cms.InputTag("electrons","elE"),
                cms.InputTag("electrons","elPt"),
                cms.InputTag("electrons","elMass"),
                cms.InputTag("electrons","elEta"),
                cms.InputTag("electrons","elPhi"),
                cms.InputTag("electrons","elCharge"),
                cms.InputTag("electrons","elD0"),
                cms.InputTag("electrons","elDz"),
                cms.InputTag("electrons","elEta"),
                cms.InputTag("electrons","elHoE"),
                cms.InputTag("electrons","elIso03"),
                cms.InputTag("electrons","elY"),
                cms.InputTag("electrons","eldEtaIn"),
                cms.InputTag("electrons","eldPhiIn"),
                cms.InputTag("electrons","elexpectedMissInHits"),
                cms.InputTag("electrons","elfull5x5siee"),
                cms.InputTag("electrons","elooEmooP"),
                cms.InputTag("electrons","elpssConVeto"),

                                       ),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            ),
        cms.PSet(
            label = mulabel,
            prefix = cms.string("mu"),
            maxInstances = leptonssize,
            variablesF = cms.VInputTag(
                cms.InputTag("muons","muE"),
                cms.InputTag("muons","muPt"),
                cms.InputTag("muons","muMass"),
                cms.InputTag("muons","muEta"),
                cms.InputTag("muons","muPhi"),
                cms.InputTag("muons","muCharge"),
                cms.InputTag("muons","muIsLooseMuon"),
                cms.InputTag("muons","muIsSoftMuon"),
                cms.InputTag("muons","muIsTightMuon"),
                cms.InputTag("muons","muD0"),
                cms.InputTag("muons","muD0err"),
                cms.InputTag("muons","muDz"),
                cms.InputTag("muons","muDzerr"),
                cms.InputTag("muons","muGenMuonCharge"),
                cms.InputTag("muons","muGenMuonEta"),
                cms.InputTag("muons","muGenMuonPt"),
                cms.InputTag("muons","muGenMuonE"),
                cms.InputTag("muons","muGenMuonPhi"),
                cms.InputTag("muons","muGenMuonY"),
                cms.InputTag("muons","muGlbTrkNormChi2"),
                cms.InputTag("muons","muHLTmuonDeltaR"),
                cms.InputTag("muons","muHLTmuonE"),
                cms.InputTag("muons","muHLTmuonEta"),
                cms.InputTag("muons","muHLTmuonPt"),
                cms.InputTag("muons","muHLTmuonPhi"),
                cms.InputTag("muons","muInTrkNormChi2"),
                cms.InputTag("muons","muIsGlobalMuon"),
                cms.InputTag("muons","muIsPFMuon"),
                cms.InputTag("muons","muIsTrackerMuon"),
                cms.InputTag("muons","muIso03"),
                cms.InputTag("muons","muNumberMatchedStations"),
                cms.InputTag("muons","muNumberOfPixelLayers"),
                cms.InputTag("muons","muNumberOfValidTrackerHits"),
                cms.InputTag("muons","muNumberTrackerLayers"),
                cms.InputTag("muons","muNumberValidMuonHits"),
                cms.InputTag("muons","muNumberValidPixelHits"),
                cms.InputTag("muons","muSumChargedHadronPt"),
                cms.InputTag("muons","muSumNeutralHadronPt"),
                cms.InputTag("muons","muSumPUPt"),
                cms.InputTag("muons","muSumPhotonPt"),
                cms.InputTag("muons","muY"),

),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            ),
        cms.PSet(
            label = jetlabel,
            prefix = cms.string("jet"),
            maxInstances = jetssize,
            variablesF = cms.VInputTag(
                cms.InputTag("jets","jetE"),
                cms.InputTag("jets","jetPt"),
                cms.InputTag("jets","jetMass"),
                cms.InputTag("jets","jetEta"),
                cms.InputTag("jets","jetPhi"),
                cms.InputTag("jets","jetPartonFlavour"),
                cms.InputTag("jets","jetPhi"),
                cms.InputTag("jets","jetCSV"),
                cms.InputTag("jets","jetCSVV1"),
                cms.InputTag("jets","jetCharge"),
                cms.InputTag("jets","jetChargeMuEnergy"),
                cms.InputTag("jets","jetChargedHadronMultiplicity"),
                cms.InputTag("jets","jetElectronEnergy"),
                cms.InputTag("jets","jetGenJetCharge"),
                cms.InputTag("jets","jetGenJetE"),
                cms.InputTag("jets","jetGenJetEta"),
                cms.InputTag("jets","jetGenJetPhi"),
                cms.InputTag("jets","jetGenJetPt"),
                cms.InputTag("jets","jetGenJetY"),
                cms.InputTag("jets","jetGenPartonCharge"),
                cms.InputTag("jets","jetGenPartonE"),
                cms.InputTag("jets","jetGenPartonEta"),
                cms.InputTag("jets","jetGenPartonPhi"),
                cms.InputTag("jets","jetGenPartonPt"),
                cms.InputTag("jets","jetGenPartonY"),
                cms.InputTag("jets","jetHFEMEnergy"),
                cms.InputTag("jets","jetHFEMMultiplicity"),
                cms.InputTag("jets","jetHFHadronEnergy"),
                cms.InputTag("jets","jetHFHadronMultiplicity"),
                cms.InputTag("jets","jetHLTjetDeltaR"),
                cms.InputTag("jets","jetHLTjetE"),
                cms.InputTag("jets","jetHLTjetEta"),
                cms.InputTag("jets","jetHLTjetPt"),
                cms.InputTag("jets","jetHLTjetPhi"),
                cms.InputTag("jets","jetHadronFlavour"),
                cms.InputTag("jets","jetIsCSVL"),
                cms.InputTag("jets","jetIsCSVM"),
                cms.InputTag("jets","jetIsCSVT"),
                cms.InputTag("jets","jetSmearedE"),
                cms.InputTag("jets","jetSmearedPt"),
                cms.InputTag("jets","jetSmearedPEta"),
                cms.InputTag("jets","jetSmearedPhi"),
                cms.InputTag("jets","jetY"),
                cms.InputTag("jets","jetelectronMultiplicity"),
                cms.InputTag("jets","jetmuonMultiplicity"),
                cms.InputTag("jets","jetneutralHadronMultiplicity"),
                cms.InputTag("jets","jetneutralMultiplicity"),
                cms.InputTag("jets","jetphotonMultiplicity"),
                                       ),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            ),
       cms.PSet(
            label = jetAK8label,
            prefix = cms.string("jetAK8"),
            maxInstances = jetssize,
            variablesF = cms.VInputTag(
                cms.InputTag("jetsAK8","jetAK8E"),
                cms.InputTag("jetsAK8","jetAK8Pt"),
                cms.InputTag("jetsAK8","jetAK8Mass"),
                cms.InputTag("jetsAK8","jetAK8Eta"),
                cms.InputTag("jetsAK8","jetAK8Phi"),
                cms.InputTag("jetsAK8","jetAK8PartonFlavour"),
                cms.InputTag("jetsAK8","jetAK8Phi"),
                cms.InputTag("jetsAK8","jetAK8CSV"),
                cms.InputTag("jetsAK8","jetAK8CSVV1"),
                cms.InputTag("jetsAK8","jetAK8Charge"),
                cms.InputTag("jetsAK8","jetAK8ChargeMuEnergy"),
                cms.InputTag("jetsAK8","jetAK8ChargedHadronMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8ElectronEnergy"),
                cms.InputTag("jetsAK8","jetAK8GenJetCharge"),
                cms.InputTag("jetsAK8","jetAK8GenJetE"),
                cms.InputTag("jetsAK8","jetAK8GenJetEta"),
                cms.InputTag("jetsAK8","jetAK8GenJetPhi"),
                cms.InputTag("jetsAK8","jetAK8GenJetPt"),
                cms.InputTag("jetsAK8","jetAK8GenJetY"),
                cms.InputTag("jetsAK8","jetAK8GenPartonCharge"),
                cms.InputTag("jetsAK8","jetAK8GenPartonE"),
                cms.InputTag("jetsAK8","jetAK8GenPartonEta"),
                cms.InputTag("jetsAK8","jetAK8GenPartonPhi"),
                cms.InputTag("jetsAK8","jetAK8GenPartonPt"),
                cms.InputTag("jetsAK8","jetAK8GenPartonY"),
                cms.InputTag("jetsAK8","jetAK8HFEMEnergy"),
                cms.InputTag("jetsAK8","jetAK8HFEMMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8HFHadronEnergy"),
                cms.InputTag("jetsAK8","jetAK8HFHadronMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8HLTjetDeltaR"),
                cms.InputTag("jetsAK8","jetAK8HLTjetE"),
                cms.InputTag("jetsAK8","jetAK8HLTjetEta"),
                cms.InputTag("jetsAK8","jetAK8HLTjetPt"),
                cms.InputTag("jetsAK8","jetAK8HLTjetPhi"),
                cms.InputTag("jetsAK8","jetAK8HadronFlavour"),
                cms.InputTag("jetsAK8","jetAK8IsCSVL"),
                cms.InputTag("jetsAK8","jetAK8IsCSVM"),
                cms.InputTag("jetsAK8","jetAK8IsCSVT"),
                cms.InputTag("jetsAK8","jetAK8SmearedE"),
                cms.InputTag("jetsAK8","jetAK8SmearedPt"),
                cms.InputTag("jetsAK8","jetAK8SmearedPEta"),
                cms.InputTag("jetsAK8","jetAK8SmearedPhi"),
                cms.InputTag("jetsAK8","jetAK8Y"),
                cms.InputTag("jetsAK8","jetAK8electronMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8muonMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8neutralHadronMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8neutralMultiplicity"),
                cms.InputTag("jetsAK8","jetAK8photonMultiplicity"),
                #cms.InputTag("jetsAK8","tau1"),
                #cms.InputTag("jetsAK8","tau2"),
                #cms.InputTag("jetsAK8","tau3"),
                #cms.InputTag("jetsAK8","trimmedMass"),
                #cms.InputTag("jetsAK8","prunedMass"),
                #cms.InputTag("jetsAK8","filteredMass"),
                # cms.InputTag("jetsAK8","minmass"),
                                       ),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            ),

        cms.PSet(
            label = metlabel,
            prefix = cms.string("met"),
            maxInstances = jetssize,
            variablesF = cms.VInputTag(
                cms.InputTag("met","metPt"),
                cms.InputTag("met","metPhi"),
                cms.InputTag("met","metPx"),
                cms.InputTag("met","metPy"),
                ),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            )
        )
    )
