#ifndef B2GTreeReader_h
#define B2GTreeReader_h

#include "TFile.h"
#include "TTree.h"

// Data Class storing ntuple variables
#include "../interface/Data.h"

class B2GTreeReader {
public:
  B2GTreeReader() {}
  ~B2GTreeReader() {}
  
  Data data;
  
  void Load_Tree(TFile &f, const char* name = "DMTreesDumper/TreeBase") {
    tree_ = (TTree*)f.Get(name);
    
    tree_->GetBranch("electrons_size")->SetAddress(&data.ele.electrons_size);
    tree_->GetBranch("elE")->SetAddress(&data.ele.elE);
    tree_->GetBranch("elPt")->SetAddress(&data.ele.elPt);
    tree_->GetBranch("elMass")->SetAddress(&data.ele.elMass);
    tree_->GetBranch("elEta")->SetAddress(&data.ele.elEta);
    tree_->GetBranch("elPhi")->SetAddress(&data.ele.elPhi);
    tree_->GetBranch("elCharge")->SetAddress(&data.ele.elCharge);
    tree_->GetBranch("elD0")->SetAddress(&data.ele.elD0);
    tree_->GetBranch("elDz")->SetAddress(&data.ele.elDz);
    tree_->GetBranch("elEta")->SetAddress(&data.ele.elEta);
    tree_->GetBranch("elHoE")->SetAddress(&data.ele.elHoE);
    tree_->GetBranch("elIso03")->SetAddress(&data.ele.elIso03);
    tree_->GetBranch("elY")->SetAddress(&data.ele.elY);
    tree_->GetBranch("eldEtaIn")->SetAddress(&data.ele.eldEtaIn);
    tree_->GetBranch("eldPhiIn")->SetAddress(&data.ele.eldPhiIn);
    tree_->GetBranch("elexpectedMissInHits")->SetAddress(&data.ele.elexpectedMissInHits);
    tree_->GetBranch("elfull5x5siee")->SetAddress(&data.ele.elfull5x5siee);
    tree_->GetBranch("elooEmooP")->SetAddress(&data.ele.elooEmooP);
    tree_->GetBranch("elpssConVeto")->SetAddress(&data.ele.elpssConVeto);
    tree_->GetBranch("electronsSFTrigger")->SetAddress(&data.ele.electronsSFTrigger);
    tree_->GetBranch("electronsSFReco")->SetAddress(&data.ele.electronsSFReco);
    tree_->GetBranch("electronsisQCD")->SetAddress(&data.ele.electronsisQCD);
    tree_->GetBranch("electronsisTightOffline")->SetAddress(&data.ele.electronsisTightOffline);
    tree_->GetBranch("electronsisLooseOffline")->SetAddress(&data.ele.electronsisLooseOffline);
    
    tree_->GetBranch("muons_size")->SetAddress(&data.mu.muons_size);
    tree_->GetBranch("muE")->SetAddress(&data.mu.muE);
    tree_->GetBranch("muPt")->SetAddress(&data.mu.muPt);
    tree_->GetBranch("muMass")->SetAddress(&data.mu.muMass);
    tree_->GetBranch("muEta")->SetAddress(&data.mu.muEta);
    tree_->GetBranch("muPhi")->SetAddress(&data.mu.muPhi);
    tree_->GetBranch("muCharge")->SetAddress(&data.mu.muCharge);
    tree_->GetBranch("muIsLooseMuon")->SetAddress(&data.mu.muIsLooseMuon);
    tree_->GetBranch("muIsSoftMuon")->SetAddress(&data.mu.muIsSoftMuon);
    tree_->GetBranch("muIsTightMuon")->SetAddress(&data.mu.muIsTightMuon);
    tree_->GetBranch("muD0")->SetAddress(&data.mu.muD0);
    tree_->GetBranch("muD0err")->SetAddress(&data.mu.muD0err);
    tree_->GetBranch("muDz")->SetAddress(&data.mu.muDz);
    tree_->GetBranch("muDzerr")->SetAddress(&data.mu.muDzerr);
    tree_->GetBranch("muGenMuonCharge")->SetAddress(&data.mu.muGenMuonCharge);
    tree_->GetBranch("muGenMuonEta")->SetAddress(&data.mu.muGenMuonEta);
    tree_->GetBranch("muGenMuonPt")->SetAddress(&data.mu.muGenMuonPt);
    tree_->GetBranch("muGenMuonE")->SetAddress(&data.mu.muGenMuonE);
    tree_->GetBranch("muGenMuonPhi")->SetAddress(&data.mu.muGenMuonPhi);
    tree_->GetBranch("muGenMuonY")->SetAddress(&data.mu.muGenMuonY);
    tree_->GetBranch("muGlbTrkNormChi2")->SetAddress(&data.mu.muGlbTrkNormChi2);
    tree_->GetBranch("muHLTmuonDeltaR")->SetAddress(&data.mu.muHLTmuonDeltaR);
    tree_->GetBranch("muHLTmuonE")->SetAddress(&data.mu.muHLTmuonE);
    tree_->GetBranch("muHLTmuonEta")->SetAddress(&data.mu.muHLTmuonEta);
    tree_->GetBranch("muHLTmuonPt")->SetAddress(&data.mu.muHLTmuonPt);
    tree_->GetBranch("muHLTmuonPhi")->SetAddress(&data.mu.muHLTmuonPhi);
    tree_->GetBranch("muInTrkNormChi2")->SetAddress(&data.mu.muInTrkNormChi2);
    tree_->GetBranch("muIsGlobalMuon")->SetAddress(&data.mu.muIsGlobalMuon);
    tree_->GetBranch("muIsPFMuon")->SetAddress(&data.mu.muIsPFMuon);
    tree_->GetBranch("muIsTrackerMuon")->SetAddress(&data.mu.muIsTrackerMuon);
    tree_->GetBranch("muIso03")->SetAddress(&data.mu.muIso03);
    tree_->GetBranch("muNumberMatchedStations")->SetAddress(&data.mu.muNumberMatchedStations);
    tree_->GetBranch("muNumberOfPixelLayers")->SetAddress(&data.mu.muNumberOfPixelLayers);
    tree_->GetBranch("muNumberOfValidTrackerHits")->SetAddress(&data.mu.muNumberOfValidTrackerHits);
    tree_->GetBranch("muNumberTrackerLayers")->SetAddress(&data.mu.muNumberTrackerLayers);
    tree_->GetBranch("muNumberValidMuonHits")->SetAddress(&data.mu.muNumberValidMuonHits);
    tree_->GetBranch("muNumberValidPixelHits")->SetAddress(&data.mu.muNumberValidPixelHits);
    tree_->GetBranch("muSumChargedHadronPt")->SetAddress(&data.mu.muSumChargedHadronPt);
    tree_->GetBranch("muSumNeutralHadronPt")->SetAddress(&data.mu.muSumNeutralHadronPt);
    tree_->GetBranch("muSumPUPt")->SetAddress(&data.mu.muSumPUPt);
    tree_->GetBranch("muSumPhotonPt")->SetAddress(&data.mu.muSumPhotonPt);
    tree_->GetBranch("muY")->SetAddress(&data.mu.muY);
    tree_->GetBranch("muonsSFTrigger")->SetAddress(&data.mu.muonsSFTrigger);
    tree_->GetBranch("muonsSFReco")->SetAddress(&data.mu.muonsSFReco);
    tree_->GetBranch("muonsisQCD")->SetAddress(&data.mu.muonsisQCD);
    tree_->GetBranch("muonsisTightOffline")->SetAddress(&data.mu.muonsisTightOffline);
    tree_->GetBranch("muonsisLooseOffline")->SetAddress(&data.mu.muonsisLooseOffline);
    
    tree_->GetBranch("jets_size")->SetAddress(&data.jet.jets_size);
    tree_->GetBranch("jetE")->SetAddress(&data.jet.jetE);
    tree_->GetBranch("jetPt")->SetAddress(&data.jet.jetPt);
    tree_->GetBranch("jetMass")->SetAddress(&data.jet.jetMass);
    tree_->GetBranch("jetEta")->SetAddress(&data.jet.jetEta);
    tree_->GetBranch("jetPhi")->SetAddress(&data.jet.jetPhi);
    tree_->GetBranch("jetPartonFlavour")->SetAddress(&data.jet.jetPartonFlavour);
    tree_->GetBranch("jetPhi")->SetAddress(&data.jet.jetPhi);
    tree_->GetBranch("jetCSV")->SetAddress(&data.jet.jetCSV);
    tree_->GetBranch("jetCSVV1")->SetAddress(&data.jet.jetCSVV1);
    tree_->GetBranch("jetCharge")->SetAddress(&data.jet.jetCharge);
    tree_->GetBranch("jetChargeMuEnergy")->SetAddress(&data.jet.jetChargeMuEnergy);
    tree_->GetBranch("jetChargedHadronMultiplicity")->SetAddress(&data.jet.jetChargedHadronMultiplicity);
    tree_->GetBranch("jetElectronEnergy")->SetAddress(&data.jet.jetElectronEnergy);
    tree_->GetBranch("jetGenJetCharge")->SetAddress(&data.jet.jetGenJetCharge);
    tree_->GetBranch("jetGenJetE")->SetAddress(&data.jet.jetGenJetE);
    tree_->GetBranch("jetGenJetEta")->SetAddress(&data.jet.jetGenJetEta);
    tree_->GetBranch("jetGenJetPhi")->SetAddress(&data.jet.jetGenJetPhi);
    tree_->GetBranch("jetGenJetPt")->SetAddress(&data.jet.jetGenJetPt);
    tree_->GetBranch("jetGenJetY")->SetAddress(&data.jet.jetGenJetY);
    tree_->GetBranch("jetGenPartonCharge")->SetAddress(&data.jet.jetGenPartonCharge);
    tree_->GetBranch("jetGenPartonE")->SetAddress(&data.jet.jetGenPartonE);
    tree_->GetBranch("jetGenPartonEta")->SetAddress(&data.jet.jetGenPartonEta);
    tree_->GetBranch("jetGenPartonPhi")->SetAddress(&data.jet.jetGenPartonPhi);
    tree_->GetBranch("jetGenPartonPt")->SetAddress(&data.jet.jetGenPartonPt);
    tree_->GetBranch("jetGenPartonY")->SetAddress(&data.jet.jetGenPartonY);
    tree_->GetBranch("jetHFEMEnergy")->SetAddress(&data.jet.jetHFEMEnergy);
    tree_->GetBranch("jetHFEMMultiplicity")->SetAddress(&data.jet.jetHFEMMultiplicity);
    tree_->GetBranch("jetHFHadronEnergy")->SetAddress(&data.jet.jetHFHadronEnergy);
    tree_->GetBranch("jetHFHadronMultiplicity")->SetAddress(&data.jet.jetHFHadronMultiplicity);
    tree_->GetBranch("jetHLTjetDeltaR")->SetAddress(&data.jet.jetHLTjetDeltaR);
    tree_->GetBranch("jetHLTjetE")->SetAddress(&data.jet.jetHLTjetE);
    tree_->GetBranch("jetHLTjetEta")->SetAddress(&data.jet.jetHLTjetEta);
    tree_->GetBranch("jetHLTjetPt")->SetAddress(&data.jet.jetHLTjetPt);
    tree_->GetBranch("jetHLTjetPhi")->SetAddress(&data.jet.jetHLTjetPhi);
    tree_->GetBranch("jetHadronFlavour")->SetAddress(&data.jet.jetHadronFlavour);
    tree_->GetBranch("jetIsCSVL")->SetAddress(&data.jet.jetIsCSVL);
    tree_->GetBranch("jetIsCSVM")->SetAddress(&data.jet.jetIsCSVM);
    tree_->GetBranch("jetIsCSVT")->SetAddress(&data.jet.jetIsCSVT);
    tree_->GetBranch("jetSmearedE")->SetAddress(&data.jet.jetSmearedE);
    tree_->GetBranch("jetSmearedPt")->SetAddress(&data.jet.jetSmearedPt);
    tree_->GetBranch("jetSmearedPEta")->SetAddress(&data.jet.jetSmearedPEta);
    tree_->GetBranch("jetSmearedPhi")->SetAddress(&data.jet.jetSmearedPhi);
    tree_->GetBranch("jetY")->SetAddress(&data.jet.jetY);
    tree_->GetBranch("jetelectronMultiplicity")->SetAddress(&data.jet.jetelectronMultiplicity);
    tree_->GetBranch("jetmuonMultiplicity")->SetAddress(&data.jet.jetmuonMultiplicity);
    tree_->GetBranch("jetneutralHadronMultiplicity")->SetAddress(&data.jet.jetneutralHadronMultiplicity);
    tree_->GetBranch("jetneutralMultiplicity")->SetAddress(&data.jet.jetneutralMultiplicity);
    tree_->GetBranch("jetphotonMultiplicity")->SetAddress(&data.jet.jetphotonMultiplicity);
    tree_->GetBranch("jetsCorrPt")->SetAddress(&data.jet.jetsCorrPt);
    tree_->GetBranch("jetsCorrEta")->SetAddress(&data.jet.jetsCorrEta);
    tree_->GetBranch("jetsCorrPhi")->SetAddress(&data.jet.jetsCorrPhi);
    tree_->GetBranch("jetsCorrE")->SetAddress(&data.jet.jetsCorrE);
    tree_->GetBranch("jetsCorrMass")->SetAddress(&data.jet.jetsCorrMass);
    tree_->GetBranch("jetsCorrNJets")->SetAddress(&data.jet.jetsCorrNJets);
    tree_->GetBranch("jetsCorrPartonFlavour")->SetAddress(&data.jet.jetsCorrPartonFlavour);
    
    tree_->GetBranch("jetsAK8_size")->SetAddress(&data.jetAK8.jetsAK8_size);
    tree_->GetBranch("jetAK8E")->SetAddress(&data.jetAK8.jetAK8E);
    tree_->GetBranch("jetAK8Pt")->SetAddress(&data.jetAK8.jetAK8Pt);
    tree_->GetBranch("jetAK8Mass")->SetAddress(&data.jetAK8.jetAK8Mass);
    tree_->GetBranch("jetAK8Eta")->SetAddress(&data.jetAK8.jetAK8Eta);
    tree_->GetBranch("jetAK8Phi")->SetAddress(&data.jetAK8.jetAK8Phi);
    tree_->GetBranch("jetAK8PartonFlavour")->SetAddress(&data.jetAK8.jetAK8PartonFlavour);
    tree_->GetBranch("jetAK8Phi")->SetAddress(&data.jetAK8.jetAK8Phi);
    tree_->GetBranch("jetAK8CSV")->SetAddress(&data.jetAK8.jetAK8CSV);
    tree_->GetBranch("jetAK8CSVV1")->SetAddress(&data.jetAK8.jetAK8CSVV1);
    tree_->GetBranch("jetAK8Charge")->SetAddress(&data.jetAK8.jetAK8Charge);
    tree_->GetBranch("jetAK8ChargeMuEnergy")->SetAddress(&data.jetAK8.jetAK8ChargeMuEnergy);
    tree_->GetBranch("jetAK8ElectronEnergy")->SetAddress(&data.jetAK8.jetAK8ElectronEnergy);
    tree_->GetBranch("jetAK8GenJetCharge")->SetAddress(&data.jetAK8.jetAK8GenJetCharge);
    tree_->GetBranch("jetAK8GenJetE")->SetAddress(&data.jetAK8.jetAK8GenJetE);
    tree_->GetBranch("jetAK8GenJetEta")->SetAddress(&data.jetAK8.jetAK8GenJetEta);
    tree_->GetBranch("jetAK8GenJetPhi")->SetAddress(&data.jetAK8.jetAK8GenJetPhi);
    tree_->GetBranch("jetAK8GenJetPt")->SetAddress(&data.jetAK8.jetAK8GenJetPt);
    tree_->GetBranch("jetAK8GenJetY")->SetAddress(&data.jetAK8.jetAK8GenJetY);
    tree_->GetBranch("jetAK8GenPartonCharge")->SetAddress(&data.jetAK8.jetAK8GenPartonCharge);
    tree_->GetBranch("jetAK8GenPartonE")->SetAddress(&data.jetAK8.jetAK8GenPartonE);
    tree_->GetBranch("jetAK8GenPartonEta")->SetAddress(&data.jetAK8.jetAK8GenPartonEta);
    tree_->GetBranch("jetAK8GenPartonPhi")->SetAddress(&data.jetAK8.jetAK8GenPartonPhi);
    tree_->GetBranch("jetAK8GenPartonPt")->SetAddress(&data.jetAK8.jetAK8GenPartonPt);
    tree_->GetBranch("jetAK8GenPartonY")->SetAddress(&data.jetAK8.jetAK8GenPartonY);
    tree_->GetBranch("jetAK8HFEMEnergy")->SetAddress(&data.jetAK8.jetAK8HFEMEnergy);
    tree_->GetBranch("jetAK8HFEMMultiplicity")->SetAddress(&data.jetAK8.jetAK8HFEMMultiplicity);
    tree_->GetBranch("jetAK8HFHadronEnergy")->SetAddress(&data.jetAK8.jetAK8HFHadronEnergy);
    tree_->GetBranch("jetAK8HFHadronMultiplicity")->SetAddress(&data.jetAK8.jetAK8HFHadronMultiplicity);
    tree_->GetBranch("jetAK8HLTjetDeltaR")->SetAddress(&data.jetAK8.jetAK8HLTjetDeltaR);
    tree_->GetBranch("jetAK8HLTjetE")->SetAddress(&data.jetAK8.jetAK8HLTjetE);
    tree_->GetBranch("jetAK8HLTjetEta")->SetAddress(&data.jetAK8.jetAK8HLTjetEta);
    tree_->GetBranch("jetAK8HLTjetPt")->SetAddress(&data.jetAK8.jetAK8HLTjetPt);
    tree_->GetBranch("jetAK8HLTjetPhi")->SetAddress(&data.jetAK8.jetAK8HLTjetPhi);
    tree_->GetBranch("jetAK8HadronFlavour")->SetAddress(&data.jetAK8.jetAK8HadronFlavour);
    tree_->GetBranch("jetAK8IsCSVL")->SetAddress(&data.jetAK8.jetAK8IsCSVL);
    tree_->GetBranch("jetAK8IsCSVM")->SetAddress(&data.jetAK8.jetAK8IsCSVM);
    tree_->GetBranch("jetAK8IsCSVT")->SetAddress(&data.jetAK8.jetAK8IsCSVT);
    tree_->GetBranch("jetAK8SmearedE")->SetAddress(&data.jetAK8.jetAK8SmearedE);
    tree_->GetBranch("jetAK8SmearedPt")->SetAddress(&data.jetAK8.jetAK8SmearedPt);
    tree_->GetBranch("jetAK8SmearedPEta")->SetAddress(&data.jetAK8.jetAK8SmearedPEta);
    tree_->GetBranch("jetAK8SmearedPhi")->SetAddress(&data.jetAK8.jetAK8SmearedPhi);
    tree_->GetBranch("jetAK8Y")->SetAddress(&data.jetAK8.jetAK8Y);
    tree_->GetBranch("jetAK8electronMultiplicity")->SetAddress(&data.jetAK8.jetAK8electronMultiplicity);
    tree_->GetBranch("jetAK8muonMultiplicity")->SetAddress(&data.jetAK8.jetAK8muonMultiplicity);
    tree_->GetBranch("jetAK8neutralMultiplicity")->SetAddress(&data.jetAK8.jetAK8neutralMultiplicity);
    tree_->GetBranch("jetAK8photonMultiplicity")->SetAddress(&data.jetAK8.jetAK8photonMultiplicity);
    tree_->GetBranch("jetsAK8CorrPt")->SetAddress(&data.jetAK8.jetsAK8CorrPt);
    tree_->GetBranch("jetsAK8CorrEta")->SetAddress(&data.jetAK8.jetsAK8CorrEta);
    tree_->GetBranch("jetsAK8CorrPhi")->SetAddress(&data.jetAK8.jetsAK8CorrPhi);
    tree_->GetBranch("jetsAK8CorrE")->SetAddress(&data.jetAK8.jetsAK8CorrE);
    tree_->GetBranch("jetsAK8CorrMass")->SetAddress(&data.jetAK8.jetsAK8CorrMass);
    tree_->GetBranch("jetsAK8CorrNJets")->SetAddress(&data.jetAK8.jetsAK8CorrNJets);
    tree_->GetBranch("jetsAK8CorrPartonFlavour")->SetAddress(&data.jetAK8.jetsAK8CorrPartonFlavour);
    
    tree_->GetBranch("met_size")->SetAddress(&data.met.met_size);
    tree_->GetBranch("metPt")->SetAddress(&data.met.metPt);
    tree_->GetBranch("metPhi")->SetAddress(&data.met.metPhi);
    tree_->GetBranch("metPx")->SetAddress(&data.met.metPx);
    tree_->GetBranch("metPy")->SetAddress(&data.met.metPy);
    tree_->GetBranch("metCorrPt")->SetAddress(&data.met.metCorrPt);
    tree_->GetBranch("metCorrPhi")->SetAddress(&data.met.metCorrPhi);
    
    tree_->GetBranch("Event_nTightMuons")->SetAddress(&data.evt.Event_nTightMuons);
    tree_->GetBranch("Event_nLooseMuons")->SetAddress(&data.evt.Event_nLooseMuons);
    tree_->GetBranch("Event_nTightElectrons")->SetAddress(&data.evt.Event_nTightElectrons);
    tree_->GetBranch("Event_nLooseElectrons")->SetAddress(&data.evt.Event_nLooseElectrons);
    tree_->GetBranch("Event_nElectronsSF")->SetAddress(&data.evt.Event_nElectronsSF);
    tree_->GetBranch("Event_nMuonsSF")->SetAddress(&data.evt.Event_nMuonsSF);
    tree_->GetBranch("Event_nCSVTJets")->SetAddress(&data.evt.Event_nCSVTJets);
    tree_->GetBranch("Event_nCSVMJets")->SetAddress(&data.evt.Event_nCSVMJets);
    tree_->GetBranch("Event_nCSVLJets")->SetAddress(&data.evt.Event_nCSVLJets);
    tree_->GetBranch("Event_nTightJets")->SetAddress(&data.evt.Event_nTightJets);
    tree_->GetBranch("Event_nLooseJets")->SetAddress(&data.evt.Event_nLooseJets);
    tree_->GetBranch("Event_bWeight1TCSVT")->SetAddress(&data.evt.Event_bWeight1TCSVT);
    tree_->GetBranch("Event_bWeight1TCSVM")->SetAddress(&data.evt.Event_bWeight1TCSVM);
    tree_->GetBranch("Event_bWeight1TCSVL")->SetAddress(&data.evt.Event_bWeight1TCSVL);
    tree_->GetBranch("Event_bWeight2TCSVT")->SetAddress(&data.evt.Event_bWeight2TCSVT);
    tree_->GetBranch("Event_bWeight2TCSVM")->SetAddress(&data.evt.Event_bWeight2TCSVM);
    tree_->GetBranch("Event_bWeight2TCSVL")->SetAddress(&data.evt.Event_bWeight2TCSVL);
    tree_->GetBranch("Event_LHEWeightSign")->SetAddress(&data.evt.Event_LHEWeightSign);
  }
  
  Long64_t GetEntries() {
    return tree_->GetEntries(); 
  }
  Int_t GetEntry(Long64_t entry, Int_t getall = 0) { 
    return tree_->GetEntry(entry,getall); 
  }
  
private: 
  TTree* tree_;
  
};

#endif
