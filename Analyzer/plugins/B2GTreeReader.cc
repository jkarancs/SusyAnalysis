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
  
  void Load_Tree(TFile &f, const char* name = "B2GTTreeMaker/B2GTree") {
    tree_ = (TTree*)f.Get(name);
    bool debug = false;
    
    set_branch_("gen_size", &data.gen.size, 1);
    set_branch_("gen_Mass", &data.gen.Mass, 1);
    set_branch_("gen_Pt", &data.gen.Pt, 1);
    set_branch_("gen_Eta", &data.gen.Eta, 1);
    set_branch_("gen_Y", &data.gen.Y, 0);
    set_branch_("gen_Phi", &data.gen.Phi, 1);
    set_branch_("gen_E", &data.gen.E, 1);
    set_branch_("gen_Charge", &data.gen.Charge, 0);
    set_branch_("gen_ID", &data.gen.ID, 1);
    set_branch_("gen_Status", &data.gen.Status, 1);
    set_branch_("gen_MomID", &data.gen.MomID, 1);
    if (debug) std::cout<<"B2GTreeReader: gen particles loaded"<<std::endl;
    
    set_branch_("el_size", &data.ele.size, 1);
    set_branch_("el_Mass", &data.ele.Mass, 1);
    set_branch_("el_Pt", &data.ele.Pt, 1);
    set_branch_("el_Eta", &data.ele.Eta, 1);
    set_branch_("el_Y", &data.ele.Y, 0);
    set_branch_("el_Phi", &data.ele.Phi, 1);
    set_branch_("el_E", &data.ele.E, 1);
    set_branch_("el_Charge", &data.ele.Charge, 0);
    set_branch_("el_Iso03", &data.ele.Iso03, 1);
    set_branch_("el_D0", &data.ele.D0, 0);
    set_branch_("el_Dz", &data.ele.Dz, 0);
    set_branch_("el_dEtaIn", &data.ele.dEtaIn, 0);
    set_branch_("el_dPhiIn", &data.ele.dPhiIn, 0);
    set_branch_("el_HoE", &data.ele.HoE, 0);
    set_branch_("el_full5x5siee", &data.ele.full5x5siee, 0);
    set_branch_("el_ooEmooP", &data.ele.ooEmooP, 0);
    set_branch_("el_missHits", &data.ele.missHits, 0);
    set_branch_("el_hasMatchedConVeto", &data.ele.hasMatchedConVeto, 0);
    set_branch_("el_isEB", &data.ele.isEB, 0);
    set_branch_("el_isVeto", &data.ele.isVeto, 1);
    set_branch_("el_isLoose", &data.ele.isLoose, 1);
    set_branch_("el_isTight", &data.ele.isTight, 1);
    set_branch_("el_isMedium", &data.ele.isMedium, 1);
    set_branch_("el_scEta", &data.ele.scEta, 0);
    if (debug) std::cout<<"B2GTreeReader: electrons loaded"<<std::endl;
    
    set_branch_("mu_size", &data.mu.size, 1);
    set_branch_("mu_Mass", &data.mu.Mass, 1);
    set_branch_("mu_Pt", &data.mu.Pt, 1);
    set_branch_("mu_Eta", &data.mu.Eta, 1);
    set_branch_("mu_Y", &data.mu.Y, 0);
    set_branch_("mu_Phi", &data.mu.Phi, 1);
    set_branch_("mu_E", &data.mu.E, 1);
    set_branch_("mu_Charge", &data.mu.Charge, 0);
    set_branch_("mu_Iso04", &data.mu.Iso04, 1);
    set_branch_("mu_D0", &data.mu.D0, 0);
    set_branch_("mu_D0err", &data.mu.D0err, 0);
    set_branch_("mu_Dxy", &data.mu.Dxy, 0);
    set_branch_("mu_Dxyerr", &data.mu.Dxyerr, 0);
    set_branch_("mu_Dz", &data.mu.Dz, 0);
    set_branch_("mu_Dzerr", &data.mu.Dzerr, 0);
    set_branch_("mu_IsLooseMuon", &data.mu.IsLooseMuon, 1);
    set_branch_("mu_IsSoftMuon", &data.mu.IsSoftMuon, 1);
    set_branch_("mu_IsTightMuon", &data.mu.IsTightMuon, 1);
    set_branch_("mu_IsPFMuon", &data.mu.IsPFMuon, 0);
    set_branch_("mu_IsGlobalMuon", &data.mu.IsGlobalMuon, 0);
    set_branch_("mu_IsTrackerMuon", &data.mu.IsTrackerMuon, 0);
    set_branch_("mu_GlbTrkNormChi2", &data.mu.GlbTrkNormChi2, 0);
    set_branch_("mu_NumberValidMuonHits", &data.mu.NumberValidMuonHits, 0);
    set_branch_("mu_NumberMatchedStations", &data.mu.NumberMatchedStations, 0);
    set_branch_("mu_NumberValidPixelHits", &data.mu.NumberValidPixelHits, 0);
    set_branch_("mu_NumberTrackerLayers", &data.mu.NumberTrackerLayers, 0);
    set_branch_("mu_NumberOfValidTrackerHits", &data.mu.NumberOfValidTrackerHits, 0);
    set_branch_("mu_NumberOfPixelLayers", &data.mu.NumberOfPixelLayers, 0);
    set_branch_("mu_InTrkNormChi2", &data.mu.InTrkNormChi2, 0);
    set_branch_("mu_SumChargedHadronPt", &data.mu.SumChargedHadronPt, 0);
    set_branch_("mu_SumNeutralHadronPt", &data.mu.SumNeutralHadronPt, 0);
    set_branch_("mu_SumPhotonPt", &data.mu.SumPhotonPt, 0);
    set_branch_("mu_SumPUPt", &data.mu.SumPUPt, 0);
    set_branch_("mu_GenMuonY", &data.mu.GenMuonY, 0);
    set_branch_("mu_GenMuonEta", &data.mu.GenMuonEta, 0);
    set_branch_("mu_GenMuonPhi", &data.mu.GenMuonPhi, 0);
    set_branch_("mu_GenMuonPt", &data.mu.GenMuonPt, 0);
    set_branch_("mu_GenMuonE", &data.mu.GenMuonE, 0);
    set_branch_("mu_GenMuonCharge", &data.mu.GenMuonCharge, 0);
    if (debug) std::cout<<"B2GTreeReader: muons loaded"<<std::endl;
    
    set_branch_("met_Pt", &data.met.Pt, 1);
    set_branch_("met_Phi", &data.met.Phi, 1);
    set_branch_("met_Px", &data.met.Px, 1);
    set_branch_("met_Py", &data.met.Py, 1);
    if (debug) std::cout<<"B2GTreeReader: met loaded"<<std::endl;
    
    set_branch_("jetAK4_size", &data.jetsAK4.size, 1);
    set_branch_("jetAK4_Mass", &data.jetsAK4.Mass, 1);			
    set_branch_("jetAK4_Pt", &data.jetsAK4.Pt, 1);
    set_branch_("jetAK4_Eta", &data.jetsAK4.Eta, 1);
    set_branch_("jetAK4_Y", &data.jetsAK4.Y, 0);
    set_branch_("jetAK4_Phi", &data.jetsAK4.Phi, 1);
    set_branch_("jetAK4_E", &data.jetsAK4.E, 1);
    set_branch_("jetAK4_Charge", &data.jetsAK4.Charge, 0);
    set_branch_("jetAK4_CSV", &data.jetsAK4.CSV, 0);
    set_branch_("jetAK4_CSVV1", &data.jetsAK4.CSVV1, 0);
    set_branch_("jetAK4_GenPartonY", &data.jetsAK4.GenPartonY, 0);
    set_branch_("jetAK4_GenPartonEta", &data.jetsAK4.GenPartonEta, 0);
    set_branch_("jetAK4_GenPartonPhi", &data.jetsAK4.GenPartonPhi, 0);
    set_branch_("jetAK4_GenPartonPt", &data.jetsAK4.GenPartonPt, 0);
    set_branch_("jetAK4_GenPartonE", &data.jetsAK4.GenPartonE, 0);
    set_branch_("jetAK4_GenPartonCharge", &data.jetsAK4.GenPartonCharge, 0);
    set_branch_("jetAK4_PartonFlavour", &data.jetsAK4.PartonFlavour, 0);
    set_branch_("jetAK4_HadronFlavour", &data.jetsAK4.HadronFlavour, 0);
    set_branch_("jetAK4_GenJetY", &data.jetsAK4.GenJetY, 0);
    set_branch_("jetAK4_GenJetEta", &data.jetsAK4.GenJetEta, 0);
    set_branch_("jetAK4_GenJetPhi", &data.jetsAK4.GenJetPhi, 0);
    set_branch_("jetAK4_GenJetPt", &data.jetsAK4.GenJetPt, 0);
    set_branch_("jetAK4_GenJetE", &data.jetsAK4.GenJetE, 0);
    set_branch_("jetAK4_GenJetCharge", &data.jetsAK4.GenJetCharge, 0);
    set_branch_("jetAK4_muonMultiplicity", &data.jetsAK4.muonMultiplicity, 0);
    set_branch_("jetAK4_PhotonEnergy", &data.jetsAK4.PhotonEnergy, 0);
    set_branch_("jetAK4_ElectronEnergy", &data.jetsAK4.ElectronEnergy, 0);
    set_branch_("jetAK4_MuonEnergy", &data.jetsAK4.MuonEnergy, 0);
    set_branch_("jetAK4_HFHadronEnergy", &data.jetsAK4.HFHadronEnergy, 0);
    set_branch_("jetAK4_HFEMEnergy", &data.jetsAK4.HFEMEnergy, 0);
    set_branch_("jetAK4_ChargedHadronMultiplicity", &data.jetsAK4.ChargedHadronMultiplicity, 0);
    set_branch_("jetAK4_numberOfDaughters", &data.jetsAK4.numberOfDaughters, 0);
    set_branch_("jetAK4_chargedMultiplicity", &data.jetsAK4.chargedMultiplicity, 0);
    set_branch_("jetAK4_neutralHadronMultiplicity", &data.jetsAK4.neutralHadronMultiplicity, 0);
    set_branch_("jetAK4_neutralHadronEnergy", &data.jetsAK4.neutralHadronEnergy, 0);
    set_branch_("jetAK4_neutralEmEnergy", &data.jetsAK4.neutralEmEnergy, 0);
    set_branch_("jetAK4_chargedEmEnergy", &data.jetsAK4.chargedEmEnergy, 0);
    set_branch_("jetAK4_chargedHadronEnergy", &data.jetsAK4.chargedHadronEnergy, 0);
    set_branch_("jetAK4_photonMultiplicity", &data.jetsAK4.photonMultiplicity, 0);
    set_branch_("jetAK4_electronMultiplicity", &data.jetsAK4.electronMultiplicity, 0);
    set_branch_("jetAK4_HFHadronMultiplicity", &data.jetsAK4.HFHadronMultiplicity, 0);
    set_branch_("jetAK4_HFEMMultiplicity", &data.jetsAK4.HFEMMultiplicity, 0);
    set_branch_("jetAK4_ChargeMuEnergy", &data.jetsAK4.ChargeMuEnergy, 0);
    set_branch_("jetAK4_neutralMultiplicity", &data.jetsAK4.neutralMultiplicity, 0);
    set_branch_("jetAK4_jecFactor0", &data.jetsAK4.jecFactor0, 0);
    set_branch_("jetAK4_jetArea", &data.jetsAK4.jetArea, 0);
    set_branch_("jetAK4_SmearedPt", &data.jetsAK4.SmearedPt, 0);
    set_branch_("jetAK4_SmearedPEta", &data.jetsAK4.SmearedPEta, 0);
    set_branch_("jetAK4_SmearedPhi", &data.jetsAK4.SmearedPhi, 0);
    set_branch_("jetAK4_SmearedE", &data.jetsAK4.SmearedE, 0);
    set_branch_("jetAK4_JERup", &data.jetsAK4.JERup, 0);
    set_branch_("jetAK4_JERdown", &data.jetsAK4.JERdown, 0);
    if (debug) std::cout<<"B2GTreeReader: AK4 jets loaded"<<std::endl;
    
    set_branch_("jetAK8_size", &data.jetsAK8.size, 1);
    set_branch_("jetAK8_Mass", &data.jetsAK8.Mass, 1);			
    set_branch_("jetAK8_Pt", &data.jetsAK8.Pt, 1);
    set_branch_("jetAK8_Eta", &data.jetsAK8.Eta, 1);
    set_branch_("jetAK8_Y", &data.jetsAK8.Y, 0);
    set_branch_("jetAK8_Phi", &data.jetsAK8.Phi, 1);
    set_branch_("jetAK8_E", &data.jetsAK8.E, 1);
    set_branch_("jetAK8_Charge", &data.jetsAK8.Charge, 0);
    set_branch_("jetAK8_CSV", &data.jetsAK8.CSV, 1);
    set_branch_("jetAK8_CSVV1", &data.jetsAK8.CSVV1, 1);
    set_branch_("jetAK8_GenPartonY", &data.jetsAK8.GenPartonY, 0);
    set_branch_("jetAK8_GenPartonEta", &data.jetsAK8.GenPartonEta, 0);
    set_branch_("jetAK8_GenPartonPhi", &data.jetsAK8.GenPartonPhi, 0);
    set_branch_("jetAK8_GenPartonPt", &data.jetsAK8.GenPartonPt, 0);
    set_branch_("jetAK8_GenPartonE", &data.jetsAK8.GenPartonE, 0);
    set_branch_("jetAK8_GenPartonCharge", &data.jetsAK8.GenPartonCharge, 0);
    set_branch_("jetAK8_PartonFlavour", &data.jetsAK8.PartonFlavour, 0);
    set_branch_("jetAK8_HadronFlavour", &data.jetsAK8.HadronFlavour, 0);
    set_branch_("jetAK8_GenJetY", &data.jetsAK8.GenJetY, 0);
    set_branch_("jetAK8_GenJetEta", &data.jetsAK8.GenJetEta, 0);
    set_branch_("jetAK8_GenJetPhi", &data.jetsAK8.GenJetPhi, 0);
    set_branch_("jetAK8_GenJetPt", &data.jetsAK8.GenJetPt, 0);
    set_branch_("jetAK8_GenJetE", &data.jetsAK8.GenJetE, 0);
    set_branch_("jetAK8_GenJetCharge", &data.jetsAK8.GenJetCharge, 0);
    set_branch_("jetAK8_muonMultiplicity", &data.jetsAK8.muonMultiplicity, 0);
    set_branch_("jetAK8_PhotonEnergy", &data.jetsAK8.PhotonEnergy, 0);
    set_branch_("jetAK8_ElectronEnergy", &data.jetsAK8.ElectronEnergy, 0);
    set_branch_("jetAK8_MuonEnergy", &data.jetsAK8.MuonEnergy, 0);
    set_branch_("jetAK8_HFHadronEnergy", &data.jetsAK8.HFHadronEnergy, 0);
    set_branch_("jetAK8_HFEMEnergy", &data.jetsAK8.HFEMEnergy, 0);
    set_branch_("jetAK8_ChargedHadronMultiplicity", &data.jetsAK8.ChargedHadronMultiplicity, 0);
    set_branch_("jetAK8_numberOfDaughters", &data.jetsAK8.numberOfDaughters, 0);
    set_branch_("jetAK8_chargedMultiplicity", &data.jetsAK8.chargedMultiplicity, 0);
    set_branch_("jetAK8_neutralHadronMultiplicity", &data.jetsAK8.neutralHadronMultiplicity, 0);
    set_branch_("jetAK8_neutralHadronEnergy", &data.jetsAK8.neutralHadronEnergy, 0);
    set_branch_("jetAK8_neutralEmEnergy", &data.jetsAK8.neutralEmEnergy, 0);
    set_branch_("jetAK8_chargedEmEnergy", &data.jetsAK8.chargedEmEnergy, 0);
    set_branch_("jetAK8_chargedHadronEnergy", &data.jetsAK8.chargedHadronEnergy, 0);
    set_branch_("jetAK8_photonMultiplicity", &data.jetsAK8.photonMultiplicity, 0);
    set_branch_("jetAK8_electronMultiplicity", &data.jetsAK8.electronMultiplicity, 0);
    set_branch_("jetAK8_HFHadronMultiplicity", &data.jetsAK8.HFHadronMultiplicity, 0);
    set_branch_("jetAK8_HFEMMultiplicity", &data.jetsAK8.HFEMMultiplicity, 0);
    set_branch_("jetAK8_ChargeMuEnergy", &data.jetsAK8.ChargeMuEnergy, 0);
    set_branch_("jetAK8_neutralMultiplicity", &data.jetsAK8.neutralMultiplicity, 0);
    set_branch_("jetAK8_jecFactor0", &data.jetsAK8.jecFactor0, 0);
    set_branch_("jetAK8_jetArea", &data.jetsAK8.jetArea, 0);
    set_branch_("jetAK8_SmearedPt", &data.jetsAK8.SmearedPt, 0);
    set_branch_("jetAK8_SmearedPEta", &data.jetsAK8.SmearedPEta, 0);
    set_branch_("jetAK8_SmearedPhi", &data.jetsAK8.SmearedPhi, 0);
    set_branch_("jetAK8_SmearedE", &data.jetsAK8.SmearedE, 0);
    set_branch_("jetAK8_JERup", &data.jetsAK8.JERup, 0);
    set_branch_("jetAK8_JERdown", &data.jetsAK8.JERdown, 0);
    set_branch_("jetAK8_vSubjetIndex0", &data.jetsAK8.vSubjetIndex0, 1);
    set_branch_("jetAK8_vSubjetIndex1", &data.jetsAK8.vSubjetIndex1, 1);
    set_branch_("jetAK8_topSubjetIndex0", &data.jetsAK8.topSubjetIndex0, 1);
    set_branch_("jetAK8_topSubjetIndex1", &data.jetsAK8.topSubjetIndex1, 1);
    set_branch_("jetAK8_topSubjetIndex2", &data.jetsAK8.topSubjetIndex2, 1);
    set_branch_("jetAK8_topSubjetIndex3", &data.jetsAK8.topSubjetIndex3, 1);
    set_branch_("jetAK8_tau1", &data.jetsAK8.tau1, 1);
    set_branch_("jetAK8_tau2", &data.jetsAK8.tau2, 1);
    set_branch_("jetAK8_tau3", &data.jetsAK8.tau3, 1);
    set_branch_("jetAK8_trimmedMass", &data.jetsAK8.trimmedMass, 1);
    set_branch_("jetAK8_prunedMass", &data.jetsAK8.prunedMass, 1);
    set_branch_("jetAK8_filteredMass", &data.jetsAK8.filteredMass, 1);
    set_branch_("jetAK8_topMass", &data.jetsAK8.topMass, 1);
    set_branch_("jetAK8_wMass", &data.jetsAK8.wMass, 1);
    set_branch_("jetAK8_nSubJets", &data.jetsAK8.nSubJets, 1);
    set_branch_("jetAK8_minmass", &data.jetsAK8.minmass, 1);
    if (debug) std::cout<<"B2GTreeReader: AK8 jets loaded"<<std::endl;
    
    set_branch_("subjetAK8_size", &data.subjetsAK8.size, 1);
    set_branch_("subjetAK8_Mass", &data.subjetsAK8.Mass, 1);			
    set_branch_("subjetAK8_Pt", &data.subjetsAK8.Pt, 1);
    set_branch_("subjetAK8_Eta", &data.subjetsAK8.Eta, 1);
    set_branch_("subjetAK8_Y", &data.subjetsAK8.Y, 0);
    set_branch_("subjetAK8_Phi", &data.subjetsAK8.Phi, 1);
    set_branch_("subjetAK8_E", &data.subjetsAK8.E, 1);
    set_branch_("subjetAK8_Charge", &data.subjetsAK8.Charge, 0);
    set_branch_("subjetAK8_CSV", &data.subjetsAK8.CSV, 1);
    set_branch_("subjetAK8_CSVV1", &data.subjetsAK8.CSVV1, 1);
    set_branch_("subjetAK8_GenPartonY", &data.subjetsAK8.GenPartonY, 0);
    set_branch_("subjetAK8_GenPartonEta", &data.subjetsAK8.GenPartonEta, 0);
    set_branch_("subjetAK8_GenPartonPhi", &data.subjetsAK8.GenPartonPhi, 0);
    set_branch_("subjetAK8_GenPartonPt", &data.subjetsAK8.GenPartonPt, 0);
    set_branch_("subjetAK8_GenPartonE", &data.subjetsAK8.GenPartonE, 0);
    set_branch_("subjetAK8_GenPartonCharge", &data.subjetsAK8.GenPartonCharge, 0);
    set_branch_("subjetAK8_PartonFlavour", &data.subjetsAK8.PartonFlavour, 0);
    set_branch_("subjetAK8_HadronFlavour", &data.subjetsAK8.HadronFlavour, 0);
    set_branch_("subjetAK8_GenJetY", &data.subjetsAK8.GenJetY, 0);
    set_branch_("subjetAK8_GenJetEta", &data.subjetsAK8.GenJetEta, 0);
    set_branch_("subjetAK8_GenJetPhi", &data.subjetsAK8.GenJetPhi, 0);
    set_branch_("subjetAK8_GenJetPt", &data.subjetsAK8.GenJetPt, 0);
    set_branch_("subjetAK8_GenJetE", &data.subjetsAK8.GenJetE, 0);
    set_branch_("subjetAK8_GenJetCharge", &data.subjetsAK8.GenJetCharge, 0);
    set_branch_("subjetAK8_muonMultiplicity", &data.subjetsAK8.muonMultiplicity, 0);
    set_branch_("subjetAK8_PhotonEnergy", &data.subjetsAK8.PhotonEnergy, 0);
    set_branch_("subjetAK8_ElectronEnergy", &data.subjetsAK8.ElectronEnergy, 0);
    set_branch_("subjetAK8_MuonEnergy", &data.subjetsAK8.MuonEnergy, 0);
    set_branch_("subjetAK8_HFHadronEnergy", &data.subjetsAK8.HFHadronEnergy, 0);
    set_branch_("subjetAK8_HFEMEnergy", &data.subjetsAK8.HFEMEnergy, 0);
    set_branch_("subjetAK8_ChargedHadronMultiplicity", &data.subjetsAK8.ChargedHadronMultiplicity, 0);
    set_branch_("subjetAK8_numberOfDaughters", &data.subjetsAK8.numberOfDaughters, 0);
    set_branch_("subjetAK8_chargedMultiplicity", &data.subjetsAK8.chargedMultiplicity, 0);
    set_branch_("subjetAK8_neutralHadronMultiplicity", &data.subjetsAK8.neutralHadronMultiplicity, 0);
    set_branch_("subjetAK8_neutralHadronEnergy", &data.subjetsAK8.neutralHadronEnergy, 0);
    set_branch_("subjetAK8_neutralEmEnergy", &data.subjetsAK8.neutralEmEnergy, 0);
    set_branch_("subjetAK8_chargedEmEnergy", &data.subjetsAK8.chargedEmEnergy, 0);
    set_branch_("subjetAK8_chargedHadronEnergy", &data.subjetsAK8.chargedHadronEnergy, 0);
    set_branch_("subjetAK8_photonMultiplicity", &data.subjetsAK8.photonMultiplicity, 0);
    set_branch_("subjetAK8_electronMultiplicity", &data.subjetsAK8.electronMultiplicity, 0);
    set_branch_("subjetAK8_HFHadronMultiplicity", &data.subjetsAK8.HFHadronMultiplicity, 0);
    set_branch_("subjetAK8_HFEMMultiplicity", &data.subjetsAK8.HFEMMultiplicity, 0);
    set_branch_("subjetAK8_ChargeMuEnergy", &data.subjetsAK8.ChargeMuEnergy, 0);
    set_branch_("subjetAK8_neutralMultiplicity", &data.subjetsAK8.neutralMultiplicity, 0);
    set_branch_("subjetAK8_jecFactor0", &data.subjetsAK8.jecFactor0, 0);
    set_branch_("subjetAK8_jetArea", &data.subjetsAK8.jetArea, 0);
    set_branch_("subjetAK8_SmearedPt", &data.subjetsAK8.SmearedPt, 0);
    set_branch_("subjetAK8_SmearedPEta", &data.subjetsAK8.SmearedPEta, 0);
    set_branch_("subjetAK8_SmearedPhi", &data.subjetsAK8.SmearedPhi, 0);
    set_branch_("subjetAK8_SmearedE", &data.subjetsAK8.SmearedE, 0);
    set_branch_("subjetAK8_JERup", &data.subjetsAK8.JERup, 0);
    set_branch_("subjetAK8_JERdown", &data.subjetsAK8.JERdown, 0);
    set_branch_("subjetAK8_subjetCSV", &data.subjetsAK8.subjetCSV, 1);
    if (debug) std::cout<<"B2GTreeReader: AK8 subjets loaded"<<std::endl;
    
    set_branch_("subjetsCmsTopTag_size", &data.subjetsCmsTopTag.size, 0);
    set_branch_("subjetsCmsTopTag_Mass", &data.subjetsCmsTopTag.Mass, 0);			
    set_branch_("subjetsCmsTopTag_Pt", &data.subjetsCmsTopTag.Pt, 0);
    set_branch_("subjetsCmsTopTag_Eta", &data.subjetsCmsTopTag.Eta, 0);
    set_branch_("subjetsCmsTopTag_Y", &data.subjetsCmsTopTag.Y, 0);
    set_branch_("subjetsCmsTopTag_Phi", &data.subjetsCmsTopTag.Phi, 0);
    set_branch_("subjetsCmsTopTag_E", &data.subjetsCmsTopTag.E, 0);
    set_branch_("subjetsCmsTopTag_Charge", &data.subjetsCmsTopTag.Charge, 0);
    set_branch_("subjetsCmsTopTag_CSV", &data.subjetsCmsTopTag.CSV, 0);
    set_branch_("subjetsCmsTopTag_CSVV1", &data.subjetsCmsTopTag.CSVV1, 0);
    set_branch_("subjetsCmsTopTag_GenPartonY", &data.subjetsCmsTopTag.GenPartonY, 0);
    set_branch_("subjetsCmsTopTag_GenPartonEta", &data.subjetsCmsTopTag.GenPartonEta, 0);
    set_branch_("subjetsCmsTopTag_GenPartonPhi", &data.subjetsCmsTopTag.GenPartonPhi, 0);
    set_branch_("subjetsCmsTopTag_GenPartonPt", &data.subjetsCmsTopTag.GenPartonPt, 0);
    set_branch_("subjetsCmsTopTag_GenPartonE", &data.subjetsCmsTopTag.GenPartonE, 0);
    set_branch_("subjetsCmsTopTag_GenPartonCharge", &data.subjetsCmsTopTag.GenPartonCharge, 0);
    set_branch_("subjetsCmsTopTag_PartonFlavour", &data.subjetsCmsTopTag.PartonFlavour, 0);
    set_branch_("subjetsCmsTopTag_HadronFlavour", &data.subjetsCmsTopTag.HadronFlavour, 0);
    set_branch_("subjetsCmsTopTag_GenJetY", &data.subjetsCmsTopTag.GenJetY, 0);
    set_branch_("subjetsCmsTopTag_GenJetEta", &data.subjetsCmsTopTag.GenJetEta, 0);
    set_branch_("subjetsCmsTopTag_GenJetPhi", &data.subjetsCmsTopTag.GenJetPhi, 0);
    set_branch_("subjetsCmsTopTag_GenJetPt", &data.subjetsCmsTopTag.GenJetPt, 0);
    set_branch_("subjetsCmsTopTag_GenJetE", &data.subjetsCmsTopTag.GenJetE, 0);
    set_branch_("subjetsCmsTopTag_GenJetCharge", &data.subjetsCmsTopTag.GenJetCharge, 0);
    set_branch_("subjetsCmsTopTag_muonMultiplicity", &data.subjetsCmsTopTag.muonMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_PhotonEnergy", &data.subjetsCmsTopTag.PhotonEnergy, 0);
    set_branch_("subjetsCmsTopTag_ElectronEnergy", &data.subjetsCmsTopTag.ElectronEnergy, 0);
    set_branch_("subjetsCmsTopTag_MuonEnergy", &data.subjetsCmsTopTag.MuonEnergy, 0);
    set_branch_("subjetsCmsTopTag_HFHadronEnergy", &data.subjetsCmsTopTag.HFHadronEnergy, 0);
    set_branch_("subjetsCmsTopTag_HFEMEnergy", &data.subjetsCmsTopTag.HFEMEnergy, 0);
    set_branch_("subjetsCmsTopTag_ChargedHadronMultiplicity", &data.subjetsCmsTopTag.ChargedHadronMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_numberOfDaughters", &data.subjetsCmsTopTag.numberOfDaughters, 0);
    set_branch_("subjetsCmsTopTag_chargedMultiplicity", &data.subjetsCmsTopTag.chargedMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_neutralHadronMultiplicity", &data.subjetsCmsTopTag.neutralHadronMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_neutralHadronEnergy", &data.subjetsCmsTopTag.neutralHadronEnergy, 0);
    set_branch_("subjetsCmsTopTag_neutralEmEnergy", &data.subjetsCmsTopTag.neutralEmEnergy, 0);
    set_branch_("subjetsCmsTopTag_chargedEmEnergy", &data.subjetsCmsTopTag.chargedEmEnergy, 0);
    set_branch_("subjetsCmsTopTag_chargedHadronEnergy", &data.subjetsCmsTopTag.chargedHadronEnergy, 0);
    set_branch_("subjetsCmsTopTag_photonMultiplicity", &data.subjetsCmsTopTag.photonMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_electronMultiplicity", &data.subjetsCmsTopTag.electronMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_HFHadronMultiplicity", &data.subjetsCmsTopTag.HFHadronMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_HFEMMultiplicity", &data.subjetsCmsTopTag.HFEMMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_ChargeMuEnergy", &data.subjetsCmsTopTag.ChargeMuEnergy, 0);
    set_branch_("subjetsCmsTopTag_neutralMultiplicity", &data.subjetsCmsTopTag.neutralMultiplicity, 0);
    set_branch_("subjetsCmsTopTag_jecFactor0", &data.subjetsCmsTopTag.jecFactor0, 0);
    set_branch_("subjetsCmsTopTag_jetArea", &data.subjetsCmsTopTag.jetArea, 0);
    set_branch_("subjetsCmsTopTag_SmearedPt", &data.subjetsCmsTopTag.SmearedPt, 0);
    set_branch_("subjetsCmsTopTag_SmearedPEta", &data.subjetsCmsTopTag.SmearedPEta, 0);
    set_branch_("subjetsCmsTopTag_SmearedPhi", &data.subjetsCmsTopTag.SmearedPhi, 0);
    set_branch_("subjetsCmsTopTag_SmearedE", &data.subjetsCmsTopTag.SmearedE, 0);
    set_branch_("subjetsCmsTopTag_JERup", &data.subjetsCmsTopTag.JERup, 0);
    set_branch_("subjetsCmsTopTag_JERdown", &data.subjetsCmsTopTag.JERdown, 0);
    if (debug) std::cout<<"B2GTreeReader: CMS top-tag subjets loaded"<<std::endl;
    
    set_branch_("evt_NLep", &data.evt.NLep, 1);
    set_branch_("evt_NTopHad", &data.evt.NTopHad, 1);
    set_branch_("evt_NTopLep", &data.evt.NTopLep, 1);
    set_branch_("evt_NTop", &data.evt.NTop, 1);
    set_branch_("evt_HtLep", &data.evt.HtLep, 1);
    set_branch_("evt_HtTop", &data.evt.HtTop, 1);
    set_branch_("evt_Ht", &data.evt.Ht, 1);
    set_branch_("evt_HtAll", &data.evt.HtAll, 1);
    set_branch_("evt_HtEx", &data.evt.HtEx, 1);
    set_branch_("evt_HtExFr", &data.evt.HtExFr, 1);
    set_branch_("evt_HtTopFr", &data.evt.HtTopFr, 1);
    set_branch_("evt_TTHadDR", &data.evt.TTHadDR, 1);
    set_branch_("evt_TTHadDPhi", &data.evt.TTHadDPhi, 1);
    set_branch_("evt_TTHadDEta", &data.evt.TTHadDEta, 1);
    set_branch_("evt_TTHadMass", &data.evt.TTHadMass, 1);
    set_branch_("evt_TTHadPz", &data.evt.TTHadPz, 1);
    set_branch_("evt_TTHadHz", &data.evt.TTHadHz, 1);
    set_branch_("evt_TTHadDPz", &data.evt.TTHadDPz, 1);
    set_branch_("evt_TTHadMR", &data.evt.TTHadMR, 1);
    set_branch_("evt_TTHadMTR", &data.evt.TTHadMTR, 1);
    set_branch_("evt_TTHadR", &data.evt.TTHadR, 1);
    set_branch_("evt_TTHadR2", &data.evt.TTHadR2, 1);
    set_branch_("evt_MR", &data.evt.MR, 1);
    set_branch_("evt_MTR", &data.evt.MTR, 1);
    set_branch_("evt_R", &data.evt.R, 1);
    set_branch_("evt_R2", &data.evt.R2, 1);
    
    if (debug) {
      std::cout<<"B2GTreeReader: Event variables loaded"<<std::endl;
      std::cout<<"B2GTreeReader: All objects loaded successfully"<<std::endl;
    }
  }
  
  Long64_t GetEntries() {
    return tree_->GetEntries(); 
  }
  Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    Int_t result = tree_->GetEntry(entry,getall);
    if (data.gen.size>500) std::cout<<"B2GReader Warning!!: Size of Data object too small for gen."<<std::endl;
    if (data.ele.size>NLEP) std::cout<<"B2GReader Warning!!: Size of Data object too small for ele. Set NLEP>="<<data.ele.size<<std::endl;
    if (data.mu.size>NLEP) std::cout<<"B2GReader Warning!!: Size of Data object too small for mu. Set NLEP>="<<data.mu.size<<std::endl;
    if (data.jetsAK4.size>NJET) std::cout<<"B2GReader Warning!!: Size of Data object too small for jetsAK4. Set NJET>="<<data.jetsAK4.size<<std::endl;
    if (data.jetsAK8.size>NJET) std::cout<<"B2GReader Warning!!: Size of Data object too small for jetsAk8. Set NJET>="<<data.jetsAK8.size<<std::endl;
    if (data.subjetsAK8.size>NJET) std::cout<<"B2GReader Warning!!: Size of Data object too small for subjetsAK8. Set NJET>="<<data.subjetsAK8.size<<std::endl;
    if (data.subjetsCmsTopTag.size>NJET) std::cout<<"B2GReader Warning!!: Size of Data object too small for subjetsCmsTopTag. Set NJET>="<<data.subjetsCmsTopTag.size<<std::endl;
    return result;
  }
  
private: 
  void set_branch_(std::string name, void* address, bool status) {
    //tree_->GetBranch(name.c_str())->SetAddress(address);
    //tree_->SetBranchStatus(name.c_str(), status);
    if (status) tree_->GetBranch(name.c_str())->SetAddress(address);
  }
  
  TTree* tree_;
  
};

#endif
