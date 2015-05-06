#ifndef Data_h
#define Data_h

#define NOVAL_I -9999
#define NOVAL_F -9999.0

#define NLEP 20
#define NJET 50
#define NGEN 500

#include <cassert>
#include <map>
#include <vector>
#include <iostream>

#include "TLorentzVector.h"

using namespace std;

class Data {
public:
  Data() {}
  ~Data() {}
  
  class GenVars {
  public:
    GenVars() { init(); }
    ~GenVars() {}
    
    // Basic
    float Mass[NGEN+1];
    float Pt[NGEN+1];
    float Eta[NGEN+1];
    float Y[NGEN+1];
    float Phi[NGEN+1];
    float E[NGEN+1];
    float Charge[NGEN+1];
    float ID[NGEN+1];
    float Status[NGEN+1];
    float MomID[NGEN+1];
    
    size_t it;
    int size;
    
    void init() {
      for (size_t i=0; i<=NGEN; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        ID[i]=NOVAL_F;
	Status[i]=NOVAL_F;
        MomID[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    void SetCurrent(size_t it) {
      Mass[NGEN]              = Mass[it];             
      Pt[NGEN]		      = Pt[it];		     
      Eta[NGEN]		      = Eta[it];		     
      Y[NGEN]		      = Y[it];		     
      Phi[NGEN]		      = Phi[it];		     
      E[NGEN]		      = E[it];		     
      Charge[NGEN]	      = Charge[it];	     
      ID[NGEN]		      = ID[it];		     
      Status[NGEN]	      = Status[it];		     
      MomID[NGEN]	      = MomID[it];		     
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } gen;
  
  class ElectronVars {
  public:
    ElectronVars() { init(); }
    ~ElectronVars() {}
    
    // Basic
    float Mass[NLEP+1];
    float Pt[NLEP+1];
    float Eta[NLEP+1];
    float Y[NLEP+1];
    float Phi[NLEP+1];
    float E[NLEP+1];
    float Charge[NLEP+1];
    // ElectronVars
    float Iso03[NLEP+1];
    float D0[NLEP+1];
    float Dz[NLEP+1];
    float dEtaIn[NLEP+1];
    float dPhiIn[NLEP+1];
    float HoE[NLEP+1];
    float full5x5siee[NLEP+1];
    float ooEmooP[NLEP+1];
    float missHits[NLEP+1];
    float hasMatchedConVeto[NLEP+1];
    float isEB[NLEP+1];
    float isVeto[NLEP+1];
    float isLoose[NLEP+1];
    float isTight[NLEP+1];
    float isMedium[NLEP+1];
    float scEta[NLEP+1];
    
    size_t it;
    int size;

    
    void init() {
      for (size_t i=0; i<=NLEP; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        Iso03[i]=NOVAL_F;
        D0[i]=NOVAL_F;
        Dz[i]=NOVAL_F;
        dEtaIn[i]=NOVAL_F;
        dPhiIn[i]=NOVAL_F;
        HoE[i]=NOVAL_F;
        full5x5siee[i]=NOVAL_F;
        ooEmooP[i]=NOVAL_F;
        missHits[i]=NOVAL_F;
        hasMatchedConVeto[i]=NOVAL_F;
        isEB[i]=NOVAL_F;
        isVeto[i]=NOVAL_F;
        isLoose[i]=NOVAL_F;
        isTight[i]=NOVAL_F;
        isMedium[i]=NOVAL_F;
        scEta[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    void SetCurrent(size_t it) {
      Mass[NLEP]              = Mass[it];             
      Pt[NLEP]		      = Pt[it];		     
      Eta[NLEP]		      = Eta[it];		     
      Y[NLEP]		      = Y[it];		     
      Phi[NLEP]		      = Phi[it];		     
      E[NLEP]		      = E[it];		     
      Charge[NLEP]	      = Charge[it];	     
      Iso03[NLEP]	      = Iso03[it];	     
      D0[NLEP]		      = D0[it];		     
      Dz[NLEP]		      = Dz[it];		     
      dEtaIn[NLEP]	      = dEtaIn[it];	     
      dPhiIn[NLEP]	      = dPhiIn[it];	     
      HoE[NLEP]		      = HoE[it];		     
      full5x5siee[NLEP]	      = full5x5siee[it];	     
      ooEmooP[NLEP]	      = ooEmooP[it];	     
      missHits[NLEP]	      = missHits[it];	     
      hasMatchedConVeto[NLEP] = hasMatchedConVeto[it];
      isEB[NLEP]	      = isEB[it];	     
      isVeto[NLEP]	      = isVeto[it];	     
      isLoose[NLEP]	      = isLoose[it];	     
      isTight[NLEP]	      = isTight[it];	     
      isMedium[NLEP]          = isMedium[it];         
      scEta[NLEP]	      = scEta[it];		     
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } ele;

  class MuonVars {
  public:
    MuonVars() { init(); }
    ~MuonVars() {}
    
    // Basic
    float Mass[NLEP+1];
    float Pt[NLEP+1];
    float Eta[NLEP+1];
    float Y[NLEP+1];
    float Phi[NLEP+1];
    float E[NLEP+1];
    float Charge[NLEP+1];
    // MuonVars
    float Iso04[NLEP+1];
    float D0[NLEP+1];
    float D0err[NLEP+1];
    float Dxy[NLEP+1];
    float Dxyerr[NLEP+1];
    float Dz[NLEP+1];
    float Dzerr[NLEP+1];
    float IsLooseMuon[NLEP+1];
    float IsSoftMuon[NLEP+1];
    float IsTightMuon[NLEP+1];
    float IsPFMuon[NLEP+1];
    float IsGlobalMuon[NLEP+1];
    float IsTrackerMuon[NLEP+1];
    float GlbTrkNormChi2[NLEP+1];
    float NumberValidMuonHits[NLEP+1];
    float NumberMatchedStations[NLEP+1];
    float NumberValidPixelHits[NLEP+1];
    float NumberTrackerLayers[NLEP+1];
    float NumberOfValidTrackerHits[NLEP+1];
    float NumberOfPixelLayers[NLEP+1];
    float InTrkNormChi2[NLEP+1];
    float SumChargedHadronPt[NLEP+1];
    float SumNeutralHadronPt[NLEP+1];
    float SumPhotonPt[NLEP+1];
    float SumPUPt[NLEP+1];
    float GenMuonY[NLEP+1];
    float GenMuonEta[NLEP+1];
    float GenMuonPhi[NLEP+1];
    float GenMuonPt[NLEP+1];
    float GenMuonE[NLEP+1];
    float GenMuonCharge[NLEP+1];
    
    size_t it;
    int size;
    
    void init() {
      for (size_t i=0; i<=NLEP; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        Iso04[i]=NOVAL_F;
        D0[i]=NOVAL_F;
        D0err[i]=NOVAL_F;
        Dxy[i]=NOVAL_F;
        Dxyerr[i]=NOVAL_F;
        Dz[i]=NOVAL_F;
        Dzerr[i]=NOVAL_F;
        IsLooseMuon[i]=NOVAL_F;
        IsSoftMuon[i]=NOVAL_F;
        IsTightMuon[i]=NOVAL_F;
        IsPFMuon[i]=NOVAL_F;
        IsGlobalMuon[i]=NOVAL_F;
        IsTrackerMuon[i]=NOVAL_F;
        GlbTrkNormChi2[i]=NOVAL_F;
        NumberValidMuonHits[i]=NOVAL_F;
        NumberMatchedStations[i]=NOVAL_F;
        NumberValidPixelHits[i]=NOVAL_F;
        NumberTrackerLayers[i]=NOVAL_F;
        NumberOfValidTrackerHits[i]=NOVAL_F;
        NumberOfPixelLayers[i]=NOVAL_F;
        InTrkNormChi2[i]=NOVAL_F;
        SumChargedHadronPt[i]=NOVAL_F;
        SumNeutralHadronPt[i]=NOVAL_F;
        SumPhotonPt[i]=NOVAL_F;
        SumPUPt[i]=NOVAL_F;
        GenMuonY[i]=NOVAL_F;
        GenMuonEta[i]=NOVAL_F;
        GenMuonPhi[i]=NOVAL_F;
        GenMuonPt[i]=NOVAL_F;
        GenMuonE[i]=NOVAL_F;
        GenMuonCharge[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    void SetCurrent(size_t it) {
      Mass[NLEP]              = Mass[it];             
      Pt[NLEP]		      = Pt[it];		     
      Eta[NLEP]		      = Eta[it];		     
      Y[NLEP]		      = Y[it];		     
      Phi[NLEP]		      = Phi[it];		     
      E[NLEP]		      = E[it];		     
      Charge[NLEP]	      = Charge[it];
      Iso04[NLEP]                    = Iso04[it];                   
      D0[NLEP] 			     = D0[it]; 			    
      D0err[NLEP] 		     = D0err[it]; 		    
      Dxy[NLEP] 		     = Dxy[it]; 			    
      Dxyerr[NLEP] 		     = Dxyerr[it]; 		    
      Dz[NLEP] 			     = Dz[it]; 			    
      Dzerr[NLEP] 		     = Dzerr[it]; 		    
      IsLooseMuon[NLEP] 	     = IsLooseMuon[it]; 	    
      IsSoftMuon[NLEP] 		     = IsSoftMuon[it]; 		    
      IsTightMuon[NLEP] 	     = IsTightMuon[it]; 	    
      IsPFMuon[NLEP] 		     = IsPFMuon[it]; 		    
      IsGlobalMuon[NLEP] 	     = IsGlobalMuon[it]; 	    
      IsTrackerMuon[NLEP] 	     = IsTrackerMuon[it]; 	    
      GlbTrkNormChi2[NLEP] 	     = GlbTrkNormChi2[it]; 	    
      NumberValidMuonHits[NLEP]      = NumberValidMuonHits[it];     
      NumberMatchedStations[NLEP]    = NumberMatchedStations[it];   
      NumberValidPixelHits[NLEP]     = NumberValidPixelHits[it];    
      NumberTrackerLayers[NLEP]      = NumberTrackerLayers[it];     
      NumberOfValidTrackerHits[NLEP] = NumberOfValidTrackerHits[it];
      NumberOfPixelLayers[NLEP]      = NumberOfPixelLayers[it];     
      InTrkNormChi2[NLEP] 	     = InTrkNormChi2[it]; 	    
      SumChargedHadronPt[NLEP] 	     = SumChargedHadronPt[it]; 	    
      SumNeutralHadronPt[NLEP] 	     = SumNeutralHadronPt[it]; 	    
      SumPhotonPt[NLEP] 	     = SumPhotonPt[it]; 	    
      SumPUPt[NLEP] 		     = SumPUPt[it]; 		    
      GenMuonY[NLEP] 		     = GenMuonY[it]; 		    
      GenMuonEta[NLEP] 		     = GenMuonEta[it];		    
      GenMuonPhi[NLEP]		     = GenMuonPhi[it]; 		    
      GenMuonPt[NLEP] 		     = GenMuonPt[it]; 		    
      GenMuonE[NLEP] 		     = GenMuonE[it]; 		    
      GenMuonCharge[NLEP] 	     = GenMuonCharge[it]; 	    
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } mu;
  
  class JetVars {
  public:
    // Basic
    float Mass[NJET+1];
    float Pt[NJET+1];
    float Eta[NJET+1];
    float Y[NJET+1];
    float Phi[NJET+1];
    float E[NJET+1];
    float Charge[NJET+1];
    // B-TAGGING
    float CSV[NJET+1];
    float CSVV1[NJET+1];
    // GEN PARTON
    float GenPartonY[NJET+1];
    float GenPartonEta[NJET+1];
    float GenPartonPhi[NJET+1];
    float GenPartonPt[NJET+1];
    float GenPartonE[NJET+1];
    float GenPartonCharge[NJET+1];
    float PartonFlavour[NJET+1];
    float HadronFlavour[NJET+1];
    // GEN JET
    float GenJetY[NJET+1];
    float GenJetEta[NJET+1];
    float GenJetPhi[NJET+1];
    float GenJetPt[NJET+1];
    float GenJetE[NJET+1];
    float GenJetCharge[NJET+1];
    // CONSTITUENTS
    float muonMultiplicity[NJET+1];
    float PhotonEnergy[NJET+1];
    float ElectronEnergy[NJET+1];
    float MuonEnergy[NJET+1];
    float HFHadronEnergy[NJET+1];
    float HFEMEnergy[NJET+1];
    float ChargedHadronMultiplicity[NJET+1];
    float numberOfDaughters[NJET+1];
    float chargedMultiplicity[NJET+1];
    float neutralHadronMultiplicity[NJET+1];
    float neutralHadronEnergy[NJET+1];
    float neutralEmEnergy[NJET+1];
    float chargedEmEnergy[NJET+1];
    float chargedHadronEnergy[NJET+1];
    float photonMultiplicity[NJET+1];
    float electronMultiplicity[NJET+1];
    float HFHadronMultiplicity[NJET+1];
    float HFEMMultiplicity[NJET+1];
    float ChargeMuEnergy[NJET+1];
    float neutralMultiplicity[NJET+1];
    //FOR JEC
    float jecFactor0[NJET+1];
    float jetArea[NJET+1];
    // FOR SYSTEMATICS
    float SmearedPt[NJET+1];
    float SmearedPEta[NJET+1];
    float SmearedPhi[NJET+1];
    float SmearedE[NJET+1];
    float JERup[NJET+1];
    float JERdown[NJET+1];
    
    void init() {
      for (size_t i=0; i<=NJET; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        CSV[i]=NOVAL_F;
        CSVV1[i]=NOVAL_F;
        GenPartonY[i]=NOVAL_F;
        GenPartonEta[i]=NOVAL_F;
        GenPartonPhi[i]=NOVAL_F;
        GenPartonPt[i]=NOVAL_F;
        GenPartonE[i]=NOVAL_F;
        GenPartonCharge[i]=NOVAL_F;
        PartonFlavour[i]=NOVAL_F;
        HadronFlavour[i]=NOVAL_F;
        GenJetY[i]=NOVAL_F;
        GenJetEta[i]=NOVAL_F;
        GenJetPhi[i]=NOVAL_F;
        GenJetPt[i]=NOVAL_F;
        GenJetE[i]=NOVAL_F;
        GenJetCharge[i]=NOVAL_F;
        muonMultiplicity[i]=NOVAL_F;
        PhotonEnergy[i]=NOVAL_F;
        ElectronEnergy[i]=NOVAL_F;
        MuonEnergy[i]=NOVAL_F;
        HFHadronEnergy[i]=NOVAL_F;
        HFEMEnergy[i]=NOVAL_F;
        ChargedHadronMultiplicity[i]=NOVAL_F;
        numberOfDaughters[i]=NOVAL_F;
        chargedMultiplicity[i]=NOVAL_F;
        neutralHadronMultiplicity[i]=NOVAL_F;
        neutralHadronEnergy[i]=NOVAL_F;
        neutralEmEnergy[i]=NOVAL_F;
        chargedEmEnergy[i]=NOVAL_F;
        chargedHadronEnergy[i]=NOVAL_F;
        photonMultiplicity[i]=NOVAL_F;
        electronMultiplicity[i]=NOVAL_F;
        HFHadronMultiplicity[i]=NOVAL_F;
        HFEMMultiplicity[i]=NOVAL_F;
        ChargeMuEnergy[i]=NOVAL_F;
        neutralMultiplicity[i]=NOVAL_F;
	jecFactor0[i]=NOVAL_F;
	jetArea[i]=NOVAL_F;
        SmearedPt[i]=NOVAL_F;
        SmearedPEta[i]=NOVAL_F;
        SmearedPhi[i]=NOVAL_F;
        SmearedE[i]=NOVAL_F;
        JERup[i]=NOVAL_F;
        JERdown[i]=NOVAL_F;
      }
    }
    
    void SetCurrent(size_t it) {
      Mass[NJET]                        = Mass[it];			      
      Pt[NJET]			        = Pt[it];			      
      Eta[NJET]			        = Eta[it];			      
      Y[NJET]			        = Y[it];			      
      Phi[NJET]			        = Phi[it];			      
      E[NJET]			        = E[it];			      
      Charge[NJET]		        = Charge[it];		      
      CSV[NJET]			        = CSV[it];			      
      CSVV1[NJET]		        = CSVV1[it];		      
      GenPartonY[NJET]		        = GenPartonY[it];		      
      GenPartonEta[NJET]		= GenPartonEta[it];		      
      GenPartonPhi[NJET]		= GenPartonPhi[it];		      
      GenPartonPt[NJET]		        = GenPartonPt[it];		      
      GenPartonE[NJET]		        = GenPartonE[it];		      
      GenPartonCharge[NJET]	        = GenPartonCharge[it];	      
      PartonFlavour[NJET]	        = PartonFlavour[it];	      
      HadronFlavour[NJET]	        = HadronFlavour[it];	      
      GenJetY[NJET]		        = GenJetY[it];		      
      GenJetEta[NJET]		        = GenJetEta[it];		      
      GenJetPhi[NJET]		        = GenJetPhi[it];		      
      GenJetPt[NJET]		        = GenJetPt[it];		      
      GenJetE[NJET]		        = GenJetE[it];		      
      GenJetCharge[NJET]		= GenJetCharge[it];		      
      muonMultiplicity[NJET]	        = muonMultiplicity[it];	      
      PhotonEnergy[NJET]		= PhotonEnergy[it];		      
      ElectronEnergy[NJET]	        = ElectronEnergy[it];	      
      MuonEnergy[NJET]		        = MuonEnergy[it];		      
      HFHadronEnergy[NJET]	        = HFHadronEnergy[it];	      
      HFEMEnergy[NJET]		        = HFEMEnergy[it];		      
      ChargedHadronMultiplicity[NJET]   = ChargedHadronMultiplicity[it];  
      numberOfDaughters[NJET]	        = numberOfDaughters[it];	      
      chargedMultiplicity[NJET]	        = chargedMultiplicity[it];	      
      neutralHadronMultiplicity[NJET]   = neutralHadronMultiplicity[it];  
      neutralHadronEnergy[NJET]         = neutralHadronEnergy[it];
      neutralEmEnergy[NJET]             = neutralEmEnergy[it];    
      chargedEmEnergy[NJET]             = chargedEmEnergy[it];    
      chargedHadronEnergy[NJET]         = chargedHadronEnergy[it];
      photonMultiplicity[NJET]	        = photonMultiplicity[it];	      
      electronMultiplicity[NJET]	= electronMultiplicity[it];	      
      HFHadronMultiplicity[NJET]	= HFHadronMultiplicity[it];	      
      HFEMMultiplicity[NJET]	        = HFEMMultiplicity[it];	      
      ChargeMuEnergy[NJET]	        = ChargeMuEnergy[it];	      
      neutralMultiplicity[NJET]	        = neutralMultiplicity[it];
      jecFactor0[NJET]                  = jecFactor0[it];
      jetArea[NJET]                     = jetArea[it];
      SmearedPt[NJET]		        = SmearedPt[it];		      
      SmearedPEta[NJET]		        = SmearedPEta[it];		      
      SmearedPhi[NJET]		        = SmearedPhi[it];		      
      SmearedE[NJET]		        = SmearedE[it];		      
      JERup[NJET]		        = JERup[it];		      
      JERdown[NJET]                     = JERdown[it];
    }
  } jet;
  
  class AK8Vars {
  public:
    float vSubjetIndex0[NJET+1];
    float vSubjetIndex1[NJET+1];
    float topSubjetIndex0[NJET+1];
    float topSubjetIndex1[NJET+1];
    float topSubjetIndex2[NJET+1];
    float topSubjetIndex3[NJET+1];
    float tau1[NJET+1];
    float tau2[NJET+1];
    float tau3[NJET+1];
    float trimmedMass[NJET+1];
    float prunedMass[NJET+1];
    float filteredMass[NJET+1];
    float topMass[NJET+1];
    float wMass[NJET+1];
    float nSubJets[NJET+1];
    float minmass[NJET+1];
    
    void init() {
      for (size_t it=0; it<=NJET; ++it) {
	vSubjetIndex0[it]=NOVAL_F;
	vSubjetIndex1[it]=NOVAL_F;
	topSubjetIndex0[it]=NOVAL_F;
	topSubjetIndex1[it]=NOVAL_F;
	topSubjetIndex2[it]=NOVAL_F;
	topSubjetIndex3[it]=NOVAL_F;
	tau1[it]=NOVAL_F;
	tau2[it]=NOVAL_F;
	tau3[it]=NOVAL_F;
	trimmedMass[it]=NOVAL_F;
	prunedMass[it]=NOVAL_F;
	filteredMass[it]=NOVAL_F;
	topMass[it]=NOVAL_F;
	wMass[it]=NOVAL_F;
	nSubJets[it]=NOVAL_F;
	minmass[it]=NOVAL_F;
      }
    }
    
    void SetCurrent(size_t it) {
      vSubjetIndex0[NJET] = vSubjetIndex0[it];
      vSubjetIndex1[NJET] = vSubjetIndex1[it];
      topSubjetIndex0[NJET] = topSubjetIndex0[it];
      topSubjetIndex1[NJET] = topSubjetIndex1[it];
      topSubjetIndex2[NJET] = topSubjetIndex2[it];
      topSubjetIndex3[NJET] = topSubjetIndex3[it];
      tau1[NJET]	 = tau1[it];
      tau2[NJET]	 = tau2[it];
      tau3[NJET]	 = tau3[it];
      trimmedMass[NJET]	 = trimmedMass[it];
      prunedMass[NJET]	 = prunedMass[it];
      filteredMass[NJET] = filteredMass[it];
      topMass[NJET]	 = topMass[it];
      wMass[NJET]	 = wMass[it];
      nSubJets[NJET]	 = nSubJets[it];
      minmass[NJET]	 = minmass[it];
    }
  } AK8;
  
  class AK4JetVars : public JetVars {
  public:
    AK4JetVars() { init(); }
    ~AK4JetVars() {}
    
    size_t it;
    int size;
    
    // Razor variables
    float MR;
    float MTR;
    float R;
    float R2;
    
    void init() {
      JetVars::init();
      
      it = -1;
      size = 0;
      
      MR = NOVAL_F;
      MTR = NOVAL_F;
      R = NOVAL_F;
      R2 = NOVAL_F;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	JetVars::SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } jetsAK4;
  
  class AK8JetVars : public JetVars, public AK8Vars {
  public:
    AK8JetVars() { init(); }
    ~AK8JetVars() {}
    
    size_t it;
    int size;
    
    // Razor variables
    float MR;
    float MTR;
    float R;
    float R2;
    
    void init() {
      JetVars::init();
      AK8Vars::init();
      
      it = -1;
      size = 0;
      
      MR = NOVAL_F;
      MTR = NOVAL_F;
      R = NOVAL_F;
      R2 = NOVAL_F;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	JetVars::SetCurrent(it);
	AK8Vars::SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } jetsAK8;
  
  class AK8SubJetVars : public JetVars {
  public:
    AK8SubJetVars() { init(); }
    ~AK8SubJetVars() {}
    
    float subjetCSV[NJET+1];
    
    size_t it;
    int size;
    
    void init() {
      JetVars::init();
      for (size_t it=0; it<=NJET; ++it) subjetCSV[it]=NOVAL_F;
      
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	JetVars::SetCurrent(it);
	subjetCSV[NJET] = subjetCSV[it];
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } subjetsAK8;
  
  class CmsTopTagSubJetVars : public JetVars {
  public:
    CmsTopTagSubJetVars() { init(); }
    ~CmsTopTagSubJetVars() {}
    
    size_t it;
    int size;
    
    void init() {
      JetVars::init();
      
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	JetVars::SetCurrent(it);
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } subjetsCmsTopTag;
  
  class MetData {
  public:
    MetData() { init(); };
    
    float Pt;
    float Phi;
    float Px;
    float Py;
    
    void init() {
      Pt=NOVAL_F;
      Phi=NOVAL_F;
      Px=NOVAL_F;
      Py=NOVAL_F;
    }
    
  } met;
  
  class EventData {
  public:
    EventData() { init(); };
    
    int NLep;
    int NTopHad;
    int NTopLep;
    int NTop;
    float HtLep;
    float HtTop;
    float Ht;
    float HtAll;
    float HtEx;
    float HtExFr;
    float HtTopFr;
    float TTHadDR;
    float TTHadDPhi;
    float TTHadDEta;
    float TTHadMass;
    float TTHadPz;
    float TTHadHz;
    float TTHadDPz;
    float TTHadMR;
    float TTHadMTR;
    float TTHadR;
    float TTHadR2;
    float MR;
    float MTR;
    float R;
    float R2;
    
    // Top tagging variables
    float nhadtops;
    float nleptops;
    float ntops;
    float tt_dR;
    float tt_dPhi;
    float tt_dEta;
    float tt_Mass;
    float tt_Pz;
    float tt_Hz;
    float tt_dPz;
    float tt_extra;
    float dHt;
    float dphi1;
    float dphi2;
    float dPhi_met_t1;
    float dPhi_met_t2;
    
    float tt_MR;
    float tt_MTR;
    float tt_R;
    float tt_R2;
    
    // Other Variables
    float HT;
    float HTall;
    float HTtt;
    float HTlep;
    float HTex;
    float HTttFraction;
    float HTexFraction;
    
    // Development 04 March
    int nmu;
    int nele;
    int neletight;
    int nmuveto;
    int neleveto;
    float DRJetLep[NJET+1];
    float EleDRJet[NLEP+1];
    float MuDRJet[NLEP+1];
    float RelPtJetLep[NJET+1];
    float EleRelPtJet[NLEP+1];
    float MuRelPtJet[NLEP+1];
    float EleJetCombMass[NLEP+1];
    float MuJetCombMass[NLEP+1];
    int nhadtoplike;

    // Development 20 April
    int JetGenTruth[NJET];
    bool JetHasMatchedGenTop[NJET];
    int JetMatchedGenTopType[NJET];
    bool JetMatchedGenTopIsMerged[NJET];
    float JetMatchedGenTopPt[NJET];
    float JetMatchedGenTopJetDR[NJET];
    float GenBJetDR[NJET];
    float GenWJetDR[NJET];
    float GenWGenBDR[NJET];
    float GenLepJetDR[NJET];
    float GenLepGenBDR[NJET];
    int NGenLepFromTop;
    bool IsGenTop[NGEN];
    int GenTopType[NGEN];
    bool GenTopHasMatchedJet[NGEN];
    bool GenTopHasMatchedTopTagJet[NGEN];
    bool JetIsHadTopTagged[NJET];
    float maxSubjetCSV[NJET];
    
    void init() {
      NLep=NOVAL_I;
      NTopHad=NOVAL_I;
      NTopLep=NOVAL_I;
      NTop=NOVAL_I;
      HtLep=NOVAL_F;
      HtTop=NOVAL_F;
      Ht=NOVAL_F;
      HtAll=NOVAL_F;
      HtEx=NOVAL_F;
      HtExFr=NOVAL_F;
      HtTopFr=NOVAL_F;
      TTHadDR=NOVAL_F;
      TTHadDPhi=NOVAL_F;
      TTHadDEta=NOVAL_F;
      TTHadMass=NOVAL_F;
      TTHadPz=NOVAL_F;
      TTHadHz=NOVAL_F;
      TTHadDPz=NOVAL_F;
      TTHadMR=NOVAL_F;
      TTHadMTR=NOVAL_F;
      TTHadR=NOVAL_F;
      TTHadR2=NOVAL_F;
      MR=NOVAL_F;
      MTR=NOVAL_F;
      R=NOVAL_F;
      R2=NOVAL_F;
      
      nhadtops=NOVAL_F;
      nleptops=NOVAL_F;
      ntops=NOVAL_F;
      tt_dR=NOVAL_F;
      tt_dPhi=NOVAL_F;
      tt_dEta=NOVAL_F;
      tt_Mass=NOVAL_F;
      tt_Pz=NOVAL_F;
      tt_Hz=NOVAL_F;
      tt_dPz=NOVAL_F;
      tt_extra=NOVAL_F;
      dHt=NOVAL_F;
      dphi1=NOVAL_F;
      dphi2=NOVAL_F;
      dPhi_met_t1=NOVAL_F;
      dPhi_met_t2=NOVAL_F;
      
      tt_MR = NOVAL_F;
      tt_MTR = NOVAL_F;
      tt_R = NOVAL_F;
      tt_R2 = NOVAL_F;
      
      // Other Variables
      HT=NOVAL_F;
      HTall=NOVAL_F;
      HTtt=NOVAL_F;
      HTex=NOVAL_F;
      HTttFraction=NOVAL_F;
      HTexFraction=NOVAL_F;

      // Development 04 March
      nmu=NOVAL_I;
      nele=NOVAL_I;
      neletight=NOVAL_I;
      nmuveto=NOVAL_I;
      neleveto=NOVAL_I;
      for (size_t i=0; i<=NJET; ++i) {
	DRJetLep[i]=NOVAL_F;
	RelPtJetLep[i]=NOVAL_F;
      }
      for (size_t i=0; i<=NLEP; ++i) {
	EleDRJet[i]=NOVAL_F;
	MuDRJet[i]=NOVAL_F;
	EleRelPtJet[i]=NOVAL_F;
	MuRelPtJet[i]=NOVAL_F;
	EleJetCombMass[i]=NOVAL_F;
	MuJetCombMass[i]=NOVAL_F;
      }
      nhadtoplike=NOVAL_I;
      
      // Development 20 April
      for (size_t i=0; i<NJET; ++i) {
	JetGenTruth[i]=NOVAL_I;
	JetHasMatchedGenTop[i]=0;
	JetMatchedGenTopType[i]=NOVAL_I;
	JetMatchedGenTopIsMerged[i]=0;
	JetMatchedGenTopPt[i]=NOVAL_F;
        JetMatchedGenTopJetDR[i]=NOVAL_F;
        GenBJetDR[i]=NOVAL_F;
        GenWJetDR[i]=NOVAL_F;
        GenWGenBDR[i]=NOVAL_F;
        GenLepJetDR[i]=NOVAL_F;
        GenLepGenBDR[i]=NOVAL_F;
	JetIsHadTopTagged[i]=0;
	maxSubjetCSV[i]=NOVAL_F;
      }
      NGenLepFromTop=NOVAL_I;
      for (size_t i=0; i<NGEN; ++i) {
	IsGenTop[i]=0;
	GenTopType[i]=NOVAL_I;
	GenTopHasMatchedJet[i]=0;
	GenTopHasMatchedTopTagJet[i]=0;
      }
    }
  } evt;
  
  void CalculateAllVariables() {
    //calcRazorAK4_();
    //calcRazorAK8_();
    //calcRazorCmsTopTag_();
    
    //--------------------------------------------------------------------------
    //                           Generator Particles
    //--------------------------------------------------------------------------
    
    // Make a list of Generator level objects and save them to vectors
    std::vector<TLorentzVector> gen_top;
    std::vector<size_t > gen_top_it;
    std::vector<int> gen_top_ID;
    std::vector<TLorentzVector> gen_W;
    std::vector<int> gen_W_ID;
    std::vector<TLorentzVector> gen_b;
    std::vector<int> gen_b_ID;
    std::vector<TLorentzVector> gen_lep;
    std::vector<TLorentzVector> gen_neu;
    while(gen.Loop()) {
      evt.IsGenTop[gen.it]=0;
      evt.GenTopType[gen.it] = NOVAL_I;
      evt.GenTopHasMatchedJet[gen.it]=0;
      evt.GenTopHasMatchedTopTagJet[gen.it]=0;
      if (gen.Pt[NGEN]>0) {
	TLorentzVector genp; genp.SetPtEtaPhiE(gen.Pt[NGEN], gen.Eta[NGEN], gen.Phi[NGEN], gen.E[NGEN]);
        if (gen.ID[NGEN]!=gen.MomID[NGEN]) {
          if (abs(gen.ID[NGEN])==6) { 
	    evt.IsGenTop[gen.it]=1; 
	    gen_top.push_back(genp); 
	    gen_top_it.push_back(gen.it); 
	    gen_top_ID.push_back(gen.ID[NGEN]);
	    evt.GenTopType[gen.it] = 0;
	  }
          if (abs(gen.ID[NGEN])==5&&abs(gen.MomID[NGEN])==6) { gen_b.push_back(genp); gen_b_ID.push_back(gen.ID[NGEN]); }
          if (abs(gen.ID[NGEN])==24&&abs(gen.MomID[NGEN])==6) { gen_W.push_back(genp); gen_W_ID.push_back(gen.ID[NGEN]); }
          if ((abs(gen.ID[NGEN])==11||abs(gen.ID[NGEN])==13)&&(abs(gen.MomID[NGEN])==24)) gen_lep.push_back(genp);
          if ((abs(gen.ID[NGEN])==12||abs(gen.ID[NGEN])==14)&&(abs(gen.MomID[NGEN])==24)) gen_neu.push_back(genp);
        } else if (gen.ID[NGEN]==gen.MomID[NGEN]) {
	  // tops emit gluons(?) and has to match consecutive tops to the original one
	  if (abs(gen.ID[NGEN])==6) {
	    size_t i=0, i_m_dR = -1, i_m_dE = -1;
	    double min_dE = 9999, min_dR = 9999;
	    while(i<gen_top.size()) {
	      if (gen_top_ID[i]==gen.MomID[NGEN]) {
		double dE = gen_top[i].E()-genp.E();
		double dR = gen_top[i].DeltaR(genp);
		if (fabs(dE)<fabs(min_dE)) {
		  min_dE = dE;
		  i_m_dE = i;
		}
		if (dR<min_dR) {
		  min_dR = dR;
		  i_m_dR = i;
		}
	      }
	      ++i;
	    }
	    //std::cout<<"match: dE: "<<gen_top[imatch].E()-genp.E()<<" Ep: "<<gen_top[imatch].E()<<" Ec: "<<genp.E()<<" dR: "<<genp.DeltaR(gen_top[imatch])<<std::endl;
	    //if (i_m_dE!=i_m_dR) printf("match: dE: %03.1f iE: %d   Ep: %04.1f Ec: %04.1f dR: %1.3f      min_dR: %1.3f Ec2: %4.1f iR: %d\n", double(min_dE), int(i_m_dE), double(gen_top[i_m_dE].E()), double(genp.E()), double(gen_top[i_m_dE].DeltaR(genp)), double(min_dR), double(gen_top[i_m_dR].E()), int(i_m_dR));
	    size_t imatch = (i_m_dE==i_m_dR) ? i_m_dE : ( (fabs(min_dE)/gen_top[i_m_dE].E()<0.1) ? i_m_dE : i_m_dR );
	    evt.IsGenTop[gen_top_it[imatch]]=0;
	    evt.IsGenTop[gen.it]=1;
	    gen_top[imatch]=genp;
	    gen_top_it[imatch]=gen.it;
	  }
	}
      }
    }
    
    // std::cout<<"Start looping on Gen\n";
    // //while(gen.Loop()) if (abs(gen.ID[NGEN])==1000021&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found ~g, it="<<gen.it<<" ID=~"<<(abs(gen.ID[NGEN])-1e6)<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[NGEN])==1000006&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found ~t, it="<<gen.it<<" ID=~"<<(abs(gen.ID[NGEN])-1e6)<<" MomID=~"<<(abs(gen.MomID[NGEN])-1e6)<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[NGEN])==1000024&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found ~chi+, it="<<gen.it<<" ID=~"<<(abs(gen.ID[NGEN])-1e6)<<" MomID=~"<<(abs(gen.MomID[NGEN])-1e6)<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[NGEN])==1000022&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found ~chi0, it="<<gen.it<<" ID=~"<<(abs(gen.ID[NGEN])-1e6)<<" MomID=~"<<(abs(gen.MomID[NGEN])-1e6)<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[NGEN])==6&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found t, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[NGEN])==5&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found b, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[NGEN])==24&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found W, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // while(gen.Loop()) if ((abs(gen.ID[NGEN])==11||abs(gen.ID[NGEN])==13)&&(abs(gen.MomID[NGEN])==24||abs(gen.MomID[NGEN])>1e6)) std::cout<<"Found mu/e, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    // std::cout<<"\n";
    //if ((gen_top[0].Pt()>400&&gen_top[1].Pt()>400)&&gen_lep.size()>0) {
    //  std::cout<<"Start looping on Gen\n";
    //  while(gen.Loop()) if (abs(gen.ID[NGEN])==6&&gen.ID[NGEN]!=gen.MomID[NGEN]) std::cout<<"Found t, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //  while(gen.Loop()) if (abs(gen.ID[NGEN])==5&&abs(gen.MomID[NGEN])==6) std::cout<<"Found b, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //  while(gen.Loop()) if (abs(gen.ID[NGEN])==24&&abs(gen.MomID[NGEN])==6) std::cout<<"Found W, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //  while(gen.Loop()) if ((abs(gen.ID[NGEN])==11||abs(gen.ID[NGEN])==13)&&abs(gen.MomID[NGEN])==24) std::cout<<"Found mu/e from W mother, it="<<gen.it<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //  if (good_W_matches&&nlep_from_top>0) for (size_t i=0; i<gen_top_matched_W_matched_lep.size(); ++i) {
    //    TLorentzVector b = gen_top_matched_b[top_parent[i]];
    //    TLorentzVector W = gen_top_matched_W[top_parent[i]];
    //    TLorentzVector lep = gen_top_matched_W_matched_lep[i];
    //    double DR_lep_to_b = lep.DeltaR(b);
    //    std::cout<<"DR(lep, b)="<<DR_lep_to_b<<std::endl;
    //  }
    //  std::cout<<"\n";
    //}
    
    // Find bs and Ws
    // Method: bs and Ws with top parent are combined
    // Best pair with lowest combined mass and DR difference is selected
    std::vector<TLorentzVector> gen_top_matched_b;
    std::vector<TLorentzVector> gen_top_matched_W;
    std::vector<bool> W_is_leptonic;
    bool good_W_matches = true;
    for (size_t i=0; i<gen_top.size(); ++i) {
      // Match b and W to t
      size_t j_b = -1, k_W = -1;
      double min_DM = 9999, min_DR = 9999;
      if (gen_b.size()<gen_top.size()||gen_W.size()<gen_top.size()) {
	//std::cout<<"Not enough b/W with top parent"<<std::endl;
	good_W_matches = false;
      } else {
        for (size_t j=0; j<gen_b.size(); ++j) {
          for (size_t k=0; k<gen_W.size(); ++k) {
            TLorentzVector bW_comb = gen_b[j]+gen_W[k];
            double DR = gen_top[i].DeltaR(bW_comb);
            double DM = fabs(gen_top[i].M()-bW_comb.M());
            if (DR<0.8) {
              if (DM<min_DM) {
                min_DM = DM;
		min_DR = DR;
                j_b = j;
                k_W = k;
              }
            }
          }
        }
	//printf("W/b to top match: %.6f %.6f\n", min_DR, min_DM);
	if (min_DR<0.8&&min_DM<1) {
	  gen_top_matched_b.push_back(gen_b[j_b]);
	  gen_top_matched_W.push_back(gen_W[k_W]);
	} else {
	  good_W_matches = false;
	}
      }
    }
    //if (!good_W_matches) {
    //  while(gen.Loop()) std::cout<<"it="<<gen.it<<" ID=   "<<gen.ID[NGEN]<<" MomID=   "<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<" E="<<gen.E[NGEN]<<std::endl;
    //  std::cout<<""<<std::endl;
    //}
    //while(gen.Loop()) if (gen.ID[NGEN]==6) std::cout<<"Found t it="<<gen.it<<" ID="<<gen.ID[NGEN]<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //while(gen.Loop()) if (gen.ID[NGEN]==-6) std::cout<<"Found t it="<<gen.it<<" ID="<<gen.ID[NGEN]<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //while(gen.Loop()) if (gen.MomID[NGEN]==6) std::cout<<"Found child it="<<gen.it<<" ID="<<gen.ID[NGEN]<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    //while(gen.Loop()) if (gen.MomID[NGEN]==-6) std::cout<<"Found child it="<<gen.it<<" ID="<<gen.ID[NGEN]<<" MomID="<<gen.MomID[NGEN]<<" Pt="<<gen.Pt[NGEN]<<" Eta="<<gen.Eta[NGEN]<<" Phi="<<gen.Phi[NGEN]<<std::endl;
    
    // If we have lepton from W, find parent
    // Do as above with tops, but use neutrino and lepton instead to find W parent
    // In the end associate with top already found
    evt.NGenLepFromTop = 0;
    std::vector<TLorentzVector> gen_top_matched_W_matched_lep;
    std::vector<TLorentzVector> gen_top_matched_W_matched_neu;
    for (size_t i=0; i<gen_top_matched_W.size(); ++i) {
      TLorentzVector lep, neu;
      // Match lep and neutrino to W
      size_t j_l = -1, k_n = -1;
      double min_DM = 9999, min_DR = 9999;
      for (size_t j=0; j<gen_lep.size(); ++j) {
	for (size_t k=0; k<gen_neu.size(); ++k) {
	  TLorentzVector ln_comb = gen_lep[j]+gen_neu[k];
	  double DR = gen_top_matched_W[i].DeltaR(ln_comb);
	  double DM = fabs(gen_top_matched_W[i].M()-ln_comb.M());
	  if (DR<0.8) {
	    if (DM<min_DM) {
	      min_DM = DM;
	      min_DR = DR;
	      j_l = j;
	      k_n = k;
	    }
	  }
	}
      }
      bool lep_found = (min_DR<0.8&&min_DM<1);
      W_is_leptonic.push_back(lep_found);
      evt.GenTopType[gen_top_it[i]] = 1;
      //printf("l/v to W match: %.6f %.6f\n", min_DR, min_DM);
      if (lep_found) {
	lep = gen_lep[j_l];
	neu = gen_neu[k_n];
	++evt.NGenLepFromTop;
      }
      gen_top_matched_W_matched_lep.push_back(lep);
      gen_top_matched_W_matched_neu.push_back(neu);
    }
    
    // Match jets to tops (find the closest jet to top, sort by distance from gen)
    // Could also do genjet matching (but this is done in B2G)
    std::vector<TLorentzVector> temp = gen_top;
    std::vector<size_t > temp_it = gen_top_it;
    std::map<size_t, size_t > jet_gentop_it;
    const bool verbose = 0;
    while (temp.size()) {
      // find gentop  - jet pair with lowest DR (associate them and remove from temp top colelction)
      double min_DR = 9999, matched_DR = 9999;
      size_t top_m_it = -1, top_closest_it = -1, jet_m_it = -1;
      for (size_t i=0; i<temp.size(); ++i) {
	TLorentzVector top = temp[i];
	while(jetsAK8.Loop()) {
	  TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
	  double DR = jet.DeltaR(top);
	  if (DR<min_DR) {
	    min_DR = DR;
	    top_closest_it = i;
	    if (!jet_gentop_it.count(jetsAK8.it)) {
	      matched_DR = DR;
	      top_m_it = i;
	      jet_m_it = jetsAK8.it;
	    }
	  }
	}
      }
      if (matched_DR<0.8) {
	if (verbose) std::cout<<"Top-jet match found, top(gen) it="<<temp_it[top_m_it]<<" jet it="<<jet_m_it<<" dR="<<matched_DR<<std::endl;
	jet_gentop_it[jet_m_it] = top_m_it;
	evt.GenTopHasMatchedJet[temp_it[top_m_it]] = 1;
	evt.GenTopHasMatchedTopTagJet[temp_it[top_m_it]] = min_DR<0.8 && jetsAK8.tau3[jet_m_it]/jetsAK8.tau2[jet_m_it]<0.75 && jetsAK8.Pt[jet_m_it]>400 && jetsAK8.prunedMass[jet_m_it]>140;
	temp.erase(temp.begin()+top_m_it);
	temp_it.erase(temp_it.begin()+top_m_it);
      } else if (jetsAK8.size) {
	if (verbose) {
	  std::cout<<"No match  found, possible pairs:"<<std::endl;
	  for (size_t i=0; i<temp.size(); ++i) {
	    TLorentzVector top = temp[i];
	    while(jetsAK8.Loop()) {
	      TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
	      double DR = jet.DeltaR(top);
	      std::cout<<"  top(gen) it="<<temp_it[i]<<" jet it="<<jetsAK8.it<<" dR="<<DR<<(jet_gentop_it.count(jetsAK8.it)?" (Already found)":"")<<std::endl;
	    }
	  }
	}
	temp.erase(temp.begin()+top_closest_it);
	temp_it.erase(temp_it.begin()+top_closest_it);
      } else {
	if (verbose) std::cout<<"No jets in event!!!!\n";
	temp.clear();
	temp_it.clear();
      }
    }
    if (verbose) std::cout<<std::endl;
    
    //--------------------------------------------------------------------------
    //                               Leptons
    //--------------------------------------------------------------------------
    
    // find good leptons (for letponic tops)
    int ngoodleptons = 0;
    std::vector<TLorentzVector> goodleps;
    evt.HTlep = 0;
    evt.nmu = 0;
    evt.nele = 0;
    evt.neletight = 0;
    evt.nmuveto = 0;
    evt.neleveto = 0;
    while(ele.Loop()) {
      TLorentzVector el; el.SetPtEtaPhiE(ele.Pt[NLEP], ele.Eta[NLEP], ele.Phi[NLEP], ele.E[NLEP]);
      evt.EleDRJet[ele.it] = 9999;
      evt.EleRelPtJet[ele.it] = 9999;
      evt.EleJetCombMass[ele.it] = 9999;
      while(jetsAK8.Loop()) {
        TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
        double dR = el.DeltaR(jet);
	if (dR<evt.EleDRJet[jetsAK8.it]) {
          evt.EleDRJet[ele.it] = dR;
          evt.EleRelPtJet[ele.it] = el.Perp(jet.Vect());
          evt.EleJetCombMass[ele.it] = (jet+el).M();
        }
      }
      if (ele.Pt[NLEP] > 35 && fabs(ele.Eta[NLEP]) < 2.5) {
	++ngoodleptons;
	++evt.nele;
	goodleps.push_back(el);
	evt.HTlep += ele.Pt[NLEP];
	if (ele.isVeto[NLEP]>0) ++evt.neleveto;
	if (ele.isTight[NLEP]>0) { // New 05 March
	  ++evt.neletight;
	}
      }
    }
    while(mu.Loop()) {
      TLorentzVector muon; muon.SetPtEtaPhiE(mu.Pt[NLEP], mu.Eta[NLEP], mu.Phi[NLEP], mu.E[NLEP]);
      evt.MuDRJet[mu.it] = 9999;
      evt.MuRelPtJet[mu.it] = 9999;
      evt.MuJetCombMass[mu.it] = 9999;
      while(jetsAK8.Loop()) {
        TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
        double dR = muon.DeltaR(jet);
        if (dR<evt.MuDRJet[jetsAK8.it]) {
          evt.MuDRJet[mu.it] = dR;
          evt.MuRelPtJet[mu.it] = muon.Perp(jet.Vect());
          evt.MuJetCombMass[mu.it] = (jet+muon).M();
        }
      }
      if (mu.Pt[NLEP] > 45 && fabs(mu.Eta[NLEP]) < 2.1) {
	if (mu.IsTightMuon[NLEP]>0) {
	  ++ngoodleptons;
	  ++evt.nmu;
	  goodleps.push_back(muon);
	  evt.HTlep += mu.Pt[NLEP];
	}
	if (mu.IsSoftMuon[NLEP]>0) ++evt.nmuveto;
      }
    }
    
    //--------------------------------------------------------------------------
    //                                 Jets
    //--------------------------------------------------------------------------
    
    std::vector<TLorentzVector> hadtop;
    std::vector<TLorentzVector> hadtoplike;
    std::vector<TLorentzVector> leptop;
    evt.ntops = 0;
    evt.nhadtops = 0;
    evt.nleptops = 0;
    evt.HT = 0;
    evt.nhadtoplike = 0;
    while(jetsAK8.Loop()) {
      bool is_top = false;
      TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
      // Gen particle truth
      evt.JetGenTruth[jetsAK8.it] = gen_top.size()>0;
      evt.JetHasMatchedGenTop[jetsAK8.it] = 0;
      evt.JetMatchedGenTopType[jetsAK8.it] = NOVAL_I;
      evt.JetMatchedGenTopIsMerged[jetsAK8.it] = false;
      evt.JetMatchedGenTopPt[jetsAK8.it] = NOVAL_F;
      evt.JetMatchedGenTopJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenBJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenWJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenWGenBDR[jetsAK8.it] = NOVAL_F;
      evt.GenLepJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenLepGenBDR[jetsAK8.it] = NOVAL_F;
      if (jet_gentop_it.count(jetsAK8.it)) {
	size_t top_it = jet_gentop_it[jetsAK8.it];
	evt.JetGenTruth[jetsAK8.it] = 2;
	evt.JetHasMatchedGenTop[jetsAK8.it] = 1;
	evt.JetMatchedGenTopType[jetsAK8.it] = 0;
	evt.JetMatchedGenTopPt[jetsAK8.it] = gen_top[top_it].Pt();
	evt.JetMatchedGenTopJetDR[jetsAK8.it] = gen_top[top_it].DeltaR(jet);
	// If W matching was successful, more information is available
	if (good_W_matches) {
	  evt.JetGenTruth[jetsAK8.it] = 2;
	  evt.JetMatchedGenTopType[jetsAK8.it] = W_is_leptonic[top_it];
	  // Both b and Whad/lepton within jet cone
	  if (jet.DeltaR(gen_top_matched_b[top_it])<0.7 && jet.DeltaR(W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it] : gen_top_matched_W[top_it])<0.7) {
	    evt.JetMatchedGenTopIsMerged[jetsAK8.it] = true;
	    evt.JetGenTruth[jetsAK8.it] = 3+W_is_leptonic[top_it];
	  }
	  evt.GenBJetDR[jetsAK8.it] = gen_top_matched_b[top_it].DeltaR(jet);
	  evt.GenWJetDR[jetsAK8.it] = gen_top_matched_W[top_it].DeltaR(jet);
	  evt.GenWGenBDR[jetsAK8.it] = gen_top_matched_W[top_it].DeltaR(gen_top_matched_b[top_it]);
	  evt.GenLepJetDR[jetsAK8.it] = W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it].DeltaR(jet) : NOVAL_F;
	  evt.GenLepGenBDR[jetsAK8.it] = W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it].DeltaR(gen_top_matched_b[top_it]) : NOVAL_F;
	} else {
	  evt.JetGenTruth[jetsAK8.it] = 5;
	}
      }
      
      // Subjets
      evt.maxSubjetCSV[jetsAK8.it] = 0;
      for (size_t i=0; i<jetsAK8.nSubJets[NJET]; ++i) if (i<2) {
	size_t subjet_it = i==0 ? jetsAK8.vSubjetIndex0[NJET] : jetsAK8.vSubjetIndex1[NJET];
	float sjCSV = subjetsAK8.subjetCSV[subjet_it];
	if (sjCSV > evt.maxSubjetCSV[jetsAK8.it]) evt.maxSubjetCSV[jetsAK8.it] = sjCSV;
      }
      //printf("Subjets of jet it: %d  pt: %3.1f  eta: %1.2f  phi: %1.2f minSjCSV: %1.2f\n", int(jetsAK8.it), jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], evt.maxSubjetCSV[jetsAK8.it]);
      //while(subjetsAK8.Loop()) {
      //  TLorentzVector subjet; subjet.SetPtEtaPhiE(subjetsAK8.Pt[NJET], subjetsAK8.Eta[NJET], subjetsAK8.Phi[NJET], subjetsAK8.E[NJET]);
      //  if (jet.DeltaR(subjet)<0.8) {
      //    printf("    subjet it: %d  pt: %3.1f  eta: %1.2f  phi: %1.2f CSV: %12f\n", int(subjetsAK8.it), subjetsAK8.Pt[NJET], subjetsAK8.Eta[NJET], subjetsAK8.Phi[NJET], subjetsAK8.subjetCSV[NJET]);
      //  }
      //}
      
      // fully hadronic tops
      evt.JetIsHadTopTagged[jetsAK8.it] = 0;
      //if (jetsAK8.tau1[NJET]>0 && jetsAK8.tau2[NJET]>0 ? jetsAK8.Pt[NJET] > 400 && jetsAK8.prunedMass[NJET] > 140 && (jetsAK8.tau2[NJET]/jetsAK8.tau1[NJET]) < 0.75 : 0) { // orig
      if ( (jetsAK8.tau2[NJET]>0 && jetsAK8.tau3[NJET]>0 ? jetsAK8.tau3[NJET]/jetsAK8.tau2[NJET] < 0.75 : 0 ) &&
           //jetsAK8.nSubJets[NJET] > 2 &&
           //jetsAK8.minmass[NJET] > 50 &&
           jetsAK8.Pt[NJET] > 400 &&
           jetsAK8.prunedMass[NJET] > 140) { // Latest
	evt.JetIsHadTopTagged[jetsAK8.it] = 1;
        ++evt.nhadtops;
        ++evt.nhadtoplike;
        is_top = true;
        hadtop.push_back(jet);
      } else if (jetsAK8.Pt[NJET] > 400 && jetsAK8.prunedMass[NJET] > 140) { // Top like jets for sideband fitting
      //} else if (jetsAK8.Pt[NJET] > 400) { // High pT jets for sideband fitting
        ++evt.nhadtoplike;
        hadtoplike.push_back(jet);
      }
      
      // semi-leptonic tops
      TLorentzVector lep;
      evt.DRJetLep[jetsAK8.it] = 9999;
      evt.RelPtJetLep[jetsAK8.it] = 9999;
      for (size_t i=0; i<goodleps.size(); ++i) {
        if (goodleps[i].DeltaR(jet)< evt.DRJetLep[jetsAK8.it]) {
          evt.DRJetLep[jetsAK8.it] = goodleps[i].DeltaR(jet);
          evt.RelPtJetLep[jetsAK8.it] = goodleps[i].Perp(jet.Vect());
          lep = goodleps[i];
        }
      }
      if (evt.DRJetLep[jetsAK8.it]<1.0) {
        evt.nleptops++;
        is_top = true;
        TLorentzVector lepjet = lep + jet;
        leptop.push_back(lepjet);
      }
      // Extra - all except above hadronic/leptonic tops
      evt.HT += jetsAK8.Pt[NJET];
      if (is_top) evt.HTtt += jetsAK8.Pt[NJET];
    }
    evt.HTall = evt.HT + met.Pt + evt.HTlep;
    //std::cout<<evt.HTall<<" "<<evt.HT<<" "<<met.Pt<<" "<<evt.HTlep<<std::endl;
    
    evt.ntops = evt.nhadtops + evt.nleptops;
    // Select 2 tops (hadronic or hadronic like)
    TLorentzVector top1;
    TLorentzVector top2;
    if (evt.nhadtoplike >= 2) {
      if (evt.nhadtops == 2) {
	top1 = hadtop[0];
	top2 = hadtop[1];
      } else if (evt.nhadtops == 1) {
        if (hadtop[0].Pt() > hadtoplike[0].Pt()) {
          top1 = hadtop[0];
          top2 = hadtoplike[0];
        } else {
          top1 = hadtoplike[0];
          top2 = hadtop[0];
        }    
      } else if (evt.nhadtops == 0) {
	top1 = hadtoplike[0];
	top2 = hadtoplike[1];
      }
    }
    
    // top pair variables
    evt.tt_dR=NOVAL_F;
    evt.tt_dPhi=NOVAL_F;
    evt.tt_dEta=NOVAL_F;
    evt.tt_Mass=NOVAL_F;
    evt.tt_MR=NOVAL_F;
    evt.tt_MTR=NOVAL_F;
    evt.tt_R=NOVAL_F;
    evt.tt_R2=NOVAL_F;
    evt.HTtt=NOVAL_F;
    evt.HTex=NOVAL_F;
    evt.HTttFraction=NOVAL_F;
    evt.HTexFraction=NOVAL_F;
    //std::cout<<evt.ntops<<std::endl;
    if (evt.nhadtoplike>=2) {
      /* python
         tt_dR[0] = top1.DeltaR(top2)
         tt_dPhi[0] = top1.DeltaPhi(top2)
         tt_dEta[0] = fabs(top1.Eta() - top2.Eta())
         tt_mtt[0] = (top1 + top2).M()
         
         tt_extra = top1 + top2
         for jet in jets:
           tt_extra += jet
         pz_tt_extra[0] = tt_extra.Pz()
         dHt[0] = math.fabs(top1.Pt() - top2.Pt())
         dphi1 = delta_phi(metphi[0], top1.Phi())�w
         dphi2 = delta_phi(metphi[0], top2.Phi())
         dPhi_met_t1[0] = max(dphi1, dphi2)
         dPhi_met_t2[0] = min(dphi1, dphi2)
      */
      evt.tt_dR = top1.DeltaR(top2);
      evt.tt_dPhi = top1.DeltaPhi(top2);
      evt.tt_dEta = fabs(top1.Eta() - top2.Eta());
      evt.tt_Mass = (top1 + top2).M();
      evt.tt_Pz = (top1 + top2).Pz();
      evt.tt_Hz = top1.Pz() + top2.Pz();
      evt.tt_dPz = fabs(top1.Pz() - top2.Pz());
      evt.HTtt = top1.Pt() + top2.Pt();
      evt.HTex = evt.HTall - evt.HTtt;
      evt.HTttFraction = evt.HTtt / evt.HTall;
      evt.HTexFraction = evt.HTex / evt.HTall;
      // Select signal region
      
      // Razor for hadronic top pair
      TVector3 metl;
      metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
      evt.tt_MR = CalcMR_(top1, top2);
      evt.tt_MTR = CalcMTR_(top1, top2, metl);
      evt.tt_R = evt.tt_MTR / evt.tt_MR;
      evt.tt_R2 = pow(evt.tt_R, 2);
      //if (evt.nhadtops == 2)
      //  std::cout<<"TT : pt1="<<top1.Pt()<<" pt2="<<top2.Pt()<<" MET="<<met.Pt<<" MR="<<evt.tt_MR<<" MTR="<<evt.tt_MTR<<" R="<<evt.tt_R<<std::endl;
	
      //std::cout<<evt.NTopHad<<" "<<evt.nhadtops<<std::endl;
      if (evt.NTopHad==2) {
      //  std::cout<<evt.NLep<<" "<<ngoodleptons<<std::endl;
      //  std::cout<<evt.NTopHad<<" "<<evt.nhadtops<<std::endl;
      //  std::cout<<evt.NTopLep<<" "<<evt.nleptops<<std::endl;
      //  std::cout<<evt.NTop<<" "<<evt.ntops<<std::endl;
	if (evt.HtLep!=evt.HTlep) 
	  std::cout<<"evt.HtLep "<<evt.HtLep<<" "<<evt.HTlep<<std::endl;
	//if (evt.HtAll!=evt.HTall) {
	//  std::cout<<"evt.HtAll "<<evt.HtAll<<" "<<evt.HTall<<std::endl;
	//  std::cout<<evt.Ht<<" "<<evt.HT<<std::endl;
	//  std::cout<<evt.HtLep<<" "<<evt.HTlep<<std::endl;
	//  std::cout<<met.Pt<<std::endl;
	//  std::cout<<evt.HtEx<<" "<<evt.HTex<<std::endl;
	//  std::cout<<evt.HtTop<<" "<<evt.HTtt<<std::endl;
	//}
      //  std::cout<<evt.HtExFr<<" "<<evt.HTexFraction<<std::endl;
      //  std::cout<<evt.HtTopFr<<" "<<evt.HTttFraction<<std::endl;
      //  std::cout<<evt.TTHadDR<<" "<<evt.tt_dR<<std::endl;
	if (evt.TTHadDPhi!=evt.tt_dPhi)
	  std::cout<<"evt.TTHadDPhi "<<evt.TTHadDPhi<<" "<<evt.tt_dPhi<<std::endl;
      //  std::cout<<evt.TTHadDEta<<" "<<evt.tt_dEta<<std::endl;
      //  std::cout<<evt.TTHadMass<<" "<<evt.tt_Mass<<std::endl;
      //  std::cout<<evt.TTHadPz<<" "<<evt.tt_Pz<<std::endl;
      //  std::cout<<evt.TTHadHz<<" "<<evt.tt_Hz<<std::endl;
      //  std::cout<<evt.TTHadDPz<<" "<<evt.tt_dPz<<std::endl;
      //  std::cout<<evt.TTHadMR<<" "<<evt.tt_MR<<std::endl;
      //  std::cout<<evt.TTHadMTR<<" "<<evt.tt_MTR<<std::endl;
        if(evt.TTHadR!=evt.tt_R)
	  std::cout<<"evt.TTHadR "<<evt.TTHadR<<" "<<evt.tt_R<<std::endl;
      //  std::cout<<evt.TTHadR2<<" "<<evt.tt_R2<<std::endl<<std::endl;
      }
    }
  }
    
  private:
  // Razor recipe taken from the RazorBoost gurus: N. Strobbe, S. Sekmen
  //   https://github.com/nstrobbe/RazorBoost/blob/master/analyzer/utils.h
  
  // Hemispheres:
  vector<TLorentzVector> CombineJets_(vector<TLorentzVector> myjets) {
    //std::cout<<"Start CombineJets with "<<myjets.size()<<" jets\n";
    vector<TLorentzVector> mynewjets;
    TLorentzVector j1, j2;
    //bool foundGood = false;
    int N_comb = 1;
    for(unsigned int i = 0; i < myjets.size(); i++){
      N_comb *= 2;
    }
    //std::cout<<"N_comb = "<<N_comb<<std::endl;
    double M_min = 9999999999.0;
    int j_count;
    for(int i = 1; i < N_comb-1; i++){
      TLorentzVector j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;
      while(j_count > 0){
        if(itemp/j_count == 1){
          j_temp1 += myjets[count];
	  //std::cout<<"  1 <- "<<count<<" M2="<<myjets[count].M2()<<std::endl;
        } else {
          j_temp2 += myjets[count];
	  //std::cout<<"  2 <- "<<count<<" M2="<<myjets[count].M2()<<std::endl;
        }
        itemp -= j_count*(itemp/j_count);
        j_count /= 2;
        count++;
      }
      double M_temp = j_temp1.M2()+j_temp2.M2();
      //std::cout<<"  --> M_temp "<<j_temp1.M2()<<" + "<<j_temp2.M2()<<" = "<<M_temp<<std::endl;
      // smallest mass
      if(M_temp < M_min){
        M_min = M_temp;
        j1 = j_temp1;
        j2 = j_temp2;
      }
      //std::cout<<" M_min = "<<M_min<<std::endl;
    }
    if(j2.Pt() > j1.Pt()){
      TLorentzVector temp = j1;
      j1 = j2;
      j2 = temp;
    }
    //std::cout<<"Result: Jet1 pT = "<<j1.Pt()<<" Jet2 pT = "<<j2.Pt()<<"\n\n";
    mynewjets.push_back(j1);
    mynewjets.push_back(j2);
    return mynewjets;
  }
  
  // MR
  double CalcMR_(TLorentzVector ja, TLorentzVector jb){
    double A = ja.P();
    double B = jb.P();
    double az = ja.Pz();
    double bz = jb.Pz();
    TVector3 jaT, jbT;
    jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
    jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
    double ATBT = (jaT+jbT).Mag2();
    double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
      		 (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
    double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
    double mygamma = 1./sqrt(1.-mybeta*mybeta);
    //gamma times MRstar
    temp *= mygamma;
    return temp;
  }
  
  // MTR
  double CalcMTR_(TLorentzVector ja, TLorentzVector jb, TVector3 met){
    double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
    temp /= 2.;
    temp = sqrt(temp);
    return temp;
  }
  
  // MT
  double CalcMT_(TLorentzVector lepton, TLorentzVector pfmet){
    return sqrt( 2 * lepton.Pt() * pfmet.Pt() * ( 1 - cos( pfmet.Phi() - lepton.Phi() ) ) );
  }
  
  // Calculate Razor variables here
  void calcRazorAK4_() {

    // Select the best pair of jets (AK4, pt>40, |eta| < 3.0)
    std::vector<TLorentzVector> sjetl;
    while(jetsAK4.Loop()) {
      if (!(jetsAK4.Pt[NJET] > 40) ) continue;
      if (!(fabs(jetsAK4.Eta[NJET]) < 3) ) continue;
      TLorentzVector jl;
      jl.SetPtEtaPhiE(jetsAK4.Pt[NJET], jetsAK4.Eta[NJET],
      		jetsAK4.Phi[NJET], jetsAK4.E[NJET]);
      sjetl.push_back(jl);
    }
    std::vector<TLorentzVector> hemis = CombineJets_(sjetl);
    
    // ---------------------
    // -- Razor variables --
    // ---------------------
    
    TVector3 metl;
    metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
    
    if (hemis.size() == 2) {
      jetsAK4.MR = CalcMR_(hemis[0], hemis[1]);
      jetsAK4.MTR = CalcMTR_(hemis[0], hemis[1], metl);
      jetsAK4.R = jetsAK4.MTR / jetsAK4.MR;
      jetsAK4.R2 = pow(jetsAK4.R, 2);
    }
  }
  
  void calcRazorAK8_() {
    
    // Select the best pair of jets (AK8, pt>40, |eta| < 3.0)
    std::vector<TLorentzVector> sjetl;
    while(jetsAK8.Loop()) {
      if (!(jetsAK8.Pt[NJET] > 40) ) continue;
      if (!(fabs(jetsAK8.Eta[NJET]) < 3) ) continue;
      TLorentzVector jl;
      jl.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET],
      		jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
      sjetl.push_back(jl);
    }
    std::vector<TLorentzVector> hemis = CombineJets_(sjetl);
    
    // ---------------------
    // -- Razor variables --
    // ---------------------
    
    TVector3 metl;
    metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
    
    if (hemis.size() == 2) {
      jetsAK8.MR = CalcMR_(hemis[0], hemis[1]);
      jetsAK8.MTR = CalcMTR_(hemis[0], hemis[1], metl);
      jetsAK8.R = jetsAK8.MTR / jetsAK8.MR;
      jetsAK8.R2 = pow(jetsAK8.R, 2);
    }
    //if (evt.NTopHad==2)
    //  std::cout<<"AK8: pt1="<<hemis[0].Pt()<<" pt2="<<hemis[1].Pt()<<" MET="<<met.Pt<<" MR="<<jetsAK8.MR<<" MTR="<<jetsAK8.MTR<<" R="<<jetsAK8.R<<std::endl;
  }
  
};

#endif
