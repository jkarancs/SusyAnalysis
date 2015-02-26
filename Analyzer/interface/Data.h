#ifndef Data_h
#define Data_h

#define NOVAL_I -9999
#define NOVAL_F -9999.0

#define NLEP 20
#define NJET 50

#include <vector>
#include "TLorentzVector.h"
#include <iostream>

using namespace std;

class Data {
public:
  Data() {}
  ~Data() {}
  
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
    float subjetIndex0[NJET+1];
    float subjetIndex1[NJET+1];
    float subjetIndex2[NJET+1];
    float subjetIndex3[NJET+1];
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
	subjetIndex0[it]=NOVAL_F;
	subjetIndex1[it]=NOVAL_F;
	subjetIndex2[it]=NOVAL_F;
	subjetIndex3[it]=NOVAL_F;
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
      subjetIndex0[NJET] = subjetIndex0[it];
      subjetIndex1[NJET] = subjetIndex1[it];
      subjetIndex2[NJET] = subjetIndex2[it];
      subjetIndex3[NJET] = subjetIndex3[it];
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
  
  class CmsTopTagJetVars : public JetVars, public AK8Vars {
  public:
    CmsTopTagJetVars() { init(); }
    ~CmsTopTagJetVars() {}
    
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
  } jetsCmsTopTag;

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
    
    // Top tagging variables
    float nhadtops;
    float nleptops;
    float ntops;
    float ncmshadtops;
    float ncmsleptops;
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
    
    void init() {
      nhadtops=NOVAL_F;
      nleptops=NOVAL_F;
      ntops=NOVAL_F;
      ncmshadtops=NOVAL_F;
      ncmsleptops=NOVAL_F;
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
    }
  } evt;
  
  void CalculateAllVariables() {
    //calcRazorAK4_();
    //calcRazorAK8_();
    //calcRazorCmsTopTag_();
    
    // find good leptons (for letponic tops)
    int ngoodleptons = 0;
    std::vector<TLorentzVector> goodleps;
    evt.HTlep = 0;
    while(ele.Loop()) if (ele.Pt[NLEP] > 35 && fabs(ele.Eta[NLEP]) < 2.5) {
    //while(ele.Loop()) if (ele.isLoose[NLEP]>0 && ele.Pt[NLEP] > 35 && fabs(ele.Eta[NLEP]) < 2.5) {
      ngoodleptons++;
      TLorentzVector goodele;
      goodele.SetPtEtaPhiE(ele.Pt[NLEP], ele.Eta[NLEP], ele.Phi[NLEP], ele.E[NLEP]);
      goodleps.push_back(goodele);
      evt.HTlep += ele.Pt[NLEP];
    }
    while(mu.Loop())  if (mu.IsTightMuon[NLEP]>0 && mu.Pt[NLEP] > 45 && fabs(mu.Eta[NLEP]) < 2.1) {
    //while(mu.Loop())  if (mu.IsLooseMuon[NLEP]>0 && mu.Pt[NLEP] > 45 && fabs(mu.Eta[NLEP]) < 2.1) {
      ngoodleptons++;
      TLorentzVector goodmu;
      goodmu.SetPtEtaPhiE(mu.Pt[NLEP], mu.Eta[NLEP], mu.Phi[NLEP], mu.E[NLEP]);
      goodleps.push_back(goodmu);
      evt.HTlep += mu.Pt[NLEP];
    }
    
    // Tag hadronic tops
    TLorentzVector hadtop1;
    TLorentzVector hadtop2;
    TLorentzVector leptop1;
    TLorentzVector leptop2;
    evt.ntops = 0;
    evt.nhadtops = 0;
    evt.nleptops = 0;
    evt.HT = 0;
    while(jetsAK8.Loop()) {
      // Find subjets
      //int nsubjet = 0;
      //std::cout<<jetsAK8.it<<" "<<jetsAK8.size<<" "<<subjetsAK8.size<<" "<<jetsAK8.subjetIndex0[NJET]<<" "<<jetsAK8.subjetIndex1[NJET]<<std::endl;
      bool is_top = false;
      // hadronic tops
      //if (jetsAK8.tau1[NJET]>0 && jetsAK8.tau2[NJET]>0 ? jetsAK8.Pt[NJET] > 400 && jetsAK8.prunedMass[NJET] > 140 && (jetsAK8.tau2[NJET]/jetsAK8.tau1[NJET]) < 0.75 : 0) { // orig
      if (jetsAK8.tau2[NJET]>0 && jetsAK8.tau3[NJET]>0 ? jetsAK8.Pt[NJET] > 400 && jetsAK8.prunedMass[NJET] > 140 && (jetsAK8.tau3[NJET]/jetsAK8.tau2[NJET]) < 0.75 : 0) { // Brandon's
      //if (jetsAK8.tau2[NJET]>0 && jetsAK8.tau3[NJET]>0 ? (jetsAK8.tau3[NJET]/jetsAK8.tau2[NJET]) < 0.75 &&
      //    jetsAK8.nSubJets[NJET] > 2 &&
      //    jetsAK8.minmass[NJET] > 50 &&
      //    jetsAK8.Pt[NJET] > 400 &&
      //    jetsAK8.prunedMass[NJET] > 140  : 0) { // New
        ++evt.nhadtops;
        is_top = true;
        if (evt.nhadtops==1) hadtop1.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
        if (evt.nhadtops==2) hadtop2.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
      }
      // leptonic tops
      TLorentzVector lep;
      TLorentzVector jet;
      jet.SetPtEtaPhiE(jetsAK8.Pt[NJET], jetsAK8.Eta[NJET], jetsAK8.Phi[NJET], jetsAK8.E[NJET]);
      double DeltaR_lep = 9999;
      for (size_t i=0; i<goodleps.size(); ++i) {
        if (goodleps[i].DeltaR(jet)< DeltaR_lep) {
	  DeltaR_lep = goodleps[i].DeltaR(jet);
	  lep = goodleps[i];
        }
      }
      if (DeltaR_lep<1.0) {
        evt.nleptops++;
        is_top = true;
        TLorentzVector leptop = lep + jet;
        if (evt.nleptops==1) leptop1 = leptop;
        if (evt.nleptops==2) leptop2 = leptop;
      }
      // Extra - all except above hadronic/leptonic tops
      evt.HT += jetsAK8.Pt[NJET];
      if (is_top) evt.HTtt += jetsAK8.Pt[NJET];
    }
    evt.HTall = evt.HT + met.Pt + evt.HTlep;
    //std::cout<<evt.HTall<<" "<<evt.HT<<" "<<met.Pt<<" "<<evt.HTlep<<std::endl;
    
    // CMS Tagged tops
    evt.ncmshadtops = 0;
    evt.ncmsleptops = 0;
    while(jetsCmsTopTag.Loop()) {
      bool is_top = false;
      // leptonic tops
      TLorentzVector lep;
      TLorentzVector cmstopjet;
      cmstopjet.SetPtEtaPhiE(jetsCmsTopTag.Pt[NJET], jetsCmsTopTag.Eta[NJET], jetsCmsTopTag.Phi[NJET], jetsCmsTopTag.E[NJET]);
      double DeltaR_lep = 9999;
      for (size_t i=0; i<goodleps.size(); ++i) if (goodleps[i].DeltaR(cmstopjet)< DeltaR_lep) {
	DeltaR_lep = goodleps[i].DeltaR(cmstopjet);
	lep = goodleps[i];
      }
      if (DeltaR_lep<1.0) {
        evt.ncmsleptops++;
        is_top = true;
        TLorentzVector leptop = lep + cmstopjet;
        //if (evt.nleptops==1) leptop1 = leptop;
        //if (evt.nleptops==2) leptop2 = leptop;
      }
      // hadronic tops
      //if (jetsCmsTopTag.Pt[NJET] > 400 && jetsCmsTopTag.Mass[NJET] > 140 ) {
      if (!is_top && jetsCmsTopTag.Pt[NJET] > 400) {
        ++evt.ncmshadtops;
        is_top = true;
        //if (evt.nhadtops==1) hadtop1.SetPtEtaPhiE(jetsCmsTopTag.Pt[NJET], jetsCmsTopTag.Eta[NJET], jetsCmsTopTag.Phi[NJET], jetsCmsTopTag.E[NJET]);
        //if (evt.nhadtops==2) hadtop2.SetPtEtaPhiE(jetsCmsTopTag.Pt[NJET], jetsCmsTopTag.Eta[NJET], jetsCmsTopTag.Phi[NJET], jetsCmsTopTag.E[NJET]);
      }
      // Extra - all except above hadronic/leptonic tops
      //evt.HT += jetsCmsTopTag.Pt[NJET];
      //if (is_top) evt.HTtt += jetsCmsTopTag.Pt[NJET];
    }
    
    // Select exactly 2 tops (hadronic or leptonic)
    // We need exactly 2 in order to calculate pair variables, eg. DeltaPhi
    TLorentzVector top1;
    TLorentzVector top2;
    evt.ntops = evt.nhadtops + evt.nleptops;
    if (evt.ntops == 2) {
      if (evt.nhadtops == 2) {
        if (hadtop1.Pt() > hadtop2.Pt()) {
          top1 = hadtop1;
          top2 = hadtop2;
        } else {
          top1 = hadtop2;
          top2 = hadtop1;
        }
      } else if (evt.nhadtops == 1) {
        if (hadtop1.Pt() > leptop1.Pt()) {
          top1 = hadtop1;
          top2 = leptop1;
        } else {
          top1 = leptop1;
          top2 = hadtop1;
        }    
      } else if (evt.nhadtops == 0) {
        if (leptop1.Pt() > leptop2.Pt()) {
          top1 = leptop1;
          top2 = leptop2;
        } else {
          top1 = leptop2;
          top2 = leptop1;
        }
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
    if (evt.ntops==2) {
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
         dphi1 = delta_phi(metphi[0], top1.Phi())íw
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
    }
  }
    
  private:
  // Razor recipe taken from the RazorBoost gurus: N. Strobbe, S. Sekmen
  //   https://github.com/nstrobbe/RazorBoost/blob/master/analyzer/utils.h
  
  // Hemispheres:
  vector<TLorentzVector> CombineJets_(vector<TLorentzVector> myjets) {
    vector<TLorentzVector> mynewjets;
    TLorentzVector j1, j2;
    //bool foundGood = false;
    int N_comb = 1;
    for(unsigned int i = 0; i < myjets.size(); i++){
      N_comb *= 2;
    }
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
        } else {
          j_temp2 += myjets[count];
        }
        itemp -= j_count*(itemp/j_count);
        j_count /= 2;
        count++;
      }
      double M_temp = j_temp1.M2()+j_temp2.M2();
      // smallest mass
      if(M_temp < M_min){
        M_min = M_temp;
        j1 = j_temp1;
        j2 = j_temp2;
      }
    }
    if(j2.Pt() > j1.Pt()){
      TLorentzVector temp = j1;
      j1 = j2;
      j2 = temp;
    }
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
  }
  
  void calcRazorCmsTopTag_() {
    
    // Select the best pair of jets (CmsTopTag, pt>40, |eta| < 3.0)
    std::vector<TLorentzVector> sjetl;
    while(jetsCmsTopTag.Loop()) {
      if (!(jetsCmsTopTag.Pt[NJET] > 40) ) continue;
      if (!(fabs(jetsCmsTopTag.Eta[NJET]) < 3) ) continue;
      TLorentzVector jl;
      jl.SetPtEtaPhiE(jetsCmsTopTag.Pt[NJET], jetsCmsTopTag.Eta[NJET],
		      jetsCmsTopTag.Phi[NJET], jetsCmsTopTag.E[NJET]);
      sjetl.push_back(jl);
    }
    std::vector<TLorentzVector> hemis = CombineJets_(sjetl);
    
    // ---------------------
    // -- Razor variables --
    // ---------------------
    
    TVector3 metl;
    metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
    
    if (hemis.size() == 2) {
      jetsCmsTopTag.MR = CalcMR_(hemis[0], hemis[1]);
      jetsCmsTopTag.MTR = CalcMTR_(hemis[0], hemis[1], metl);
      jetsCmsTopTag.R = jetsCmsTopTag.MTR / jetsCmsTopTag.MR;
      jetsCmsTopTag.R2 = pow(jetsCmsTopTag.R, 2);
    }
  }
  
};

#endif
