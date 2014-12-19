#ifndef Data_h
#define Data_h

#define NOVAL_F -9999.0

class Data {
public:
  Data() {};
  ~Data() {};
  
  class ElectronData {
  public:
    ElectronData() { init(); };
    
    int electrons_size;
    float elE[5];
    float elPt[5];
    float elMass[5];
    float elEta[5];
    float elPhi[5];
    float elCharge[5];
    float elD0[5];
    float elDz[5];
    float elHoE[5];
    float elIso03[5];
    float elY[5];
    float eldEtaIn[5];
    float eldPhiIn[5];
    float elexpectedMissInHits[5];
    float elfull5x5siee[5];
    float elooEmooP[5];
    float elpssConVeto[5];
    float electronsSFTrigger[5];
    float electronsSFReco[5];
    float electronsisQCD[5];
    float electronsisTightOffline[5];
    float electronsisLooseOffline[5];
    
    float E;
    float Pt;
    float Mass;
    float Eta;
    float Phi;
    float Charge;
    float D0;
    float Dz;
    float HoE;
    float Iso03;
    float Y;
    float dEtaIn;
    float dPhiIn;
    float expectedMissInHits;
    float full5x5siee;
    float ooEmooP;
    float pssConVeto;
    float SFTrigger;
    float SFReco;
    float isQCD;
    float isTightOffline;
    float isLooseOffline;

    void init() {
      electrons_size = 0;
      for (int i=0; i<5; ++i) {
	elE[i]=NOVAL_F;
        elPt[i]=NOVAL_F;
        elMass[i]=NOVAL_F;
        elEta[i]=NOVAL_F;
        elPhi[i]=NOVAL_F;
        elCharge[i]=NOVAL_F;
        elD0[i]=NOVAL_F;
        elDz[i]=NOVAL_F;
        elHoE[i]=NOVAL_F;
        elIso03[i]=NOVAL_F;
        elY[i]=NOVAL_F;
        eldEtaIn[i]=NOVAL_F;
        eldPhiIn[i]=NOVAL_F;
        elexpectedMissInHits[i]=NOVAL_F;
        elfull5x5siee[i]=NOVAL_F;
        elooEmooP[i]=NOVAL_F;
        elpssConVeto[i]=NOVAL_F;
        electronsSFTrigger[i]=NOVAL_F;
        electronsSFReco[i]=NOVAL_F;
        electronsisQCD[i]=NOVAL_F;
        electronsisTightOffline[i]=NOVAL_F;
        electronsisLooseOffline[i]=NOVAL_F;
      }
      E=NOVAL_F;
      Pt=NOVAL_F;
      Mass=NOVAL_F;
      Eta=NOVAL_F;
      Phi=NOVAL_F;
      Charge=NOVAL_F;
      D0=NOVAL_F;
      Dz=NOVAL_F;
      HoE=NOVAL_F;
      Iso03=NOVAL_F;
      Y=NOVAL_F;
      dEtaIn=NOVAL_F;
      dPhiIn=NOVAL_F;
      expectedMissInHits=NOVAL_F;
      full5x5siee=NOVAL_F;
      ooEmooP=NOVAL_F;
      pssConVeto=NOVAL_F;
      SFTrigger=NOVAL_F;
      SFReco=NOVAL_F;
      isQCD=NOVAL_F;
      isTightOffline=NOVAL_F;
      isLooseOffline=NOVAL_F;
      
      it_ = 0;
    }

    bool Loop() {
      if (it_<electrons_size) {
	E                  = elE[it_];			   
	Pt		   = elPt[it_];		   
        Mass		   = elMass[it_];		   
        Eta		   = elEta[it_];		   
        Phi		   = elPhi[it_];		   
        Charge		   = elCharge[it_];		   
        D0		   = elD0[it_];		   
        Dz		   = elDz[it_];		   
        HoE		   = elHoE[it_];		   
        Iso03		   = elIso03[it_];		   
        Y		   = elY[it_];			   
        dEtaIn		   = eldEtaIn[it_];		   
        dPhiIn		   = eldPhiIn[it_];		   
        expectedMissInHits = elexpectedMissInHits[it_];   
        full5x5siee	   = elfull5x5siee[it_];	   
        ooEmooP		   = elooEmooP[it_];		   
        pssConVeto	   = elpssConVeto[it_];	   
        SFTrigger	   = electronsSFTrigger[it_];	   
        SFReco		   = electronsSFReco[it_];	   
        isQCD		   = electronsisQCD[it_];	   
        isTightOffline	   = electronsisTightOffline[it_];
        isLooseOffline	   = electronsisLooseOffline[it_];
	
	++it_;
	return 1;
      } else {
	it_=0;
	return 0;
      }
    }

  private:
    int it_;
    
  } ele;
  
  class MuonData {
  public:
    MuonData() { init(); };
    
    int muons_size;
    float muE[5];
    float muPt[5];
    float muMass[5];
    float muEta[5];
    float muPhi[5];
    float muCharge[5];
    float muIsLooseMuon[5];
    float muIsSoftMuon[5];
    float muIsTightMuon[5];
    float muD0[5];
    float muD0err[5];
    float muDz[5];
    float muDzerr[5];
    float muGenMuonCharge[5];
    float muGenMuonEta[5];
    float muGenMuonPt[5];
    float muGenMuonE[5];
    float muGenMuonPhi[5];
    float muGenMuonY[5];
    float muGlbTrkNormChi2[5];
    float muHLTmuonDeltaR[5];
    float muHLTmuonE[5];
    float muHLTmuonEta[5];
    float muHLTmuonPt[5];
    float muHLTmuonPhi[5];
    float muInTrkNormChi2[5];
    float muIsGlobalMuon[5];
    float muIsPFMuon[5];
    float muIsTrackerMuon[5];
    float muIso03[5];
    float muNumberMatchedStations[5];
    float muNumberOfPixelLayers[5];
    float muNumberOfValidTrackerHits[5];
    float muNumberTrackerLayers[5];
    float muNumberValidMuonHits[5];
    float muNumberValidPixelHits[5];
    float muSumChargedHadronPt[5];
    float muSumNeutralHadronPt[5];
    float muSumPUPt[5];
    float muSumPhotonPt[5];
    float muY[5];
    float muonsSFTrigger[5];
    float muonsSFReco[5];
    float muonsisQCD[5];
    float muonsisTightOffline[5];
    float muonsisLooseOffline[5];
    
    float E;
    float Pt;
    float Mass;
    float Eta;
    float Phi;
    float Charge;
    float IsLooseMuon;
    float IsSoftMuon;
    float IsTightMuon;
    float D0;
    float D0err;
    float Dz;
    float Dzerr;
    float GenMuonCharge;
    float GenMuonEta;
    float GenMuonPt;
    float GenMuonE;
    float GenMuonPhi;
    float GenMuonY;
    float GlbTrkNormChi2;
    float HLTmuonDeltaR;
    float HLTmuonE;
    float HLTmuonEta;
    float HLTmuonPt;
    float HLTmuonPhi;
    float InTrkNormChi2;
    float IsGlobalMuon;
    float IsPFMuon;
    float IsTrackerMuon;
    float Iso03;
    float NumberMatchedStations;
    float NumberOfPixelLayers;
    float NumberOfValidTrackerHits;
    float NumberTrackerLayers;
    float NumberValidMuonHits;
    float NumberValidPixelHits;
    float SumChargedHadronPt;
    float SumNeutralHadronPt;
    float SumPUPt;
    float SumPhotonPt;
    float Y;
    float SFTrigger;
    float SFReco;
    float isQCD;
    float isTightOffline;
    float isLooseOffline;
    
    void init() {
      muons_size = 0;
      for (int i=0; i<5; ++i) {
        muE[i]=NOVAL_F;
        muPt[i]=NOVAL_F;
        muMass[i]=NOVAL_F;
        muEta[i]=NOVAL_F;
        muPhi[i]=NOVAL_F;
        muCharge[i]=NOVAL_F;
        muIsLooseMuon[i]=NOVAL_F;
        muIsSoftMuon[i]=NOVAL_F;
        muIsTightMuon[i]=NOVAL_F;
        muD0[i]=NOVAL_F;
        muD0err[i]=NOVAL_F;
        muDz[i]=NOVAL_F;
        muDzerr[i]=NOVAL_F;
        muGenMuonCharge[i]=NOVAL_F;
        muGenMuonEta[i]=NOVAL_F;
        muGenMuonPt[i]=NOVAL_F;
        muGenMuonE[i]=NOVAL_F;
        muGenMuonPhi[i]=NOVAL_F;
        muGenMuonY[i]=NOVAL_F;
        muGlbTrkNormChi2[i]=NOVAL_F;
        muHLTmuonDeltaR[i]=NOVAL_F;
        muHLTmuonE[i]=NOVAL_F;
        muHLTmuonEta[i]=NOVAL_F;
        muHLTmuonPt[i]=NOVAL_F;
        muHLTmuonPhi[i]=NOVAL_F;
        muInTrkNormChi2[i]=NOVAL_F;
        muIsGlobalMuon[i]=NOVAL_F;
        muIsPFMuon[i]=NOVAL_F;
        muIsTrackerMuon[i]=NOVAL_F;
        muIso03[i]=NOVAL_F;
        muNumberMatchedStations[i]=NOVAL_F;
        muNumberOfPixelLayers[i]=NOVAL_F;
        muNumberOfValidTrackerHits[i]=NOVAL_F;
        muNumberTrackerLayers[i]=NOVAL_F;
        muNumberValidMuonHits[i]=NOVAL_F;
        muNumberValidPixelHits[i]=NOVAL_F;
        muSumChargedHadronPt[i]=NOVAL_F;
        muSumNeutralHadronPt[i]=NOVAL_F;
        muSumPUPt[i]=NOVAL_F;
        muSumPhotonPt[i]=NOVAL_F;
        muY[i]=NOVAL_F;
        muonsSFTrigger[i]=NOVAL_F;
        muonsSFReco[i]=NOVAL_F;
        muonsisQCD[i]=NOVAL_F;
        muonsisTightOffline[i]=NOVAL_F;
        muonsisLooseOffline[i]=NOVAL_F;
      }
      E=NOVAL_F;
      Pt=NOVAL_F;
      Mass=NOVAL_F;
      Eta=NOVAL_F;
      Phi=NOVAL_F;
      Charge=NOVAL_F;
      IsLooseMuon=NOVAL_F;
      IsSoftMuon=NOVAL_F;
      IsTightMuon=NOVAL_F;
      D0=NOVAL_F;
      D0err=NOVAL_F;
      Dz=NOVAL_F;
      Dzerr=NOVAL_F;
      GenMuonCharge=NOVAL_F;
      GenMuonEta=NOVAL_F;
      GenMuonPt=NOVAL_F;
      GenMuonE=NOVAL_F;
      GenMuonPhi=NOVAL_F;
      GenMuonY=NOVAL_F;
      GlbTrkNormChi2=NOVAL_F;
      HLTmuonDeltaR=NOVAL_F;
      HLTmuonE=NOVAL_F;
      HLTmuonEta=NOVAL_F;
      HLTmuonPt=NOVAL_F;
      HLTmuonPhi=NOVAL_F;
      InTrkNormChi2=NOVAL_F;
      IsGlobalMuon=NOVAL_F;
      IsPFMuon=NOVAL_F;
      IsTrackerMuon=NOVAL_F;
      Iso03=NOVAL_F;
      NumberMatchedStations=NOVAL_F;
      NumberOfPixelLayers=NOVAL_F;
      NumberOfValidTrackerHits=NOVAL_F;
      NumberTrackerLayers=NOVAL_F;
      NumberValidMuonHits=NOVAL_F;
      NumberValidPixelHits=NOVAL_F;
      SumChargedHadronPt=NOVAL_F;
      SumNeutralHadronPt=NOVAL_F;
      SumPUPt=NOVAL_F;
      SumPhotonPt=NOVAL_F;
      Y=NOVAL_F;
      SFTrigger=NOVAL_F;
      SFReco=NOVAL_F;
      isQCD=NOVAL_F;
      isTightOffline=NOVAL_F;
      isLooseOffline=NOVAL_F;
      
      it_ = 0;
    }

    bool Loop() {
      if (it_<muons_size) {
	E                        = muE[it_];			  
        Pt			 = muPt[it_];			  
        Mass			 = muMass[it_];			  
        Eta			 = muEta[it_];			  
        Phi			 = muPhi[it_];			  
        Charge			 = muCharge[it_];			  
        IsLooseMuon		 = muIsLooseMuon[it_];		  
        IsSoftMuon		 = muIsSoftMuon[it_];		  
        IsTightMuon		 = muIsTightMuon[it_];		  
        D0			 = muD0[it_];			  
        D0err			 = muD0err[it_];			  
        Dz			 = muDz[it_];			  
        Dzerr			 = muDzerr[it_];			  
        GenMuonCharge		 = muGenMuonCharge[it_];		  
        GenMuonEta		 = muGenMuonEta[it_];		  
        GenMuonPt		 = muGenMuonPt[it_];		  
        GenMuonE		 = muGenMuonE[it_];		  
        GenMuonPhi		 = muGenMuonPhi[it_];		  
        GenMuonY		 = muGenMuonY[it_];		  
        GlbTrkNormChi2		 = muGlbTrkNormChi2[it_];		  
        HLTmuonDeltaR		 = muHLTmuonDeltaR[it_];		  
        HLTmuonE		 = muHLTmuonE[it_];		  
        HLTmuonEta		 = muHLTmuonEta[it_];		  
        HLTmuonPt		 = muHLTmuonPt[it_];		  
        HLTmuonPhi		 = muHLTmuonPhi[it_];		  
        InTrkNormChi2		 = muInTrkNormChi2[it_];		  
        IsGlobalMuon		 = muIsGlobalMuon[it_];		  
        IsPFMuon		 = muIsPFMuon[it_];		  
        IsTrackerMuon		 = muIsTrackerMuon[it_];		  
        Iso03			 = muIso03[it_];			  
        NumberMatchedStations	 = muNumberMatchedStations[it_];	  
        NumberOfPixelLayers	 = muNumberOfPixelLayers[it_];	  
        NumberOfValidTrackerHits = muNumberOfValidTrackerHits[it_]; 
        NumberTrackerLayers	 = muNumberTrackerLayers[it_];	  
        NumberValidMuonHits	 = muNumberValidMuonHits[it_];	  
        NumberValidPixelHits	 = muNumberValidPixelHits[it_];	  
        SumChargedHadronPt	 = muSumChargedHadronPt[it_];	  
        SumNeutralHadronPt	 = muSumNeutralHadronPt[it_];	  
        SumPUPt			 = muSumPUPt[it_];		  
        SumPhotonPt		 = muSumPhotonPt[it_];		  
        Y			 = muY[it_];			  
        SFTrigger		 = muonsSFTrigger[it_];		  
        SFReco			 = muonsSFReco[it_];		  
        isQCD			 = muonsisQCD[it_];		  
        isTightOffline		 = muonsisTightOffline[it_];	  
        isLooseOffline		 = muonsisLooseOffline[it_];        
	
	++it_;
	return 1;
      } else {
	it_=0;
	return 0;
      }
    }

  private:
    int it_;
    
  } mu;
  
  class JetData {
  public:
    JetData() { init(); };
    
    int jets_size;
    float jetE[20];
    float jetPt[20];
    float jetMass[20];
    float jetEta[20];
    float jetPhi[20];
    float jetPartonFlavour[20];
    float jetCSV[20];
    float jetCSVV1[20];
    float jetCharge[20];
    float jetChargeMuEnergy[20];
    float jetChargedHadronMultiplicity[20];
    float jetElectronEnergy[20];
    float jetGenJetCharge[20];
    float jetGenJetE[20];
    float jetGenJetEta[20];
    float jetGenJetPhi[20];
    float jetGenJetPt[20];
    float jetGenJetY[20];
    float jetGenPartonCharge[20];
    float jetGenPartonE[20];
    float jetGenPartonEta[20];
    float jetGenPartonPhi[20];
    float jetGenPartonPt[20];
    float jetGenPartonY[20];
    float jetHFEMEnergy[20];
    float jetHFEMMultiplicity[20];
    float jetHFHadronEnergy[20];
    float jetHFHadronMultiplicity[20];
    float jetHLTjetDeltaR[20];
    float jetHLTjetE[20];
    float jetHLTjetEta[20];
    float jetHLTjetPt[20];
    float jetHLTjetPhi[20];
    float jetHadronFlavour[20];
    float jetIsCSVL[20];
    float jetIsCSVM[20];
    float jetIsCSVT[20];
    float jetSmearedE[20];
    float jetSmearedPt[20];
    float jetSmearedPEta[20];
    float jetSmearedPhi[20];
    float jetY[20];
    float jetelectronMultiplicity[20];
    float jetmuonMultiplicity[20];
    float jetneutralHadronMultiplicity[20];
    float jetneutralMultiplicity[20];
    float jetphotonMultiplicity[20];
    float jetsCorrPt[20];
    float jetsCorrEta[20];
    float jetsCorrPhi[20];
    float jetsCorrE[20];
    float jetsCorrMass[20];
    float jetsCorrNJets[20];
    float jetsCorrPartonFlavour[20];
    
    float E;
    float Pt;
    float Mass;
    float Eta;
    float Phi;
    float PartonFlavour;
    float CSV;
    float CSVV1;
    float Charge;
    float ChargeMuEnergy;
    float ChargedHadronMultiplicity;
    float ElectronEnergy;
    float GenJetCharge;
    float GenJetE;
    float GenJetEta;
    float GenJetPhi;
    float GenJetPt;
    float GenJetY;
    float GenPartonCharge;
    float GenPartonE;
    float GenPartonEta;
    float GenPartonPhi;
    float GenPartonPt;
    float GenPartonY;
    float HFEMEnergy;
    float HFEMMultiplicity;
    float HFHadronEnergy;
    float HFHadronMultiplicity;
    float HLTjetDeltaR;
    float HLTjetE;
    float HLTjetEta;
    float HLTjetPt;
    float HLTjetPhi;
    float HadronFlavour;
    float IsCSVL;
    float IsCSVM;
    float IsCSVT;
    float SmearedE;
    float SmearedPt;
    float SmearedPEta;
    float SmearedPhi;
    float Y;
    float electronMultiplicity;
    float muonMultiplicity;
    float neutralHadronMultiplicity;
    float neutralMultiplicity;
    float photonMultiplicity;
    float CorrPt;
    float CorrEta;
    float CorrPhi;
    float CorrE;
    float CorrMass;
    float CorrNJets;
    float CorrPartonFlavour;
    
    void init() {
      jets_size=0;
      for (int i=0; i<20; ++i) {
        jetE[i]=NOVAL_F;
        jetPt[i]=NOVAL_F;
        jetMass[i]=NOVAL_F;
        jetEta[i]=NOVAL_F;
        jetPhi[i]=NOVAL_F;
        jetPartonFlavour[i]=NOVAL_F;
        jetCSV[i]=NOVAL_F;
        jetCSVV1[i]=NOVAL_F;
        jetCharge[i]=NOVAL_F;
        jetChargeMuEnergy[i]=NOVAL_F;
        jetChargedHadronMultiplicity[i]=NOVAL_F;
        jetElectronEnergy[i]=NOVAL_F;
        jetGenJetCharge[i]=NOVAL_F;
        jetGenJetE[i]=NOVAL_F;
        jetGenJetEta[i]=NOVAL_F;
        jetGenJetPhi[i]=NOVAL_F;
        jetGenJetPt[i]=NOVAL_F;
        jetGenJetY[i]=NOVAL_F;
        jetGenPartonCharge[i]=NOVAL_F;
        jetGenPartonE[i]=NOVAL_F;
        jetGenPartonEta[i]=NOVAL_F;
        jetGenPartonPhi[i]=NOVAL_F;
        jetGenPartonPt[i]=NOVAL_F;
        jetGenPartonY[i]=NOVAL_F;
        jetHFEMEnergy[i]=NOVAL_F;
        jetHFEMMultiplicity[i]=NOVAL_F;
        jetHFHadronEnergy[i]=NOVAL_F;
        jetHFHadronMultiplicity[i]=NOVAL_F;
        jetHLTjetDeltaR[i]=NOVAL_F;
        jetHLTjetE[i]=NOVAL_F;
        jetHLTjetEta[i]=NOVAL_F;
        jetHLTjetPt[i]=NOVAL_F;
        jetHLTjetPhi[i]=NOVAL_F;
        jetHadronFlavour[i]=NOVAL_F;
        jetIsCSVL[i]=NOVAL_F;
        jetIsCSVM[i]=NOVAL_F;
        jetIsCSVT[i]=NOVAL_F;
        jetSmearedE[i]=NOVAL_F;
        jetSmearedPt[i]=NOVAL_F;
        jetSmearedPEta[i]=NOVAL_F;
        jetSmearedPhi[i]=NOVAL_F;
        jetY[i]=NOVAL_F;
        jetelectronMultiplicity[i]=NOVAL_F;
        jetmuonMultiplicity[i]=NOVAL_F;
        jetneutralHadronMultiplicity[i]=NOVAL_F;
        jetneutralMultiplicity[i]=NOVAL_F;
        jetphotonMultiplicity[i]=NOVAL_F;
        jetsCorrPt[i]=NOVAL_F;
        jetsCorrEta[i]=NOVAL_F;
        jetsCorrPhi[i]=NOVAL_F;
        jetsCorrE[i]=NOVAL_F;
        jetsCorrMass[i]=NOVAL_F;
        jetsCorrNJets[i]=NOVAL_F;
        jetsCorrPartonFlavour[i]=NOVAL_F;
      }
      E=NOVAL_F;
      Pt=NOVAL_F;
      Mass=NOVAL_F;
      Eta=NOVAL_F;
      PartonFlavour=NOVAL_F;
      Phi=NOVAL_F;
      CSV=NOVAL_F;
      CSVV1=NOVAL_F;
      Charge=NOVAL_F;
      ChargeMuEnergy=NOVAL_F;
      ChargedHadronMultiplicity=NOVAL_F;
      ElectronEnergy=NOVAL_F;
      GenJetCharge=NOVAL_F;
      GenJetE=NOVAL_F;
      GenJetEta=NOVAL_F;
      GenJetPhi=NOVAL_F;
      GenJetPt=NOVAL_F;
      GenJetY=NOVAL_F;
      GenPartonCharge=NOVAL_F;
      GenPartonE=NOVAL_F;
      GenPartonEta=NOVAL_F;
      GenPartonPhi=NOVAL_F;
      GenPartonPt=NOVAL_F;
      GenPartonY=NOVAL_F;
      HFEMEnergy=NOVAL_F;
      HFEMMultiplicity=NOVAL_F;
      HFHadronEnergy=NOVAL_F;
      HFHadronMultiplicity=NOVAL_F;
      HLTjetDeltaR=NOVAL_F;
      HLTjetE=NOVAL_F;
      HLTjetEta=NOVAL_F;
      HLTjetPt=NOVAL_F;
      HLTjetPhi=NOVAL_F;
      HadronFlavour=NOVAL_F;
      IsCSVL=NOVAL_F;
      IsCSVM=NOVAL_F;
      IsCSVT=NOVAL_F;
      SmearedE=NOVAL_F;
      SmearedPt=NOVAL_F;
      SmearedPEta=NOVAL_F;
      SmearedPhi=NOVAL_F;
      Y=NOVAL_F;
      electronMultiplicity=NOVAL_F;
      muonMultiplicity=NOVAL_F;
      neutralHadronMultiplicity=NOVAL_F;
      neutralMultiplicity=NOVAL_F;
      photonMultiplicity=NOVAL_F;
      CorrPt=NOVAL_F;
      CorrEta=NOVAL_F;
      CorrPhi=NOVAL_F;
      CorrE=NOVAL_F;
      CorrMass=NOVAL_F;
      CorrNJets=NOVAL_F;
      CorrPartonFlavour=NOVAL_F;
      
      it_ = 0;
    }
    
    bool Loop() {
      if (it_<jets_size) {
	E                         = jetE[it_];			 
        Pt			  = jetPt[it_];			 
        Mass			  = jetMass[it_];			 
        Eta			  = jetEta[it_];			 
        PartonFlavour		  = jetPhi[it_];			 
        Phi			  = jetPartonFlavour[it_];		 
        CSV			  = jetCSV[it_];			 
        CSVV1			  = jetCSVV1[it_];			 
        Charge			  = jetCharge[it_];			 
        ChargeMuEnergy		  = jetChargeMuEnergy[it_];		 
        ChargedHadronMultiplicity = jetChargedHadronMultiplicity[it_]; 
        ElectronEnergy		  = jetElectronEnergy[it_];		 
        GenJetCharge		  = jetGenJetCharge[it_];		 
        GenJetE			  = jetGenJetE[it_];			 
        GenJetEta		  = jetGenJetEta[it_];		 
        GenJetPhi		  = jetGenJetPhi[it_];		 
        GenJetPt		  = jetGenJetPt[it_];			 
        GenJetY			  = jetGenJetY[it_];			 
        GenPartonCharge		  = jetGenPartonCharge[it_];		 
        GenPartonE		  = jetGenPartonE[it_];		 
        GenPartonEta		  = jetGenPartonEta[it_];		 
        GenPartonPhi		  = jetGenPartonPhi[it_];		 
        GenPartonPt		  = jetGenPartonPt[it_];		 
        GenPartonY		  = jetGenPartonY[it_];		 
        HFEMEnergy		  = jetHFEMEnergy[it_];		 
        HFEMMultiplicity	  = jetHFEMMultiplicity[it_];		 
        HFHadronEnergy		  = jetHFHadronEnergy[it_];		 
        HFHadronMultiplicity	  = jetHFHadronMultiplicity[it_];	 
        HLTjetDeltaR		  = jetHLTjetDeltaR[it_];		 
        HLTjetE			  = jetHLTjetE[it_];			 
        HLTjetEta		  = jetHLTjetEta[it_];		 
        HLTjetPt		  = jetHLTjetPt[it_];			 
        HLTjetPhi		  = jetHLTjetPhi[it_];		 
        HadronFlavour		  = jetHadronFlavour[it_];		 
        IsCSVL			  = jetIsCSVL[it_];			 
        IsCSVM			  = jetIsCSVM[it_];			 
        IsCSVT			  = jetIsCSVT[it_];			 
        SmearedE		  = jetSmearedE[it_];			 
        SmearedPt		  = jetSmearedPt[it_];		 
        SmearedPEta		  = jetSmearedPEta[it_];		 
        SmearedPhi		  = jetSmearedPhi[it_];		 
        Y			  = jetY[it_];			 
        electronMultiplicity	  = jetelectronMultiplicity[it_];	 
        muonMultiplicity	  = jetmuonMultiplicity[it_];		 
        neutralHadronMultiplicity = jetneutralHadronMultiplicity[it_]; 
        neutralMultiplicity	  = jetneutralMultiplicity[it_];	 
        photonMultiplicity	  = jetphotonMultiplicity[it_];	 
        CorrPt			  = jetsCorrPt[it_];			 
        CorrEta			  = jetsCorrEta[it_];			 
        CorrPhi			  = jetsCorrPhi[it_];			 
        CorrE			  = jetsCorrE[it_];			 
        CorrMass		  = jetsCorrMass[it_];		 
        CorrNJets		  = jetsCorrNJets[it_];		 
        CorrPartonFlavour	  = jetsCorrPartonFlavour[it_];        
	
	++it_;
	return 1;
      } else {
	it_=0;
	return 0;
      }
    }

  private:
    int it_;
    
  } jet;
  
  class JetAK8Data {
  public:
    JetAK8Data() { init(); };
    
    int jetsAK8_size;
    float jetAK8E[20];
    float jetAK8Pt[20];
    float jetAK8Mass[20];
    float jetAK8Eta[20];
    float jetAK8Phi[20];
    float jetAK8PartonFlavour[20];
    float jetAK8CSV[20];
    float jetAK8CSVV1[20];
    float jetAK8Charge[20];
    float jetAK8ChargeMuEnergy[20];
    float jetAK8ElectronEnergy[20];
    float jetAK8GenJetCharge[20];
    float jetAK8GenJetE[20];
    float jetAK8GenJetEta[20];
    float jetAK8GenJetPhi[20];
    float jetAK8GenJetPt[20];
    float jetAK8GenJetY[20];
    float jetAK8GenPartonCharge[20];
    float jetAK8GenPartonE[20];
    float jetAK8GenPartonEta[20];
    float jetAK8GenPartonPhi[20];
    float jetAK8GenPartonPt[20];
    float jetAK8GenPartonY[20];
    float jetAK8HFEMEnergy[20];
    float jetAK8HFEMMultiplicity[20];
    float jetAK8HFHadronEnergy[20];
    float jetAK8HFHadronMultiplicity[20];
    float jetAK8HLTjetDeltaR[20];
    float jetAK8HLTjetE[20];
    float jetAK8HLTjetEta[20];
    float jetAK8HLTjetPt[20];
    float jetAK8HLTjetPhi[20];
    float jetAK8HadronFlavour[20];
    float jetAK8IsCSVL[20];
    float jetAK8IsCSVM[20];
    float jetAK8IsCSVT[20];
    float jetAK8SmearedE[20];
    float jetAK8SmearedPt[20];
    float jetAK8SmearedPEta[20];
    float jetAK8SmearedPhi[20];
    float jetAK8Y[20];
    float jetAK8electronMultiplicity[20];
    float jetAK8muonMultiplicity[20];
    float jetAK8neutralMultiplicity[20];
    float jetAK8photonMultiplicity[20];
    float jetsAK8CorrPt[20];
    float jetsAK8CorrEta[20];
    float jetsAK8CorrPhi[20];
    float jetsAK8CorrE[20];
    float jetsAK8CorrMass[20];
    float jetsAK8CorrNJets[20];
    float jetsAK8CorrPartonFlavour[20];
    
    float E;
    float Pt;
    float Mass;
    float Eta;
    float Phi;
    float PartonFlavour;
    float CSV;
    float CSVV1;
    float Charge;
    float ChargeMuEnergy;
    float ElectronEnergy;
    float GenJetCharge;
    float GenJetE;
    float GenJetEta;
    float GenJetPhi;
    float GenJetPt;
    float GenJetY;
    float GenPartonCharge;
    float GenPartonE;
    float GenPartonEta;
    float GenPartonPhi;
    float GenPartonPt;
    float GenPartonY;
    float HFEMEnergy;
    float HFEMMultiplicity;
    float HFHadronEnergy;
    float HFHadronMultiplicity;
    float HLTjetDeltaR;
    float HLTjetE;
    float HLTjetEta;
    float HLTjetPt;
    float HLTjetPhi;
    float HadronFlavour;
    float IsCSVL;
    float IsCSVM;
    float IsCSVT;
    float SmearedE;
    float SmearedPt;
    float SmearedPEta;
    float SmearedPhi;
    float Y;
    float electronMultiplicity;
    float muonMultiplicity;
    float neutralMultiplicity;
    float photonMultiplicity;
    float CorrPt;
    float CorrEta;
    float CorrPhi;
    float CorrE;
    float CorrMass;
    float CorrNJets;
    float CorrPartonFlavour;
    
    void init() {
      jetsAK8_size=0;
      for (int i=0; i<20; ++i) {
        jetAK8E[i]=NOVAL_F;
        jetAK8Pt[i]=NOVAL_F;
        jetAK8Mass[i]=NOVAL_F;
        jetAK8Eta[i]=NOVAL_F;
        jetAK8Phi[i]=NOVAL_F;
        jetAK8PartonFlavour[i]=NOVAL_F;
        jetAK8CSV[i]=NOVAL_F;
        jetAK8CSVV1[i]=NOVAL_F;
        jetAK8Charge[i]=NOVAL_F;
        jetAK8ChargeMuEnergy[i]=NOVAL_F;
        jetAK8ElectronEnergy[i]=NOVAL_F;
        jetAK8GenJetCharge[i]=NOVAL_F;
        jetAK8GenJetE[i]=NOVAL_F;
        jetAK8GenJetEta[i]=NOVAL_F;
        jetAK8GenJetPhi[i]=NOVAL_F;
        jetAK8GenJetPt[i]=NOVAL_F;
        jetAK8GenJetY[i]=NOVAL_F;
        jetAK8GenPartonCharge[i]=NOVAL_F;
        jetAK8GenPartonE[i]=NOVAL_F;
        jetAK8GenPartonEta[i]=NOVAL_F;
        jetAK8GenPartonPhi[i]=NOVAL_F;
        jetAK8GenPartonPt[i]=NOVAL_F;
        jetAK8GenPartonY[i]=NOVAL_F;
        jetAK8HFEMEnergy[i]=NOVAL_F;
        jetAK8HFEMMultiplicity[i]=NOVAL_F;
        jetAK8HFHadronEnergy[i]=NOVAL_F;
        jetAK8HFHadronMultiplicity[i]=NOVAL_F;
        jetAK8HLTjetDeltaR[i]=NOVAL_F;
        jetAK8HLTjetE[i]=NOVAL_F;
        jetAK8HLTjetEta[i]=NOVAL_F;
        jetAK8HLTjetPt[i]=NOVAL_F;
        jetAK8HLTjetPhi[i]=NOVAL_F;
        jetAK8HadronFlavour[i]=NOVAL_F;
        jetAK8IsCSVL[i]=NOVAL_F;
        jetAK8IsCSVM[i]=NOVAL_F;
        jetAK8IsCSVT[i]=NOVAL_F;
        jetAK8SmearedE[i]=NOVAL_F;
        jetAK8SmearedPt[i]=NOVAL_F;
        jetAK8SmearedPEta[i]=NOVAL_F;
        jetAK8SmearedPhi[i]=NOVAL_F;
        jetAK8Y[i]=NOVAL_F;
        jetAK8electronMultiplicity[i]=NOVAL_F;
        jetAK8muonMultiplicity[i]=NOVAL_F;
        jetAK8neutralMultiplicity[i]=NOVAL_F;
        jetAK8photonMultiplicity[i]=NOVAL_F;
        jetsAK8CorrPt[i]=NOVAL_F;
        jetsAK8CorrEta[i]=NOVAL_F;
        jetsAK8CorrPhi[i]=NOVAL_F;
        jetsAK8CorrE[i]=NOVAL_F;
        jetsAK8CorrMass[i]=NOVAL_F;
        jetsAK8CorrNJets[i]=NOVAL_F;
        jetsAK8CorrPartonFlavour[i]=NOVAL_F;
      }
      E=NOVAL_F;
      Pt=NOVAL_F;
      Mass=NOVAL_F;
      Eta=NOVAL_F;
      Phi=NOVAL_F;
      PartonFlavour=NOVAL_F;
      CSV=NOVAL_F;
      CSVV1=NOVAL_F;
      Charge=NOVAL_F;
      ChargeMuEnergy=NOVAL_F;
      ElectronEnergy=NOVAL_F;
      GenJetCharge=NOVAL_F;
      GenJetE=NOVAL_F;
      GenJetEta=NOVAL_F;
      GenJetPhi=NOVAL_F;
      GenJetPt=NOVAL_F;
      GenJetY=NOVAL_F;
      GenPartonCharge=NOVAL_F;
      GenPartonE=NOVAL_F;
      GenPartonEta=NOVAL_F;
      GenPartonPhi=NOVAL_F;
      GenPartonPt=NOVAL_F;
      GenPartonY=NOVAL_F;
      HFEMEnergy=NOVAL_F;
      HFEMMultiplicity=NOVAL_F;
      HFHadronEnergy=NOVAL_F;
      HFHadronMultiplicity=NOVAL_F;
      HLTjetDeltaR=NOVAL_F;
      HLTjetE=NOVAL_F;
      HLTjetEta=NOVAL_F;
      HLTjetPt=NOVAL_F;
      HLTjetPhi=NOVAL_F;
      HadronFlavour=NOVAL_F;
      IsCSVL=NOVAL_F;
      IsCSVM=NOVAL_F;
      IsCSVT=NOVAL_F;
      SmearedE=NOVAL_F;
      SmearedPt=NOVAL_F;
      SmearedPEta=NOVAL_F;
      SmearedPhi=NOVAL_F;
      Y=NOVAL_F;
      electronMultiplicity=NOVAL_F;
      muonMultiplicity=NOVAL_F;
      neutralMultiplicity=NOVAL_F;
      photonMultiplicity=NOVAL_F;
      CorrPt=NOVAL_F;
      CorrEta=NOVAL_F;
      CorrPhi=NOVAL_F;
      CorrE=NOVAL_F;
      CorrMass=NOVAL_F;
      CorrNJets=NOVAL_F;
      CorrPartonFlavour=NOVAL_F;
      
      it_ = 0;
    }

    bool Loop() {
      if (it_<jetsAK8_size) {
	E                    = jetAK8E[it_];		       
        Pt		     = jetAK8Pt[it_];		       
        Mass		     = jetAK8Mass[it_];		       
        Eta		     = jetAK8Eta[it_];		       
        Phi		     = jetAK8Phi[it_];		       
        PartonFlavour	     = jetAK8PartonFlavour[it_];	       
        CSV		     = jetAK8CSV[it_];		       
        CSVV1		     = jetAK8CSVV1[it_];		       
        Charge		     = jetAK8Charge[it_];	       
        ChargeMuEnergy	     = jetAK8ChargeMuEnergy[it_];       
        ElectronEnergy	     = jetAK8ElectronEnergy[it_];       
        GenJetCharge	     = jetAK8GenJetCharge[it_];	       
        GenJetE		     = jetAK8GenJetE[it_];	       
        GenJetEta	     = jetAK8GenJetEta[it_];	       
        GenJetPhi	     = jetAK8GenJetPhi[it_];	       
        GenJetPt	     = jetAK8GenJetPt[it_];	       
        GenJetY		     = jetAK8GenJetY[it_];	       
        GenPartonCharge	     = jetAK8GenPartonCharge[it_];      
        GenPartonE	     = jetAK8GenPartonE[it_];	       
        GenPartonEta	     = jetAK8GenPartonEta[it_];	       
        GenPartonPhi	     = jetAK8GenPartonPhi[it_];	       
        GenPartonPt	     = jetAK8GenPartonPt[it_];	       
        GenPartonY	     = jetAK8GenPartonY[it_];	       
        HFEMEnergy	     = jetAK8HFEMEnergy[it_];	       
        HFEMMultiplicity     = jetAK8HFEMMultiplicity[it_];     
        HFHadronEnergy	     = jetAK8HFHadronEnergy[it_];       
        HFHadronMultiplicity = jetAK8HFHadronMultiplicity[it_]; 
        HLTjetDeltaR	     = jetAK8HLTjetDeltaR[it_];	       
        HLTjetE		     = jetAK8HLTjetE[it_];	       
        HLTjetEta	     = jetAK8HLTjetEta[it_];	       
        HLTjetPt	     = jetAK8HLTjetPt[it_];	       
        HLTjetPhi	     = jetAK8HLTjetPhi[it_];	       
        HadronFlavour	     = jetAK8HadronFlavour[it_];	       
        IsCSVL		     = jetAK8IsCSVL[it_];	       
        IsCSVM		     = jetAK8IsCSVM[it_];	       
        IsCSVT		     = jetAK8IsCSVT[it_];	       
        SmearedE	     = jetAK8SmearedE[it_];	       
        SmearedPt	     = jetAK8SmearedPt[it_];	       
        SmearedPEta	     = jetAK8SmearedPEta[it_];	       
        SmearedPhi	     = jetAK8SmearedPhi[it_];	       
        Y		     = jetAK8Y[it_];		       
        electronMultiplicity = jetAK8electronMultiplicity[it_]; 
        muonMultiplicity     = jetAK8muonMultiplicity[it_];     
        neutralMultiplicity  = jetAK8neutralMultiplicity[it_];  
        photonMultiplicity   = jetAK8photonMultiplicity[it_];   
        CorrPt		     = jetsAK8CorrPt[it_];	       
        CorrEta		     = jetsAK8CorrEta[it_];	       
        CorrPhi		     = jetsAK8CorrPhi[it_];	       
        CorrE		     = jetsAK8CorrE[it_];	       
        CorrMass	     = jetsAK8CorrMass[it_];	       
        CorrNJets	     = jetsAK8CorrNJets[it_];	       
        CorrPartonFlavour    = jetsAK8CorrPartonFlavour[it_];   
	
	++it_;
	return 1;
      } else {
	it_=0;
	return 0;
      }
    }
    
  private:
    int it_;
    
  } jetAK8;
  
  class MetData {
  public:
    MetData() { init(); };
    
    int met_size;
    float metPt[20];
    float metPhi[20];
    float metPx[20];
    float metPy[20];
    float metCorrPt[20];
    float metCorrPhi[20];
    
    float Pt;
    float Phi;
    float Px;
    float Py;
    float CorrPt;
    float CorrPhi;
    
    void init() {
      met_size=0;
      for (int i=0; i<20; ++i) {
        metPt[i]=NOVAL_F;
        metPhi[i]=NOVAL_F;
        metPx[i]=NOVAL_F;
        metPy[i]=NOVAL_F;
        metCorrPt[i]=NOVAL_F;
        metCorrPhi[i]=NOVAL_F;
      }
      Pt=NOVAL_F;
      Phi=NOVAL_F;
      Px=NOVAL_F;
      Py=NOVAL_F;
      CorrPt=NOVAL_F;
      CorrPhi=NOVAL_F;
      
      it_ = 0;
    }

    bool Loop() {
      if (it_<met_size) {
	Pt      = metPt[it_];     
	Phi	= metPhi[it_];    
	Px	= metPx[it_];     
	Py	= metPy[it_];     
	CorrPt	= metCorrPt[it_]; 
	CorrPhi	= metCorrPhi[it_];
	
	++it_;
	return 1;
      } else {
	it_=0;
	return 0;
      }
    }

  private:
    int it_;
    
  } met;
  
  class EventData {
  public:
    EventData() { init(); };
    
    float Event_nTightMuons;
    float Event_nLooseMuons;
    float Event_nTightElectrons;
    float Event_nLooseElectrons;
    float Event_nElectronsSF;
    float Event_nMuonsSF;
    float Event_nCSVTJets;
    float Event_nCSVMJets;
    float Event_nCSVLJets;
    float Event_nTightJets;
    float Event_nLooseJets;
    float Event_bWeight1TCSVT;
    float Event_bWeight1TCSVM;
    float Event_bWeight1TCSVL;
    float Event_bWeight2TCSVT;
    float Event_bWeight2TCSVM;
    float Event_bWeight2TCSVL;
    float Event_LHEWeightSign;
    
    void init() {
      Event_nTightMuons=NOVAL_F;
      Event_nLooseMuons=NOVAL_F;
      Event_nTightElectrons=NOVAL_F;
      Event_nLooseElectrons=NOVAL_F;
      Event_nElectronsSF=NOVAL_F;
      Event_nMuonsSF=NOVAL_F;
      Event_nCSVTJets=NOVAL_F;
      Event_nCSVMJets=NOVAL_F;
      Event_nCSVLJets=NOVAL_F;
      Event_nTightJets=NOVAL_F;
      Event_nLooseJets=NOVAL_F;
      Event_bWeight1TCSVT=NOVAL_F;
      Event_bWeight1TCSVM=NOVAL_F;
      Event_bWeight1TCSVL=NOVAL_F;
      Event_bWeight2TCSVT=NOVAL_F;
      Event_bWeight2TCSVM=NOVAL_F;
      Event_bWeight2TCSVL=NOVAL_F;
      Event_LHEWeightSign=NOVAL_F;
    }
    
  } evt;
  
};

#endif
