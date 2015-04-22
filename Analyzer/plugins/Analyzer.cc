// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Physics objects
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Filter Results
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// User classes
#include "SusyAnalysis/Analyzer/interface/SmartHistos.h"

//
// class declaration
//

class Analyzer : public edm::EDAnalyzer {
public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void beginJob() override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  //edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;

  SmartHistos* sh_;
  // objects for SmartHistos
  pat::Muon muon_;
  pat::Electron ele_;
  pat::Photon pho_;
  pat::Tau tau_;
  pat::Jet jet_;
  pat::MET met_;
};

Analyzer::Analyzer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  //fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
  sh_ = new SmartHistos();
}

Analyzer::~Analyzer() {}

void Analyzer::beginJob() {
  
  sh_->AddHistoType("muons");
  sh_->AddHistoType("electrons");
  sh_->AddHistoType("met");
  
  // Define Postfixes here:
  //sh_->AddNewPostfix("Layers",     &v.pf_layers,   "Lay[1to3]", "Layer [1to3]", "2-4");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh_->AddNewFillParam("MuPt",     { .nbin=200, .bins={0, 200}, .fill=[this](){ return muon_.pt(); }, .axis_title=";Muon p_{T} (GeV/c)"});
  sh_->AddNewFillParam("ElePt",    { .nbin=200, .bins={0, 200}, .fill=[this](){ return ele_.pt(); },  .axis_title=";Electron p_{T} (GeV/c)"});
  sh_->AddNewFillParam("MetPt",    { .nbin=200, .bins={0, 200}, .fill=[this](){ return met_.pt(); },  .axis_title=";MET p_{T} (GeV/c)"});
  sh_->AddNewFillParam("GenMetPt", { .nbin=200, .bins={0, 200}, .fill=[this](){ return met_.genMET()->pt(); },  .axis_title=";Gen-MET p_{T} (GeV/c)"});
  
  // Define Cuts here:
  sh_->AddNewCut("GenMetPt>1", [this](){ return met_.genMET()->pt()>1; });
  
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  // Validation  
  //sh_->SetHistoWeights({&v.corr_pu_weight, &v.norm_factor});
  //sh_->AddHistos("evt", { .fill="NVertices", .pfs={"Validation"}, .cuts={"ZeroBias","Nvtx"}, .draw="", .opt="", .ranges={} });
  
  sh_->AddHistos("muons",     { .fill="MuPt",     .pfs={}, .cuts={}, .draw="", .opt="", .ranges={} });
  sh_->AddHistos("electrons", { .fill="ElePt",    .pfs={}, .cuts={}, .draw="", .opt="", .ranges={} });
  sh_->AddHistos("met",       { .fill="MetPt",    .pfs={}, .cuts={"GenMetPt>1"}, .draw="", .opt="", .ranges={} });
  sh_->AddHistos("met",       { .fill="GenMetPt", .pfs={}, .cuts={"GenMetPt>1"}, .draw="", .opt="", .ranges={} });
  std::cout<<"-----------------------------------------------------------------\n";
  std::cout<<"SmartHistos: Creating the following plots:\n"; sh_->PrintNames();
  std::cout<<"-----------------------------------------------------------------\n";
  
  //sh_->Add("PHM_out/Validation2_VER5.root");
   
}

void Analyzer::endJob() {
  // Write histograms and canvases in an output file
  TFile* f_out = new TFile("test.root","recreate");
  std::cout<<"Writing plots to file: "<<f_out->GetName()<<" done.\n";
  //sh_->DrawPlots();
  sh_->Write();
  f_out->Close();
}


void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bool print = false;
  
  //___________________________________________________________________
  //                          Vertices  

  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();
  /*
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin();
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // The "good vertex" selection is borrowed from Giovanni Zevi Della Porta
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if ( //!vtx->isFake() &&
	!(vtx->chi2()==0 && vtx->ndof()==0)
	&& vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	&& fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }
  
  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs
  */

  //___________________________________________________________________
  //                           Muons  
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    muon_ = mu;
    if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
    sh_->Fill("muons");
    if (print) printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
		      mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
  }
  
  //___________________________________________________________________
  //                          Electrons
  
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &el : *electrons) {
    ele_ = el;
    // Keep only electrons above 5 GeV.
    // NOTE: miniAOD does not store some of the info for electrons <5 GeV at all!
    if (el.pt() < 5) continue;
    sh_->Fill("electrons");
    if (print) printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), pass conv veto %d\n",
		      el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.passConversionVeto());
    
    // lines below taken from:
    // https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/ElectronNtupler.cc
    // Kinematics

    /*
    double pt = el.pt();
    double etaSC = el.superCluster()->eta();
    
    // ID and matching
    double dEtaIn = el.deltaEtaSuperClusterTrackAtVtx();
    double dPhiIn = el.deltaPhiSuperClusterTrackAtVtx();
    //double sigmaIetaIeta = el.sigmaIetaIeta();
    double full5x5_sigmaIetaIeta = el.full5x5_sigmaIetaIeta();
    double hOverE = el.hcalOverEcal();
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero (always
    // the case for electrons below 5 GeV in miniAOD)
    double ooEmooP = 1e30;
    if( el.ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP = 1e30;
    }else if( !std::isfinite(el.ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP = 1e30;
    }else{
      ooEmooP = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
    }



    //
    double cut_abs_dEtaIn[2][2][4] = { { { 0.012, 0.012, 0.009, 0.009 },
					 { 0.022, 0.021, 0.019, 0.017 } },
				       { { 0.012, 0.012, 0.009, 0.009 },
					 { 0.015, 0.014, 0.012, 0.010 } } };
    double cut_abs_dPhiIn[2][2][4] = { { { 0.8 	0.15 	0.06 	0.03  },
					 { 0.7 	0.10 	0.03 	0.02  } },
				       { { 0.8 	0.15 	0.06 	0.03  },
					 { 0.7 	0.10 	0.03 	0.02  } } };
    double cut_full5x5_sigmaIetaIeta[2][2][4] = { { { 0.012 	0.012 	0.01 	0.01 },
						    { 0.033 	0.033 	0.031 	0.031 } },
						  { { 0.01 	0.01 	0.01 	0.01 },
						    { 0.033 	0.033 	0.031 	0.031 } } };
    double cut_hOverE[2][2][4] = { { {  	0.15 	0.12 	0.12 	0.12  },
				     {  	N/A 	0.12 	0.12 	0.12 } },
				   { {  	0.15 	0.12 	0.12 	0.12  },
				     {  	N/A 	0.12 	0.12 	0.12 } } };
    double cut_abs_d0[2][2][4] = { { {  	0.04 	0.02 	0.02 	0.02 },
				   { 0.04 	0.02 	0.02 	0.02 } },
				 { {  	0.04 	0.02 	0.02 	0.02  },
				   {  	0.04 	0.02 	0.02 	0.02  } } };
    double cut_abs_dz[2][2][4] = { { {  	0.2 	0.2 	0.1 	0.1  },
				   {  	0.2 	0.2 	0.1 	0.1 } },
				 { {  	0.2 	0.2 	0.1 	0.1  },
				   { 0.2 	0.2 	0.1 	0.1  } } };
    double cut_ooEmOOP[2][2][4] = { { {  	N/A 	0.05 	0.05 	0.05  },
				   { N/A 	0.05 	0.05 	0.05 } },
				 { {  	N/A 	0.05 	0.05 	0.05  },
				   {  	N/A 	0.05 	0.05 	0.05  } } };
    double cut_abs_[2][2][4] = { { {  },
				   {  } },
				 { {  },
				   {  } } };
    double cut_abs_[2][2][4] = { { {  },
				   {  } },
				 { {  },
				   {  } } };
    double cut_abs_[2][2][4] = { { {  },
				   {  } },
				 { {  },
				   {  } } };
    double cut_abs_[2][2][4] = { { {  },
				   {  } },
				 { {  },
				   {  } } };
    double cut_abs_[2][2][4] = { { {  },
				   {  } },
				 { {  },
				   {  } } };
    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
    // Compute isolation with delta beta correction for PU
    double absiso = pfIso.sumChargedHadronPt
      + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    relIsoWithDBeta = absiso/pt_;
    
    // Impact parameter
    double d0 = (-1) * el.gsfTrack()->dxy(firstGoodVertex->position() );
    double dz = el.gsfTrack()->dz( firstGoodVertex->position() );
    
    // Conversion rejection
    int expectedMissingInnerHits = el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    bool passConversionVeto = el.passConversionVeto();
    
    // Match to generator level truth
    //
    // Explicit loop over gen candidates method
    // 
    // Functions implemented in:
    //   https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/ElectronNtupler.cc
    // bool isTrueElectron = matchToTruth( el, prunedGenParticles);
    // bool isTrueElectronAlternative = matchToTruthAlternative( el );
    
    // Electron Cut Based ID:
    bool is_alignment_good = true;
    */
    
  }
  
  //___________________________________________________________________
  //                          Photons
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);
  for (const pat::Photon &pho : *photons) {
    pho_ = pho;
    if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    if (print) printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f\n",
		      pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta());
  }
  
  //___________________________________________________________________
  //                           Taus
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  for (const pat::Tau &tau : *taus) {
    tau_ = tau;
    if (tau.pt() < 20) continue;
    if (print) printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
		      tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
  }
  
  //___________________________________________________________________
  //                            Jets
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  int ijet = 0;
  for (const pat::Jet &j : *jets) {
    jet_ = j;
    if (j.pt() < 20) continue;
    if (print) printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n", 
		      j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
    if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
      std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
      for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
        if (print) printf("         constituent %3d: pt %6.2f, pdgId %+3d\n", i2,(*daus[i2]).pt(),(*daus[i2]).pdgId());
	// for MiniAOD (PackedCandidate) only:
        //const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
        //printf("         constituent %3d: pt %6.2f, dz(PV) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
      }
    }
  }
  
  // MiniAOD:
  //edm::Handle<pat::JetCollection> fatjets;
  //iEvent.getByToken(fatjetToken_, fatjets);
  //for (const pat::Jet &j : *fatjets) {
  //    printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
  //        j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSPrunedLinks"), j.userFloat("ak8PFJetsCHSTrimmedLinks"), j.userFloat("ak8PFJetsCHSFilteredLinks"), j.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
  //}
  
  //___________________________________________________________________
  //                             MET
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  met_ = met;
  sh_->Fill("met");
  // MiniAOD Only:
  //printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
  //       met.pt(), met.phi(), met.sumEt(),
  //       met.genMET()->pt(),
  //       met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
  if (print) printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f.\n",
		    met.pt(), met.phi(), met.sumEt(),
		    met.genMET()->pt());
  if (print) printf("\n");
  

  //___________________________________________________________________
  //                     Filter decisions

  edm::Handle<edm::TriggerResults> filterBits;
  edm::InputTag labfilterBits("TriggerResults");
  iEvent.getByLabel(labfilterBits,filterBits);
  const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
  for (unsigned int i = 0, n = filterBits->size(); i < n; ++i) {
    std::string filterName = fnames.triggerName(i);
    int filterdecision = filterBits->accept(i);
    if (print) { // all lines = 1 for MiniAOD
      if (filterName=="Flag_trackingFailureFilter") std::cout<<"trackingfailturefilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_goodVertices") std::cout<<"goodVerticesfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_CSCTightHaloFilter") std::cout<<"cschalofilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_trkPOGFilters") std::cout<<"trkPOGfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_trkPOG_logErrorTooManyClusters") std::cout<<"trkPOG_logErrorTooManyClustersfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_EcalDeadCellTriggerPrimitiveFilter") std::cout<<"EcalDeadCellTriggerPrimitivefilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_ecalLaserCorrFilter") std::cout<<"ecallaserfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_trkPOG_manystripclus53X") std::cout<<"trkPOG_manystripclus53XfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_eeBadScFilter") std::cout<<"eebadscfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_METFilters") std::cout<<"METFiltersfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_HBHENoiseFilter") std::cout<<"HBHENoisefilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_trkPOG_toomanystripclus53X") std::cout<<"trkPOG_toomanystripclus53XfilterResult = "<<filterdecision<<"\n";
      else if (filterName=="Flag_hcalLaserEventFilter") std::cout<<"hcallaserfilterResult = "<<filterdecision<<"\n";
    }
  }
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
