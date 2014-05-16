//

/*
  \class    	  MuonDeltaPfRecoPt
  \description:
  	  	 Produces valuemaps for the absolute pt difference between PF and Reco Muons
  	  	 To be used in T&P in the following way:
		 
			 process.deltaPfRecoPt = cms.EDProducer("MuonDeltaPfRecoPt",
				probes = cms.InputTag("patMuonCollection"),
			 )
			 
			 // in TagProbeFitTreeProducer
			 variables = cms.PSet(
			 	...
				absdeltapt = cms.InputTag("deltaPfRecoPt","absdeltapt"),
				...
			 )
			 
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

class MuonDeltaPfRecoPt : public edm::EDProducer {
public:
  explicit MuonDeltaPfRecoPt(const edm::ParameterSet & iConfig);
  virtual ~MuonDeltaPfRecoPt() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;

};

MuonDeltaPfRecoPt::MuonDeltaPfRecoPt(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes"))
{
  produces<edm::ValueMap<float> >("absdeltapt");
}

MuonDeltaPfRecoPt::~MuonDeltaPfRecoPt()
{
}

void MuonDeltaPfRecoPt::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Read input
  edm::Handle<pat::MuonCollection> probes;
  bool found = iEvent.getByLabel(probes_,probes);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get PatMuons: "
       <<probes_<<std::endl;
    edm::LogError("MuonDeltaPfRecoPt")<<err.str();
    throw cms::Exception( "MissingProduct", err.str());
  }
 
  // First loop on Reco Muons to get their references and pt
  std::vector<float> deltapt;

  unsigned nmuon=probes->size();
  for (unsigned imuon=0; imuon<nmuon; ++imuon) {
    pat::MuonRef muon(probes, imuon);
    deltapt.push_back(muon->isPFMuon() ? (muon->pfP4().pt()>10 ? fabs(muon->pt() - muon->pfP4().pt()) : -9999.0) : -9999.0);
    //if (muon->isPFMuon()) std::cout<<deltapt[imuon]<<std::endl;
  }
  
  /* Check Maching
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel("particleFlow",pfCandidates);
  reco::PFCandidateCollection pfmuons = *pfCandidates;
  for (reco::PFCandidateCollection::const_iterator ipfMu = pfmuons.begin();
       ipfMu != pfmuons.end(); ++ipfMu) {
    if ((*ipfMu).particleId() != reco::PFCandidate::mu) continue;

    for (unsigned imuon=0; imuon<nmuon; ++imuon) {
      pat::MuonRef muon(probes, imuon);
      if (fabs(muon->pt()-(*ipfMu).muonRef()->pt())<0.001)
	std::cout<<">> "<<imuon<<" "<<(*ipfMu).muonRef()->pt()<<" "<<(*ipfMu).pt()<<std::endl;
    }
  }
  */

  // Convert into ValueMap and store
  std::auto_ptr<edm::ValueMap<float> > vm_deltapt(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fill_deltapt(*vm_deltapt);
  fill_deltapt.insert(probes, deltapt.begin(), deltapt.end());
  fill_deltapt.fill();
  iEvent.put(vm_deltapt,"absdeltapt");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonDeltaPfRecoPt);
