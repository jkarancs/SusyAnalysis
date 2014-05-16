//

/*
  \class    	  ElectronDeltaPfRecoPt
  \description:
  	  	 Produces valuemaps for the absolute pt difference between PF and Reco Electrons
  	  	 To be used in T&P in the following way:
		 
			 process.deltaPfRecoPt = cms.EDProducer("ElectronDeltaPfRecoPt",
				probes = cms.InputTag("patElectronCollection"),
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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

class ElectronDeltaPfRecoPt : public edm::EDProducer {
public:
  explicit ElectronDeltaPfRecoPt(const edm::ParameterSet & iConfig);
  virtual ~ElectronDeltaPfRecoPt() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;

};

ElectronDeltaPfRecoPt::ElectronDeltaPfRecoPt(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes"))
{
  produces<edm::ValueMap<float> >("absdeltapt");
}

ElectronDeltaPfRecoPt::~ElectronDeltaPfRecoPt()
{
}

void ElectronDeltaPfRecoPt::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Read input
  edm::Handle<pat::ElectronCollection> probes;
  bool found = iEvent.getByLabel(probes_,probes);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get PatElectrons: "
       <<probes_<<std::endl;
    edm::LogError("ElectronDeltaPfRecoPt")<<err.str();
    throw cms::Exception( "MissingProduct", err.str());
  }
  
  // Prepare vector for output
  std::vector<float> deltapt;
  
  // First loop to fill non-pf Electrons
  unsigned nele=probes->size();
  for (unsigned iele=0; iele<nele; ++iele) {
    pat::ElectronRef ele(probes, iele);
    deltapt.push_back(ele->corrections().pflowP4.pt()>10 ? fabs(ele->pt() - ele->corrections().pflowP4.pt()) : -9999.0);
    //std::cout<<ele->corrections().pflowP4.pt()<<" "<<deltapt[iele]<<std::endl;
  }
  
  /*
  // Matching PF Electrons to Pat electrons by minimising DeltaR
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel("particleFlow",pfCandidates);
  reco::PFCandidateCollection pfeles = *pfCandidates;
  for (reco::PFCandidateCollection::const_iterator ipfEl = pfeles.begin();
       ipfEl != pfeles.end(); ++ipfEl) {
    if ((*ipfEl).particleId() != reco::PFCandidate::e) continue;
    
    unsigned imin = -1;
    double dr2min = 9999.0;
    for (unsigned iele=0; iele<nele; ++iele) {
      pat::ElectronRef ele(probes, iele);
      double dr2 = deltaR2(probes->at(iele), (*ipfEl));
      if (dr2 < dr2min) {
	dr2min = dr2; 
	imin = iele; 
      }
    }
    pat::ElectronRef ele(probes, imin);
    deltapt[imin] = (*ipfEl).pt()>10 ? fabs(ele->pt() - (*ipfEl).pt()) : -9999.0;
    //std::cout<<"DeltaR: "<<sqrt(dr2min)<<" delta pt: "<<fabs(ele->pt() - (*ipfEl).pt())<<std::endl;
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
DEFINE_FWK_MODULE(ElectronDeltaPfRecoPt);
