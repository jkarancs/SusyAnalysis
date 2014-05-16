//

/*
  \class    	  GsfElectronDeltaPfRecoPt
  \description:
  	  	 Produces valuemaps for the absolute pt difference between PF and Reco Electrons
  	  	 To be used in T&P in the following way:
		 
			 process.deltaPfRecoPt = cms.EDProducer("GsfElectronDeltaPfRecoPt",
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

class GsfElectronDeltaPfRecoPt : public edm::EDProducer {
public:
  explicit GsfElectronDeltaPfRecoPt(const edm::ParameterSet & iConfig);
  virtual ~GsfElectronDeltaPfRecoPt() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;

};

GsfElectronDeltaPfRecoPt::GsfElectronDeltaPfRecoPt(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes"))
{
  produces<edm::ValueMap<float> >("absdeltapt");
}

GsfElectronDeltaPfRecoPt::~GsfElectronDeltaPfRecoPt()
{
}

void GsfElectronDeltaPfRecoPt::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Read input
  edm::Handle<reco::GsfElectronCollection> probes;
  bool found = iEvent.getByLabel(probes_,probes);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get GsfElectrons: "
       <<probes_<<std::endl;
    edm::LogError("GsfElectronDeltaPfRecoPt")<<err.str();
    throw cms::Exception("MissingProduct", err.str());
  }
  
  // Prepare vector for output
  std::vector<float> deltapt;
  
  // First loop to fill non-pf Electrons
  unsigned nele=probes->size();
  for (unsigned iele=0; iele<nele; ++iele) {
    reco::GsfElectronRef ele(probes, iele);
    deltapt.push_back(ele->corrections().pflowP4.pt()>10 ? fabs(ele->pt() - ele->corrections().pflowP4.pt()) : -9999.0);
    //std::cout<<ele->corrections().pflowP4.pt()<<" "<<deltapt[iele]<<std::endl;
  }
  
  // Convert into ValueMap and store
  std::auto_ptr<edm::ValueMap<float> > vm_deltapt(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fill_deltapt(*vm_deltapt);
  fill_deltapt.insert(probes, deltapt.begin(), deltapt.end());
  fill_deltapt.fill();
  iEvent.put(vm_deltapt,"absdeltapt");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GsfElectronDeltaPfRecoPt);
