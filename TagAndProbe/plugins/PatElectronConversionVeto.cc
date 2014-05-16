//

/*
  \class    	  PatElectronConversionVeto
  \description:
  	  	 Produces valuemaps for Pat Electron conversion veto
  	  	 To be used in T&P in the following way:
		 
			 process.patEleConvVeto = cms.EDProducer("PatElectronConversionVeto",
				probes = cms.InputTag("patElectronCollection"),
			 )
			 
			 // in TagProbeFitTreeProducer
			 variables = cms.PSet(
			 	...
				noconv = cms.InputTag("patEleConvVeto","noconv"),
				...
			 )
			 
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

class PatElectronConversionVeto : public edm::EDProducer {
public:
  explicit PatElectronConversionVeto(const edm::ParameterSet & iConfig);
  virtual ~PatElectronConversionVeto() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;

};

PatElectronConversionVeto::PatElectronConversionVeto(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes"))
{
  produces<edm::ValueMap<bool> >("noconv");
}

PatElectronConversionVeto::~PatElectronConversionVeto()
{
}

void PatElectronConversionVeto::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Read input
  edm::Handle<pat::ElectronCollection> probes;
  bool found = iEvent.getByLabel(probes_,probes);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get PatElectrons: "
       <<probes_<<std::endl;
    edm::LogError("PatElectronConversionVeto")<<err.str();
    throw cms::Exception( "MissingProduct", err.str());
  }
  
  // Conversions
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions",hConversions);

  // BeamSpot
  math::XYZPoint beamSpotPosition;
  beamSpotPosition.SetCoordinates(NAN,NAN,NAN);
  edm::Handle<reco::BeamSpot> bsHandle;
  try {
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    if (!bsHandle.isValid() || bsHandle.failedToGet()) {
      std::cout << "[PatElectronConversionVeto] BeamSpot not valid!." << std::endl;
    } else {
      beamSpotPosition = bsHandle->position();
    }
  } catch (cms::Exception & e) {
    std::cout  << "[PatElectronConversionVeto] error (BeamSpot): " << e.what() << std::endl;
  }

  // Prepare vector for output
  std::vector<bool> noconv;

  // First loop to fill non-pf Electrons
  unsigned nele=probes->size();
  for (unsigned iele=0; iele<nele; ++iele) {
    pat::ElectronRef ele(probes, iele);
    // Conversions check
    bool convveto = !ConversionTools::hasMatchedConversion(ele->core(),hConversions,beamSpotPosition);
    if (convveto!=ele->passConversionVeto()) std::cout<<"Conversion disagreement "<<convveto<<" "<<ele->passConversionVeto()<<std::endl;
    noconv.push_back(convveto);
  }
  
  // Convert into ValueMap and store
  std::auto_ptr<edm::ValueMap<bool> > vm_noconv(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler fill_noconv(*vm_noconv);
  fill_noconv.insert(probes, noconv.begin(), noconv.end());
  fill_noconv.fill();
  iEvent.put(vm_noconv,"noconv");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatElectronConversionVeto);
