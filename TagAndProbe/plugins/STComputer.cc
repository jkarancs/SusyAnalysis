#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

template<typename T>
class STComputer : public edm::EDProducer {
public:
  explicit STComputer(const edm::ParameterSet & iConfig);
  virtual ~STComputer() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;
  edm::InputTag metTag_;
};

template<typename T>
STComputer<T>::STComputer(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes")),
    metTag_(iConfig.getParameter<edm::InputTag>("metTag")) {
  produces<edm::ValueMap<float> >();
}


template<typename T>
STComputer<T>::~STComputer()
{
}

template<typename T>
void 
STComputer<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  
  // read input
  Handle<View<reco::Candidate> > probes;
  Handle<View<T> > metHandle;
  iEvent.getByLabel(probes_,  probes);
  iEvent.getByLabel(metTag_, metHandle);
  
  // met
  float met = metHandle.isValid() ? metHandle->front().pt() : 0.0;
  
  // prepare and fill vector for output
  std::vector<float> v_stlep;
  View<reco::Candidate>::const_iterator probe, endprobes = probes->end();
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    v_stlep.push_back(probe->pt()+met);
  }
  
  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > vm_stlep(new ValueMap<float>());
  ValueMap<float>::Filler filler_stlep(*vm_stlep);
  filler_stlep.insert(probes, v_stlep.begin(), v_stlep.end());
  filler_stlep.fill();
  iEvent.put(vm_stlep);
}


typedef STComputer<pat::MET> PatMetSTComputer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatMetSTComputer);
