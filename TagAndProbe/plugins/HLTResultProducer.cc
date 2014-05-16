#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>


class HLTResultProducer : public edm::EDProducer {
 public:
  explicit HLTResultProducer(const edm::ParameterSet & iConfig);
  virtual ~HLTResultProducer();
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

 private:
  edm::InputTag probes_;
  edm::InputTag triggerTag_;
  std::vector<std::string > HLTPaths_;
  bool andOr_;
  bool first;
};

HLTResultProducer::HLTResultProducer(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes")),
  triggerTag_(iConfig.getParameter<edm::InputTag>("TriggerResultsTag")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string > >("HLTPaths")),
  andOr_(iConfig.getParameter<bool>("andOr")) {

  first=false;
  produces<edm::ValueMap<float> >();
}


HLTResultProducer::~HLTResultProducer() {
}

void HLTResultProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // get HLT decision
    int hltPassed = 0;

    if (HLTPaths_.size()>0) {

      edm::Handle<edm::TriggerResults> triggerResults;
      iEvent.getByLabel(triggerTag_, triggerResults);

      if (triggerResults.isValid()) {

	edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerResults);
	std::vector<unsigned int> indices(HLTPaths_.size(), triggerNames.size());
	size_t nMatches=0;

	for (size_t j=0; j<triggerNames.size(); j++) {
	  //if (!first) std::cout<<j<<" "<<triggerNames.triggerNames()[j]<<std::endl;
	  for (size_t i=0; i<HLTPaths_.size(); i++) {
	    if (indices[i]!=triggerNames.size()) continue;
	    boost::regex e(HLTPaths_[i]);
	    boost::match_results<std::string::const_iterator> m;
	    if (boost::regex_match(triggerNames.triggerName(j), m, e, boost::match_default)) {
	      indices[i]=j;
	      nMatches++;
	      if (!first) std::cout<<HLTPaths_[i]<<" matched "<<triggerNames.triggerName(j)<<std::endl;
	    }
	  }
	  if (nMatches==HLTPaths_.size()) break;
	}

	for (size_t i=0; i<indices.size(); i++) {
	  if (indices[i]>=triggerResults->size()) {
	    if (andOr_) continue;
	    hltPassed=0;
	    break;
	  }
	  bool pass=triggerResults->accept(indices[i]);
	  if (andOr_) {
	    if (!pass) continue;
	    hltPassed=1;
	    break;
	  }
	  // in case of AND between triggers (which have valid results)
	  if (!pass) {
	    hltPassed=0;
	    break;
	  }
	  hltPassed=1;
	}
      } // triggerResults.isValid()

    } // HLTPaths_.size()

    if (!first) std::cout<<"Trigger result is: "<<hltPassed<<std::endl;
    first=true;

    // read input
    edm::Handle<edm::View<reco::Candidate> > probes;
    iEvent.getByLabel(probes_,  probes);
    
    // prepare vector for output    
    edm::View<reco::Candidate>::const_iterator probe, endprobes = probes->end(); 
    std::vector<float> values(probes->size(), float(hltPassed));

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*valMap);
    filler.insert(probes, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTResultProducer);
