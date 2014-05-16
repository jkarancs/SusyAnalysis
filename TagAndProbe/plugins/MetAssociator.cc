#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"

template<typename T>
class MetAssociator : public edm::EDProducer {
    public:
        explicit MetAssociator(const edm::ParameterSet & iConfig);
        virtual ~MetAssociator() ;

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::InputTag probes_;
        edm::InputTag metTag_;
};

template<typename T>
MetAssociator<T>::MetAssociator(const edm::ParameterSet & iConfig) :
    probes_(iConfig.getParameter<edm::InputTag>("probes")),
    metTag_(iConfig.getParameter<edm::InputTag>("metTag")) {
    produces<edm::ValueMap<float> >();
}


template<typename T>
MetAssociator<T>::~MetAssociator()
{
}

template<typename T>
void 
MetAssociator<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > probes;
    Handle<View<T> > metHandle;
    iEvent.getByLabel(probes_,  probes);
    iEvent.getByLabel(metTag_, metHandle);
    
    // fill
    float met = 0.0;
    if (metHandle.isValid()) met = metHandle->front().pt();

    // prepare vector for output    
    std::vector<float> values(probes->size(), met);

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(probes, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}


typedef MetAssociator<pat::MET> PatMetAssociator;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatMetAssociator);
