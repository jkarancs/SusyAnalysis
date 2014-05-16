/*****************************************************************************
 * modification of /PhysicsTools/TagAndProbe/plugins/DeltaRNearestObjectComputer.cc
 * considers only objects > minDeltaR_ to mimic RA4 lepton-jet cross-cleaning
 *****************************************************************************/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

template<typename T>
class minCutDeltaRNearestObjectComputer : public edm::EDProducer {
    public:
        explicit minCutDeltaRNearestObjectComputer(const edm::ParameterSet & iConfig);
        virtual ~minCutDeltaRNearestObjectComputer() ;

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::InputTag probes_;            
        edm::InputTag objects_;
        StringCutObjectSelector<T,true> objCut_; // lazy parsing, to allow cutting on variables not in reco::Candidate class
        double minDeltaR_;
};

template<typename T>
minCutDeltaRNearestObjectComputer<T>::minCutDeltaRNearestObjectComputer(const edm::ParameterSet & iConfig) :
    probes_(iConfig.getParameter<edm::InputTag>("probes")),
    objects_(iConfig.getParameter<edm::InputTag>("objects")),
    objCut_(iConfig.existsAs<std::string>("objectSelection") ? iConfig.getParameter<std::string>("objectSelection") : "", true),
    minDeltaR_ (iConfig.getParameter<double>("minDeltaR"))
{
    produces<edm::ValueMap<float> >();
}


template<typename T>
minCutDeltaRNearestObjectComputer<T>::~minCutDeltaRNearestObjectComputer()
{
}

template<typename T>
void 
minCutDeltaRNearestObjectComputer<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > probes;
    iEvent.getByLabel(probes_,  probes);

    Handle<View<T> > objects;
    iEvent.getByLabel(objects_, objects);

    // prepare vector for output    
    std::vector<float> values;
    
    // fill
    View<reco::Candidate>::const_iterator probe, endprobes = probes->end();
    for (probe = probes->begin(); probe != endprobes; ++probe) {
        double dr2min = 10000;
        for (unsigned int iObj=0; iObj<objects->size(); iObj++) {
	  const T& obj = objects->at(iObj);
	  if (!objCut_(obj)) continue;
            double dr2 = deltaR2(*probe, obj);
            if (dr2 < dr2min && dr2 > minDeltaR_*minDeltaR_) { dr2min = dr2; }
        }
        values.push_back(sqrt(dr2min));
    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(probes, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

typedef minCutDeltaRNearestObjectComputer<reco::Candidate>     minCutDeltaRNearestCandidateComputer;
typedef minCutDeltaRNearestObjectComputer<reco::Muon>          minCutDeltaRNearestMuonComputer;
typedef minCutDeltaRNearestObjectComputer<reco::Electron>      minCutDeltaRNearestElectronComputer;
typedef minCutDeltaRNearestObjectComputer<reco::GsfElectron>   minCutDeltaRNearestGsfElectronComputer;
typedef minCutDeltaRNearestObjectComputer<reco::Photon>        minCutDeltaRNearestPhotonComputer;
typedef minCutDeltaRNearestObjectComputer<reco::Jet>           minCutDeltaRNearestJetComputer;
typedef minCutDeltaRNearestObjectComputer<pat::Jet>            minCutDeltaRNearestPatJetComputer;


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(minCutDeltaRNearestCandidateComputer);          
DEFINE_FWK_MODULE(minCutDeltaRNearestMuonComputer);          
DEFINE_FWK_MODULE(minCutDeltaRNearestElectronComputer);          
DEFINE_FWK_MODULE(minCutDeltaRNearestGsfElectronComputer);          
DEFINE_FWK_MODULE(minCutDeltaRNearestPhotonComputer);          
DEFINE_FWK_MODULE(minCutDeltaRNearestJetComputer);         
DEFINE_FWK_MODULE(minCutDeltaRNearestPatJetComputer);          
