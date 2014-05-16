#include <cmath>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class MatchGsfElectronsToPAT : public edm::EDProducer {

    public:
        /// constructor
        MatchGsfElectronsToPAT(const edm::ParameterSet &iConfig) ;
        /// destructor
        ~MatchGsfElectronsToPAT() ;
    
        /// method to be called at each event
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;

    private:
        /// Input collection of electrons and of partice flow
        edm::InputTag gsfElectrons_, pat_;

        StringCutObjectSelector<pat::Electron> patCut_;

        /// Perform matching by reference (works only if gsfElectrons_ = "electrons")
        bool matchByReference_;

}; // C++ note: you need a ';' at the end of the class declaration.


MatchGsfElectronsToPAT::MatchGsfElectronsToPAT(const edm::ParameterSet &iConfig) :
    gsfElectrons_(iConfig.getParameter<edm::InputTag>("electrons")),
    pat_(iConfig.getParameter<edm::InputTag>("pat")),
    patCut_(iConfig.getParameter<std::string>("patCut")),
    matchByReference_(iConfig.getParameter<bool>("matchByReference"))
{
    produces<edm::RefToBaseVector<reco::GsfElectron> >();
}

MatchGsfElectronsToPAT::~MatchGsfElectronsToPAT()
{
}

void 
MatchGsfElectronsToPAT::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    using namespace edm; 

    Handle<View<reco::GsfElectron> > electrons;
    iEvent.getByLabel(gsfElectrons_, electrons);

    Handle<View<pat::Electron> > pat;
    iEvent.getByLabel(pat_, pat);
    
    /// prepare the vector for the output
    std::auto_ptr<RefToBaseVector<reco::GsfElectron> > out(new RefToBaseVector<reco::GsfElectron>());

    View<pat::Electron>::const_iterator patit, patbegin = pat->begin(), patend = pat->end();

    /// Now loop
    for (size_t i = 0, n = electrons->size(); i < n; ++i) {
        // read the reference to the gsfelectron track
    	reco::GsfTrackRef trackRef =  electrons->at(i).gsfTrack();

        /// Loop on pat
        for (patit = patbegin; patit != patend; ++patit) {
        	reco::GsfTrackRef  patTrackRef = patit->gsfTrack();
            if (patTrackRef.isNull()) continue;
            if (!patCut_(*patit)) continue;

            // Perform the matching
            if (matchByReference_) {
                if (patTrackRef.id() != trackRef.id()) throw cms::Exception("Configuration")
                    << "Cannot match by reference the electrons from " << gsfElectrons_.encode() << " (id " << trackRef.id() << ")"
                    << " with the ones referenced by " << pat_.encode() << " (id " << patTrackRef.id() << ")\n";
                if (patTrackRef.key() == trackRef.key()) {
                    out->push_back(electrons->refAt(i));
                    break;
                }
            } else if (std::abs(trackRef->eta() - patTrackRef->eta()) < 1e-5 &&
                       std::abs(trackRef->phi() - patTrackRef->phi()) < 1e-5 &&
                       std::abs(trackRef->pt()  - patTrackRef->pt() ) < 1e-5) {
                out->push_back(electrons->refAt(i));
                break;
            }
        }
    }

    // Write the output to the event
    iEvent.put(out);
}

/// Register this as a CMSSW Plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MatchGsfElectronsToPAT);
