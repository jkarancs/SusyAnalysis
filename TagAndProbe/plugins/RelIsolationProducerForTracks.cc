#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

class RelIsolationProducerForTracks : public edm::EDProducer {
public:
  RelIsolationProducerForTracks(const edm::ParameterSet & );
private:
  void produce(edm::Event& event, const edm::EventSetup& setup);

  edm::InputTag tracks_;
  edm::InputTag highPtTracks_;
  edm::InputTag trkisoDeps_;
  edm::InputTag ecalisoDeps_;
  edm::InputTag hcalisoDeps_;
  double trackPtMin_;
  double coneSize_;

};


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include <iostream>
#include <iterator>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

typedef edm::ValueMap<float> TkIsoMap;

RelIsolationProducerForTracks::RelIsolationProducerForTracks(const ParameterSet & pset) :
  tracks_( pset.getParameter<InputTag>( "tracks" ) ),
  highPtTracks_( pset.getParameter<InputTag>( "highPtTracks" ) ),
  trkisoDeps_( pset.getParameter<InputTag>( "trkisoDeps" ) ),
  ecalisoDeps_( pset.getParameter<InputTag>( "ecalisoDeps" ) ),
  hcalisoDeps_( pset.getParameter<InputTag>( "hcalisoDeps" ) ),
  trackPtMin_( pset.getParameter<double>( "trackPtMin" ) ),
  coneSize_( pset.getParameter<double>( "coneSize" ) )
 {
  produces<TkIsoMap>();
 }

void RelIsolationProducerForTracks::produce(Event & event, const EventSetup & setup) {
  std::auto_ptr<TkIsoMap> relIsolations(new TkIsoMap);
  TkIsoMap::Filler filler(*relIsolations);
  {
    Handle<CandidateView> tracks;
    event.getByLabel(tracks_, tracks);

    Handle<CandidateView> highPtTracks;
    event.getByLabel(highPtTracks_, highPtTracks);

    Handle<IsoDepositMap> trkisoDeps;
    event.getByLabel(trkisoDeps_, trkisoDeps);

    Handle<IsoDepositMap> ecalisoDeps;
    event.getByLabel(ecalisoDeps_, ecalisoDeps);

    Handle<IsoDepositMap> hcalisoDeps;
    event.getByLabel(hcalisoDeps_, hcalisoDeps);

    int nTracks = tracks->size();
    int nHighPtTracks = highPtTracks->size();
    std::vector<double> iso(nTracks);

    OverlapChecker overlap;

    for(int i = 0; i < nTracks; ++i ) {
      const Candidate & tkCand = (*tracks)[ i ];
      double relIso = - 1.0;
      if( tkCand.pt() > trackPtMin_) {
        for(int j = 0; j < nHighPtTracks; ++j ) {
          const Candidate & highPtTkCand = (*highPtTracks)[ j ];
          if(overlap(tkCand, highPtTkCand) ) {
            CandidateBaseRef tkRef = highPtTracks->refAt(j);
            const IsoDeposit &trkisoDep = (*trkisoDeps)[tkRef];
            const IsoDeposit &ecalisoDep = (*ecalisoDeps)[tkRef];
            const IsoDeposit &hcalisoDep = (*hcalisoDeps)[tkRef];
            relIso = trkisoDep.depositWithin(coneSize_) + ecalisoDep.depositWithin(coneSize_) + hcalisoDep.depositWithin(coneSize_);
            relIso /= tkCand.pt();
            break;
          }
        }
      }
      iso[i] = relIso;
    }
    filler.insert(tracks, iso.begin(), iso.end());
  }

  // really fill the association map
  filler.fill();
  event.put(relIsolations);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RelIsolationProducerForTracks);

