
#ifndef PhysicsTools_TagAndProbe_PFJetSelector_h
#define PhysicsTools_TagAndProbe_PFJetSelector_h

#include "FWCore/Framework/interface/EDFilter.h"

#include "DataFormats/Common/interface/RefVector.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"


#include <vector>


class PFJetSelector : public edm::EDFilter {
public:


	PFJetSelector( edm::ParameterSet const & params ) :
		edm::EDFilter( ),
		src_( params.getParameter<edm::InputTag>("src") ),
		cut_( params.getParameter<std::string>("cut") ),
		filter_(false),
		selector_( cut_ )
	{
		produces< std::vector<reco::PFJet> >();
		produces<reco::PFCandidateCollection > ("pfCandidates");

		if ( params.exists("filter") ) {
			filter_ = params.getParameter<bool>("filter");
		}
	}

	virtual ~PFJetSelector() {}

	virtual void beginJob() {}
	virtual void endJob() {}

	virtual bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

		std::auto_ptr< std::vector<reco::PFJet> > pfJets ( new std::vector<reco::PFJet>() );

		edm::Handle< edm::View<reco::PFJet> > h_jets;
		iEvent.getByLabel( src_, h_jets );

		for ( edm::View<reco::PFJet>::const_iterator ibegin = h_jets->begin(),
				iend = h_jets->end(), ijet = ibegin;
				ijet != iend; ++ijet ) {

			// Check the selection
			if ( selector_(*ijet) ) {
				// Add the jets that pass to the output collection
				pfJets->push_back( *ijet );

			}
		}


		// put genEvt  in Event
		bool pass = pfJets->size() > 0;
		iEvent.put(pfJets);

		if ( filter_ )
			return pass;
		else
			return true;
	}

protected:
	edm::InputTag                  src_;
	std::string                    cut_;
	bool                           filter_;
	StringCutObjectSelector<reco::PFJet> selector_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFJetSelector);


#endif
