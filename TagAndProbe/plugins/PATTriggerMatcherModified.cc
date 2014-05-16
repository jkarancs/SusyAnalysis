#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "CommonTools/UtilAlgos/interface/MasterCollectionHelper.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "CommonTools/UtilAlgos/interface/MatchByDR.h"
#include "PhysicsTools/PatAlgos/plugins/PATTriggerMatchSelector.h"
#include "CommonTools/UtilAlgos/interface/MatchByDR.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include "CommonTools/UtilAlgos/interface/MatchLessByDPt.h"
#include "CommonTools/UtilAlgos/interface/MatchByDEta.h"
#include "CommonTools/UtilAlgos/interface/MatchLessByDEta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

/// Default class for ranking matches: sorting by smaller distance.
class LessByMatchDistance {
public:
  LessByMatchDistance (const edm::ParameterSet& cfg,
		       const reco::CandidateView& c1, const pat::TriggerObjectStandAloneCollection& c2) :
    distance_(reco::modules::make<reco::MatchByDRDPt< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type >>(cfg)), c1_(c1), c2_(c2) {}
  bool operator() (const std::pair<size_t,size_t>& p1,
		   const std::pair<size_t,size_t>& p2) const {
    return 
      distance_(c1_[p1.first],c2_[p1.second])<
      distance_(c1_[p2.first],c2_[p2.second]);
  }
private:
  reco::MatchByDRDPt< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type > distance_;
  const reco::CandidateView& c1_;
  const pat::TriggerObjectStandAloneCollection& c2_;
};

class PATTriggerMatcherModified : public edm::EDProducer {
public:
  PATTriggerMatcherModified(const edm::ParameterSet & cfg);
  ~PATTriggerMatcherModified();
private:
  typedef typename reco::CandidateView::value_type T1;
  typedef typename pat::TriggerObjectStandAloneCollection::value_type T2;
  typedef edm::Association<pat::TriggerObjectStandAloneCollection> MatchMap;
  typedef std::pair<size_t, size_t> IndexPair;
  typedef std::vector<IndexPair> MatchContainer;
  void produce(edm::Event&, const edm::EventSetup&);
  edm::ParameterSet config_;
  edm::InputTag src_;
  edm::InputTag matched_;
  bool resolveAmbiguities_;            // resolve ambiguities after
                                       //   first pass?
  bool resolveByMatchQuality_;         // resolve by (global) quality
                                       //   of match (otherwise: by order
                                       //   of test candidates)
  bool select(const T1 & c1, const T2 & c2) const { 
    return select_(c1, c2); 
  }
  pat::PATTriggerMatchSelector< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type > select_;
  reco::MatchByDRDPt< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type > distance_;
};

PATTriggerMatcherModified::PATTriggerMatcherModified(const edm::ParameterSet & cfg) :
  config_(cfg),
  src_(cfg.getParameter<edm::InputTag>("src")),
  matched_(cfg.getParameter<edm::InputTag>("matched")), 
  resolveAmbiguities_(cfg.getParameter<bool>("resolveAmbiguities")),
  resolveByMatchQuality_(cfg.getParameter<bool>("resolveByMatchQuality")),
  select_(reco::modules::make<pat::PATTriggerMatchSelector< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type>>(cfg)), 
  distance_(reco::modules::make<reco::MatchByDRDPt< reco::CandidateView::value_type, pat::TriggerObjectStandAloneCollection::value_type >>(cfg)) {
  // definition of the product
  produces<MatchMap>();
  // set resolveByMatchQuality only if ambiguities are to be resolved
  resolveByMatchQuality_ = resolveByMatchQuality_ && resolveAmbiguities_;
}

PATTriggerMatcherModified::~PATTriggerMatcherModified() { }

void PATTriggerMatcherModified::produce(edm::Event& evt, const edm::EventSetup&) {
  using namespace edm;
  using namespace std;
  typedef std::pair<size_t, size_t> IndexPair;
  typedef std::vector<IndexPair> MatchContainer;
  // get collections from event
  Handle<pat::TriggerObjectStandAloneCollection> matched;  
  evt.getByLabel(matched_, matched);
  Handle<reco::CandidateView> cands;
  evt.getByLabel(src_, cands);
  // create product
  auto_ptr<MatchMap> matchMap(new MatchMap(matched));
  size_t size = cands->size();
  //printf( "matched.isValid(): %d matched.size(): %d cands.isValid(): %d cands->size(): %d\n", (int)matched.isValid(), (int)matched->size(), (int)cands.isValid(), (int)cands->size());

  if( size != 0 ) {
    //
    // create helpers
    //
    LessByMatchDistance comparator(config_,*cands,*matched);
    typename MatchMap::Filler filler(*matchMap);
    ::helper::MasterCollection<reco::CandidateView> master(cands);
    vector<int> indices(master.size(), -1);      // result: indices in target collection
    vector<bool> mLock(matched->size(),false);   // locks in target collection
    MatchContainer matchPairs;                   // container of matched pairs
    // loop over candidates
    for(size_t c = 0; c != size; ++ c) {
      const T1 & cand = (*cands)[c];
      // no global comparison of match quality -> reset the container for each candidate
      if ( !resolveByMatchQuality_ )  matchPairs.clear();
      // loop over target collection
      for(size_t m = 0; m != matched->size(); ++m) {
        const T2 & match = (* matched)[m];
        // check lock and preselection
        if ( !mLock[m] && select(cand, match)) {
	  std::cout<<"Match found muon["<<c<<"] trigger["<<m<<"] paths:"<<std::endl;
	  if (match.path("HLT_Mu5_v*",0,0))
	    std::cout<<"  path(HLT_Mu5_v*,g0,0))\n";
	  if (match.path("HLT_Mu8_v*",0,0))
	    std::cout<<"  path(HLT_Mu8_v*,0,0))\n";
	  if (match.path("HLT_Mu12_v*",0,0))
	    std::cout<<"  path(HLT_Mu12_v*,0,0))\n";
	  if (match.path("HLT_Mu17_v*",0,0))
	    std::cout<<"  path(HLT_Mu17_v*,0,0))\n";
	  if (match.path("HLT_Mu15_eta2p1_v*",0,0))
	    std::cout<<"  path(HLT_Mu15_eta2p1_v*,0,0))\n";
	  if (match.path("HLT_Mu24_eta2p1_v*",0,0))
	    std::cout<<"  path(HLT_Mu24_eta2p1_v*,0,0))\n";
	  if (match.path("HLT_IsoMu24_v*",0,0))
	    std::cout<<"  path(HLT_IsoMu24_v*,0,0))\n";
	  /*
	  std::cout<<"Match found electron["<<c<<"] trigger["<<m<<"] paths:"<<std::endl;
	  if (match.path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",1,0))
	    std::cout<<"  path(HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,1,0)"<<std::endl;
	  if (match.path("HLT_Ele27_WP80_v*",1,0))
	    std::cout<<"  path(HLT_Ele27_WP80_*,1,0)"<<std::endl;
	  
	  if (match.path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"     ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"     ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*" ,0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*" ,0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"    ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"    ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*"                              ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*"                              ,0,1))
	    std::cout<<"  path(HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*"                          ,0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*,0,1)"<<std::endl;
	  if (match.path("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*"                          ,0,1))
	    std::cout<<"  path(HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*,0,1)"<<std::endl;
	  */
	  // matching requirement fulfilled -> store pair of indices
          if ( distance_(cand,match) )  matchPairs.push_back(make_pair(c,m));
        }
      }
      // if match(es) found and no global ambiguity resolution requested
      if ( matchPairs.size()>0 && !resolveByMatchQuality_ ) {
	std::cout<<"Ambiguity found, order by reco candidate"<<std::endl;
        // look for and store best match
        size_t idx = master.index(c);
        assert(idx < indices.size());
        size_t index = min_element(matchPairs.begin(), matchPairs.end(), comparator)->second;
        indices[idx] = index;
        // if ambiguity resolution by order of (reco) candidates:
        //   lock element in target collection
        if ( resolveAmbiguities_ )  mLock[index] = true;
      }
    }
    // ambiguity resolution by global match quality (if requested)
    if ( resolveByMatchQuality_ ) {
      // sort container of all matches by quality
      sort(matchPairs.begin(),matchPairs.end(),comparator);
      vector<bool> cLock(master.size(),false);
      // loop over sorted container
      for ( MatchContainer::const_iterator i=matchPairs.begin();
            i!=matchPairs.end(); ++i ) {
        size_t c = (*i).first;
        size_t m = (*i).second;
        // accept only pairs without any lock
        if ( mLock[m] || cLock[c] )  continue;
        // store index to target collection and lock the two items
        size_t idx = master.index(c);
        assert(idx < indices.size());
        indices[idx] = m;
	std::cout<<"----------------------------------------------------------------------------"<<std::endl;
	std::cout<<"Resolution found muon["<<c<<"] trigger["<<m<<"] paths:"<<std::endl;
	if ((*matched)[m].path("HLT_Mu5_v*",0,0))
	  std::cout<<"  path(HLT_Mu5_v*,0,0))\n";
	if ((*matched)[m].path("HLT_Mu8_v*",0,0))
	  std::cout<<"  path(HLT_Mu8_v*,0,0))\n";
	if ((*matched)[m].path("HLT_Mu12_v*",0,0))
	  std::cout<<"  path(HLT_Mu12_v*,0,0))\n";
	if ((*matched)[m].path("HLT_Mu17_v*",0,0))
	  std::cout<<"  path(HLT_Mu17_v*,0,0))\n";
	if ((*matched)[m].path("HLT_Mu15_eta2p1_v*",0,0))
	  std::cout<<"  path(HLT_Mu15_eta2p1_v*,0,0))\n";
	if ((*matched)[m].path("HLT_Mu24_eta2p1_v*",0,0))
	  std::cout<<"  path(HLT_Mu24_eta2p1_v*,0,0))\n";
	if ((*matched)[m].path("HLT_IsoMu24_v*",0,0))
	  std::cout<<"  path(HLT_IsoMu24_v*,0,0))\n";
	/*
	std::cout<<"  Resolution found for electron["<<c<<"] trigger["<<m<<"] paths:"<<std::endl;
	if ((*matched)[m].path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",1,0))
	  std::cout<<"    path(HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,1,0)"<<std::endl;
	if ((*matched)[m].path("HLT_Ele27_WP80_v*",1,0))
	  std::cout<<"    path(HLT_Ele27_WP80_*,1,0)"<<std::endl;
	
	if ((*matched)[m].path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"     ,0,1))
	  std::cout<<"    path(HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"     ,0,1))
	  std::cout<<"    path(HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*" ,0,1))
	  std::cout<<"    path(HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*" ,0,1))
	  std::cout<<"    path(HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"    ,0,1))
	  std::cout<<"    path(HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"    ,0,1))
	  std::cout<<"    path(HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",0,1))
	  std::cout<<"    path(HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",0,1))
	  std::cout<<"    path(HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*,0,1)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*"                              ,0,0))
	  std::cout<<"    path(HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*,0,0)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*"                              ,0,0))
	  std::cout<<"    path(HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*,0,0)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*"                          ,0,0))
	  std::cout<<"    path(HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*,0,0)"<<std::endl;
	if ((*matched)[m].path("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*"                          ,0,0))
	  std::cout<<"    path(HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*,0,0)"<<std::endl;
	*/
        mLock[m] = true;
        cLock[c] = true;
      }
    }
    filler.insert(master.get(), indices.begin(), indices.end());
    filler.fill();
  }
  evt.put(matchMap);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATTriggerMatcherModified);
