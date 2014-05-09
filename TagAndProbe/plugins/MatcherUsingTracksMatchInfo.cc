// modification of MatcherUsingTracks (needs MuonAnalysis/MuonAssociators/ to run)

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "MuonAnalysis/MuonAssociators/interface/MatcherUsingTracksAlgorithm.h"




// template-related workaround for bug in OwnVector+Ptr
namespace edm { using std::advance; }

namespace pat {
  
  class MatcherUsingTracksMatchInfo : public edm::EDProducer {
  public:
    explicit MatcherUsingTracksMatchInfo(const edm::ParameterSet & iConfig);
    virtual ~MatcherUsingTracksMatchInfo() { }
    
    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
    
  private:
    /// Labels for input collections
    edm::InputTag src_, matched_;
    
    /// The real workhorse
    MatcherUsingTracksAlgorithm algo_;
    
    /// Some extra configurables
    bool dontFailOnMissingInput_;
    bool writeExtraPATOutput_;
    
    /// Store extra information in a ValueMap
    template<typename T>
    void storeValueMap(edm::Event &iEvent, 
		       const edm::Handle<edm::View<reco::Candidate> > & handle,
		       const std::vector<T> & values,
		       const std::string    & label) const ;
    
  };
  
} // namespace

pat::MatcherUsingTracksMatchInfo::MatcherUsingTracksMatchInfo(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")),
  matched_(iConfig.getParameter<edm::InputTag>("matched")),
  algo_(iConfig),
  dontFailOnMissingInput_(iConfig.existsAs<bool>("dontFailOnMissingInput") ? iConfig.getParameter<bool>("dontFailOnMissingInput") : false),
  writeExtraPATOutput_(iConfig.existsAs<bool>("writeExtraPATOutput") ? iConfig.getParameter<bool>("writeExtraPATOutput") : false)
{
  // this is the basic output (edm::Association is not generic)
  produces<edm::ValueMap<reco::CandidatePtr> >(); 
  if (writeExtraPATOutput_) {
    // this is the crazy stuff to get the same with UserData
    produces<edm::OwnVector<pat::UserData> >();
    produces<edm::ValueMap<edm::Ptr<pat::UserData> > >();
  }
  // this is the extra stuff
  if (algo_.hasMetrics()) {
    produces<edm::ValueMap<float> >("deltaR");
    produces<edm::ValueMap<float> >("deltaEta");
    produces<edm::ValueMap<float> >("deltaPhi");
    produces<edm::ValueMap<float> >("deltaLocalPos");
    produces<edm::ValueMap<float> >("deltaPtRel");
    if (algo_.hasChi2()) {
      produces<edm::ValueMap<float> >("chi2");
    }
  } else {
    produces<edm::ValueMap<int> >("matched");
  }
  produces<edm::ValueMap<float> >("ptErrByPt2");
  produces<edm::ValueMap<float> >("deltaPtRecoPF");
  produces<edm::ValueMap<float> >("absdeltaPtRecoPF");
  produces<edm::ValueMap<float> >("pfreliso");
}

void 
pat::MatcherUsingTracksMatchInfo::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;
  
  algo_.init(iSetup);
  
  Handle<View<reco::Candidate> > src, matched;
  
  iEvent.getByLabel(src_, src);
  iEvent.getByLabel(matched_, matched);
  
  // declare loop variables and some intermediate stuff
  View<reco::Candidate>::const_iterator itsrc, edsrc;
  int isrc, nsrc = src->size();
  
  // working and output variables
  vector<int>   match(nsrc, -1);
  vector<float> deltaRs(nsrc, 999);
  vector<float> deltaEtas(nsrc, 999);
  vector<float> deltaPhis(nsrc, 999);
  vector<float> deltaPtRel(nsrc, 999);
  vector<float> deltaLocalPos(nsrc, 999);
  vector<float> chi2(nsrc, 999999);
  vector<float> ptErrByPt2(nsrc, 999999);
  vector<float> deltaPtRecoPF(nsrc, 999999);
  vector<float> absdeltaPtRecoPF(nsrc, 999999);
  vector<float> pfreliso(nsrc, 999999);
  
  // don't try matching if the input collection is missing and the module is configured to fail silently
  if (!(matched.failedToGet() && dontFailOnMissingInput_)) {
    // loop on the source collection, and request for the match
    for (itsrc = src->begin(), edsrc = src->end(), isrc = 0; itsrc != edsrc; ++itsrc, ++isrc) {
      match[isrc] = algo_.match(*itsrc, *matched, deltaRs[isrc], deltaEtas[isrc], deltaPhis[isrc], deltaLocalPos[isrc], deltaPtRel[isrc], chi2[isrc]); 
    }
  }
  
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  try {
    iEvent.getByLabel("particleFlow",pfCandidates);
  } catch ( cms::Exception & e ) {
    cout <<"[MatcherUsingTracksMatchInfo] Error: " << e.what() << endl;
  }
  reco::PFCandidateCollection::const_iterator pfit, pfbegin = pfCandidates->begin(), pfend = pfCandidates->end(), curr_pf = pfend;
  
  std::vector<reco::CandidatePtr> ptrs(nsrc);
  for (isrc = 0; isrc < nsrc; ++isrc) {
    if (match[isrc] != -1) {
      ptrs[isrc] = matched->ptrAt(match[isrc]);
      const pat::Muon* pMu = dynamic_cast< const pat::Muon* >(ptrs[isrc].get());
      if (pMu != 0) {
	ptErrByPt2[isrc] = pMu->globalTrack()->ptError()/pMu->pt()/pMu->pt();
        
	double mindR = 999999;
	double curr_dR = 999999;
	for (pfit = pfbegin;  pfit != pfend; ++pfit) {
	  if (pfit->particleId() != reco::PFCandidate::mu) continue;
	  if (pfit->pt() < 1.) continue;
	  curr_dR = reco::deltaR(pfit->p4(), pMu->p4());
	  if (curr_dR < mindR){
	    mindR = curr_dR;
	    curr_pf = pfit;
	  }
	}
	if (curr_pf !=  pfend) {
	  deltaPtRecoPF[isrc] = curr_pf->pt() > 10 ? fabs((pMu->pt()-curr_pf->pt())/pMu->pt()) : -9999;
	  absdeltaPtRecoPF[isrc] = curr_pf->pt() > 10 ? fabs((pMu->pt()-curr_pf->pt())) : -9999;
	  pfreliso[isrc] = (pMu->pfIsolationR03().sumChargedHadronPt + max(0., pMu->pfIsolationR03().sumNeutralHadronEt + pMu->pfIsolationR03().sumPhotonEt - 0.5*pMu->pfIsolationR03().sumPUPt ) ) / pMu->pt();
	}
      }
    }
  }
  storeValueMap<float>(iEvent, src, ptErrByPt2,      "ptErrByPt2");
  storeValueMap<float>(iEvent, src, deltaPtRecoPF,   "deltaPtRecoPF");
  storeValueMap<float>(iEvent, src, absdeltaPtRecoPF,"absdeltaPtRecoPF");
  storeValueMap<float>(iEvent, src, pfreliso,        "pfreliso");
  
  storeValueMap<reco::CandidatePtr>(iEvent, src, ptrs, "");
  
  if (writeExtraPATOutput_) {
    std::auto_ptr<edm::OwnVector<pat::UserData> > outUDVect(new edm::OwnVector<pat::UserData>());
    std::vector<int>                              idxUD(nsrc, -1);
    for (isrc = 0; isrc < nsrc; ++isrc) {
      if (match[isrc] != -1) {
	outUDVect->push_back(pat::UserData::make(ptrs[isrc]));
	idxUD[isrc] = outUDVect->size() - 1; 
      }
    }
    edm::OrphanHandle<edm::OwnVector<pat::UserData> > doneUDVect = iEvent.put(outUDVect);
    std::vector<edm::Ptr<pat::UserData> > ptrUD(nsrc);
    for (isrc = 0; isrc < nsrc; ++isrc) {
      if (idxUD[isrc] != -1) ptrUD[isrc] = edm::Ptr<pat::UserData>(doneUDVect, idxUD[isrc]);
    }
    storeValueMap<edm::Ptr<pat::UserData> >(iEvent, src, ptrUD, "");
  }
  
  if (algo_.hasMetrics()) {
    storeValueMap<float>(iEvent, src, deltaRs,   "deltaR");
    storeValueMap<float>(iEvent, src, deltaEtas, "deltaEta");
    storeValueMap<float>(iEvent, src, deltaPhis, "deltaPhi");
    storeValueMap<float>(iEvent, src, deltaLocalPos, "deltaLocalPos");
    storeValueMap<float>(iEvent, src, deltaPtRel,    "deltaPtRel");
    if (algo_.hasChi2()) {
      storeValueMap<float>(iEvent, src, chi2, "chi2");
    }
  } else  {
    std::vector<int> ismatched(nsrc, 0);
    for (isrc = 0; isrc < nsrc; ++isrc) {
      ismatched[isrc] = (match[isrc] != -1);
    }
    storeValueMap<int>(iEvent, src, ismatched, "matched");
  }
}

template<typename T>
void
pat::MatcherUsingTracksMatchInfo::storeValueMap(edm::Event &iEvent,
						const edm::Handle<edm::View<reco::Candidate> > & handle,
						const std::vector<T> & values,
						const std::string    & label) const {
  using namespace edm; using namespace std;
  auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(valMap, label);
}


#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(MatcherUsingTracksMatchInfo);
