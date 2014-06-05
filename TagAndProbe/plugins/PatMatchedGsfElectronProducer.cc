#include "DataFormats/PatCandidates/interface/Electron.h"
#include "PhysicsTools/TagAndProbe/plugins/ObjectViewMatcher.cc"

typedef ObjectViewMatcher<reco::GsfElectron, pat::Electron> PatMatchedGsfElectronProducer;
DEFINE_FWK_MODULE(PatMatchedGsfElectronProducer);
