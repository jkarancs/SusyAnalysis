#ifndef SusyAnalysis_TagAndProbe_TreeFillerModified_h
#define SusyAnalysis_TagAndProbe_TreeFillerModified_h

#include "PhysicsTools/TagAndProbe/interface/BaseTreeFiller.h"

namespace tnp {
class TPTreeFillerModified : public BaseTreeFiller {
    public:
        TPTreeFillerModified(const edm::ParameterSet config);
        ~TPTreeFillerModified();

        // We declare 'const' the methods which don't change the configuration
        void init(const edm::Event &iEvent) const ;
        void fill(const reco::CandidateBaseRef &probe, double mass, bool mcTrue=false, float mcMass=0.0) const ;

    protected:
        /// extra branch for the mass
        mutable float  mass_;
        /// extra branch for the mc truth
        mutable int32_t mcTrue_;
        /// extra branch for the mc-truth mass
        mutable float  mcMass_;
};
}

#endif
