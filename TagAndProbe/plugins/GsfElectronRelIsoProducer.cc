//

/*
  \class    	  GsfElectronRelIsoProducer
  \description:
  	  	 produces valuemaps for the Rho corrected gsfElectron Isolation
  	  	 To be used in T&P in the following way:

			 process.GSFRelIso = cms.EDProducer("GsfElectronRelIsoProducer",
				ElectronProbes  = cms.InputTag("gsfElectrons"),
				rhoIsoInputTag  = cms.InputTag("kt6PFJetsForIsolation", "rho"),
				isoValInputTags = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
				                                cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
								cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
				isMC            = cms.bool(False), # or True for MC
			 )

			 // in TagProbeFitTreeProducer
			 variables = cms.PSet(
			 	...
				reliso = cms.InputTag("GSFRelIso","reliso"),
				...
			 )

*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/Common/interface/ValueMap.h"

class GsfElectronRelIsoProducer : public edm::EDProducer {
 public:
  explicit GsfElectronRelIsoProducer(const edm::ParameterSet & iConfig);
  virtual ~GsfElectronRelIsoProducer() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
  
 private:
  edm::InputTag inputTagGsfElectrons_;
  edm::InputTag rhoIsoInputTag_;
  std::vector<edm::InputTag>  isoValInputTags_;
  bool isMC_;
  
};

GsfElectronRelIsoProducer::GsfElectronRelIsoProducer(const edm::ParameterSet & iConfig) :
  inputTagGsfElectrons_(iConfig.getParameter<edm::InputTag>("ElectronProbes")),
  rhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("rhoIsoInputTag")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  isoValInputTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");
  produces<edm::ValueMap<float> >("reliso");
}

GsfElectronRelIsoProducer::~GsfElectronRelIsoProducer()
{
}

void GsfElectronRelIsoProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  
  // Electron Collection (from Reco)
  edm::Handle<reco::GsfElectronCollection> h_gsfElectron;
  bool found = iEvent.getByLabel(inputTagGsfElectrons_,h_gsfElectron);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get GsfElectrons: "
       <<inputTagGsfElectrons_<<std::endl;
    edm::LogError("GsfElectronRelIsoProducer")<<err.str();
    throw cms::Exception( "MissingProduct", err.str());
  }
  
  // Iso Deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
    iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }

  // Rho for Isolation
  edm::Handle<double> rhoIso;
  iEvent.getByLabel(rhoIsoInputTag_, rhoIso);
  double rho = *(rhoIso.product());
  
  // prepare vector for output
  std::vector<double> v_reliso;

  // Loop on Electrons
  unsigned nele=h_gsfElectron->size();
  for(unsigned iele=0; iele<nele; ++iele) {
    reco::GsfElectronRef gsfel(h_gsfElectron, iele);
    
    double charged = (*(isoVals)[0])[gsfel];
    double photon = (*(isoVals)[1])[gsfel];
    double neutral = (*(isoVals)[2])[gsfel];
  
    // Effective Area
    double AEff03 = 0.0;
    if (isMC_) {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, fabs(gsfel->superCluster()->eta()), ElectronEffectiveArea::kEleEAFall11MC);
    } else {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, fabs(gsfel->superCluster()->eta()), ElectronEffectiveArea::kEleEAData2012);
    }
    
    float reliso = (charged + std::max(0., photon + neutral - rho * AEff03)) / gsfel->pt();
    v_reliso.push_back(reliso);
  }

  // convert into ValueMap and store
  std::auto_ptr<edm::ValueMap<float> > vm_reliso(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fill_reliso(*vm_reliso);
  fill_reliso.insert(h_gsfElectron, v_reliso.begin(), v_reliso.end());
  fill_reliso.fill();
  iEvent.put(vm_reliso,"reliso");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GsfElectronRelIsoProducer);
