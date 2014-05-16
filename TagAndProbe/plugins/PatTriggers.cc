/*
  \class        ElectronDeltaPfRecoPt
  \description: Print out variables of PatTriggers

*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

class PatTriggers : public edm::EDProducer {
public:
  explicit PatTriggers(const edm::ParameterSet & iConfig);
  virtual ~PatTriggers() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

private:
  edm::InputTag src_;

};

PatTriggers::PatTriggers(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")) {}

PatTriggers::~PatTriggers() {}

void PatTriggers::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Read input
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggers;
  iEvent.getByLabel(src_,triggers);

  // Loop Trigger Objects
  unsigned ntrig=triggers->size();
  //int nele=0;
  
  std::cout<<"Trigger paths in event: "<<std::endl;
  for (unsigned itrig=0; itrig<ntrig; ++itrig) {
    pat::TriggerObjectStandAloneRef trig(triggers, itrig);


    std::vector<std::string> paths=trig->pathNames(0,1);
    for (size_t i=0; i<paths.size(); i++) {
      std::cout<<" "<<itrig<<" "<<i<<" "<<paths[i]<<std::endl;
    }
  }

  /*
    if (trig->type(82)) { // Get TriggerElectrons
      nele++;
      std::cout<<"Electron "<<nele<<std::endl;
      // 2012 SingleElectron Triggers
      if (trig->path("HLT_Ele22_CaloIdL_CaloIsoVL_v*"                                         ,1,0))
	std::cout<<"    HLT_Ele22_CaloIdL_CaloIsoVL_v*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_v*"                                         ,1,0))
	std::cout<<"    HLT_Ele27_CaloIdL_CaloIsoVL_v*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                        ,1,0))
	std::cout<<"    HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_WP80_v*"                                                      ,1,0))
	std::cout<<"    HLT_Ele27_WP80_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_WP80_PFMET_MT50_v*"                                           ,1,0))
	std::cout<<"    HLT_Ele27_WP80_PFMET_MT50_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_WP80_CentralPFJet80_v*"                                       ,1,0))
	std::cout<<"    HLT_Ele27_WP80_CentralPFJet80_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele27_WP80_WCandPt80_v*"                                            ,1,0))
	std::cout<<"    HLT_Ele27_WP80_WCandPt80_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele30_CaloIdVT_TrkIdT_v*"                                           ,1,0))
	std::cout<<"    HLT_Ele30_CaloIdVT_TrkIdT_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                        ,1,0))
	std::cout<<"    HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele65_CaloIdVT_TrkIdT_v*"                                           ,1,0))
	std::cout<<"    HLT_Ele65_CaloIdVT_TrkIdT_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele80_CaloIdVT_TrkIdT_v*"                                           ,1,0))
	std::cout<<"    HLT_Ele80_CaloIdVT_TrkIdT_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele80_CaloIdVT_GsfTrkIdT_v*"                                        ,1,0))
	std::cout<<"    HLT_Ele80_CaloIdVT_GsfTrkIdT_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele90_CaloIdVT_GsfTrkIdT_v*"                                        ,1,0))
	std::cout<<"    HLT_Ele90_CaloIdVT_GsfTrkIdT_*,1,0"<<std::endl;
      if (trig->path("HLT_Ele100_CaloIdVT_TrkIdT_v*"                                          ,1,0))
	std::cout<<"    HLT_Ele100_CaloIdVT_TrkIdT_*,1,0"<<std::endl;
      // 2012 EleHad Triggers
      if (trig->path("HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*"                                      ,1,0))
	std::cout<<"    HLT_Ele8_CaloIdT_TrkIdVL_Jet30_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"     ,0,0))
	std::cout<<"    HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"     ,0,0))
	std::cout<<"    HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*" ,0,0))
	std::cout<<"    HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*" ,0,0))
	std::cout<<"    HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*"    ,0,0))
	std::cout<<"    HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*"    ,0,0))
	std::cout<<"    HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",0,0))
	std::cout<<"    HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",0,0))
	std::cout<<"    HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*"                              ,0,0))
	std::cout<<"    HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*"                              ,0,0))
	std::cout<<"    HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*"                          ,0,0))
	std::cout<<"    HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_*,0,0"<<std::endl;
      if (trig->path("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*"                          ,0,0))
	std::cout<<"    HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_*,0,0"<<std::endl;

      // Cross trigger filters
      if (trig->filter("hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET45"))
	std::cout<<"      hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET45"<<std::endl;
      if (trig->filter("hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET50"))
	std::cout<<"      hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET50"<<std::endl;
      if (trig->filter("hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET45"))
	std::cout<<"      hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET45"<<std::endl;
      if (trig->filter("hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET50"))
	std::cout<<"      hltEle5CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET50"<<std::endl;
      if (trig->filter("hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET45"))
	std::cout<<"      hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET45"<<std::endl;
      if (trig->filter("hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET50"))
	std::cout<<"      hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350PFMET50"<<std::endl;
      if (trig->filter("hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET45"))
	std::cout<<"      hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET45"<<std::endl;
      if (trig->filter("hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET50"))
	std::cout<<"      hltElectron15CaloIdTCaloIsoVLTrkIdTTrkIsoVLCleanedPFHT350NoPUPFMET50"<<std::endl;
      if (trig->filter("hltElectron40CaloIdTTrkIdTCleanedPFHT300"))
	std::cout<<"      hltElectron40CaloIdTTrkIdTCleanedPFHT300"<<std::endl;
      if (trig->filter("hltElectron60CaloIdTTrkIdTCleanedPFHT300"))
	std::cout<<"      hltElectron60CaloIdTTrkIdTCleanedPFHT300"<<std::endl;
      if (trig->filter("hltElectron40CaloIdTTrkIdTCleanedPFHT300NoPU"))
	std::cout<<"      hltElectron40CaloIdTTrkIdTCleanedPFHT300NoPU"<<std::endl;
      if (trig->filter("hltElectron60CaloIdTTrkIdTCleanedPFHT300NoPU"))
	std::cout<<"      hltElectron60CaloIdTTrkIdTCleanedPFHT300NoPU"<<std::endl;

      nele++;
      if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",0,0) || trig->path("HLT_Ele27_WP80_v*",0,0)) {
        std::cout<<std::endl;
	if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                        ,0,0))
	  std::cout<<"  HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,0,0"<<std::endl;
	if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                        ,1,0))
	  std::cout<<"  HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,1,0"<<std::endl;
	if (trig->path("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                        ,0,1))
	  std::cout<<"  HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_*,0,1"<<std::endl;
	if (trig->path("HLT_Ele27_WP80_v*"                                                      ,0,0))
	  std::cout<<"    HLT_Ele27_WP80_*,0,0"<<std::endl;
	if (trig->path("HLT_Ele27_WP80_v*"                                                      ,1,0))
	  std::cout<<"    HLT_Ele27_WP80_*,1,0"<<std::endl;
	if (trig->path("HLT_Ele27_WP80_v*"                                                      ,0,1))
	  std::cout<<"    HLT_Ele27_WP80_*,0,1"<<std::endl;
	
	// std::cout<<"Filters:"<<std::endl;
	// size_t nlabel = trig->filterLabels().size();
	// for (size_t i=0; i<nlabel; i++)
	//   std::cout<<"  "<<trig->filterLabels()[i]<<std::endl;
        std::cout<<"--------------------------------------------------------------"<<std::endl;
      }
    }
  }
  */
  
  /*
      /// Gets all HLT filter labels
      std::vector< std::string > filterLabels() const { return filtersOrConditions(); };
      /// Gets all L1 condition names
      std::vector< std::string > conditionNames() const { return filtersOrConditions(); };
      /// Gets all HLT path names
      std::vector< std::string > pathNames( bool pathLastFilterAccepted = false, bool pathL3FilterAccepted = true ) const { return pathsOrAlgorithms( pathLastFilterAccepted, pathL3FilterAccepted ); };
      /// Gets all L1 algorithm names
      std::vector< std::string > algorithmNames( bool algoCondAccepted = true ) const { return pathsOrAlgorithms( algoCondAccepted, false ); };
      /// Gets the pat::TriggerObject (parent class)
      TriggerObject triggerObject();
      /// Checks, if a certain HLT filter label is assigned
      bool hasFilterLabel( const std::string & filterLabel ) const { return hasFilterOrCondition( filterLabel ); };
      /// Checks, if a certain L1 condition name is assigned
      bool hasConditionName( const std::string & conditionName ) const { return hasFilterOrCondition( conditionName ); };
      /// Checks, if a certain HLT path name is assigned
      bool hasPathName( const std::string & pathName, bool pathLastFilterAccepted = false, bool pathL3FilterAccepted = true ) const { return hasPathOrAlgorithm( pathName, pathLastFilterAccepted, pathL3FilterAccepted ); };
      /// Checks, if a certain L1 algorithm name is assigned
      bool hasAlgorithmName( const std::string & algorithmName, bool algoCondAccepted = true ) const { return hasPathOrAlgorithm( algorithmName, algoCondAccepted, false ); };
      /// Checks, if a certain label of original collection is assigned (method overrides)
      virtual bool hasCollection( const std::string & collName ) const;
      virtual bool hasCollection( const edm::InputTag & collName ) const { return hasCollection( collName.encode() ); };
      /// Checks, if the usage indicator vector has been filled
      bool hasPathLastFilterAccepted() const { return hasLastFilter(); };
      bool hasAlgoCondAccepted() const { return hasLastFilter(); };
      bool hasPathL3FilterAccepted() const { return hasL3Filter(); };

      /// Special methods for the cut string parser
      /// - argument types usable in the cut string parser
      /// - short names for readable configuration files

      /// Calls 'hasFilterLabel(...)'
      bool filter( const std::string & filterLabel ) const { return hasFilterLabel( filterLabel ); };
      /// Calls 'hasConditionName(...)'
      bool cond( const std::string & conditionName ) const { return hasConditionName( conditionName ); };
      /// Calls 'hasPathName(...)'
      bool path( const std::string & pathName, unsigned pathLastFilterAccepted = 0, unsigned pathL3FilterAccepted = 1 ) const { return hasPathName( pathName, bool( pathLastFilterAccepted ), bool( pathL3FilterAccepted ) ); };
      /// Calls 'hasAlgorithmName(...)'
      bool algo( const std::string & algorithmName, unsigned algoCondAccepted = 1 ) const { return hasAlgorithmName( algorithmName, bool( algoCondAccepted ) ); };
      /// Calls 'hasCollection(...)' (method override)
      virtual bool coll( const std::string & collName ) const { return hasCollection( collName ); };





      /// Get the label of the collection the trigger object originates from
      std::string collection() const { return collection_; };
      /// Special methods for 'l1extra' particles

      /// Special methods for the cut string parser
      /// - argument types usable in the cut string parser
      /// - short names for readable configuration files

      /// Calls 'hasCollection(...)'
      virtual bool coll( const std::string & collName ) const { return hasCollection( collName ); };
      /// Calls 'hasTriggerObjectType(...)'
      bool type( trigger::TriggerObjectType triggerObjectType ) const { return hasTriggerObjectType( triggerObjectType ); };
      bool type( int                        triggerObjectType ) const { return hasTriggerObjectType( trigger::TriggerObjectType ( triggerObjectType ) ); };


  */
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatTriggers);
