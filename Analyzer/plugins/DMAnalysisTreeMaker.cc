#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "TFile.h"
#include "TTree.h"
#include "SusyAnalysis/Analyzer/interface/Data.h"

#include<vector>

//using namespace reco;
using namespace edm;
using namespace std;


class  DMAnalysisTreeMaker : public edm::EDAnalyzer 
{
public:
  explicit DMAnalysisTreeMaker( const edm::ParameterSet & );   

private:
  virtual void analyze(const edm::Event &, const edm::EventSetup & );
  vector<string> additionalVariables(string);
  string makeName(string label,string pref,string var);
  std::vector<edm::ParameterSet > physObjects;
  std::vector<edm::InputTag > variablesFloat, variablesInt, singleFloat, singleInt;
  edm::InputTag lhes_;
  TTree * treesBase;
  map<string, TTree * > treesByCategory;
  std::vector<string> names;
  map< string , float[50] > vfloats_values;
  map< string , int[50] > vints_values;
  map< string , string > obs_to_obj;
  map< string , string > obj_to_pref;
  
  map< string , float > float_values;
  map< string , int > int_values;
  map< string , int > sizes;

  map< string , bool > got_label; 
  map< string , int > max_instances; 

  map<string, edm::Handle<std::vector<float> > > h_floats;
  map<string, edm::Handle<std::vector<int> > > h_ints;
  map<string, edm::Handle<float> > h_float;
  map<string, edm::Handle<int> >h_int;

  string mu_label, ele_label, jetsAK4_label, jetsAK8_label, jetsCmsTopTag_label, subjetsAK8_label, subjetsCmsTopTag_label, met_label ; 

};


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig){
  
  mu_label = iConfig.getParameter<std::string >("muLabel");
  ele_label = iConfig.getParameter<std::string >("eleLabel");
  jetsAK4_label = iConfig.getParameter<std::string >("jetsAK4Label");
  jetsAK8_label = iConfig.getParameter<std::string >("jetsAK8Label");
  jetsCmsTopTag_label = iConfig.getParameter<std::string >("jetsCmsTopTagLabel");
  subjetsAK8_label = iConfig.getParameter<std::string >("subjetsAK8Label");
  subjetsCmsTopTag_label = iConfig.getParameter<std::string >("subjetsCmsTopTagLabel");
  met_label = iConfig.getParameter<std::string >("metLabel");

  physObjects = iConfig.template getParameter<std::vector<edm::ParameterSet> >("physicsObjects");
  
  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();
  
  Service<TFileService> fs;
  TFileDirectory DMTrees = fs->mkdir( "systematics_trees" );

  treesBase = new TTree("TreeBase", "TreeBase");
  
  lhes_ = iConfig.getParameter<edm::InputTag>( "lhes" );
  

  for (;itPsets!=physObjects.end();++itPsets){ 
    int maxI = itPsets->getUntrackedParameter< int >("maxInstances",10);
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI"); 
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    string namelabel = itPsets->getParameter< string >("label");
    string nameprefix = itPsets->getParameter< string >("prefix");
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();
    stringstream max_instance_str;
    max_instance_str<<maxI;
    max_instances[namelabel]=maxI;
    string nameobs = namelabel;
    string prefix = nameprefix;
    for (;itF != variablesFloat.end();++itF){
      
      string name=itF->instance()+"_"+itF->label();
      string nameinstance=itF->instance();
      string nameshort=itF->instance();
      //if(nameobs!=itF->label())cout<<" warning! label not matching the module! check if members of pset are consistent. If intentional , ignore this warning";
      //cout << "nameobs "<< nameobs<< " name " << name <<" nameshort " <<nameshort << " strsizw "<< (nameshort+"["+max_instance_str.str()+"]/F").c_str()<<endl;
      treesBase->Branch(nameshort.c_str(), &vfloats_values[name],(nameshort+"["+max_instance_str.str()+"]/F").c_str());
      names.push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;
      cout << " branching name "<< name<< " for obs "<< nameobs << " instance "<< nameinstance << endl;
    }
    
    for (;itI != variablesInt.end();++itI){
      string name=itI->instance()+"_"+itI->label();
      string nameshort=itF->instance();
      treesBase->Branch(nameshort.c_str(), &vints_values[name],(nameshort+"["+max_instance_str.str()+"]/I").c_str());
      names.push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;
    }  
    
    if (variablesFloat.size()>0){
      string nameshortv = namelabel;
      vector<string> extravars = additionalVariables(nameshortv);
      for(size_t addv = 0; addv < extravars.size();++addv){
	string name = nameshortv+extravars.at(addv);
	treesBase->Branch(name.c_str(), &vfloats_values[name],(name+"["+max_instance_str.str()+"]/F").c_str());
	obs_to_obj[name] = nameobs;
	obj_to_pref[nameobs] = prefix;
      }
    }
    names.push_back(nameobs);
    treesBase->Branch((nameobs+"_size").c_str(), &sizes[nameobs]);

    //Initialize single pset objects
    for (;itsF != singleFloat.end();++itsF){
      string name=itsF->instance()+itsF->label();
      string nameshort=itsF->instance();
      treesBase->Branch(nameshort.c_str(), &float_values[name]);
    }
    for (;itsI != singleInt.end();++itsI){
      string name=itsI->instance()+itsI->label();
      string nameshort=itsI->instance();
      treesBase->Branch(nameshort.c_str(), &int_values[name]);
    }

  }
  vector<string> extravars = additionalVariables("Event");
  for(size_t addv = 0; addv < extravars.size();++addv){
    string name = extravars.at(addv);
    treesBase->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
  }
}

void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();

  edm::Handle<LHEEventProduct > lhes;
  //iEvent.getByLabel(lhes_, lhes);

  
  //Part 1 taking the obs values from the edm file
  for (;itPsets!=physObjects.end();++itPsets){ 
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI"); 
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();
    string namelabel = itPsets->getParameter< string >("label");
    size_t maxInstance=(size_t)max_instances[namelabel];
    
    //Vectors of floats
    for (;itF != variablesFloat.end();++itF){
      string name=itF->instance()+"_"+itF->label();
      //string namelabel;
      float tmp =1.0;
      iEvent.getByLabel(*(itF),h_floats[name]);
      //      cout << "name "<< name <<endl;
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_floats[name]->size()){tmp = h_floats[name]->at(fi);}
	else { tmp = -9999.; }
	//cout << " setting name "<< name<< " at instance "<< fi <<" to value "<< tmp <<endl;
	vfloats_values[name][fi]=tmp;
      }
      sizes[namelabel]=h_floats[name]->size();
      //      cout<< " size for "<< namelabel <<" is then "<< sizes[namelabel]<<endl; 
    }
    
    //Vectors of ints
    for (;itI != variablesInt.end();++itI){
      string name=itI->instance()+"_"+itI->label();
      int tmp = 1;
      iEvent.getByLabel(*(itI),h_ints[name]);
      for (size_t fi = 0;fi < maxInstance;++fi){
	if(fi <h_ints[name]->size()){tmp = h_ints[name]->at(fi);}
	else { tmp = -9999.; }
	vints_values[name][fi]=tmp;
      }
    }  
    
    //Single floats/ints
    for (;itsF != singleFloat.end();++itsF){
      string name=itsF->instance()+itsF->label();
      iEvent.getByLabel(*(itsF),h_float[name]);
      float_values[name]=*h_float[name];
    }
    for (;itsI != singleInt.end();++itsI){
      string name=itsI->instance()+itsI->label();
      iEvent.getByLabel(*(itsI),h_int[name]);
      int_values[name]=*h_int[name];
    }
    
  }
  
  // Code added by Janos Karancsi
  // Objectives: Read all variables in B2G edm tuple, use same data structure as that of B2GPlotter
  // The data structure class can calculate all Analysis variables by itself
  // After calculation the new variables are saved also in the TTree, so plotting by hand is possible
  Data d;
  
  //Part 2: selection and analysis-level changes
  //Electrons:
  string pref = obj_to_pref[ele_label];
  d.ele.size = sizes[ele_label];
  for(int it = 0;it < max_instances[ele_label] ;++it){
    d.ele.Mass[it] = vfloats_values[makeName(ele_label,pref,"Mass")][it];
    d.ele.Pt[it] = vfloats_values[makeName(ele_label,pref,"Pt")][it];
    d.ele.Eta[it] = vfloats_values[makeName(ele_label,pref,"Eta")][it];
    d.ele.Y[it] = vfloats_values[makeName(ele_label,pref,"Y")][it];
    d.ele.Phi[it] = vfloats_values[makeName(ele_label,pref,"Phi")][it];
    d.ele.E[it] = vfloats_values[makeName(ele_label,pref,"E")][it];
    d.ele.Charge[it] = vfloats_values[makeName(ele_label,pref,"Charge")][it];
    d.ele.Iso03[it] = vfloats_values[makeName(ele_label,pref,"Iso03")][it];
    d.ele.D0[it] = vfloats_values[makeName(ele_label,pref,"D0")][it];
    d.ele.Dz[it] = vfloats_values[makeName(ele_label,pref,"Dz")][it];
    d.ele.dEtaIn[it] = vfloats_values[makeName(ele_label,pref,"dEtaIn")][it];
    d.ele.dPhiIn[it] = vfloats_values[makeName(ele_label,pref,"dPhiIn")][it];
    d.ele.HoE[it] = vfloats_values[makeName(ele_label,pref,"HoE")][it];
    d.ele.full5x5siee[it] = vfloats_values[makeName(ele_label,pref,"full5x5siee")][it];
    d.ele.ooEmooP[it] = vfloats_values[makeName(ele_label,pref,"ooEmooP")][it];
    d.ele.missHits[it] = vfloats_values[makeName(ele_label,pref,"missHits")][it];
    d.ele.hasMatchedConVeto[it] = vfloats_values[makeName(ele_label,pref,"hasMatchedConVeto")][it];
    d.ele.isEB[it] = vfloats_values[makeName(ele_label,pref,"isEB")][it];
    d.ele.isVeto[it] = vfloats_values[makeName(ele_label,pref,"isVeto")][it];
    d.ele.isLoose[it] = vfloats_values[makeName(ele_label,pref,"isLoose")][it];
    d.ele.isTight[it] = vfloats_values[makeName(ele_label,pref,"isTight")][it];
    d.ele.isMedium[it] = vfloats_values[makeName(ele_label,pref,"isMedium")][it];
    
    float pt = vfloats_values[makeName(ele_label,pref,"Pt")][it];
    float eta = vfloats_values[makeName(ele_label,pref,"Eta")][it];
    bool isTightElectron = pt>30 && fabs(eta)<2.1;
    bool isLooseElectron= pt>20;
    //    cout << " ele_label "<< ele_label <<  makeName(ele_label,pref,"Pt") << " eta is "<< endl;
    //cout <<" makename with pt "<<   makeName(ele_label,pref,"Pt")<<endl;
    //    pt = pt+1.;
    if(isTightElectron){
      ++float_values["nTightElectrons"];
    }
    if(isLooseElectron){
      ++float_values["nLooseElectrons"]; 
    }
  } 
  
  pref = obj_to_pref[mu_label];
  d.mu.size = sizes[mu_label];
  for(int it = 0;it < max_instances[mu_label] ;++it){
    d.mu.Mass[it] = vfloats_values[makeName(mu_label,pref,"Mass")][it];
    d.mu.Pt[it] = vfloats_values[makeName(mu_label,pref,"Pt")][it];
    d.mu.Eta[it] = vfloats_values[makeName(mu_label,pref,"Eta")][it];
    d.mu.Y[it] = vfloats_values[makeName(mu_label,pref,"Y")][it];
    d.mu.Phi[it] = vfloats_values[makeName(mu_label,pref,"Phi")][it];
    d.mu.E[it] = vfloats_values[makeName(mu_label,pref,"E")][it];
    d.mu.Charge[it] = vfloats_values[makeName(mu_label,pref,"Charge")][it];
    d.mu.Iso03[it] = vfloats_values[makeName(mu_label,pref,"Iso03")][it];
    d.mu.D0[it] = vfloats_values[makeName(mu_label,pref,"D0")][it];
    d.mu.D0err[it] = vfloats_values[makeName(mu_label,pref,"D0err")][it];
    d.mu.Dz[it] = vfloats_values[makeName(mu_label,pref,"Dz")][it];
    d.mu.Dzerr[it] = vfloats_values[makeName(mu_label,pref,"Dzerr")][it];
    d.mu.IsLooseMuon[it] = vfloats_values[makeName(mu_label,pref,"IsLooseMuon")][it];
    d.mu.IsSoftMuon[it] = vfloats_values[makeName(mu_label,pref,"IsSoftMuon")][it];
    d.mu.IsTightMuon[it] = vfloats_values[makeName(mu_label,pref,"IsTightMuon")][it];
    d.mu.IsPFMuon[it] = vfloats_values[makeName(mu_label,pref,"IsPFMuon")][it];
    d.mu.IsGlobalMuon[it] = vfloats_values[makeName(mu_label,pref,"IsGlobalMuon")][it];
    d.mu.IsTrackerMuon[it] = vfloats_values[makeName(mu_label,pref,"IsTrackerMuon")][it];
    d.mu.GlbTrkNormChi2[it] = vfloats_values[makeName(mu_label,pref,"GlbTrkNormChi2")][it];
    d.mu.NumberValidMuonHits[it] = vfloats_values[makeName(mu_label,pref,"NumberValidMuonHits")][it];
    d.mu.NumberMatchedStations[it] = vfloats_values[makeName(mu_label,pref,"NumberMatchedStations")][it];
    d.mu.NumberValidPixelHits[it] = vfloats_values[makeName(mu_label,pref,"NumberValidPixelHits")][it];
    d.mu.NumberTrackerLayers[it] = vfloats_values[makeName(mu_label,pref,"NumberTrackerLayers")][it];
    d.mu.NumberOfValidTrackerHits[it] = vfloats_values[makeName(mu_label,pref,"NumberOfValidTrackerHits")][it];
    d.mu.NumberOfPixelLayers[it] = vfloats_values[makeName(mu_label,pref,"NumberOfPixelLayers")][it];
    d.mu.InTrkNormChi2[it] = vfloats_values[makeName(mu_label,pref,"InTrkNormChi2")][it];
    d.mu.SumChargedHadronPt[it] = vfloats_values[makeName(mu_label,pref,"SumChargedHadronPt")][it];
    d.mu.SumNeutralHadronPt[it] = vfloats_values[makeName(mu_label,pref,"SumNeutralHadronPt")][it];
    d.mu.SumPhotonPt[it] = vfloats_values[makeName(mu_label,pref,"SumPhotonPt")][it];
    d.mu.SumPUPt[it] = vfloats_values[makeName(mu_label,pref,"SumPUPt")][it];
    d.mu.GenMuonY[it] = vfloats_values[makeName(mu_label,pref,"GenMuonY")][it];
    d.mu.GenMuonEta[it] = vfloats_values[makeName(mu_label,pref,"GenMuonEta")][it];
    d.mu.GenMuonPhi[it] = vfloats_values[makeName(mu_label,pref,"GenMuonPhi")][it];
    d.mu.GenMuonPt[it] = vfloats_values[makeName(mu_label,pref,"GenMuonPt")][it];
    d.mu.GenMuonE[it] = vfloats_values[makeName(mu_label,pref,"GenMuonE")][it];
    d.mu.GenMuonCharge[it] = vfloats_values[makeName(mu_label,pref,"GenMuonCharge")][it];
    d.mu.HLTmuonDeltaR[it] = vfloats_values[makeName(mu_label,pref,"HLTmuonDeltaR")][it];
    d.mu.HLTmuonPt[it] = vfloats_values[makeName(mu_label,pref,"HLTmuonPt")][it];
    d.mu.HLTmuonEta[it] = vfloats_values[makeName(mu_label,pref,"HLTmuonEta")][it];
    d.mu.HLTmuonPhi[it] = vfloats_values[makeName(mu_label,pref,"HLTmuonPhi")][it];
    d.mu.HLTmuonE[it] = vfloats_values[makeName(mu_label,pref,"HLTmuonE")][it];
    
    float pt = vfloats_values[makeName(mu_label,pref,"Pt")][it];
    float eta = vfloats_values[makeName(mu_label,pref,"Eta")][it];
    bool isTightMuon = pt>30 && fabs(eta)<2.1;
    bool isLooseMuon= pt>20;
    //    cout << " ele_label "<< mu_label <<  makeName(mu_label,pref,"Pt") << " eta is "<< endl;
    //cout <<" makename with pt "<<   makeName(ele_label,pref,"Pt")<<endl;
    //    pt = pt+1.;
    if(isTightMuon){
      ++float_values["nTightMuons"];
    }
    if(isLooseMuon){
      ++float_values["nLooseMuons"]; 
    }
  }
  
  //MET:
  pref = obj_to_pref[met_label];
  d.met.Pt = vfloats_values[makeName(met_label,pref,"Pt")][0];
  d.met.Phi = vfloats_values[makeName(met_label,pref,"Phi")][0];
  d.met.Px = vfloats_values[makeName(met_label,pref,"Px")][0];
  d.met.Py = vfloats_values[makeName(met_label,pref,"Py")][0];
  
  //Jets:
  pref = obj_to_pref[jetsAK4_label];
  d.jetsAK4.size = sizes[jetsAK4_label];
  for(int it = 0;it < max_instances[jetsAK4_label] ;++it){
    d.jetsAK4.Mass[it] = vfloats_values[makeName(jetsAK4_label,pref,"Mass")][it];
    d.jetsAK4.Pt[it] = vfloats_values[makeName(jetsAK4_label,pref,"Pt")][it];
    d.jetsAK4.Eta[it] = vfloats_values[makeName(jetsAK4_label,pref,"Eta")][it];
    d.jetsAK4.Y[it] = vfloats_values[makeName(jetsAK4_label,pref,"Y")][it];
    d.jetsAK4.Phi[it] = vfloats_values[makeName(jetsAK4_label,pref,"Phi")][it];
    d.jetsAK4.E[it] = vfloats_values[makeName(jetsAK4_label,pref,"E")][it];
    d.jetsAK4.Charge[it] = vfloats_values[makeName(jetsAK4_label,pref,"Charge")][it];
    d.jetsAK4.IsCSVL[it] = vfloats_values[makeName(jetsAK4_label,pref,"IsCSVL")][it];
    d.jetsAK4.IsCSVM[it] = vfloats_values[makeName(jetsAK4_label,pref,"IsCSVM")][it];
    d.jetsAK4.IsCSVT[it] = vfloats_values[makeName(jetsAK4_label,pref,"IsCSVT")][it];
    d.jetsAK4.CSV[it] = vfloats_values[makeName(jetsAK4_label,pref,"CSV")][it];
    d.jetsAK4.CSVV1[it] = vfloats_values[makeName(jetsAK4_label,pref,"CSVV1")][it];
    d.jetsAK4.GenPartonY[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonY")][it];
    d.jetsAK4.GenPartonEta[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonEta")][it];
    d.jetsAK4.GenPartonPhi[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonPhi")][it];
    d.jetsAK4.GenPartonPt[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonPt")][it];
    d.jetsAK4.GenPartonE[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonE")][it];
    d.jetsAK4.GenPartonCharge[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenPartonCharge")][it];
    d.jetsAK4.PartonFlavour[it] = vfloats_values[makeName(jetsAK4_label,pref,"PartonFlavour")][it];
    d.jetsAK4.HadronFlavour[it] = vfloats_values[makeName(jetsAK4_label,pref,"HadronFlavour")][it];
    d.jetsAK4.GenJetY[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetY")][it];
    d.jetsAK4.GenJetEta[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetEta")][it];
    d.jetsAK4.GenJetPhi[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetPhi")][it];
    d.jetsAK4.GenJetPt[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetPt")][it];
    d.jetsAK4.GenJetE[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetE")][it];
    d.jetsAK4.GenJetCharge[it] = vfloats_values[makeName(jetsAK4_label,pref,"GenJetCharge")][it];
    d.jetsAK4.HLTjetEta[it] = vfloats_values[makeName(jetsAK4_label,pref,"HLTjetEta")][it];
    d.jetsAK4.HLTjetPhi[it] = vfloats_values[makeName(jetsAK4_label,pref,"HLTjetPhi")][it];
    d.jetsAK4.HLTjetPt[it] = vfloats_values[makeName(jetsAK4_label,pref,"HLTjetPt")][it];
    d.jetsAK4.HLTjetE[it] = vfloats_values[makeName(jetsAK4_label,pref,"HLTjetE")][it];
    d.jetsAK4.HLTjetDeltaR[it] = vfloats_values[makeName(jetsAK4_label,pref,"HLTjetDeltaR")][it];
    d.jetsAK4.muonMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"muonMultiplicity")][it];
    d.jetsAK4.PhotonEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"PhotonEnergy")][it];
    d.jetsAK4.ElectronEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"ElectronEnergy")][it];
    d.jetsAK4.MuonEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"MuonEnergy")][it];
    d.jetsAK4.HFHadronEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"HFHadronEnergy")][it];
    d.jetsAK4.HFEMEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"HFEMEnergy")][it];
    d.jetsAK4.ChargedHadronMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"ChargedHadronMultiplicity")][it];
    d.jetsAK4.numberOfDaughters[it] = vfloats_values[makeName(jetsAK4_label,pref,"numberOfDaughters")][it];
    d.jetsAK4.chargedMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"chargedMultiplicity")][it];
    d.jetsAK4.neutralHadronMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"neutralHadronMultiplicity")][it];
    d.jetsAK4.neutralHadronEnergyFraction[it] = vfloats_values[makeName(jetsAK4_label,pref,"neutralHadronEnergyFraction")][it];
    d.jetsAK4.neutralEmEnergyFraction[it] = vfloats_values[makeName(jetsAK4_label,pref,"neutralEmEnergyFraction")][it];
    d.jetsAK4.chargedEmEnergyFraction[it] = vfloats_values[makeName(jetsAK4_label,pref,"chargedEmEnergyFraction")][it];
    d.jetsAK4.chargedHadronEnergyFraction[it] = vfloats_values[makeName(jetsAK4_label,pref,"chargedHadronEnergyFraction")][it];
    d.jetsAK4.photonMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"photonMultiplicity")][it];
    d.jetsAK4.electronMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"electronMultiplicity")][it];
    d.jetsAK4.HFHadronMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"HFHadronMultiplicity")][it];
    d.jetsAK4.HFEMMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"HFEMMultiplicity")][it];
    d.jetsAK4.ChargeMuEnergy[it] = vfloats_values[makeName(jetsAK4_label,pref,"ChargeMuEnergy")][it];
    d.jetsAK4.neutralMultiplicity[it] = vfloats_values[makeName(jetsAK4_label,pref,"neutralMultiplicity")][it];
    d.jetsAK4.SmearedPt[it] = vfloats_values[makeName(jetsAK4_label,pref,"SmearedPt")][it];
    d.jetsAK4.SmearedPEta[it] = vfloats_values[makeName(jetsAK4_label,pref,"SmearedPEta")][it];
    d.jetsAK4.SmearedPhi[it] = vfloats_values[makeName(jetsAK4_label,pref,"SmearedPhi")][it];
    d.jetsAK4.SmearedE[it] = vfloats_values[makeName(jetsAK4_label,pref,"SmearedE")][it];
    d.jetsAK4.JERup[it] = vfloats_values[makeName(jetsAK4_label,pref,"JERup")][it];
    d.jetsAK4.JERdown[it] = vfloats_values[makeName(jetsAK4_label,pref,"JERdown")][it];
    
    float pt = vfloats_values[makeName(jetsAK4_label,pref,"Pt")][it];
    float eta = vfloats_values[makeName(jetsAK4_label,pref,"Eta")][it];
    bool isTightJet = pt>40 && fabs(eta)<4.7;
    bool isTightJetCSVT = isTightJet && (vfloats_values[makeName(jetsAK4_label,pref,"CSV")][it]>0.7);
    //    bool isTightJetCSVT = isTightJet && (vfloats_values[makeName(jetsAK4_label,pref,"IsCSVT")][it]>0.);
    bool isLooseJet= pt>20;
    //    cout << " jetsAK4_label "<< jetsAK4_label <<  makeName(jetsAK4_label,pref,"Pt") << " eta is "<< endl;
    //cout <<" makename with pt "<<   makeName(jetsAK4_label,pref,"Pt")<<endl;
    //    pt = pt+1.;
    if(isTightJet){
      ++float_values["nTightJets"];
      if(isTightJetCSVT)++float_values["nCSVTJets"];
    }
    if(isLooseJet){
      ++float_values["nLooseJets"]; 
    }
  } 
  //float weightsign = lhes->hepeup().XWGTUP;
  //float LHEWeightSign = weightsign/fabs(weightsign);
  //float_values["LHEWeightSign"]=LHEWeightSign;

  pref = obj_to_pref[jetsAK8_label];
  d.jetsAK8.size = sizes[jetsAK8_label];
  for(int it = 0;it < max_instances[jetsAK8_label] ;++it){
    d.jetsAK8.Mass[it] = vfloats_values[makeName(jetsAK8_label,pref,"Mass")][it];
    d.jetsAK8.Pt[it] = vfloats_values[makeName(jetsAK8_label,pref,"Pt")][it];
    d.jetsAK8.Eta[it] = vfloats_values[makeName(jetsAK8_label,pref,"Eta")][it];
    d.jetsAK8.Y[it] = vfloats_values[makeName(jetsAK8_label,pref,"Y")][it];
    d.jetsAK8.Phi[it] = vfloats_values[makeName(jetsAK8_label,pref,"Phi")][it];
    d.jetsAK8.E[it] = vfloats_values[makeName(jetsAK8_label,pref,"E")][it];
    d.jetsAK8.Charge[it] = vfloats_values[makeName(jetsAK8_label,pref,"Charge")][it];
    d.jetsAK8.IsCSVL[it] = vfloats_values[makeName(jetsAK8_label,pref,"IsCSVL")][it];
    d.jetsAK8.IsCSVM[it] = vfloats_values[makeName(jetsAK8_label,pref,"IsCSVM")][it];
    d.jetsAK8.IsCSVT[it] = vfloats_values[makeName(jetsAK8_label,pref,"IsCSVT")][it];
    d.jetsAK8.CSV[it] = vfloats_values[makeName(jetsAK8_label,pref,"CSV")][it];
    d.jetsAK8.CSVV1[it] = vfloats_values[makeName(jetsAK8_label,pref,"CSVV1")][it];
    d.jetsAK8.GenPartonY[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonY")][it];
    d.jetsAK8.GenPartonEta[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonEta")][it];
    d.jetsAK8.GenPartonPhi[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonPhi")][it];
    d.jetsAK8.GenPartonPt[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonPt")][it];
    d.jetsAK8.GenPartonE[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonE")][it];
    d.jetsAK8.GenPartonCharge[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenPartonCharge")][it];
    d.jetsAK8.PartonFlavour[it] = vfloats_values[makeName(jetsAK8_label,pref,"PartonFlavour")][it];
    d.jetsAK8.HadronFlavour[it] = vfloats_values[makeName(jetsAK8_label,pref,"HadronFlavour")][it];
    d.jetsAK8.GenJetY[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetY")][it];
    d.jetsAK8.GenJetEta[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetEta")][it];
    d.jetsAK8.GenJetPhi[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetPhi")][it];
    d.jetsAK8.GenJetPt[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetPt")][it];
    d.jetsAK8.GenJetE[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetE")][it];
    d.jetsAK8.GenJetCharge[it] = vfloats_values[makeName(jetsAK8_label,pref,"GenJetCharge")][it];
    d.jetsAK8.HLTjetEta[it] = vfloats_values[makeName(jetsAK8_label,pref,"HLTjetEta")][it];
    d.jetsAK8.HLTjetPhi[it] = vfloats_values[makeName(jetsAK8_label,pref,"HLTjetPhi")][it];
    d.jetsAK8.HLTjetPt[it] = vfloats_values[makeName(jetsAK8_label,pref,"HLTjetPt")][it];
    d.jetsAK8.HLTjetE[it] = vfloats_values[makeName(jetsAK8_label,pref,"HLTjetE")][it];
    d.jetsAK8.HLTjetDeltaR[it] = vfloats_values[makeName(jetsAK8_label,pref,"HLTjetDeltaR")][it];
    d.jetsAK8.muonMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"muonMultiplicity")][it];
    d.jetsAK8.PhotonEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"PhotonEnergy")][it];
    d.jetsAK8.ElectronEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"ElectronEnergy")][it];
    d.jetsAK8.MuonEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"MuonEnergy")][it];
    d.jetsAK8.HFHadronEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"HFHadronEnergy")][it];
    d.jetsAK8.HFEMEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"HFEMEnergy")][it];
    d.jetsAK8.ChargedHadronMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"ChargedHadronMultiplicity")][it];
    d.jetsAK8.numberOfDaughters[it] = vfloats_values[makeName(jetsAK8_label,pref,"numberOfDaughters")][it];
    d.jetsAK8.chargedMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"chargedMultiplicity")][it];
    d.jetsAK8.neutralHadronMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"neutralHadronMultiplicity")][it];
    d.jetsAK8.neutralHadronEnergyFraction[it] = vfloats_values[makeName(jetsAK8_label,pref,"neutralHadronEnergyFraction")][it];
    d.jetsAK8.neutralEmEnergyFraction[it] = vfloats_values[makeName(jetsAK8_label,pref,"neutralEmEnergyFraction")][it];
    d.jetsAK8.chargedEmEnergyFraction[it] = vfloats_values[makeName(jetsAK8_label,pref,"chargedEmEnergyFraction")][it];
    d.jetsAK8.chargedHadronEnergyFraction[it] = vfloats_values[makeName(jetsAK8_label,pref,"chargedHadronEnergyFraction")][it];
    d.jetsAK8.photonMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"photonMultiplicity")][it];
    d.jetsAK8.electronMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"electronMultiplicity")][it];
    d.jetsAK8.HFHadronMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"HFHadronMultiplicity")][it];
    d.jetsAK8.HFEMMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"HFEMMultiplicity")][it];
    d.jetsAK8.ChargeMuEnergy[it] = vfloats_values[makeName(jetsAK8_label,pref,"ChargeMuEnergy")][it];
    d.jetsAK8.neutralMultiplicity[it] = vfloats_values[makeName(jetsAK8_label,pref,"neutralMultiplicity")][it];
    d.jetsAK8.SmearedPt[it] = vfloats_values[makeName(jetsAK8_label,pref,"SmearedPt")][it];
    d.jetsAK8.SmearedPEta[it] = vfloats_values[makeName(jetsAK8_label,pref,"SmearedPEta")][it];
    d.jetsAK8.SmearedPhi[it] = vfloats_values[makeName(jetsAK8_label,pref,"SmearedPhi")][it];
    d.jetsAK8.SmearedE[it] = vfloats_values[makeName(jetsAK8_label,pref,"SmearedE")][it];
    d.jetsAK8.JERup[it] = vfloats_values[makeName(jetsAK8_label,pref,"JERup")][it];
    d.jetsAK8.JERdown[it] = vfloats_values[makeName(jetsAK8_label,pref,"JERdown")][it];
    d.jetsAK8.subjetIndex0[it] = vfloats_values[makeName(jetsAK8_label,pref,"subjetIndex0")][it];
    d.jetsAK8.subjetIndex1[it] = vfloats_values[makeName(jetsAK8_label,pref,"subjetIndex1")][it];
    d.jetsAK8.tau1[it] = vfloats_values[makeName(jetsAK8_label,pref,"tau1")][it];
    d.jetsAK8.tau2[it] = vfloats_values[makeName(jetsAK8_label,pref,"tau2")][it];
    d.jetsAK8.tau3[it] = vfloats_values[makeName(jetsAK8_label,pref,"tau3")][it];
    d.jetsAK8.trimmedMass[it] = vfloats_values[makeName(jetsAK8_label,pref,"trimmedMass")][it];
    d.jetsAK8.prunedMass[it] = vfloats_values[makeName(jetsAK8_label,pref,"prunedMass")][it];
    d.jetsAK8.filteredMass[it] = vfloats_values[makeName(jetsAK8_label,pref,"filteredMass")][it];
  }
  d.CalculateAllVariables();
  // Razor variables
  float_values[makeName(jetsAK4_label,obj_to_pref[jetsAK4_label],"MR")] = d.jetsAK4.MR;
  float_values[makeName(jetsAK4_label,obj_to_pref[jetsAK4_label],"MTR")] = d.jetsAK4.MTR;
  float_values[makeName(jetsAK4_label,obj_to_pref[jetsAK4_label],"R")] = d.jetsAK4.R;
  float_values[makeName(jetsAK4_label,obj_to_pref[jetsAK4_label],"R2")] = d.jetsAK4.R2;
  float_values[makeName(jetsAK8_label,obj_to_pref[jetsAK8_label],"MR")] = d.jetsAK8.MR;
  float_values[makeName(jetsAK8_label,obj_to_pref[jetsAK8_label],"MTR")] = d.jetsAK8.MTR;
  float_values[makeName(jetsAK8_label,obj_to_pref[jetsAK8_label],"R")] = d.jetsAK8.R;
  float_values[makeName(jetsAK8_label,obj_to_pref[jetsAK8_label],"R2")] = d.jetsAK8.R2;
  // Top tagging variables
  float_values["nhadtops"] = d.evt.nhadtops;    
  float_values["nleptops"] = d.evt.nleptops;    
  float_values["ntops"] = d.evt.ntops;	     
  float_values["tt_dR"] = d.evt.tt_dR;	     
  float_values["tt_dPhi"] = d.evt.tt_dPhi;     
  float_values["tt_dEta"] = d.evt.tt_dEta;     
  float_values["tt_Mass"] = d.evt.tt_Mass;     
  float_values["tt_Pz"] = d.evt.tt_Pz;	     
  float_values["tt_Hz"] = d.evt.tt_Hz;	     
  float_values["tt_dPz"] = d.evt.tt_dPz;	     
  float_values["tt_extra"] = d.evt.tt_extra;    
  float_values["dHt"] = d.evt.dHt;	     
  float_values["dphi1"] = d.evt.dphi1;	     
  float_values["dphi2"] = d.evt.dphi2;	     
  float_values["dPhi_met_t1"] = d.evt.dPhi_met_t1; 
  float_values["dPhi_met_t2"] = d.evt.dPhi_met_t2; 
  float_values["tt_MR"] = d.evt.tt_MR;	     
  float_values["tt_MTR"] = d.evt.tt_MTR;	     
  float_values["tt_R"] = d.evt.tt_R;	     
  float_values["tt_R2"] = d.evt.tt_R2;	     
  float_values["HT"] = d.evt.HT;	     
  float_values["HTall"] = d.evt.HTall;	     
  float_values["HTtt"] = d.evt.HTtt;	     
  float_values["HTlep"] = d.evt.HTlep;	     
  float_values["HTex"] = d.evt.HTex;	     
  float_values["HTttFraction"] = d.evt.HTttFraction;
  float_values["HTexFraction"] = d.evt.HTexFraction;
  // JK End

  //Part 3: filling the additional variables
  treesBase->Fill(); 
  //Reset event weights/#objects
  vector<string> extravars = additionalVariables("Event");
  for(size_t addv = 0; addv < extravars.size();++addv){
    float_values[extravars.at(addv)]=0;
  }
}
string DMAnalysisTreeMaker::makeName(string label,string pref,string var){
  return pref+var+"_"+label;
}

vector<string> DMAnalysisTreeMaker::additionalVariables(string object){
  vector<string> addvar;
  bool ismuon=object.find("muon")!=std::string::npos;
  bool iselectron=object.find("electron")!=std::string::npos;
  bool ismet=object.find("met")!=std::string::npos;
  bool isAK4jet=object.find("jetAK4")==0;
  bool isAK8jet=object.find("jetsAK8")==0;
  //bool isCmsTopTagjet=object.find("jetsCmsTopTag")==0;
  //bool isAK8subjet=object.find("subjetsAK8")==0;
  //bool isCmsTopTagjet=object.find("subjetsCmsTopTag")==0;
  bool isevent=object.find("Event")!=std::string::npos;
  if(ismuon || iselectron){
    //addvar.push_back("SFTrigger");
    //addvar.push_back("SFReco");
    //addvar.push_back("isQCD");
    //addvar.push_back("isTightOffline");
    //addvar.push_back("isLooseOffline");
  }
  if(ismet){
    //addvar.push_back("CorrPt");
    //addvar.push_back("CorrPhi");
  }
  if(isAK4jet){
    addvar.push_back("MR");
    addvar.push_back("MTR");
    addvar.push_back("R");
    addvar.push_back("R2");
    
    //addvar.push_back("CorrPt");
    //addvar.push_back("CorrEta");
    //addvar.push_back("CorrPhi");
    //addvar.push_back("CorrE");
    //addvar.push_back("CorrMass");
    //addvar.push_back("CorrNJets");
    //addvar.push_back("CorrPartonFlavour");
  }
  if(isAK8jet){
    addvar.push_back("MR");
    addvar.push_back("MTR");
    addvar.push_back("R");
    addvar.push_back("R2");
    
    //addvar.push_back("CorrPt");
    //addvar.push_back("CorrEta");
    //addvar.push_back("CorrPhi");
    //addvar.push_back("CorrE");
    //addvar.push_back("CorrMass");
    //addvar.push_back("CorrNJets");
    //addvar.push_back("CorrPartonFlavour");
  }
  if(isevent){
    addvar.push_back("nhadtops");
    addvar.push_back("nleptops");
    addvar.push_back("ntops");
    addvar.push_back("tt_dR");
    addvar.push_back("tt_dPhi");
    addvar.push_back("tt_dEta");
    addvar.push_back("tt_Mass");
    addvar.push_back("tt_Pz");
    addvar.push_back("tt_Hz");
    addvar.push_back("tt_dPz");
    addvar.push_back("tt_extra");
    addvar.push_back("dHt");
    addvar.push_back("dphi1");
    addvar.push_back("dphi2");
    addvar.push_back("dPhi_met_t1");
    addvar.push_back("dPhi_met_t2");
    addvar.push_back("tt_MR");
    addvar.push_back("tt_MTR");
    addvar.push_back("tt_R");
    addvar.push_back("tt_R2");
    addvar.push_back("HT");
    addvar.push_back("HTall");
    addvar.push_back("HTtt");
    addvar.push_back("HTlep");
    addvar.push_back("HTex");
    addvar.push_back("HTttFraction");
    addvar.push_back("HTexFraction");
    
    // DMTT
    addvar.push_back("nTightMuons");
    addvar.push_back("nLooseMuons");
    addvar.push_back("nTightElectrons");
    addvar.push_back("nLooseElectrons");
    addvar.push_back("nTightJets");
    addvar.push_back("nLooseJets");
    
    //addvar.push_back("nElectronsSF");
    //addvar.push_back("nMuonsSF");
    //addvar.push_back("nCSVTJets");
    //addvar.push_back("nCSVMJets");
    //addvar.push_back("nCSVLJets");
    //addvar.push_back("bWeight1TCSVT");
    //addvar.push_back("bWeight1TCSVM");
    //addvar.push_back("bWeight1TCSVL");
    //addvar.push_back("bWeight2TCSVT");
    //addvar.push_back("bWeight2TCSVM");
    //addvar.push_back("bWeight2TCSVL");
    //addvar.push_back("LHEWeightSign");
  }
  return addvar;
}
//DMAnalysisTreeMaker::~DMAnalysisTreeMaker(const edm::ParameterSet& iConfig)
// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(DMAnalysisTreeMaker);
