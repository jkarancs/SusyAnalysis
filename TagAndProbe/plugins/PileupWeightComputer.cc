/*
  \class        PileupWeightComputer
  \author       Janos Karancsi
  \description: Add In Time PU info and Event weight from an input root file
                Containing Data and MC Pileup Distributions
		
		Reweighting (Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData):
		- First Check out the LumiDB package:
                  Official Twiki page you can check for latest info
                  https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc

		cmsenv
		git cms-addpkg RecoLuminosity/LumiDB
		scram b

		- Get the Good run/ls JSON file corresponding to the processed dataset
		- Get the latest official 2012 pileup JSON file                              
		+ "--minBiasXsec" option to set it to the approved value of 68000 (for 2011) or 69400 (for 2012)
		- Get Pileup distribution from MC
		- Use output root file 
		
		Quick test:
		- Generate Pileup Histo for data using pileupCalc.py script
		Use the "--minBiasXsec" option to set the MinBias cross section to
		the approved value of 68000 (for 2011) or 69400 (for 2012))
		You will need the following package:
		git cms-addpkg RecoLuminosity/LumiDB
		
		dbs search --query "find run, lumi where file=/store/Run2012C/MinimumBias/RECO/22Jan2013-v1/20006/04E396DF-9172-E211-B250-003048FEAF90.root"
		echo '{"202504": [[127, 128], [249, 253]]}' > test_input_JSON.txt
		cp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt .
		RecoLuminosity/LumiDB/scripts/pileupCalc.py -i test_input_JSON.txt --inputLumiJSON pileup_latest.txt \
		--calcMode true --minBiasXsec 69400 --maxPileupBin 50 --numPileupBins 50 \
		PileupHistogram_test.root
		
		- Get MC Pileup distribution from an input file and write it to previous file
		
		root -l MC_input_file.root
		TH1D* mcpileup = new TH1D("mcpileup","MC Pileup", 50, 0, 50);
		Events->Draw("(TrueNumInteractions_)>>mcpileup");
		TFile *f= new TFile("PileupHistogram_test.root","update");
		mcpileup->Write();
		
		- Use this root file as input for the PileupWeightComputer
		process.pileupWeightComputer = cms.EDProducer("PileupWeightComputer",
		    isMC = cms.bool(True),
		    dataPileupFile = cms.string("PileupHistogram_test.root"),
		    mcPileupFile   = cms.string("PileupHistogram_test.root"),
		    dataPileupHistoName = cms.string("pileup"),
		    mcPileupHistoName = cms.string("mcpileup")
		)

		- dataInstlumiInputFile
		  How to download and create a text file with
                  the offical instlumi/pileup values for 2012 Data ordered by Run, Ls
		
		  Run lumiCalc2.py script and reformat output to:
		  [Run] [LumiSection] [InstLumi (delivered - nb-1/LS)] [Pileup]
		
		RecoLuminosity/LumiDB/scripts/lumiCalc2.py --begin "01/01/10 00:00:00" --end "05/09/13 12:00:00" --nowarning -b stable overview | grep -v WARNING | grep -v n/a | tail -n+13 | head -n-6 | sed "s;:; ;;s;(/mb);0.001;;s;(/ub);1;;s;(/nb);1000;;s;(/pb);1000000;" | awk '{ printf "%d %f\n", $2, $7*$8 }' > ! run_ls_instlumi_pileup_2012.txt
		
		

*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class PileupWeightComputer : public edm::EDProducer {
public:
  explicit PileupWeightComputer(const edm::ParameterSet & iConfig);
  virtual ~PileupWeightComputer() ;
  
  virtual void beginJob();
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  edm::InputTag probes_;

  bool isMC_;

  bool calcWeights_;
  
  std::string mcPileupFile_;
  std::string mcPileupHistoName_;
  std::string dataPileupFile_;
  std::string dataPileupHistoName_;

  double mcLumiScale_;
  
  edm::LumiReWeighting LumiWeights_;

  std::string dataPileupInputFile_;
  std::map<unsigned long int, float> runls_pileup_;
  std::map<unsigned long int, float> runls_instlumi_;
  
};

PileupWeightComputer::PileupWeightComputer(const edm::ParameterSet & iConfig) :
  probes_(iConfig.getParameter<edm::InputTag>("probes")),
  isMC_(iConfig.getParameter<bool>("isMC")) {

  calcWeights_ = false;

  produces<edm::ValueMap<float> >("pileup");
  produces<edm::ValueMap<float> >("instlumi");
  produces<edm::ValueMap<float> >("weight");

  if (isMC_) {
    if (iConfig.exists("mcPileupFile")&&iConfig.exists("mcPileupHistoName")&&
	iConfig.exists("dataPileupFile")&&iConfig.exists("dataPileupHistoName")) {
      calcWeights_ = true;
      mcPileupFile_=iConfig.getParameter<std::string>("mcPileupFile");
      mcPileupHistoName_=iConfig.getParameter<std::string>("mcPileupHistoName");
      dataPileupFile_=iConfig.getParameter<std::string>("dataPileupFile");
      dataPileupHistoName_=iConfig.getParameter<std::string>("dataPileupHistoName");
    }
    if (iConfig.exists("mcLumiScale"))
      mcLumiScale_=iConfig.getParameter<double>("mcLumiScale");
    else mcLumiScale_ = 221.95;
  } else {
    if (iConfig.exists("dataPileupInputFile"))
      dataPileupInputFile_=iConfig.getParameter<std::string>("dataPileupInputFile");
    else {
      std::cout<<"PileupWeightComputer WARNING!: dataPileupInputFile is not specified."<<std::endl;
    }
  }
}

void PileupWeightComputer::beginJob() {
  // Initialize LumiReWeighting
  if (calcWeights_) LumiWeights_ = edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupHistoName_, dataPileupHistoName_);

  if (!isMC_) {
    // Read InstLumi (delivered) and Pileup from text file generated with lumiCalc2.py
    int run = 0;
    int ls = 0;
    unsigned long int runls = 0;
    float instlumi = 0.0;
    float pileup = 0.0;
    FILE * input = fopen (dataPileupInputFile_.c_str(),"r");
    if (input) std::cout<<"PileupWeightComputer: Reading data pileup info from file: "<<dataPileupInputFile_<<std::endl;
    else std::cout<<"PileupWeightComputer ERROR!: dataPileupInputFile ="<<dataPileupInputFile_<<" is not found."<<std::endl;

    
    // LumiCalc integrates lumi in the whole LumiSection
    float ls_length = 23.3104;

    int a = (input!=0);
    while (a==1) {
      a = fscanf (input, "%d", &run);
      a = fscanf (input, "%d", &ls);
      a = fscanf (input, "%f", &instlumi);
      a = fscanf (input, "%f", &pileup);
      if (a==1) {
        runls = run * 100000 + ls;
	runls_instlumi_[runls] = instlumi * 1000 / ls_length; // Convert nb-1/LS to ub-1-s
        runls_pileup_[runls] = pileup;
      }
    }
    if (input) fclose (input);
  }

}

PileupWeightComputer::~PileupWeightComputer() {}

void PileupWeightComputer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  float pileup = -9999.0;
  float instlumi = -9999.0;
  float weight = -9999.0;

  if (isMC_) {
    edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
    iEvent.getByLabel("addPileupInfo", puInfo);
    
    if (puInfo.isValid()) {
      // look for the intime PileupSummaryInfo
      std::vector<PileupSummaryInfo>::const_iterator pu;
      std::vector<PileupSummaryInfo>::const_iterator pu0=puInfo->end();
      
      for(pu=puInfo->begin(); pu!=puInfo->end(); ++pu) if (pu->getBunchCrossing()==0) pu0=pu;
      
      if(pu0!=puInfo->end()) {
	pileup = pu0->getTrueNumInteractions();
	instlumi = pileup * mcLumiScale_; // Units: ub-1-s
	if (calcWeights_) {
	  weight = LumiWeights_.weight(pileup);
	}
      } else {
	std::cout<<"** ERROR (PileupWeightComputer::produce): Cannot find the in-time pileup info\n";
      }

    } else std::cout<<"** WARNING (PileupWeightComputer::produce): PileupInfo invalid\n";
  } else {
    
    unsigned long int runls = iEvent.id().run()*100000 + iEvent.luminosityBlock();
    if (runls_pileup_.count(runls)) {
      pileup = runls_pileup_[runls];
      instlumi = runls_instlumi_[runls];
    }
    
    weight = 1;
  }
  //std::cout<<"Pileup: "<<pileup<<" Instlumi: "<<instlumi<<" Weight: "<<weight<<std::endl;

  // read input
  edm::Handle<edm::View<reco::Candidate> > probes;
  iEvent.getByLabel(probes_,  probes);
  
  // prepare and fill vector for output
  std::vector<float> v_pileup(probes->size(),pileup);
  std::vector<float> v_instlumi(probes->size(),instlumi);
  std::vector<float> v_weight(probes->size(),weight);
  
  // convert into edm::ValueMap and store
  std::auto_ptr<edm::ValueMap<float> > vm_pileup(new edm::ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > vm_instlumi(new edm::ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > vm_weight(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_pileup(*vm_pileup);
  edm::ValueMap<float>::Filler filler_instlumi(*vm_instlumi);
  edm::ValueMap<float>::Filler filler_weight(*vm_weight);
  filler_pileup.insert(probes, v_pileup.begin(), v_pileup.end());
  filler_instlumi.insert(probes, v_instlumi.begin(), v_instlumi.end());
  filler_weight.insert(probes, v_weight.begin(), v_weight.end());
  filler_pileup.fill();
  filler_instlumi.fill();
  filler_weight.fill();
  iEvent.put(vm_pileup,"pileup");
  iEvent.put(vm_instlumi,"instlumi");
  iEvent.put(vm_weight,"weight");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PileupWeightComputer);
