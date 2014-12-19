#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"

int main() {
  const char* Input = "b2gTree_T2tt_2J_mStop-650_mLSP-325.root/DMTreesDumper/TreeBase";
  const char* Output = "plots.root";
  bool Run = 1;
  
  // Data variable
  Data d;
  
  // Histogram storage class
  SmartHistos sh;
  sh.AddHistoType("evt");
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jet");
  sh.AddHistoType("jetAK8");
  sh.AddHistoType("met");
  
  // Define Postfixes here:
  sh.AddNewPostfix("AK4JetsPtOrdered", &d.jet.it, "Jet[1to10]", "Jet [1to10]", "1-10");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh.AddNewFillParam("MuonEnergy",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.mu.E;        }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuonPt",            { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.mu.Pt;       }, .axis_title="Muon p_{T} (GeV/c)"});
  sh.AddNewFillParam("EleEnergy",         { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.ele.E;       }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.ele.Pt;      }, .axis_title="Electron p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jet.E;       }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jet.Pt;      }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jet.Mass;    }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetEnergy",      { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.E;    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("AK8JetPt",          { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.Pt;   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMass",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.Mass; }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("MetPt",             { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.met.Pt;      }, .axis_title="MET p_{T} (GeV/c)"});
  
  // Define Cuts here:
  sh.AddNewCut("AK4Highest2Jet", [&d](){ return d.jet.jets_size>=2 && d.jet.it<2; });
  sh.AddNewCut("AK8Highest2Jet", [&d](){ return d.jetAK8.jetsAK8_size>=2 && d.jetAK8.it<2; });
  
  // Set Histogram weight (empty = 1)
  sh.SetHistoWeights({});
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  sh.AddHistos("mu",     { .fill="MuonEnergy",    .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("mu",     { .fill="MuonPt",        .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("ele",    { .fill="EleEnergy",     .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="ElePt",         .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("jet",    { .fill="AK4JetEnergy",  .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jet",    { .fill="AK4JetPt",      .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,0, 0,1000} });
  sh.AddHistos("jet",    { .fill="AK4JetMass",    .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,250, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetEnergy",  .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetPt",      .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetMass",    .pfs={}, .cuts={}, .draw="", .ranges={0,250, 0,0} });
  
  sh.AddHistos("met",    { .fill="MetPt",         .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  TFile *file;
  if (Run) {
    // Initialize TreeReader
    B2GTreeReader reader;
    
    // Loop on files and read the Trees
    B2GTreeLooper looper(1,1);
    looper.AddFiles(Input);
    while (looper.LoopOnFiles()) {
      TFile *curr_file = looper.CurrentFile();
      reader.Load_Tree(*curr_file,looper.TreeName());
      while(looper.LoopOnEntries()) {
        reader.GetEntry(looper.CurrentEntry());
	d = reader.data;
	
	sh.Fill("evt");
	// loop on objects and fill their histos
	while(d.mu.Loop())     sh.Fill("mu");
	while(d.ele.Loop())    sh.Fill("ele");
	while(d.jet.Loop()) {
	  //std::cout<<d.jet.it<<" "<<d.jet.Pt<<std::endl;
	  sh.Fill("jet");
	}
	while(d.jetAK8.Loop()) sh.Fill("jetAK8");
	while(d.met.Loop())    sh.Fill("met");
      }
      curr_file->Close();
    }
    
    std::cout<<std::endl;
    std::cout<<"Finished ..."<<std::endl;
    std::cout<<"Writing Histograms to File: "<<Output<<std::endl;
    file = new TFile(Output,"recreate");
    sh.DrawPlots();
    sh.Write();
    file->Close();
  } else {
    std::cout<<"Loading Histos from file: "<<Output<<std::endl;
    sh.Load(Output);
  }
  
  return 0;
}

