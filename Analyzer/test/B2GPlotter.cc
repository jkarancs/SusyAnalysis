#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"
#include "../interface/Razor.h"

#define NTHFILE 1

int main(int argc, char* argv[]) {
  // Get arguments from shell
  std::vector<std::string> filelist;
  std::string outputfile="plots.root";
  // -o <output file> option:
  // Specify the output root filename
  // Rest of the arguments are treated as files added
  // Note:
  // If using postfixes with the v.pf_fila_add variable
  // each added file will increase this variable so when using *
  // add ""-s so instead of the shell TChain will parse the argument
  bool is_o = false;
  for(int i=1; i<argc; i++) {
    std::string arg = argv[i];
    if (arg[0]=='-'&&arg.size()==2) { is_o = (arg[1]=='o'); }
    else if (is_o) { outputfile=arg; is_o=0; }
    else filelist.push_back(arg);
  }
  
  bool Run = 1;
  
  // Data variable
  Data d;
  
  // Histogram storage class
  SmartHistos sh;
  sh.AddHistoType("evt");
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jetAK4");
  sh.AddHistoType("jetAK8");
  sh.AddHistoType("met");
  
  // Define Postfixes here:
  sh.AddNewPostfix("AK4JetsPtOrdered", &d.jetAK4.it, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("AK8JetsPtOrdered", &d.jetAK8.it, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh.AddNewFillParam("MuonEnergy",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.mu.E;              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuonPt",            { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.mu.Pt;             }, .axis_title="Muon p_{T} (GeV/c)"});
  sh.AddNewFillParam("EleEnergy",         { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.ele.E;             }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.ele.Pt;            }, .axis_title="Electron p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK4.E;          

}, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK4.Pt;                 }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK4.Mass;               }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetEnergy",      { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.E;                  }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("AK8JetPt",          { .nbin= 200, .low= 0,   .high=1000, .fill=[&d](){ return d.jetAK8.Pt;                 }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMass",        { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.Mass;               }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetPrunedMass",  { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.jetAK8.prunedMass;         }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetTau1",        { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau1;               }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau2",        { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau2;               }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("AK8JetTau3",        { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3;               }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("AK8JetTau21",       { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau2/d.jetAK8.tau1; }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau31",       { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3/d.jetAK8.tau1; }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau32",       { .nbin= 100, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3/d.jetAK8.tau2; }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("AK8JetMR",          { .nbin=  50, .low= 0,   .high=5000, .fill=[&d](){ return d.jetAK8.MR;                 }, .axis_title="M_{R} (GeV)"});
  sh.AddNewFillParam("AK8JetMTR",         { .nbin=  50, .low= 0,   .high=5000, .fill=[&d](){ return d.jetAK8.MTR;                }, .axis_title="M_{T}^{R} (GeV)"});
  sh.AddNewFillParam("AK8JetR",           { .nbin=  50, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.R;                  }, .axis_title="R"});
  sh.AddNewFillParam("AK8JetR2",          { .nbin=  50, .low= 0,   .high=   1, .fill=[&d](){ return d.jetAK8.R2;                 }, .axis_title="R^2"});
  sh.AddNewFillParam("MetPt",             { .nbin= 100, .low= 0,   .high= 500, .fill=[&d](){ return d.met.Pt;                    }, .axis_title="MET p_{T} (GeV/c)"});
  
  // Define Cuts here:
  sh.AddNewCut("AK4Highest2Jet", [&d](){ return d.jetAK4.jets_size>=2 && d.jetAK4.it<2; });
  sh.AddNewCut("AK8Highest2Jet", [&d](){ return d.jetAK8.jetsAK8_size>=2 && d.jetAK8.it<2; });
  sh.AddNewCut("AK4Highest3Jet", [&d](){ return d.jetAK4.jets_size>=3 && d.jetAK4.it<3; });
  sh.AddNewCut("AK8Highest3Jet", [&d](){ return d.jetAK8.jetsAK8_size>=3 && d.jetAK8.it<3; });
  
  // Set Histogram weight (empty = 1)
  sh.SetHistoWeights({});
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  sh.AddHistos("mu",     { .fill="MuonEnergy",       .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("mu",     { .fill="MuonPt",           .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("ele",    { .fill="EleEnergy",        .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="ElePt",            .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("jetAK4", { .fill="AK4JetEnergy",     .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK4", { .fill="AK4JetPt",         .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,0, 0,1000} });
  sh.AddHistos("jetAK4", { .fill="AK4JetMass",       .pfs={"AK4JetsPtOrdered"}, .cuts={"AK4Highest2Jet"}, .draw="", .ranges={0,250, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetEnergy",     .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetPt",         .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,300} });
  sh.AddHistos("jetAK8", { .fill="AK8JetMass",       .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,250, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetPrunedMass", .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau1",       .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,100} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau2",       .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,150} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau3",       .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,200} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau21",      .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest3Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau31",      .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest3Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau32",      .pfs={"AK8JetsPtOrdered"}, .cuts={"AK8Highest3Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetMR",         .pfs={}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetMTR",        .pfs={}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetR",          .pfs={}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetR2",         .pfs={}, .cuts={"AK8Highest2Jet"}, .draw="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("met",    { .fill="MetPt",         .pfs={}, .cuts={}, .draw="", .ranges={0,0, 0,0} });
  
  std::cout<<"-----------------------------------------------------------------\n";
  std::cout<<"Creating the following plots:\n"; sh.PrintNames();
  std::cout<<"-----------------------------------------------------------------\n";
  
  TFile *file;
  if (Run) {
    // Initialize TreeReader
    B2GTreeReader reader;
    
    // Class to Loop on files and read the Trees
    B2GTreeLooper looper(NTHFILE,1);
    if (filelist.size()) {
      std::cout<<"Adding "<<filelist.size()<<" files from the shell arguments.\n";
      for (size_t i=0; i<filelist.size(); ++i) looper.AddFile(filelist[i]);
    } else {
      looper.AddFile("b2gTree_T2tt_2J_mStop-650_mLSP-325.root");
    }
    
    while (looper.LoopOnSamples()) {
      while (looper.LoopOnFiles()) {
        TFile *curr_file = looper.CurrentFile();
        reader.Load_Tree(*curr_file,looper.TreeName());
        while(looper.LoopOnEntries()) {
          reader.GetEntry(looper.CurrentEntry());
          d = reader.data;
	  calcRazorAK4(d);
	  calcRazorAK8(d);
	  
          sh.Fill("evt");
          // loop on objects and fill their histos
          while(d.mu.Loop())     sh.Fill("mu");
          while(d.ele.Loop())    sh.Fill("ele");
          while(d.jetAK4.Loop()) {
            //std::cout<<d.jetAK4.it<<" "<<d.jetAK4.Pt<<std::endl;
            sh.Fill("jetAK4");
          }
          while(d.jetAK8.Loop()) sh.Fill("jetAK8");
          while(d.met.Loop())    sh.Fill("met");
        }
        curr_file->Close();
      }
    }
    
    std::cout<<std::endl;
    std::cout<<"Finished ..."<<std::endl;
    std::cout<<"Writing Histograms to File: "<<outputfile<<std::endl;
    file = new TFile(outputfile.c_str(),"recreate");
    sh.DrawPlots();
    sh.Write();
    file->Close();
  } else {
    std::cout<<"Loading Histos from file: "<<outputfile<<std::endl;
    sh.Load(outputfile.c_str());
  }
  
  return 0;
}

