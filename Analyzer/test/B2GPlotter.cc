#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/Samples.h"
#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"

#define TEST 0
#if TEST == 1
#define NTHSTAT 40
#else
#define NTHSTAT 1
#endif

#define IntLumi_invfb 10.0

int main(int argc, char* argv[]) {
  // Get arguments from shell
  std::vector<std::string> filelist_fromshell;
  std::string inputfile="";
  std::string outputfile="plots.root";
  // -o <output file> option:
  // Specify the output root filename
  // Rest of the arguments are treated as files added
  // Note:
  // If using postfixes with the v.pf_fila_add variable
  // each added file will increase this variable so when using *
  // add ""-s so instead of the shell TChain will parse the argument
  bool is_i = false;
  bool is_o = false;
  Samples samples;
  
  for(int i=1; i<argc; i++) {
    std::string arg = argv[i];
    if (arg[0]=='-'&&arg.size()==2) { 
      is_i = (arg[1]=='i'); 
      is_o = (arg[1]=='o'); 
    } else if (is_i) { 
      inputfile=arg; is_i=0; 
    } else if (is_o) {
      outputfile=arg; is_o=0; 
    } else {
      filelist_fromshell.push_back(arg);
    }
  }
  bool Run = inputfile.size()==0;
  if (filelist_fromshell.size()) {
    samples.AddSample("test", "test", "1", { { .dir="", .xsec_pb=1/(IntLumi_invfb*1000) } });
  } else {
    bool test = 0;
    //if (test) samples.AddSample("test", "T5ttttDeg (#tilde{g} #rightarrow t + #tilde{t}_{4-body decay} )", { { .dir="../../../B2GTTreeNtupleExtra_susy.root", .xsec_pb=0.0460525 } });
    if (test) samples.AddSample("test", "T5ttttDeg (#tilde{g}#rightarrowt(#tilde{t}#rightarrow#tilde{#chi}^{0}_{1}b l#nu/q#bar{q}) )", "1", { { .dir="../../../B2GTTreeNtupleExtra_susy.root", .xsec_pb=0.0460525 } });
    else {
      //"T5tttt (#tilde{t} 2/3body decay - M_{#tilde{g}}=1.3TeV, M_{#tilde{t}}=, M_{#tilde{#chi^{#pm}_{1}}}, M_{#tilde{#chi^{0}_{1}}})"
      //"T5ttttDeg (#tilde{g} #rightarrow t + (#tilde{t} #rightarrow b + #tilde{#chi^{0}_{1}} + l#nu/q#bar{q} ) )"
      //"T5ttttDeg (#tilde{g} #rightarrow t + (#tilde{t} #rightarrow b + ( #tilde{#chi^{#pm}_{1}} #rightarrow #tilde{#chi^{0}_{1}} + l#nu/q#bar{q} ) ) )"
      //#tilde{g}#rightarrowt(#tilde{t}#rightarrowb#tilde{#chi}^{0}_{1} l#nu/q#bar{q}) 
      //#tilde{g}#rightarrowt(#tilde{t}#rightarrowb(#tilde{#chi^{#pm}_{1}}#rightarrow#tilde{#chi}^{0}_{1} l#nu/q#bar{q}))
      
      //std::string sample_dir = "/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb18_edm_Feb13/";
      std::string sample_dir = "/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb28_edm_Feb20/";
      samples.AddSample("T5ttttDeg_4bodydec_mGo1300", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.3TeV)", "1",
                        { { .dir=sample_dir+"T5ttttDeg_mGo1300_4bodydec/*.root", .xsec_pb=0.0460525 } });
      samples.AddSample("T5ttttDeg_3bodydec_mGo1300", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.3TeV)", "12",
                        { { .dir=sample_dir+"T5ttttDeg_mGo1300_23bodydec/*.root", .xsec_pb=0.0460525 } });
      samples.AddSample("T5ttttDeg_4bodydec_mGo1000", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.0TeV)", "14",
                        { { .dir=sample_dir+"T5ttttDeg_mGo1000_4bodydec/*.root", .xsec_pb=0.325388 } });
      samples.AddSample("T5ttttDeg_3bodydec_mGo1000", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.0TeV)", "16",
                        { { .dir=sample_dir+"T5ttttDeg_mGo1000_23bodydec/*.root", .xsec_pb=0.325388 } });
      samples.AddSample("QCD_HT", "QCD (HT bins)", "4",
        		{ { .dir=sample_dir+"QCD_HT-100To250/*.root",  .xsec_pb=28730000 },
                          { .dir=sample_dir+"QCD_HT_250To500/*.root",  .xsec_pb=670500 },
            	          { .dir=sample_dir+"QCD_HT-500To1000/*.root", .xsec_pb=26740 },
            	          { .dir=sample_dir+"QCD_HT_1000ToInf/*.root", .xsec_pb=769.7 } });
      //samples.AddSample("QCD_Pt_bcToE", "QCD (Pt bins, b/c#rightarrowe)", "38",
      //                  { { .dir=sample_dir+"QCD_Pt_20to30_bcToE/*.root",   .xsec_pb=675900000 },
      //                    { .dir=sample_dir+"QCD_Pt_30to80_bcToE/*.root",   .xsec_pb=185900000 },
      //                    { .dir=sample_dir+"QCD_Pt_80to170_bcToE/*.root",  .xsec_pb=3495000 },
      //                    { .dir=sample_dir+"QCD_Pt_170toInf_bcToE/*.root", .xsec_pb=128500 } });
      samples.AddSample("TTBar", "t#bar{t}", "2", { { .dir=sample_dir+"TT/*.root", .xsec_pb=806.1 } });
      samples.AddSample("GJets_HT", "G+Jets (HT bins)", "804",
                        { { .dir=sample_dir+"GJets_HT-100to200/*.root",    .xsec_pb=1534 },
                          { .dir=sample_dir+"GJets_HT-200to400/*.root",    .xsec_pb=489.9 },
                          { .dir=sample_dir+"GJets_HT-400to600/*.root",    .xsec_pb=62.05 },
                          { .dir=sample_dir+"GJets_HT-600toInf/*.root",    .xsec_pb=20.87 } });
      //samples.AddSample("WJets", "W+Jets #rightarrow l+#nu", "3", { { .dir=sample_dir+"WJetsToLNu/*.root", .xsec_pb=61526.7 } });
      samples.AddSample("WJets_HT", "W+Jets #rightarrow l+#nu (HT bins)", "3",
        		{ { .dir=sample_dir+"WJetsToLNu_HT-100to200/*.root", .xsec_pb=2234.9 },
        		  { .dir=sample_dir+"WJetsToLNu_HT-200to400/*.root", .xsec_pb=580.1 }, 
        		  { .dir=sample_dir+"WJetsToLNu_HT-400to600/*.root", .xsec_pb=68.4 },
            	          { .dir=sample_dir+"WJetsToLNu_HT-600toInf/*.root", .xsec_pb=23.14 } });
      samples.AddSample("T_tW", "single t/#bar{t} (tW channel)", "617",
                        { { .dir=sample_dir+"T_tW-channel/*.root", .xsec_pb=35 },
                          { .dir=sample_dir+"Tbar_tW-channel/*.root", .xsec_pb=35 } });
      //samples.AddSample("TToLep_s_t", "single t/#bar{t}#rightarrowl (s,t channel)", "882",
      //                  { { .dir=sample_dir+"TToLeptons_s-channel/*.root",    .xsec_pb=2 },
      //                    { .dir=sample_dir+"TBarToLeptons_s-channel/*.root", .xsec_pb=1 },
      //                    { .dir=sample_dir+"TToLeptons_t-channel/*.root",    .xsec_pb=103.4 },
      //                    { .dir=sample_dir+"TBarToLeptons_t-channel/*.root", .xsec_pb=61.6 } });
      samples.AddSample("ZJets_HT", "Z+Jets #rightarrow #nu#nu (HT bins)", "401",
        		{ { .dir=sample_dir+"ZJetsToNuNu_HT-100to200/*.root", .xsec_pb=372.6 },
        		  { .dir=sample_dir+"ZJetsToNuNu_HT-200to400/*.root", .xsec_pb=100.8 }, 
        		  { .dir=sample_dir+"ZJetsToNuNu_HT-400to600/*.root", .xsec_pb=11.99 },
            	          { .dir=sample_dir+"ZJetsToNuNu_HT-600toInf/*.root", .xsec_pb=4.113 } });
      //samples.AddSample("DYJets", "DY+Jets #rightarrow l+l", "797", { { .dir=sample_dir+"DYJetsToLL/*.root", .xsec_pb=4746 } });
      //samples.AddSample("DYJets_HT", "DY+Jets #rightarrow l+l (HT bins)", "797",
      //                  { { .dir=sample_dir+"DYJetsToLL_M-50_HT-100to200/*.root", .xsec_pb=194.3 },
      //                    { .dir=sample_dir+"DYJetsToLL_M-50_HT-200to400/*.root", .xsec_pb=52.24 },
      //                    { .dir=sample_dir+"DYJetsToLL_M-50_HT-400to600/*.root", .xsec_pb=6.546 },
      //                    { .dir=sample_dir+"DYJetsToLL_M-50_HT-600toInf/*.root", .xsec_pb=2.179 } });
      //samples.AddSample("GGJets_M", "GG+Jets (M bins)", "803",
      //  		{ { .dir=sample_dir+"GGJets_M-200To500/*.root",    .xsec_pb=2.43383 },
      //  		  { .dir=sample_dir+"GGJets_M-500To1000/*.root",   .xsec_pb=0.172872 },
      //  		  { .dir=sample_dir+"GGJets_M-1000To2000/*.root",  .xsec_pb=0.0104901 },
      //      	          { .dir=sample_dir+"GGJets_M-2000To4000/*.root",  .xsec_pb=0.000439813 },
      //      	          { .dir=sample_dir+"GGJets_M-4000To8000/*.root",  .xsec_pb=0.00000219697 },
      //      	          { .dir=sample_dir+"GGJets_M-8000To13000/*.root", .xsec_pb=7.05314e-11 } });
    }
  }
  std::vector<size_t > dir_to_index = samples.GetDirToIndex();
  std::vector<double> sample_xsec_pb = samples.GetCrossSections();
  
  // Initialize TreeReader
  B2GTreeReader reader;
  
  // Class to Loop on files and read the Trees
  B2GTreeLooper looper(NTHSTAT,1);
  
  // Data variable
  Data d;
  
  // Histogram storage class
  SmartHistos sh;
  sh.AddHistoType("evt");
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jetsAK4");
  sh.AddHistoType("jetsAK8");
  sh.AddHistoType("met");
  
  // Define Postfixes here:
  sh.AddNewPostfix("ttbar,qcd",                 [&looper](){ return looper.it_sample; }, "ttbar;qcd", "t#bar{t};QCD", "2,6");
  sh.AddNewPostfix("ttbar,qcd,Susy3,Susy4",     [&looper](){ return looper.it_sample; }, "ttbar;qcd;susy3body;susy4body", "t#bar{t};QCD;T5tttt - 3body;T5tttt - 4body", "2,6,4,3");
  sh.AddNewPostfix("AllSamples",                [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, samples.GetPFNames(), samples.GetLatexNames(), samples.GetColors());
  //const char* Samples = "ttbar,qcd,Susy3,Susy4";
  const char* Samples = "AllSamples";
  
  sh.AddNewPostfix("SideBand,Signal", [&d](){ return d.evt.HTall > 1500; }, "SideBand;Signal", "H_{T,all} < 1500 GeV/c;H_{T,all} > 1500 GeV/c", "4,2");
  sh.AddNewPostfix("R0p25",           [&d](){ return d.evt.R > 0.25; },  "RBelow0p25;RAbove0p25", "R < 0.25;R > 0.25", "4,2");
  sh.AddNewPostfix("R0p3",            [&d](){ return d.evt.R > 0.3; },   "RBelow0p3;RAbove0p3", "R < 0.3;R > 0.3", "4,2");
  sh.AddNewPostfix("R0p32",           [&d](){ return d.evt.R > 0.32; },  "RBelow0p32;RAbove0p32", "R < 0.32;R > 0.32", "4,2"); // Best cut
  sh.AddNewPostfix("R0p35",           [&d](){ return d.evt.R > 0.35; },  "RBelow0p35;RAbove0p35", "R < 0.35;R > 0.35", "4,2");
  sh.AddNewPostfix("R0p4",            [&d](){ return d.evt.R > 0.4; },   "RBelow0p4;RAbove0p4", "R < 0.4;R > 0.4", "4,2");
  sh.AddNewPostfix("DPhi2p8",         [&d](){ return fabs(d.evt.TTHadDPhi) > 2.8; },  "DPhiBelow2p8;DPhiAbove2p8", "#Delta#phi_{t#bar{t}} < 2.8;#Delta#phi_{t#bar{t}} > 2.8", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1450",       [&d](){ return d.evt.HTall > 1450; },  "HtAllBelow1450;HtAllAbove1450", "H_{T,all} < 1450;H_{T,all} > 1450", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1500",       [&d](){ return d.evt.HTall > 1500; },  "HtAllBelow1500;HtAllAbove1500", "H_{T,all} < 1500;H_{T,all} > 1500", "4,2");
  
  //sh.AddNewPostfix("ttbar,Susy3,Susy4", &looper.it_sample, "ttbar;susy3body;susy4body", "SM t#bar{t};T5tttt - 3body;T5tttt - 4body", "2,4,3");
  sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetsAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("JetsPtOrdered",    [&d](){ return d.jetsAK8.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("NSubJet",          [&d](){ return (size_t)d.jetsAK8.nSubJets[NJET]; }, "NSubJet[0to4]", "N_{subjet}=[0to4]", "1-5");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh.AddNewFillParam("MuEnergy",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.E[NLEP];              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuPt",              { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.Pt[NLEP];             }, .axis_title="Muon p_{T} (GeV/c)"});
  sh.AddNewFillParam("MuDRJet",           { .nbin=  50, .low=   0,   .high=   5, .fill=[&d](){ return d.evt.MuDRJet[d.mu.it];    }, .axis_title="#DeltaR (#mu, jet)"});
  sh.AddNewFillParam("MuRelPtJet",        { .nbin= 100, .low=   0,   .high= 500, .fill=[&d](){ return d.evt.MuRelPtJet[d.mu.it]; }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV/c)"});
  sh.AddNewFillParam("MuJetCombMass",     { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.evt.MuJetCombMass[d.mu.it];    }, .axis_title="Mass_{#mu+jet comb.} (GeV/c^{2})"});
  
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.E[NLEP];              }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.Pt[NLEP];             }, .axis_title="Electron p_{T} (GeV/c)"});
  sh.AddNewFillParam("EleDRJet",          { .nbin=  50, .low=   0,   .high=   5, .fill=[&d](){ return d.evt.EleDRJet[d.ele.it];   }, .axis_title="#DeltaR (e, jet)"});
  sh.AddNewFillParam("EleRelPtJet",       { .nbin= 100, .low=   0,   .high= 500, .fill=[&d](){ return d.evt.EleRelPtJet[d.ele.it];}, .axis_title="p_{T}^{rel} (e, jet) (GeV/c)"});
  sh.AddNewFillParam("EleJetCombMass",    { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.evt.EleJetCombMass[d.ele.it];    }, .axis_title="Mass_{e+jet comb.} (GeV/c^{2})"});
  
  sh.AddNewFillParam("MetPt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.met.Pt;                   }, .axis_title="MET p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.E[NJET];          }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.Pt[NJET];         }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.Mass[NJET];       }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  
  // Jets
  // AK8
  sh.AddNewFillParam("NJet",           { .nbin=  21, .low=-0.5,   .high=20.5, .fill=[&d](){ return d.jetsAK8.size;                       }, .axis_title="N_{AK8-jet}"});
  sh.AddNewFillParam("NJetSelected",   { .nbin=  21, .low=-0.5,   .high=20.5, .fill=[&d](){ return d.evt.NTopHad+d.evt.NTopLep;          }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetHadronic",   { .nbin=  21, .low=-0.5,   .high=20.5, .fill=[&d](){ return d.evt.NTopHad;                        }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetLeptonic",   { .nbin=  21, .low=-0.5,   .high=20.5, .fill=[&d](){ return d.evt.NTopLep;                        }, .axis_title="N_{Leptonic AK8-jet}"});
  sh.AddNewFillParam("JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.E[NJET];                    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.Pt[NJET];                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.Mass[NJET];                 }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMass",  { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMassCoarse",{ .nbin= 200, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetFilteredMass",{ .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.filteredMass[NJET];         }, .axis_title="AK8-jet Filtered Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTrimmedMass", { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.trimmedMass[NJET];          }, .axis_title="AK8-jet Trimmed Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTopMass",     { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.topMass[NJET];              }, .axis_title="AK8-jet Top Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetMinMass",     { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.minmass[NJET];              }, .axis_title="AK8-jet Min. Subjet-pair Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetNSubJets",    { .nbin=  11, .low=-0.5,   .high=10.5, .fill=[&d](){ return d.jetsAK8.nSubJets[NJET];             }, .axis_title="AK8-jet N_{subjet}"});
  sh.AddNewFillParam("JetTau1",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau1[NJET];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("JetTau2",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("JetTau3",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("JetTau21",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau31",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau32",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("JetTau21Coarse", { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetTau31Coarse", { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau32Coarse", { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("JetDRLep",       { .nbin=  50, .low=   0,   .high=   5, .fill=[&d](){ return d.evt.DRJetLep[d.jetsAK8.it];         }, .axis_title="#DeltaR (lepton, jet)"});
  sh.AddNewFillParam("JetRelPtLep",    { .nbin= 100, .low=   0,   .high= 500, .fill=[&d](){ return d.evt.RelPtJetLep[d.jetsAK8.it];      }, .axis_title="p_{T}^{rel} (lepton, jet) [GeV/c]"});
  
  // Event variables
  sh.AddNewFillParam("NLep",             { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.NLep;                      }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NLepTight",        { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.neletight+d.evt.nmu;       }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMu",              { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.nmu;                       }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEle",             { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.nele;                      }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NEleTight",        { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.neletight;                 }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NLepVeto",         { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.neleveto+d.evt.nmuveto;    }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMuVeto",          { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.nmuveto;                   }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEleVeto",         { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.neleveto;                  }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NTopLep",          { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{top, leptonic}"});
  sh.AddNewFillParam("NTopHad",          { .nbin=   6, .low=-0.5,   .high=  5.5, .fill=[&d](){ return d.evt.NTopHad;                   }, .axis_title="N_{top, hadronic}"});
  sh.AddNewFillParam("TTHadMR",          { .nbin=  50, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMRFine",      { .nbin= 200, .low=   0,   .high=10000, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMTR",         { .nbin=  50, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadMTR;                  }, .axis_title="M_{T,t#bar{t}}^{R} (GeV/c)"});
  sh.AddNewFillParam("TTHadDPhi",        { .nbin=  16, .low=   0,   .high=  3.2, .fill=[&d](){ return fabs(d.evt.TTHadDPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("TTHadDPhiFine",    { .nbin=  64, .low=   0,   .high=  3.2, .fill=[&d](){ return fabs(d.evt.TTHadDPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("TTHadDEta",        { .nbin=  50, .low=   0,   .high=    5, .fill=[&d](){ return d.evt.TTHadDEta;                 }, .axis_title="#Delta#eta_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadDR",          { .nbin=  50, .low=   0,   .high=    5, .fill=[&d](){ return d.evt.TTHadDR;                   }, .axis_title="#DeltaR_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadPz",          { .nbin= 100, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadPz;                   }, .axis_title="P_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadDPz",         { .nbin= 100, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadDPz;                  }, .axis_title="#DeltaP_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadHz",          { .nbin= 100, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadHz;                   }, .axis_title="H_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMass",        { .nbin= 100, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.TTHadMass;                 }, .axis_title="M_{t#bar{t}} (GeV/c^{2})"});
  sh.AddNewFillParam("TTHadR",           { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadRFine",       { .nbin= 100, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadR2",          { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.TTHadR2;                   }, .axis_title="R_{t#bar{t}}^{2}"});
  sh.AddNewFillParam("R",                { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("RFine",            { .nbin= 100, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("R2",               { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.R2;                        }, .axis_title="R^{2}"});
  sh.AddNewFillParam("RFine",            { .nbin= 100, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("MR",               { .nbin=  50, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.MR;                        }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("MTR",              { .nbin=  50, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.MTR;                       }, .axis_title="M_{T}^{R} (GeV/c)"});
  sh.AddNewFillParam("HtTopFraction",    { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.HtTopFr;                   }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("HtExFraction",     { .nbin=  20, .low=   0,   .high=    1, .fill=[&d](){ return d.evt.HtExFr;                    }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("Ht",               { .nbin=  50, .low=   0,   .high=10000, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAllCoarse",      { .nbin=  20, .low=   0,   .high= 6000, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAll",            { .nbin=  50, .low=   0,   .high=10000, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAllFine",        { .nbin= 200, .low=   0,   .high=10000, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtTop",            { .nbin=  25, .low=   0,   .high= 5000, .fill=[&d](){ return d.evt.HtTop;                     }, .axis_title="H_{T,tops} (GeV/c)"});
  sh.AddNewFillParam("HtEx",             { .nbin=  50, .low=   0,   .high=10000, .fill=[&d](){ return d.evt.HtEx;                      }, .axis_title="H_{T,extra} (GeV/c)"});
  
  // Define Cuts here:
  sh.AddNewCut("GoodEle",       [&d](){ return d.ele.Pt[NLEP] > 35 && fabs(d.ele.Eta[NLEP]) < 2.5 && d.ele.isTight[NLEP] > 0; });
  sh.AddNewCut("EleJetCombMass>90", [&d](){ return d.evt.EleJetCombMass[d.ele.it] > 90; });
  
  sh.AddNewCut("GoodMu",        [&d](){ return d.mu.Pt[NLEP] > 45 && fabs(d.mu.Eta[NLEP]) < 2.1 && d.mu.IsTightMuon[NLEP] > 0; });
  sh.AddNewCut("MuJetCombMass>90",  [&d](){ return d.evt.MuJetCombMass[d.mu.it] > 90; });
  
  sh.AddNewCut("AK4Highest2Jet",   [&d](){ return d.jetsAK4.size>=2 && d.jetsAK4.it<2; });
  sh.AddNewCut("AK4Highest3Jet",   [&d](){ return d.jetsAK4.size>=3 && d.jetsAK4.it<3; });
  sh.AddNewCut("Highest2Jet",   [&d](){ return d.jetsAK8.size>=2 && d.jetsAK8.it<2; });
  sh.AddNewCut("Highest3Jet",   [&d](){ return d.jetsAK8.size>=3 && d.jetsAK8.it<3; });
  
  sh.AddNewCut("HadTop",           [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoPtCut",    [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoTauCut",   [&d](){ return d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140; });
  sh.AddNewCut("HadTopNoMassCut",  [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.Pt[NJET] > 400 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoMassNoTauCut", [&d](){ return d.jetsAK8.Pt[NJET] > 400; });
  
  sh.AddNewCut("NTop==2",          [&d](){ return d.evt.NTopHad+d.evt.NTopLep==2; });
  //sh.AddNewCut("NTopHad==2",       [&d](){ return d.evt.NTopHad==2&&d.evt.NTopLep==0; }); // This cut was enforced previously
  sh.AddNewCut("NTopHad==2",       [&d](){ return d.evt.NTopHad==2; });
  sh.AddNewCut("NTopLep==1",       [&d](){ return d.evt.NTopLep==1; });
  sh.AddNewCut("NLepTight==0",     [&d](){ return (d.evt.nmu+d.evt.neletight)==0; });
  sh.AddNewCut("NLepVeto==0",      [&d](){ return (d.evt.nmuveto+d.evt.neleveto)==0; });
  sh.AddNewCut("|DPhi|<2.8",       [&d](){ return fabs(d.evt.TTHadDPhi)<2.8; });
  sh.AddNewCut("ttbar",            [&looper](){ return looper.it_sample==0; });
  sh.AddNewCut("ttbar,qcd",        [&looper](){ return looper.it_sample<2; });
  sh.AddNewCut("noqcd",            [&looper](){ return looper.it_sample!=1; });
  sh.AddNewCut("NLepTight==1",     [&d](){ return (d.evt.neletight+d.evt.nmu)==1; });
  sh.AddNewCut("NEleTight==1",     [&d](){ return d.evt.neletight==1; });
  sh.AddNewCut("NMuTight==1",      [&d](){ return d.evt.nmu==1; });
  
  // Set Histogram weight (empty = 1)
#if TEST == 1
  sh.SetHistoWeights({});
#else
  sh.SetHistoWeights({[&looper,sample_xsec_pb](){ return IntLumi_invfb /*IntLumi (fb-1)*/ * 1000 * sample_xsec_pb[looper.it_sample] / looper.nevents[looper.it_sample]; }});
#endif
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  //sh.AddHistos("mu",     { .fill="MuonEnergy",       .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("mu",     { .fill="MuonPt",           .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("ele",    { .fill="EleEnergy",        .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("ele",    { .fill="ElePt",            .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,1000, 0,0} });
  sh.AddHistos("met",    { .fill="MetPt",            .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  
  //sh.AddHistos("jetsAK4", { .fill="AK4JetEnergy",     .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK4", { .fill="AK4JetPt",         .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK4", { .fill="AK4JetMass",       .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,500, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetEnergy",     .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  //sh.AddHistos("jetsAK8", { .fill="JetPt",         .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetMass",       .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau1",       .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau2",       .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau3",       .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau21",      .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau31",      .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau32",      .pfs={"JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });

  //sh.AddHistos("evt",   { .fill="NTopHad",          .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,4, 0,0} });
  //sh.AddHistos("evt",   { .fill="MR",         .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="MTR",        .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,3000, 0,0} });
  //sh.AddHistos("evt",   { .fill="R",          .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadMR",            .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadMTR",           .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadR",             .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,2500, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhiFine",      .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDEta",      .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,4, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDR",        .pfs={Samples}, .cuts={"NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadMass",          .pfs={Samples}, .cuts={"noqcd","NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadPz",            .pfs={Samples}, .cuts={"noqcd","NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadHz",            .pfs={Samples}, .cuts={"noqcd","NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPz",           .pfs={Samples}, .cuts={"noqcd","NTopHad==2"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  
  // 01/16
  //// Signal cut variables
  //sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTop",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadR",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadR2",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R2",         .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi",   .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="", .ranges={0,0, 0,0} });
  //
  //// 2D Correlation plots
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadDPhi",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_R",               .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadR",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadDPhi", .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_R",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadR",      .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R_vs_TTHadDPhi",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadR_vs_TTHadDPhi",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  
  //// 01/17
  //// Signal cut variables
  //sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples},                   .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="HtTop",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="TTHadR",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="TTHadR2",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="R",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="R2",         .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi",   .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //
  //// 2D Correlation plots
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadDPhi",        .pfs={Samples},                   .cuts={"NTopHad==2"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_R",               .pfs={Samples},                   .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadR",             .pfs={Samples},                   .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadDPhi", .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_R",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadR",      .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2","|DPhi|<2.8"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R_vs_TTHadDPhi",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadR_vs_TTHadDPhi",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NTopHad==2"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //
  //// 01/20
  //// Signal cut variables
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p25,RAbove0p25"},                       .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p25,RAbove0p25","SideBand,Signal"},     .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p3,RAbove0p3"},                         .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p3,RAbove0p3","SideBand,Signal"},       .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p35,RAbove0p35"},                       .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={"RBelow0p35,RAbove0p35","SideBand,Signal"},     .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R",        .pfs={"DPhiBelow2p8,DPhiAbove2p8"},                   .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R",        .pfs={"DPhiBelow2p8,DPhiAbove2p8","SideBand,Signal"}, .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R",        .pfs={"SideBand,Signal","DPhiBelow2p8,DPhiAbove2p8"}, .cuts={"ttbar","NTopHad==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  
  // latest
  // Top Tagging (CMS-PAS-JME-13-007)
  // 3 >= subjets, Mmin > 50, tau32 < 0.55, 250 > Mjet > 140
  //sh.AddHistos("jetsAK8", { .fill="JetPt",         .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetMass",       .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetPrunedMass", .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau1",       .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau2",       .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau3",       .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau21",      .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau31",      .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetTau32",      .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0,0} });
  //  N-1 plots
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32Coarse", .pfs={Samples}, .cuts={"HadTopNoMassNoTauCut"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetTau21Coarse_vs_JetTau32Coarse",      .pfs={Samples}, .cuts={"HadTopNoMassNoTauCut"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetTau21",        .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31",        .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",        .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",   .pfs={Samples}, .cuts={"HadTopNoMassCut"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass", .pfs={Samples}, .cuts={"HadTopNoMassCut"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass ", .pfs={Samples}, .cuts={"HadTopNoMassCut"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetPt",           .pfs={Samples}, .cuts={"HadTopNoPtCut"},   .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",      .pfs={Samples}, .cuts={"HadTop"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",      .pfs={Samples}, .cuts={"HadTop"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("jetsAK8", { .fill="JetNSubJets",     .pfs={Samples}, .cuts={"HadTop"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  
  // Number of tagged tops
  sh.AddHistos("jetsAK8", { .fill="NJet",          .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("jetsAK8", { .fill="NJetSelected",  .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("jetsAK8", { .fill="NJetHadronic",  .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  //sh.AddHistos("jetsAK8", { .fill="NJetLeptonic",  .pfs={Samples}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  
  const char* cut1 = "R0p32";
  const char* cut2 = "DPhi2p8";
  const char* cut3 = "HtAll1450";
  
  // 02/24 - Razor presentation
  // Signal cut variables
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut1,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadR",          .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="Ht",              .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",           .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtEx" ,           .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtExFraction",    .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi",       .pfs={Samples,cut1,cut2,cut3}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  
  // 2D Correlation plots
  sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadDPhi",   .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_R",           .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtAllCoarse_vs_TTHadR",      .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadDPhi", .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtTopFraction_vs_R",         .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtTopFraction_vs_TTHadR",    .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="R_vs_TTHadDPhi",             .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TTHadR_vs_TTHadDPhi",        .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  
  // Latest
  // Try MR instead of HtAll
  sh.AddHistos("evt",   { .fill="TTHadMR",   .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="MR",        .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadMR",   .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="MR",        .pfs={Samples,cut1}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadMR",   .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="MR",        .pfs={Samples,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadMR",   .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="MR",        .pfs={Samples,cut1,cut2}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  // Plots for Background estimation
  sh.AddHistos("evt",   { .fill="R",         .pfs={cut2,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={cut3,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={cut2,cut3,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={cut3,cut2,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={cut1,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={cut3,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={cut1,cut3,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="TTHadDPhi", .pfs={cut3,cut1,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",     .pfs={cut1,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",     .pfs={cut2,Samples},      .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",     .pfs={cut1,cut2,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="HtAll",     .pfs={cut2,cut1,Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="Log", .ranges={0,0, 0.001,100000} });
  // 3D Plots to get best signal cuts (Maximize Smin)
  sh.AddHistos("evt",   { .fill="HtAllFine_vs_TTHadDPhiFine_vs_RFine",        .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtAllFine_vs_TTHadDPhiFine_vs_TTHadRFine",   .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TTHadMRFine_vs_TTHadDPhiFine_vs_TTHadRFine", .pfs={Samples}, .cuts={"NTopHad==2","NLepVeto==0"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  
  // 04 March
  sh.AddHistos("evt",   { .fill="NTopHad",   .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NTopLep",   .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NLep",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NEle",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NMu",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NLepTight", .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NEleTight", .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NLepVeto",  .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NEleVeto",  .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("evt",   { .fill="NMuVeto",   .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,1.1} });
  sh.AddHistos("jetsAK8", { .fill="JetDRLep",                .pfs={Samples},                 .cuts={"NLepTight==1"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetRelPtLep",             .pfs={Samples},                 .cuts={"NLepTight==1"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={Samples},                 .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={Samples,"NSubJet"},       .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={Samples,"JetsPtOrdered"}, .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("ele",   { .fill="EleJetCombMass",              .pfs={Samples},               .cuts={"NEleTight==1","GoodEle"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("ele",   { .fill="EleDRJet",                    .pfs={Samples},               .cuts={"NEleTight==1","GoodEle","EleJetCombMass>90"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("ele",   { .fill="EleRelPtJet",                 .pfs={Samples},               .cuts={"NEleTight==1","GoodEle","EleJetCombMass>90"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("mu",    { .fill="MuJetCombMass",               .pfs={Samples},               .cuts={"NMuTight==1","GoodMu"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("mu",    { .fill="MuDRJet",                     .pfs={Samples},               .cuts={"NMuTight==1","GoodMu","MuJetCombMass>90"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("mu",    { .fill="MuRelPtJet",                  .pfs={Samples},               .cuts={"NMuTight==1","GoodMu","MuJetCombMass>90"}, .draw="", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  std::cout<<"-----------------------------------------------------------------\n";
  std::cout<<"Creating the following plots:\n"; sh.PrintNames();
  std::cout<<"-----------------------------------------------------------------\n";
  
  TFile *file;
  if (Run) {
    if (filelist_fromshell.size()) {
      std::cout<<"Adding "<<filelist_fromshell.size()<<" files from the shell arguments.\n";
      for (size_t i=0; i<filelist_fromshell.size(); ++i) looper.AddFile(filelist_fromshell[i], !i);
    } else {
      std::vector<std::string> dirs = samples.GetListOfDirectories();
      for ( std::string dir : dirs ) looper.AddFile(dir);
    }
    bool debug = 0;
    if (debug) cout<<"Start ok\n";
    while (looper.LoopOnSamples()) {
      if (debug) cout<<"Sample ok\n";
      while (looper.LoopOnFiles()) {
        if (debug) cout<<"File ok\n";
        TFile *curr_file = looper.CurrentFile();
        reader.Load_Tree(*curr_file,looper.TreeName());
        if (debug) cout<<"TreeReader ok\n";
        while(looper.LoopOnEntries()) {
          reader.GetEntry(looper.CurrentEntry());
          d = reader.data;
          if (debug) cout<<"GetEntry ok\n";
          d.CalculateAllVariables();
          if (debug) cout<<"CalcVar ok\n";
          
          // loop on objects and fill their histos
          while(d.mu.Loop())     sh.Fill("mu");  if (debug) cout<<"Fill Muons ok\n";
          while(d.ele.Loop())    sh.Fill("ele"); if (debug) cout<<"Fill Electrons ok\n";
          while(d.jetsAK4.Loop()) {
            //std::cout<<d.jetsAK4.it<<" "<<d.jetsAK4.Pt[NJET]<<std::endl;
            sh.Fill("jetsAK4");
          } if (debug) cout<<"Fill AK4Jets ok\n";
          while(d.jetsAK8.Loop()) sh.Fill("jetsAK8"); if (debug) cout<<"Fill AK8Jets ok\n";
          sh.Fill("met"); if (debug) cout<<"Fill MET ok\n";
          sh.Fill("evt"); if (debug) cout<<"Fill Evt ok\n";
        }
        curr_file->Close();
      }
    }
    std::cout<<std::endl;
  } else {
    std::cout<<"Loading Histos from file: "<<inputfile<<std::endl;
    sh.Load(inputfile.c_str());
  }
  std::cout<<"Finished ..."<<std::endl;
  std::cout<<"Writing Histograms to File: "<<outputfile<<std::endl;
  file = new TFile(outputfile.c_str(),"recreate");
  sh.DrawPlots();
  sh.Write();
  file->Close();
  
  return 0;
}

