#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"

#define NTHSTAT 1

int main(int argc, char* argv[]) {
#include "Samples_Feb04.txt"
  
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
  sh.AddHistoType("jetsCmsTopTag");
  sh.AddHistoType("met");
  
  
  // Define Postfixes here:
  sh.AddNewPostfix("ttbar,qcd",                 [&looper](){ return looper.it_sample; }, "ttbar;qcd", "t#bar{t};QCD", "2,6");
  sh.AddNewPostfix("ttbar,qcd,Susy3,Susy4",     [&looper](){ return looper.it_sample; }, "ttbar;qcd;susy3body;susy4body", "t#bar{t};QCD;T5tttt - 3body;T5tttt - 4body", "2,6,4,3");
  //std::vector<double> Sample_xsec_pb;
  sh.AddNewPostfix("AllSamples",                [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, PF_names, Latex_names, "1-9");
  //const char* Samples = "ttbar,qcd,Susy3,Susy4";
  const char* Samples = "AllSamples";
  
  sh.AddNewPostfix("SideBand,Signal",           [&d](){ return d.evt.HTall > 1500; }, "SideBand;Signal", "H_{T,all} < 1500 GeV/c;H_{T,all} > 1500 GeV/c", "4,2");
  sh.AddNewPostfix("RBelow0p25,RAbove0p25",     [&d](){ return d.jetsAK8.R > 0.25; },  "RBelow0p25;RAbove0p25", "R < 0.25;R > 0.25", "4,2");
  sh.AddNewPostfix("DPhiBelow2p8,DPhiAbove2p8", [&d](){ return fabs(d.evt.tt_dPhi) > 2.8; },  "DPhiBelow2p8;DPhiAbove2p8", "#Delta#phi_{t#bar{t}} < 2.8;#Delta#phi_{t#bar{t}} > 2.8", "4,2");
  sh.AddNewPostfix("NCmsLepTop",                [&d](){ return d.evt.ncmsleptops; },  "[0to8]", "[0to8] Cms Leptonic Top", "1-8");
  
  //sh.AddNewPostfix("ttbar,Susy3,Susy4", &looper.it_sample, "ttbar;susy3body;susy4body", "SM t#bar{t};T5tttt - 3body;T5tttt - 4body", "2,4,3");
  sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetsAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("AK8JetsPtOrdered", [&d](){ return d.jetsAK8.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh.AddNewFillParam("MuonEnergy",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.E[NLEP];              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuonPt",            { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.Pt[NLEP];             }, .axis_title="Muon p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.E[NLEP];             }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.Pt[NLEP];            }, .axis_title="Electron p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("MetPt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.met.Pt;                  }, .axis_title="MET p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.E[NJET];                    }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.Pt[NJET];                   }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK4.Mass[NJET];                 }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  
  // Jets
  // AK8
  sh.AddNewFillParam("NAK8Jet",           { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.jetsAK8.size;                       }, .axis_title="N_{AK8-jet}"});
  sh.AddNewFillParam("NAK8JetSelected",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.nhadtops+d.evt.nleptops;        }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NAK8JetHadronic",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.nhadtops;                       }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NAK8JetLeptonic",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.nleptops;                       }, .axis_title="N_{Leptonic AK8-jet}"});
  sh.AddNewFillParam("AK8JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.E[NJET];                    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("AK8JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.Pt[NJET];                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.Mass[NJET];                 }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetPrunedMass",  { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("2dAK8JetPrunedMass",{ .nbin= 200, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetTau1",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau1[NJET];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau2",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("AK8JetTau3",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("AK8JetTau21",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau31",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau32",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("2dAK8JetTau21",     { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau2[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("2dAK8JetTau31",     { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("2dAK8JetTau32",     { .nbin=  50, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("AK8JetMR",          { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetsAK8.MR;                   }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMTR",         { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetsAK8.MTR;                  }, .axis_title="M_{T}^{R} (GeV/c)"});
  // CMS Top tagged jets
  sh.AddNewFillParam("NCmsTopTagJet",           { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.jetsCmsTopTag.size;                       }, .axis_title="N_{CmsTopTag-jet}"});
  sh.AddNewFillParam("NCmsTopTagJetSelected",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.ncmshadtops+d.evt.ncmsleptops;        }, .axis_title="N_{Hadronic CmsTopTag-jet}"});
  sh.AddNewFillParam("NCmsTopTagJetHadronic",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.ncmshadtops;                          }, .axis_title="N_{Hadronic CmsTopTag-jet}"});
  sh.AddNewFillParam("NCmsTopTagJetLeptonic",   { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.ncmsleptops;                          }, .axis_title="N_{Leptonic CmsTopTag-jet}"});
  sh.AddNewFillParam("CmsTopTagJetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsCmsTopTag.E[NJET];                    }, .axis_title="CmsTopTag-jet Energy (GeV)"});
  sh.AddNewFillParam("CmsTopTagJetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsCmsTopTag.Pt[NJET];                   }, .axis_title="CmsTopTag-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("CmsTopTagJetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsCmsTopTag.Mass[NJET];                 }, .axis_title="CmsTopTag-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("CmsTopTagJetPrunedMass",  { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetsCmsTopTag.prunedMass[NJET];           }, .axis_title="CmsTopTag-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("CmsTopTagJetTau1",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau1[NJET];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("CmsTopTagJetTau2",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau2[NJET];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("CmsTopTagJetTau3",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau3[NJET];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("CmsTopTagJetTau21",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau2[NJET]/d.jetsCmsTopTag.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("CmsTopTagJetTau31",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau3[NJET]/d.jetsCmsTopTag.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("CmsTopTagJetTau32",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsCmsTopTag.tau3[NJET]/d.jetsCmsTopTag.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("CmsTopTagJetMR",          { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetsCmsTopTag.MR;                   }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("CmsTopTagJetMTR",         { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetsCmsTopTag.MTR;                  }, .axis_title="M_{T}^{R} (GeV/c)"});
  
  // Event variables
  sh.AddNewFillParam("NHadTop",           { .nbin=   6, .low=-0.5,   .high= 5.5, .fill=[&d](){ return d.evt.nhadtops;                }, .axis_title="N_{hadronic top}"});
  sh.AddNewFillParam("TT_MR",             { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_MR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TT_MTR",            { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_MTR;                  }, .axis_title="M_{T,t#bar{t}}^{R} (GeV/c)"});
  sh.AddNewFillParam("TT_DeltaR",         { .nbin=  50, .low=   0,   .high=   5, .fill=[&d](){ return d.evt.tt_dR;                   }, .axis_title="#DeltaR_{t#bar{t}}"});
  sh.AddNewFillParam("TT_DeltaPhi",       { .nbin=  64, .low=-3.2,   .high= 3.2, .fill=[&d](){ return d.evt.tt_dPhi;                 }, .axis_title="#Delta#phi_{t#bar{t}}"});
  sh.AddNewFillParam("TT_DeltaEta",       { .nbin=  50, .low=   0,   .high=   5, .fill=[&d](){ return d.evt.tt_dEta;                 }, .axis_title="#Delta#eta_{t#bar{t}}"});
  sh.AddNewFillParam("TT_Mass",           { .nbin= 100, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_Mass;                 }, .axis_title="M_{t#bar{t}} (GeV/c^{2})"});
  sh.AddNewFillParam("TT_Pz",             { .nbin= 100, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_Pz;                   }, .axis_title="P_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TT_Hz",             { .nbin= 100, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_Hz;                   }, .axis_title="H_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TT_dPz",            { .nbin= 100, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.tt_dPz;                  }, .axis_title="#DeltaP_{Z,t#bar{t}} (GeV/c)"});
  
  sh.AddNewFillParam("TT_AbsDeltaPhi",    { .nbin=  16, .low=   0,   .high= 3.2, .fill=[&d](){ return fabs(d.evt.tt_dPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("TT_R",              { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.tt_R;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TT_R2",             { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.tt_R2;                   }, .axis_title="R_{t#bar{t}}^{2}"});
  sh.AddNewFillParam("AK8JetR",           { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.R;                    }, .axis_title="R"});
  sh.AddNewFillParam("AK8JetR2",          { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.R2;                   }, .axis_title="R^{2}"});
  sh.AddNewFillParam("H_TttFraction",     { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTttFraction;            }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("H_TexFraction",     { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTexFraction;            }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("H_T",               { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HT;                      }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("H_Tevt",            { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTall;                   }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("H_Ttt",             { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTtt;                    }, .axis_title="H_{T,tops} (GeV/c)"});
  sh.AddNewFillParam("H_Tex",             { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTex;                    }, .axis_title="H_{T,extra} (GeV/c)"});
  
  // 2D binning of Signal cut variables
  sh.AddNewFillParam("DeltaPhi",          { .nbin=  16, .low=   0,   .high= 3.2, .fill=[&d](){ return fabs(d.evt.tt_dPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("Rtt",               { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.tt_R;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("R",                 { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetsAK8.R;                   }, .axis_title="R"});
  sh.AddNewFillParam("HT",                { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HT;                      }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("HTall",             { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HTall;                   }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HTex",              { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HTex;                    }, .axis_title="H_{T,extra} (GeV/c)"});
  sh.AddNewFillParam("HTttFraction",      { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTttFraction;            }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  
  // Define Cuts here:
  sh.AddNewCut("AK4Highest2Jet",   [&d](){ return d.jetsAK4.size>=2 && d.jetsAK4.it<2; });
  sh.AddNewCut("AK4Highest3Jet",   [&d](){ return d.jetsAK4.size>=3 && d.jetsAK4.it<3; });
  sh.AddNewCut("AK8Highest2Jet",   [&d](){ return d.jetsAK8.size>=2 && d.jetsAK8.it<2; });
  sh.AddNewCut("AK8Highest3Jet",   [&d](){ return d.jetsAK8.size>=3 && d.jetsAK8.it<3; });
  
  sh.AddNewCut("HadTop",           [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoPtCut",    [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoTauCut",   [&d](){ return d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140; });
  sh.AddNewCut("HadTopNoMassCut",  [&d](){ return d.jetsAK8.tau2[NJET]>0 && d.jetsAK8.tau3[NJET]>0 ? d.jetsAK8.Pt[NJET] > 400 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoMassNoTauCut", [&d](){ return d.jetsAK8.Pt[NJET] > 400; });
  sh.AddNewCut("CmsHadTop",           [&d](){ return d.jetsCmsTopTag.tau2[NJET]>0 && d.jetsCmsTopTag.tau3[NJET]>0 ? d.jetsCmsTopTag.Pt[NJET] > 400 && d.jetsCmsTopTag.prunedMass[NJET] > 140 && (d.jetsCmsTopTag.tau3[NJET]/d.jetsCmsTopTag.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("CmsHadTopNoPtCut",    [&d](){ return d.jetsCmsTopTag.tau2[NJET]>0 && d.jetsCmsTopTag.tau3[NJET]>0 ? d.jetsCmsTopTag.prunedMass[NJET] > 140 && (d.jetsCmsTopTag.tau3[NJET]/d.jetsCmsTopTag.tau2[NJET]) < 0.75 : 0; });
  sh.AddNewCut("CmsHadTopNoTauCut",   [&d](){ return d.jetsCmsTopTag.Pt[NJET] > 400 && d.jetsCmsTopTag.prunedMass[NJET] > 140; });
  sh.AddNewCut("CmsHadTopNoMassCut",  [&d](){ return d.jetsCmsTopTag.tau2[NJET]>0 && d.jetsCmsTopTag.tau3[NJET]>0 ? d.jetsCmsTopTag.Pt[NJET] > 400 && (d.jetsCmsTopTag.tau3[NJET]/d.jetsCmsTopTag.tau2[NJET]) < 0.75 : 0; });
  
  sh.AddNewCut("NTop==2",          [&d](){ return d.evt.nhadtops+d.evt.nleptops==2; });
  sh.AddNewCut("NHadTop==2",       [&d](){ return d.evt.nhadtops==2; });
  sh.AddNewCut("|DPhi|<2.8",       [&d](){ return fabs(d.evt.tt_dPhi)<2.8; });
  sh.AddNewCut("NCmsTop==2",       [&d](){ return d.jetsCmsTopTag.size==2; });
  sh.AddNewCut("NCmsSelTop==2",    [&d](){ return d.evt.ncmshadtops+d.evt.ncmsleptops==2; });
  sh.AddNewCut("NCmsHadTop==2",    [&d](){ return d.evt.ncmshadtops==2; });
  sh.AddNewCut("NCmsLepTop==2",    [&d](){ return d.evt.ncmsleptops==2; });
  sh.AddNewCut("ttbar",            [&looper](){ return looper.it_sample==0; });
  sh.AddNewCut("ttbar,qcd",        [&looper](){ return looper.it_sample<2; });
  sh.AddNewCut("noqcd",            [&looper](){ return looper.it_sample!=1; });
  
  // Set Histogram weight (empty = 1)
  sh.SetHistoWeights({[&looper](){ return 10000 * (looper.it_sample==0 ? 806.1 /*ttbar*/ : looper.it_sample==1 ? 1722 /*qcd*/ : looper.it_sample==2 ? 0.0460525 /*susy3body*/ : looper.it_sample==3 ? 0.0460525 /*susy4body*/ : 0) / looper.nevents[looper.it_sample]; }});
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  sh.AddHistos("mu",     { .fill="MuonEnergy",       .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("mu",     { .fill="MuonPt",           .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="EleEnergy",        .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="ElePt",            .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,1000, 0,0} });
  sh.AddHistos("met",    { .fill="MetPt",            .pfs={Samples}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  
  //sh.AddHistos("jetsAK4", { .fill="AK4JetEnergy",     .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK4", { .fill="AK4JetPt",         .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK4", { .fill="AK4JetMass",       .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,500, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetEnergy",     .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  //sh.AddHistos("jetsAK8", { .fill="AK8JetPt",         .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetMass",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau1",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau2",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau3",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau21",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau31",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau32",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });

  sh.AddHistos("evt",   { .fill="NHadTop",          .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,4, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetMR",         .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetMTR",        .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,3000, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_MR",            .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_MTR",           .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R",             .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,2500, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaPhi",      .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaEta",      .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,4, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaR",        .pfs={Samples}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Mass",          .pfs={Samples}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Pz",            .pfs={Samples}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Hz",            .pfs={Samples}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_dPz",           .pfs={Samples}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // 01/16
  //// Signal cut variables
  //sh.AddHistos("evt",   { .fill="H_T",              .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Tevt",           .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Ttt",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Tex" ,           .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_TttFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_TexFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_R",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_R2",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="AK8JetR2",         .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi",   .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //
  //// 2D Correlation plots
  //sh.AddHistos("evt",   { .fill="HTall_vs_DeltaPhi",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTall_vs_R",               .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTall_vs_Rtt",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_DeltaPhi", .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_R",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_Rtt",      .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R_vs_DeltaPhi",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="Rtt_vs_DeltaPhi",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  // 01/17
  // Signal cut variables
  sh.AddHistos("evt",   { .fill="H_T",              .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Tevt",           .pfs={Samples},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Ttt",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Tex" ,           .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_TttFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_TexFraction",    .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R",             .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R2",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR2",         .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi",   .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // 2D Correlation plots
  sh.AddHistos("evt",   { .fill="HTall_vs_DeltaPhi",        .pfs={Samples},                   .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTall_vs_R",               .pfs={Samples},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTall_vs_Rtt",             .pfs={Samples},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_DeltaPhi", .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_R",        .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_Rtt",      .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="R_vs_DeltaPhi",            .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="Rtt_vs_DeltaPhi",          .pfs={Samples,"SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  // 01/20
  // Signal cut variables
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi", .pfs={"RBelow0p25,RAbove0p25"},                       .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi", .pfs={"RBelow0p25,RAbove0p25","SideBand,Signal"},     .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"DPhiBelow2p8,DPhiAbove2p8"},                   .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"DPhiBelow2p8,DPhiAbove2p8","SideBand,Signal"}, .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"SideBand,Signal","DPhiBelow2p8,DPhiAbove2p8"}, .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  
  // Warning: Currently Pruned mass for AK8Jets is not available
  // For above plot, I had to switch back to normal Mass for the NHadTop cut
  // Below I try the CMS top tagging jets to better tune our top selection
  
  // latest
  // CMS top tagging plots (tau, prunedMass etc. doesn't work yet)
  //sh.AddHistos("jetsAK8", { .fill="AK8JetPt",         .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetMass",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetPrunedMass", .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau1",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau2",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau3",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau21",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau31",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="AK8JetTau32",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetPt",         .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetMass",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetPrunedMass", .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau1",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau2",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau3",       .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau21",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau31",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau32",      .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // Top Tagging (CMS-PAS-JME-13-007)
  // 3 >= subjets, Mmin > 50, tau32 < 0.55, 250 > Mjet > 140
  //  N-1 plots
  sh.AddHistos("jetsAK8", { .fill="AK8JetTau21",                     .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AK8JetTau31",                     .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AK8JetTau32",                     .pfs={Samples}, .cuts={"HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AK8JetPrunedMass",                .pfs={Samples}, .cuts={"HadTopNoMassCut"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AK8JetPt",                        .pfs={Samples}, .cuts={"HadTopNoPtCut"},   .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="2dAK8JetPrunedMass_vs_2dAK8JetTau32", .pfs={Samples}, .cuts={"HadTopNoMassNoTauCut"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="2dAK8JetTau21_vs_2dAK8JetTau32",  .pfs={Samples}, .cuts={"HadTopNoMassNoTauCut"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,500} });
  
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetMass",       .pfs={Samples}, .cuts={"CmsHadTopNoMassCut"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetPrunedMass", .pfs={Samples}, .cuts={"CmsHadTopNoMassCut"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetPt",         .pfs={Samples}, .cuts={"CmsHadTopNoPtCut"},   .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau21",      .pfs={Samples}, .cuts={"CmsHadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau31",      .pfs={Samples}, .cuts={"CmsHadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsCmsTopTag", { .fill="CmsTopTagJetTau32",      .pfs={Samples}, .cuts={"CmsHadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // Number of tagged tops
  sh.AddHistos("jetsAK8", { .fill="NAK8Jet",          .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="NAK8JetSelected",  .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="NAK8JetHadronic",  .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="NAK8JetLeptonic",  .pfs={Samples}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  
  std::cout<<"-----------------------------------------------------------------\n";
  std::cout<<"Creating the following plots:\n"; sh.PrintNames();
  std::cout<<"-----------------------------------------------------------------\n";
  
  TFile *file;
  if (Run) {
    if (filelist.size()) {
      std::cout<<"Adding "<<filelist.size()<<" files from the shell arguments.\n";
      for (size_t i=0; i<filelist.size(); ++i) looper.AddFile(filelist[i]);
    } else {
      //looper.AddFile("../../../b2gTree_T2tt_2J_mStop-650_mLSP-325.root");
      //  looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/ttbar/*.root");
      //  looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/qcd/*.root");
      //  looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/susy3body/*.root");
      //  looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/susy4body/*.root");
      
      // test
      looper.AddFile("../../../B2GTTreeNtupleExtra_susy.root");
      
      // Feb04/feb12 background samples
      //for (std::string dirname : Sample_dirnames) looper.AddFile(std::string("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/")+dirname+"/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/susy4body/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/susy3body/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/TT/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/T_tW-channel/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/Tbar_tW-channel/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/DYJetsToLL/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/WJetsToLNu/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/WJetsToLNu_HT-100to200/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/WJetsToLNu_HT-200to400/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/WJetsToLNu_HT-400to600/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/WJetsToLNu_HT-600toInf/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_Pt_20to30_bcToE/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_Pt_30to80_bcToE/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_Pt_80to170_bcToE/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_Pt_170toInf_bcToE/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_HT-100To250/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_HT_250To500/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_HT-500To1000/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/QCD_HT_1000ToInf/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-200To500/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-500To1000/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-1000To2000/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-2000To4000/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-4000To8000/*.root");
      // looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb12_edm_Feb04/GGJets_M-8000To13000/*.root");
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
	  
          sh.Fill("evt");
	  if (debug) cout<<"Fill Evt ok\n";
          // loop on objects and fill their histos
          while(d.mu.Loop())     sh.Fill("mu");
	  if (debug) cout<<"Fill Muons ok\n";
          while(d.ele.Loop())    sh.Fill("ele");
	  if (debug) cout<<"Fill Electrons ok\n";
          while(d.jetsAK4.Loop()) {
            //std::cout<<d.jetsAK4.it<<" "<<d.jetsAK4.Pt[NJET]<<std::endl;
            sh.Fill("jetsAK4");
          }
	  if (debug) cout<<"Fill AK4Jets ok\n";
          while(d.jetsAK8.Loop()) sh.Fill("jetsAK8");
	  if (debug) cout<<"Fill AK8Jets ok\n";
          while(d.jetsCmsTopTag.Loop()) sh.Fill("jetsCmsTopTag");
	  if (debug) cout<<"Fill Tops ok\n";
	  sh.Fill("met");
	  if (debug) cout<<"Fill Met ok\n";
        }
        curr_file->Close();
      }
    }
    
    std::cout<<std::endl;
  } else {
    std::cout<<"Loading Histos from file: "<<outputfile<<std::endl;
    sh.Load(outputfile.c_str());
  }
  std::cout<<"Finished ..."<<std::endl;
  std::cout<<"Writing Histograms to File: "<<outputfile<<std::endl;
  file = new TFile(outputfile.c_str(),"recreate");
  sh.DrawPlots();
  sh.Write();
  file->Close();
  
  return 0;
}

