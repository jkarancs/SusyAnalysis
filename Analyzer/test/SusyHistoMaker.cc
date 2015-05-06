// Simplified version of B2GPlotter
// Latest FillParameters, Cuts and Postfixes are kept to see how the SmartHistos framework works
// Class conventions:
// 
//  Samples:
//    AddSample ( "sample name postfix", "title on legend", "TColor", Array of { { .dir = "sample_dir/filnames.root", .xsec_pb = Sample XSEC (pb) } } )
//

#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/Samples.h"
#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"

#define TEST 0 // If you want to switch to some sample and plot quickly
#if TEST == 0
#define NTHSTAT 100 // You can reduce the number of entries to run on to 1/NTHSTAT (again for quick plotting)
#else
#define NTHSTAT 1
#endif

#define IntLumi_invfb 10.0 // What is the Integrated luminosity you want to normalise to

int main(int argc, char* argv[]) {
  // Get arguments from shell
  std::vector<std::string> filelist_fromshell;
  std::string inputfile="";
  std::string outputfile="plots.root";
  // -o <output file> option:
  // Specify the output root filename
  // Rest of the arguments are treated as files added
  // Note:
  // If using postfixes with the v.pf_file_add variable
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
  bool test = TEST;

  // these numbers are needed so some of the samples can be ignored but still added to some plots specifically
  const size_t nsignal = 4; // You can show all types of signals or just one you pick specifically
  const size_t nttbar = 3; // Only TTJets is used, but for some plots you can compare different TT MCs
  if (filelist_fromshell.size()) {
    samples.AddSample("test", "test", "1", { { .dir="", .xsec_pb=100 } });
  } else {
    if (test) {
      samples.AddSample("TTJets", "t#bar{t}+jets", "633", { { .dir="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Apr02_edm_Apr01/TTJets/*.root", .xsec_pb=831.76 } }); //
    } else {
      std::string sample_dir = "/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Apr02_edm_Apr01/";
      if (nsignal>0)
	samples.AddSample("T5ttttDeg_4bodydec_mGo1300", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.3TeV)", "1",
			  { { .dir=sample_dir+"T5ttttDeg_mGo1300_4bodydec/*.root", .xsec_pb=0.0460525 } });
      if (nsignal>1)
	samples.AddSample("T5ttttDeg_3bodydec_mGo1300", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.3TeV)", "12",
			  { { .dir=sample_dir+"T5ttttDeg_mGo1300_23bodydec/*.root", .xsec_pb=0.0460525 } });
      if (nsignal>2)
	samples.AddSample("T5ttttDeg_4bodydec_mGo1000", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.0TeV)", "14",
			  { { .dir=sample_dir+"T5ttttDeg_mGo1000_4bodydec/*.root", .xsec_pb=0.325388 } });
      if (nsignal>3)
	samples.AddSample("T5ttttDeg_3bodydec_mGo1000", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.0TeV)", "16",
			  { { .dir=sample_dir+"T5ttttDeg_mGo1000_23bodydec/*.root", .xsec_pb=0.325388 } });
      if (nttbar>2)
	samples.AddSample("TT_FastSim_73X", "t#bar{t} (FastSim)", "6", { { .dir=sample_dir+"TT_FastSim_73X/*.root", .xsec_pb=831.76 } }); //
      if (nttbar>1)
	samples.AddSample("TT", "t#bar{t}", "2", { { .dir=sample_dir+"TT/*.root", .xsec_pb=831.76 } }); //
      samples.AddSample("TTJets", "t#bar{t}+jets", "634", { { .dir=sample_dir+"TTJets/*.root", .xsec_pb=831.76 } }); //
      samples.AddSample("WJets", "W+Jets", "3",
        		{ { .dir=sample_dir+"WJetsToLNu_HT-100to200/*.root", .xsec_pb=2234.9 },
        		  { .dir=sample_dir+"WJetsToLNu_HT-200to400/*.root", .xsec_pb=580.1 }, 
        		  { .dir=sample_dir+"WJetsToLNu_HT-400to600/*.root", .xsec_pb=68.4 },
            	          { .dir=sample_dir+"WJetsToLNu_HT-600toInf/*.root", .xsec_pb=23.14 } }); //
      samples.AddSample("ZJets", "Z+Jets", "401",
        		{ { .dir=sample_dir+"ZJetsToNuNu_HT-100to200/*.root", .xsec_pb=372.6 },
        		  { .dir=sample_dir+"ZJetsToNuNu_HT-200to400/*.root", .xsec_pb=100.8 }, 
        		  { .dir=sample_dir+"ZJetsToNuNu_HT-400to600/*.root", .xsec_pb=11.99 },
            	          { .dir=sample_dir+"ZJetsToNuNu_HT-600toInf/*.root", .xsec_pb=4.113 } }); //
      samples.AddSample("T_tW", "single t/#bar{t} (tW channel)", "617",
                        { { .dir=sample_dir+"T_tW-channel/*.root", .xsec_pb=35.6 }, //
                          { .dir=sample_dir+"Tbar_tW-channel/*.root", .xsec_pb=35.6 } }); //
      samples.AddSample("QCD", "QCD", "4",
        		{ { .dir=sample_dir+"QCD_HT-100To250/*.root",  .xsec_pb=28730000 },
                          { .dir=sample_dir+"QCD_HT_250To500/*.root",  .xsec_pb=670500 },
        		  { .dir=sample_dir+"QCD_HT-500To1000/*.root", .xsec_pb=26740 }, //
            	          { .dir=sample_dir+"QCD_HT_1000ToInf/*.root", .xsec_pb=769.7 } }); //
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
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jetsAK4");
  sh.AddHistoType("jetsAK8");
  sh.AddHistoType("met");
  sh.AddHistoType("evt");
  sh.AddHistoType("gen");
  
  // Histogram weights
  // here you can specify how you want to normalise to XSec for
  // a given integrated luminosity
#if TEST == 1
  // (empty = 1)
  sh.SetHistoWeights({});
#else
  sh.SetHistoWeights({[&looper,sample_xsec_pb](){ return IntLumi_invfb /*IntLumi (fb-1)*/ * 1000 * sample_xsec_pb[looper.it_sample] / looper.nevents[looper.it_sample]; }});
#endif
  
  // --------------------------------------------------------------------
  //                             CUTS
  // --------------------------------------------------------------------
  
  // Sample
  sh.AddNewCut("Background",         [&looper](){ return looper.it_sample>4; });
  sh.AddNewCut("NoDoubleTT",         [&looper,nsignal,nttbar](){ return !((looper.it_sample>=nsignal)&&(looper.it_sample<nsignal+nttbar-1)); });
  sh.AddNewCut("TTJets",             [&looper,nsignal,nttbar](){ return looper.it_sample==nsignal+nttbar-1; });
  sh.AddNewCut("NoMultiJet",         [&looper,nsignal,nttbar](){ return looper.it_sample<nsignal+nttbar; });
  
  // Ele
  sh.AddNewCut("GoodEle",       [&d](){ return d.ele.Pt[NLEP] > 35 && fabs(d.ele.Eta[NLEP]) < 2.5 && d.ele.isTight[NLEP] > 0; });
  sh.AddNewCut("EleJetCombMass>90", [&d](){ return d.evt.EleJetCombMass[d.ele.it] > 90; });
  // Mu
  sh.AddNewCut("GoodMu",        [&d](){ return d.mu.Pt[NLEP] > 45 && fabs(d.mu.Eta[NLEP]) < 2.1 && d.mu.IsTightMuon[NLEP] > 0; });
  sh.AddNewCut("MuJetCombMass>90",  [&d](){ return d.evt.MuJetCombMass[d.mu.it] > 90; });
  
  // Jets
  sh.AddNewCut("AK4Highest2Jet",   [&d](){ return d.jetsAK4.size>=2 && d.jetsAK4.it<2; });
  sh.AddNewCut("AK4Highest3Jet",   [&d](){ return d.jetsAK4.size>=3 && d.jetsAK4.it<3; });
  sh.AddNewCut("Highest2Jet",   [&d](){ return d.jetsAK8.size>=2 && d.jetsAK8.it<2; });
  sh.AddNewCut("Highest3Jet",   [&d](){ return d.jetsAK8.size>=3 && d.jetsAK8.it<3; });
  sh.AddNewCut("HadTop",           [&d](){ return d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75; });
  sh.AddNewCut("HadTopNoPtCut",    [&d](){ return d.jetsAK8.prunedMass[NJET] > 140 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75; });
  sh.AddNewCut("HadTopNoTauCut",   [&d](){ return d.jetsAK8.Pt[NJET] > 400 && d.jetsAK8.prunedMass[NJET] > 140; });
  sh.AddNewCut("HadTopNoMassCut",  [&d](){ return d.jetsAK8.Pt[NJET] > 400 && (d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET]) < 0.75; });
  sh.AddNewCut("HadTopNoMassNoTauCut", [&d](){ return d.jetsAK8.Pt[NJET] > 400; });
  sh.AddNewCut("JetHasMatchedGenTop",  [&d](){ return d.evt.JetHasMatchedGenTop[d.jetsAK8.it]; });
  sh.AddNewCut("GenLepTop",            [&d](){ return d.evt.JetMatchedGenTopType[d.jetsAK8.it]==1; });
  sh.AddNewCut("JetIsHadTopTagged",    [&d](){ return d.evt.JetIsHadTopTagged[d.jetsAK8.it]; });
  
  // Events
  sh.AddNewCut("NTop==2",          [&d](){ return d.evt.NTopHad+d.evt.NTopLep==2; });
  sh.AddNewCut("NTopHad<3",        [&d](){ return d.evt.NTopHad<3; });
  sh.AddNewCut("NTopHad==2",       [&d](){ return d.evt.NTopHad==2; });
  sh.AddNewCut("NTopLep==1",       [&d](){ return d.evt.NTopLep==1; });
  sh.AddNewCut("NLepTight==0",     [&d](){ return (d.evt.nmu+d.evt.neletight)==0; });
  sh.AddNewCut("NLepVeto==0",      [&d](){ return (d.evt.nmuveto+d.evt.neleveto)==0; });
  sh.AddNewCut("|DPhi|<2.8",       [&d](){ return fabs(d.evt.tt_dPhi)<2.8; });
  sh.AddNewCut("ttbar",            [&looper](){ return looper.it_sample==0; });
  sh.AddNewCut("ttbar,qcd",        [&looper](){ return looper.it_sample<2; });
  sh.AddNewCut("noqcd",            [&looper](){ return looper.it_sample!=1; });
  sh.AddNewCut("NLepTight==1",     [&d](){ return (d.evt.neletight+d.evt.nmu)==1; });
  sh.AddNewCut("NEleTight==1",     [&d](){ return d.evt.neletight==1; });
  sh.AddNewCut("NMuTight==1",      [&d](){ return d.evt.nmu==1; });
  sh.AddNewCut("NTopLikeHad>=2",   [&d](){ return d.evt.nhadtoplike>=2; });
  sh.AddNewCut("NGenLepFromTop==0",[&d](){ return d.evt.NGenLepFromTop==0; });
  
  // Gen Particles
  sh.AddNewCut("IsGenTop",         [&d](){ return d.evt.IsGenTop[d.gen.it]; });
  
  // --------------------------------------------------------------------
  //                         Postfixes
  //                  (Alternatives to cuts)
  // --------------------------------------------------------------------
  
  // Samples
  std::string Samples_PFs = samples.GetPFNames();
  std::string Samples_Lat = samples.GetLatexNames();
  std::string Samples_Cols = samples.GetColors();
  
  sh.AddNewPostfix("Directories",                 [&looper](){ return looper.it_sample; }, "[0to50]", "Sample [0to50]", "1-51");
  sh.AddNewPostfix("AllSamples",                  [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
  if (test) {
    sh.AddNewPostfix("Signals,TT,MultiJet", [](){ return 0; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Signals,Background",  [](){ return 0; }, Samples_PFs, Samples_Lat, Samples_Cols);
  } else {
    sh.AddNewPostfix("Signals,TT,MultiJet",         [&looper,&dir_to_index,nsignal,nttbar](){ return dir_to_index[looper.it_sample]<nsignal+nttbar ? dir_to_index[looper.it_sample] : nsignal+nttbar; }, 
		     std::string(Samples_PFs).replace(Samples_PFs.find("WJets"),5,"Multijet"), std::string(Samples_Lat).replace(Samples_Lat.find("W+Jets"),6,"Multijet"), Samples_Cols);
    sh.AddNewPostfix("Signals,Background",          [&looper,&dir_to_index,nsignal,nttbar](){ return dir_to_index[looper.it_sample]<nsignal+nttbar-1 ? dir_to_index[looper.it_sample] : nsignal+nttbar-1; }, 
		     std::string(Samples_PFs).replace(Samples_PFs.find("TTJets"),6,"SMBkg"), std::string(Samples_Lat).replace(Samples_Lat.find("t#bar{t}+jets"),13,"SM Background"), Samples_Cols);
  }
  const char* Samples = "AllSamples";
  
  // Jets
  sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetsAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("JetsPtOrdered",    [&d](){ return d.jetsAK8.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("NSubJet",          [&d](){ return (size_t)d.jetsAK8.nSubJets[NJET]; }, "NSubJet[0to4]", "N_{subjet}=[0to4]", "1-5");
  sh.AddNewPostfix("JetMassCut",       [&d](){ return (size_t)(d.jetsAK8.prunedMass[NJET] > 140); }, "MassBelow140;MassAbove140", "M_{pruned} < 140;M_{pruned} > 140", "2,3");
  sh.AddNewPostfix("JetTau32Cut",      [&d](){ return (size_t)(d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET] < 0.75); }, "Tau32Above0p75;Tau32Below0p75", "#tau_{3}/#tau_{2} > 0.75;#tau_{3}/#tau_{2} < 0.75", "2,3");
  sh.AddNewPostfix("JetPtCut",         [&d](){ return (size_t)(d.jetsAK8.Pt[NJET] > 400); }, "PtBelow400;PtAbove400", "p_{T} < 400;p_{T} > 400", "2,3");
  
  sh.AddNewPostfix("JetGenTruth",              [&d](){ return (size_t)(d.evt.JetGenTruth[d.jetsAK8.it]); }, "NoTopInEvt;NoTopMatch;NonMergedTopJet;MergedHadTop;MergedLepTop;BadWMatch", "No top in event;No top matched;Matched non-merged top;Merged (R<0.7) top - hadronic;Merged (R<0.7) top - leptonic;Bad W Match (unknown)", "1,13,2,3,4,5");
  sh.AddNewPostfix("JetMatchedGenTopType",     [&d](){ return (size_t)(d.evt.JetMatchedGenTopType[d.jetsAK8.it]!=-9999 ? d.evt.JetMatchedGenTopType[d.jetsAK8.it] : -1); }, "MatchedGenTopHad;MatchedGenTopLep", "Hadronic top;Semi-leptonic top", "2,4");
  sh.AddNewPostfix("IsMerged",     [&d](){ return (size_t)d.evt.JetMatchedGenTopIsMerged[d.jetsAK8.it]; }, "NonMergedTop;MergedTop", "Non-merged top;Merged top", "2,4");
  
  // Event
  sh.AddNewPostfix("R0p25",           [&d](){ return d.evt.R > 0.25; },  "RBelow0p25;RAbove0p25", "R < 0.25;R > 0.25", "4,2");
  sh.AddNewPostfix("R0p3",            [&d](){ return d.evt.R > 0.3; },   "RBelow0p3;RAbove0p3", "R < 0.3;R > 0.3", "4,2");
  sh.AddNewPostfix("R0p32",           [&d](){ return d.evt.R > 0.32; },  "RBelow0p32;RAbove0p32", "R < 0.32;R > 0.32", "4,2"); // Best cut
  sh.AddNewPostfix("R0p35",           [&d](){ return d.evt.R > 0.35; },  "RBelow0p35;RAbove0p35", "R < 0.35;R > 0.35", "4,2");
  sh.AddNewPostfix("CutR",            [&d](){ return d.evt.R > 0.4; },   "RBelow0p4;RAbove0p4", "R < 0.4;R > 0.4", "4,2");
  sh.AddNewPostfix("CutDPhi",         [&d](){ return fabs(d.evt.tt_dPhi) > 2.8; },  "DPhiBelow2p8;DPhiAbove2p8", "#Delta#phi_{t#bar{t}} < 2.8;#Delta#phi_{t#bar{t}} > 2.8", "2,4");
  sh.AddNewPostfix("CutHtAll",       [&d](){ return d.evt.HTall > 1200; },  "HtAllBelow1200;HtAllAbove1200", "H_{T,all} < 1200;H_{T,all} > 1200", "4,2"); // Best cut
  sh.AddNewPostfix("NTopHad",         [&d](){ return (d.evt.NTopHad>1)+(d.evt.NTopHad>2); }, "0To1HadTop;2HadTop;3HadTop", "N_{top,hadronic}<2;N_{top,hadronic}=2;N_{top,hadronic}=3", "4,2,3");
  sh.AddNewPostfix("HtAll1450",       [&d](){ return d.evt.HTall > 1450; },  "HtAllBelow1450;HtAllAbove1450", "H_{T,all} < 1450;H_{T,all} > 1450", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1500",       [&d](){ return d.evt.HTall > 1500; },  "HtAllBelow1500;HtAllAbove1500", "H_{T,all} < 1500;H_{T,all} > 1500", "4,2");
  sh.AddNewPostfix("NTopHad",         [&d](){ return (d.evt.NTopHad>1)+(d.evt.NTopHad>2); }, "0To1HadTop;2HadTop", "N_{top,hadronic}<2;N_{top,hadronic}=2", "4,2");
  sh.AddNewPostfix("NGenLepFromTop",   [&d](){ return (size_t)d.evt.NGenLepFromTop; }, "FullHad;[1to4]LepTop", "[0to4]l (e/#mu, from top)", "1-5");
  
  // Gen Particles
  sh.AddNewPostfix("GenTopType",      [&d](){ return (size_t)(d.evt.GenTopType[d.gen.it]!=-9999 ? d.evt.GenTopType[d.gen.it] : -1); }, "GenTopHad;GenTopLep", "Hadronic top;Semi-leptonic top", "2,4");
  
  // --------------------------------------------------------------------
  //                         Fill Parameters
  // --------------------------------------------------------------------
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:

  // Samples
  sh.AddNewFillParam("Sample",            { .nbin= 50,  .bins={   0,     50}, .fill=[&looper](){ return looper.it_sample;     }, .axis_title="iSample"});
  
  // Muons
  sh.AddNewFillParam("MuEnergy",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.E[NLEP];              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuPt",              { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.Pt[NLEP];             }, .axis_title="Muon p_{T} (GeV/c)"});
  sh.AddNewFillParam("MuDRJet",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.MuDRJet[d.mu.it];    }, .axis_title="#DeltaR (#mu, jet)"});
  sh.AddNewFillParam("MuRelPtJet",        { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.MuRelPtJet[d.mu.it]; }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV/c)"});
  sh.AddNewFillParam("MuJetCombMass",     { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.MuJetCombMass[d.mu.it]; }, .axis_title="Mass_{#mu+jet comb.} (GeV/c^{2})"});
  
  // Electrons
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.E[NLEP];              }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.Pt[NLEP];             }, .axis_title="Electron p_{T} (GeV/c)"});
  sh.AddNewFillParam("EleDRJet",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.EleDRJet[d.ele.it];   }, .axis_title="#DeltaR (e, jet)"});
  sh.AddNewFillParam("EleRelPtJet",       { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.EleRelPtJet[d.ele.it];}, .axis_title="p_{T}^{rel} (e, jet) (GeV/c)"});
  sh.AddNewFillParam("EleJetCombMass",    { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.EleJetCombMass[d.ele.it]; }, .axis_title="Mass_{e+jet comb.} (GeV/c^{2})"});
  
  // MET
  sh.AddNewFillParam("MetPt",             { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.met.Pt;                   }, .axis_title="MET p_{T} (GeV/c)"});
  
  // AK4 Jets
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.E[NJET];          }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Pt[NJET];         }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Mass[NJET];       }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  
  // Jets (AK8)
  sh.AddNewFillParam("JetEnergy",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.E[NJET];                    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("JetPt",              { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.Pt[NJET];                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtCoarse",        { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.Pt[NJET];                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtBins",          { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.jetsAK8.Pt[NJET]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtFewBins",       { .nbin=   5, .bins={0, 300, 400, 600, 1000, 2000}, .fill=[&d](){ return d.jetsAK8.Pt[NJET]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtOneBin",        { .nbin=   1, .bins={400, 5000}, .fill=[&d](){ return d.jetsAK8.Pt[NJET]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("JetMass",            { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.Mass[NJET];                 }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMass",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMassCoarse",{ .nbin= 200, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.prunedMass[NJET];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetFilteredMass",    { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.filteredMass[NJET];         }, .axis_title="AK8-jet Filtered Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTrimmedMass",     { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.trimmedMass[NJET];          }, .axis_title="AK8-jet Trimmed Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTopMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.topMass[NJET];              }, .axis_title="AK8-jet Top Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetMinMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.minmass[NJET];              }, .axis_title="AK8-jet Min. Subjet-pair Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetNSubJets",        { .nbin=  11, .bins={-0.5,   10.5}, .fill=[&d](){ return d.jetsAK8.nSubJets[NJET];             }, .axis_title="AK8-jet N_{subjet}"});
  sh.AddNewFillParam("JetTau1",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau1[NJET];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("JetTau2",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau2[NJET];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("JetTau3",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[NJET];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("JetTau21",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau2[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetTau31",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau1[NJET];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau32",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[NJET]/d.jetsAK8.tau2[NJET];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("JetDRLep",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.DRJetLep[d.jetsAK8.it];         }, .axis_title="#DeltaR (lepton, jet)"});
  sh.AddNewFillParam("JetRelPtLep",        { .nbin= 100, .bins={   0,    500}, .fill=[&d](){ return d.evt.RelPtJetLep[d.jetsAK8.it];      }, .axis_title="p_{T}^{rel} (lepton, jet) [GeV/c]"});
  
  sh.AddNewFillParam("JetMatchedGenTopPt",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.evt.JetMatchedGenTopPt[d.jetsAK8.it];         }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopPtCoarse",  { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.evt.JetMatchedGenTopPt[d.jetsAK8.it];         }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopPtBins",    { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.evt.JetMatchedGenTopPt[d.jetsAK8.it]; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDR",     { .nbin=  50, .bins={   0,      6}, .fill=[&d](){ return d.evt.JetMatchedGenTopJetDR[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDRFine", { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.JetMatchedGenTopJetDR[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("GenBJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenBJetDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenBJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenBJetDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenWJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWJetDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenWJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWJetDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenWGenBDR",         { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWGenBDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, Gen. b)"});
  sh.AddNewFillParam("GenWGenBDRFine",     { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWGenBDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, Gen. b)"});
  sh.AddNewFillParam("GenLepJetDR",        { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepJetDR[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. lep, jet)"});
  sh.AddNewFillParam("GenLepJetDRFine",    { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepJetDR[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. lep, jet)"});
  sh.AddNewFillParam("GenLepGenBDR",       { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepGenBDR[d.jetsAK8.it];     }, .axis_title="#DeltaR (Gen. lep, Gen. b)"});
  sh.AddNewFillParam("GenLepGenBDRFine",   { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepGenBDR[d.jetsAK8.it];     }, .axis_title="#DeltaR (Gen. lep, Gen. b)"});
  
  sh.AddNewFillParam("MaxSubJetCSV",       { .nbin=  9, .bins={ 0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 }, .fill=[&d](){ return d.evt.maxSubjetCSV[d.jetsAK8.it]; }, .axis_title="Maximum Subjet CSV"});
  
  // Special Y/Z axis parameters:
  // You can define how to calculate them in SmartHistos
  // Rates/Efficiencies/Averages etc can be easily calculated from 2D plots (eg. Profile)
  // SmartHistos::init_() for extra histo name and axis title
  // SmartHistos::calc_spec_[1,2]d() for the method
  sh.AddNewFillParam("MergedTopFraction",  { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.JetMatchedGenTopIsMerged[d.jetsAK8.it]; }, .axis_title="Fraction of Merged Tops" });   
  sh.AddNewFillParam("TopTagEfficiency",   { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.JetIsHadTopTagged[d.jetsAK8.it]; }, .axis_title="Top-tagging Efficiency" });
  sh.AddNewFillParam("MisTagRate",         { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return (size_t)(!d.evt.JetHasMatchedGenTop[d.jetsAK8.it]); }, .axis_title="Mis-tag Rate" });
  
  // Event variables
  sh.AddNewFillParam("NJet",             { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.jetsAK8.size;                  }, .axis_title="N_{AK8-jet}"});
  sh.AddNewFillParam("NJetSelected",     { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopHad+d.evt.NTopLep;     }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetHadronic",     { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopHad;                   }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetLeptonic",     { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{Leptonic AK8-jet}"});
  sh.AddNewFillParam("NLep",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NLep;                      }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NLepTight",        { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight+d.evt.nmu;       }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMu",              { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmu;                       }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEle",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nele;                      }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NEleTight",        { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight;                 }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NLepVeto",         { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto+d.evt.nmuveto;    }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMuVeto",          { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmuveto;                   }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEleVeto",         { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto;                  }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NTopLikeHad",      { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nhadtoplike;               }, .axis_title="N_{jet} (p_{T}>400, M>100)"});
  sh.AddNewFillParam("NTopLep",          { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{top, leptonic}"});
  sh.AddNewFillParam("NTopHad",          { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopHad;                   }, .axis_title="N_{top, hadronic}"});
  sh.AddNewFillParam("TTHadMR",          { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMRFine",      { .nbin= 200, .bins={   0,   10000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMTR",         { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMTR;                  }, .axis_title="M_{T,t#bar{t}}^{R} (GeV/c)"});
  sh.AddNewFillParam("DPhi",             { .nbin=  16, .bins={   0,     3.2}, .fill=[&d](){ return fabs(d.evt.tt_dPhi);             }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhiFine",         { .nbin=  64, .bins={   0,     3.2}, .fill=[&d](){ return fabs(d.evt.tt_dPhi);             }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhiBins",         { .nbin=  11, .bins={ 0, 0.5, 1.0, 1.5, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.15 }, .fill=[&d](){ return fabs(d.evt.tt_dPhi); }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("TTHadDEta",        { .nbin=  50, .bins={   0,       5}, .fill=[&d](){ return d.evt.TTHadDEta;                 }, .axis_title="#Delta#eta_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadDR",          { .nbin=  60, .bins={   0,       6}, .fill=[&d](){ return d.evt.TTHadDR;                   }, .axis_title="#DeltaR_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadPz",          { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadPz;                   }, .axis_title="P_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadDPz",         { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadDPz;                  }, .axis_title="#DeltaP_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadHz",          { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadHz;                   }, .axis_title="H_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMass",        { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMass;                 }, .axis_title="M_{t#bar{t}} (GeV/c^{2})"});
  sh.AddNewFillParam("TTHadR",           { .nbin=  24, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadRFine",       { .nbin= 240, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadR2",          { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.TTHadR2;                   }, .axis_title="R_{t#bar{t}}^{2}"});
  sh.AddNewFillParam("R",                { .nbin=  24, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("RFine",            { .nbin= 240, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("RBins",            { .nbin=  14, .bins={ 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2 }, .fill=[&d](){ return d.evt.R; }, .axis_title="R"});
  sh.AddNewFillParam("R2",               { .nbin=  30, .bins={   0,     1.5}, .fill=[&d](){ return d.evt.R2;                        }, .axis_title="R^{2}"});
  sh.AddNewFillParam("MR",               { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.MR;                        }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("MRFine",           { .nbin= 200, .bins={   0,   10000}, .fill=[&d](){ return d.evt.MR;                        }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("MTR",              { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.MTR;                       }, .axis_title="M_{T}^{R} (GeV/c)"});
  sh.AddNewFillParam("HtTopFraction",    { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtTopFr;                   }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("HtExFraction",     { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtExFr;                    }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("Ht",               { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAllCoarse",      { .nbin=  20, .bins={   0,    6000}, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAll",            { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAllFine",        { .nbin= 200, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HTall;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtTop",            { .nbin=  25, .bins={   0,    5000}, .fill=[&d](){ return d.evt.HtTop;                     }, .axis_title="H_{T,tops} (GeV/c)"});
  sh.AddNewFillParam("HtEx",             { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtEx;                      }, .axis_title="H_{T,extra} (GeV/c)"});
  
  // Gen Particles
  sh.AddNewFillParam("GenTopPtBins",      { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.evt.IsGenTop[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("GenTopPtFewBins",   { .nbin=  4, .bins={0, 300, 400, 800, 2000}, .fill=[&d](){ return d.evt.IsGenTop[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  
  // Special Y/Z axis parameters:
  sh.AddNewFillParam("JetFindingEfficiency",  { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.GenTopHasMatchedJet[d.gen.it]; }, .axis_title="Jet finding Efficiency" });
  sh.AddNewFillParam("TopFindingEfficiency",  { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.GenTopHasMatchedTopTagJet[d.gen.it]; }, .axis_title="Top finding Efficiency" });
  
  // --------------------------------------------------------------------
  //                       Histogram Definitions
  // --------------------------------------------------------------------
  
  // MET Plots
  sh.AddHistos("met",    { .fill="MetPt",            .pfs={Samples}, .cuts={}, .draw="HISTE1", .opt="", .ranges={0,0, 0,0} });
  
  // Events
  sh.AddHistos("evt",   { .fill="RFine",     .pfs={"NTopHad","CutDPhi",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="RFine",     .pfs={"CutDPhi","NTopHad",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"NTopHad","CutDPhi",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"CutDPhi","NTopHad",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhiFine",  .pfs={"NTopHad","CutR",   Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhiFine",  .pfs={"CutR",   "NTopHad",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"NTopHad","CutR",   Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"CutR",   "NTopHad",Samples}, .cuts={"NTopHad<3","NTopLikeHad>=2"/*,lepveto*/}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 0.001,100000} });
  
  sh.AddHistos("evt",   { .fill="TTHadMR", .pfs={Samples},           .cuts={}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,1000000} });
  sh.AddHistos("evt",   { .fill="TTHadMR", .pfs={Samples,"CutR"},      .cuts={}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,1000000} });
  sh.AddHistos("evt",   { .fill="TTHadMR", .pfs={Samples,"CutDPhi"},      .cuts={}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,1000000} });
  sh.AddHistos("evt",   { .fill="TTHadMR", .pfs={Samples,"CutR","CutDPhi"}, .cuts={}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,1000000} });
  
  // Distributions and N-1 plots
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples},           .cuts={"NTopHad==2"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhi",            .pfs={Samples},           .cuts={"NTopHad==2"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="R",               .pfs={Samples,"CutDPhi"}, .cuts={"NTopHad==2"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,100000} });
  sh.AddHistos("evt",   { .fill="DPhi",            .pfs={Samples,"CutR"},    .cuts={"NTopHad==2"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 0.001,100000} });
  
  // Jets
  
  // 1D distributions
  // No Cut
  sh.AddHistos("jetsAK8", { .fill="JetPtCoarse",                           .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau1",                               .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau2",                               .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau3",                               .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau21",                              .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31",                              .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass ",                       .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetNSubJets",                           .pfs={"Signals,TT,MultiJet"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // Apply 2 Cuts (N-1)
  sh.AddHistos("jetsAK8", { .fill="JetPtCoarse",                           .pfs={"Signals,TT,MultiJet","JetMassCut","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,MultiJet","JetPtCut","JetMassCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // Apply 3 Cuts
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,MultiJet","JetTau32Cut","JetPtCut","JetMassCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  
  // 2D Plots
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32",       .pfs={"Signals,TT,MultiJet","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  
  // Jets with GenTruth information (N-1 Cuts)
  sh.AddHistos("jetsAK8", { .fill="JetPtCoarse",                           .pfs={"JetGenTruth","Signals,Background","JetMassCut","JetTau32Cut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // 2 Cuts
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},  .cuts={"NoDoubleTT"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},  .cuts={"NoDoubleTT"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},  .cuts={"NoDoubleTT"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut","JetMassCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} }); // 3 Cuts
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut","JetMassCut"}, .cuts={"NoDoubleTT"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  
  // Jets with Matched GenTop (or constituents) are found in cone
  // --> Switch to Top to Jet matching (calc Efficiency)
  sh.AddHistos("jetsAK8", { .fill="JetPtCoarse_vs_JetMatchedGenTopPtCoarse",      .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  sh.AddHistos("jetsAK8", { .fill="JetPtCoarse_vs_JetMatchedGenTopPtCoarse",      .pfs={"Signals,Background","IsMerged","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"Signals,Background","IsMerged","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"IsMerged","JetMatchedGenTopType","Signals,Background"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"JetMatchedGenTopType","IsMerged","Signals,Background"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  
  // Distances - Merged/Type
  sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR",                                  .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenBJetDR",                                              .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenWJetDR",                                              .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenWGenBDR",                                             .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenLepJetDR",                                            .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR",                                           .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR_vs_JetMatchedGenTopPtBins",        .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenBJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenWJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenWGenBDR_vs_JetMatchedGenTopPtBins",                   .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, } });
  sh.AddHistos("jetsAK8", { .fill="GenLepJetDR_vs_JetMatchedGenTopPtBins",                  .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR_vs_JetMatchedGenTopPtBins",                 .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgJetMatchedGenTopJetDRFine_vs_JetMatchedGenTopPtBins", .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgGenBJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgGenWJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgGenWGenBDRFine_vs_JetMatchedGenTopPtBins",            .pfs={"IsMerged","Signals,Background","JetMatchedGenTopType"}, .cuts={"NoDoubleTT","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgGenLepJetDRFine_vs_JetMatchedGenTopPtBins",           .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="AvgGenLepGenBDRFine_vs_JetMatchedGenTopPtBins",          .pfs={"IsMerged","Signals,Background"},                        .cuts={"NoDoubleTT","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  
  // Fraction of jets that have all constituents merged in the cone
  // Since it's gen particles we consider tops to be merged if both
  // b and W is found R<0.7 from the Jet
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins",  .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType"},            .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetPtBins",               .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType"},            .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  
  // Top tagging efficiency
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins",   .pfs={"Signals,TT,MultiJet"},                                   .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins",   .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType"},            .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins",   .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType","IsMerged"}, .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetPtBins",                .pfs={"Signals,TT,MultiJet"},                                   .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetPtBins",                .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType"},            .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetPtBins",                .pfs={"Signals,TT,MultiJet","JetMatchedGenTopType","IsMerged"}, .cuts={"NoMultiJet","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });  
  
  // Probability of not having a gen top within the jet cone
  // for tagged hadronic top jets
  sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtBins",     .pfs={Samples}, .cuts={"TTJets","JetIsHadTopTagged"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtOneBin",   .pfs={Samples}, .cuts={"TTJets","JetIsHadTopTagged"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",  .pfs={Samples}, .cuts={"TTJets","JetIsHadTopTagged"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  // GenParticles
  
  // Probability to find a jet close to the top
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT,MultiJet"},              .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT,MultiJet"},              .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT,MultiJet","GenTopType"}, .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT,MultiJet","GenTopType"}, .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  // Full probability to find a jet that is also top tagged
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT,MultiJet"},              .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT,MultiJet"},              .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT,MultiJet","GenTopType"}, .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT,MultiJet","GenTopType"}, .cuts={"NoMultiJet","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  
  
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
          while(d.gen.Loop()) sh.Fill("gen"); if (debug) cout<<"Fill Gen ok\n";
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

