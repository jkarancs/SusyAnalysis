#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"
#include "../interface/Razor.h"

#define NTHFILE 1

void calcVariables(Data& d) {
  
  // find good leptons (for letponic tops)
  int ngoodleptons = 0;
  std::vector<TLorentzVector> goodleps;
  d.evt.HTlep = 0;
  while(d.ele.Loop()) if (d.ele.Pt > 35 && fabs(d.ele.Eta) < 2.5) {
  //while(d.ele.Loop()) if (d.ele.isLoose>0 && d.ele.Pt > 35 && fabs(d.ele.Eta) < 2.5) {
    ngoodleptons++;
    TLorentzVector goodele;
    goodele.SetPtEtaPhiE(d.ele.Pt, d.ele.Eta, d.ele.Phi, d.ele.E);
    goodleps.push_back(goodele);
    d.evt.HTlep += d.ele.Pt;
  }
  while(d.mu.Loop())  if (d.mu.IsTightMuon>0 && d.mu.Pt > 45 && fabs(d.mu.Eta) < 2.1) {
  //while(d.mu.Loop())  if (d.mu.IsLooseMuon>0 && d.mu.Pt > 45 && fabs(d.mu.Eta) < 2.1) {
    ngoodleptons++;
    TLorentzVector goodmu;
    goodmu.SetPtEtaPhiE(d.mu.Pt, d.mu.Eta, d.mu.Phi, d.mu.E);
    goodleps.push_back(goodmu);
    d.evt.HTlep += d.mu.Pt;
  }
  
  // Tag hadronic tops
  TLorentzVector hadtop1;
  TLorentzVector hadtop2;
  TLorentzVector leptop1;
  TLorentzVector leptop2;
  d.evt.ntops = 0;
  d.evt.nhadtops = 0;
  d.evt.nleptops = 0;
  d.evt.HT = 0;
  while(d.jetAK8.Loop()) {
    bool is_top = false;
    // hadronic tops
    if (d.jetAK8.tau1>0 && d.jetAK8.tau2>0 ? d.jetAK8.Pt > 400 && d.jetAK8.prunedMass > 140 && (d.jetAK8.tau2/d.jetAK8.tau1) < 0.75 : 0) {
      ++d.evt.nhadtops;
      is_top = true;
      if (d.evt.nhadtops==1) hadtop1.SetPtEtaPhiE(d.jetAK8.Pt, d.jetAK8.Eta, d.jetAK8.Phi, d.jetAK8.E);
      if (d.evt.nhadtops==2) hadtop2.SetPtEtaPhiE(d.jetAK8.Pt, d.jetAK8.Eta, d.jetAK8.Phi, d.jetAK8.E);
    }
    // leptonic tops
    TLorentzVector lep;
    TLorentzVector jet;
    jet.SetPtEtaPhiE(d.jetAK8.Pt, d.jetAK8.Eta, d.jetAK8.Phi, d.jetAK8.E);
    double DeltaR_lep = 9999;
    for (size_t i=0; i<goodleps.size(); ++i) {
      if (goodleps[i].DeltaR(jet)< DeltaR_lep) {
	DeltaR_lep = goodleps[i].DeltaR(jet);
	lep = goodleps[i];
      }
    }
    if (DeltaR_lep<1.0) {
      d.evt.nleptops++;
      is_top = true;
      TLorentzVector leptop = lep + jet;
      if (d.evt.nleptops==1) leptop1 = leptop;
      if (d.evt.nleptops==2) leptop2 = leptop;
    }
    // Extra - all except above hadronic/leptonic tops
    d.evt.HT += d.jetAK8.Pt;
    if (is_top) d.evt.HTtt += d.jetAK8.Pt;
  }
  d.evt.HTall = d.evt.HT + d.met.metPt[0] + d.evt.HTlep;

  // Select exactly 2 tops (hadronic or leptonic)
  // We need exactly 2 in order to calculate pair variables, eg. DeltaPhi
  TLorentzVector top1;
  TLorentzVector top2;
  d.evt.ntops = d.evt.nhadtops + d.evt.nleptops;
  if (d.evt.ntops == 2) {
    if (d.evt.nhadtops == 2) {
      if (hadtop1.Pt() > hadtop2.Pt()) {
        top1 = hadtop1;
        top2 = hadtop2;
      } else {
        top1 = hadtop2;
        top2 = hadtop1;
      }
    } else if (d.evt.nhadtops == 1) {
      if (hadtop1.Pt() > leptop1.Pt()) {
        top1 = hadtop1;
        top2 = leptop1;
      } else {
        top1 = leptop1;
        top2 = hadtop1;
      }    
    } else if (d.evt.nhadtops == 0) {
      if (leptop1.Pt() > leptop2.Pt()) {
        top1 = leptop1;
        top2 = leptop2;
      } else {
        top1 = leptop2;
        top2 = leptop1;
      }
    }
  }
  
  // top pair variables
  d.evt.tt_dR=NOVAL_F;
  d.evt.tt_dPhi=NOVAL_F;
  d.evt.tt_dEta=NOVAL_F;
  d.evt.tt_Mass=NOVAL_F;
  d.evt.tt_MR=NOVAL_F;
  d.evt.tt_MTR=NOVAL_F;
  d.evt.tt_R=NOVAL_F;
  d.evt.tt_R2=NOVAL_F;
  d.evt.HTtt=NOVAL_F;
  d.evt.HTex=NOVAL_F;
  d.evt.HTttFraction=NOVAL_F;
  d.evt.HTexFraction=NOVAL_F;
  if (d.evt.ntops==2) {
    /* python
       tt_dR[0] = top1.DeltaR(top2)
       tt_dPhi[0] = top1.DeltaPhi(top2)
       tt_dEta[0] = fabs(top1.Eta() - top2.Eta())
       tt_mtt[0] = (top1 + top2).M()
       
       tt_extra = top1 + top2
       for jet in jets:
         tt_extra += jet
       pz_tt_extra[0] = tt_extra.Pz()
       dHt[0] = math.fabs(top1.Pt() - top2.Pt())
       dphi1 = delta_phi(metphi[0], top1.Phi())
       dphi2 = delta_phi(metphi[0], top2.Phi())
       dPhi_met_t1[0] = max(dphi1, dphi2)
       dPhi_met_t2[0] = min(dphi1, dphi2)
    */
    d.evt.tt_dR = top1.DeltaR(top2);
    d.evt.tt_dPhi = top1.DeltaPhi(top2);
    d.evt.tt_dEta = fabs(top1.Eta() - top2.Eta());
    d.evt.tt_Mass = (top1 + top2).M();
    d.evt.tt_Pz = (top1 + top2).Pz();
    d.evt.tt_Hz = top1.Pz() + top2.Pz();
    d.evt.tt_dPz = fabs(top1.Pz() - top2.Pz());
    d.evt.HTtt = top1.Pt() + top2.Pt();
    d.evt.HTex = d.evt.HTall - d.evt.HTtt;
    d.evt.HTttFraction = d.evt.HTtt / d.evt.HTall;
    d.evt.HTexFraction = d.evt.HTex / d.evt.HTall;
    // Select signal region

    // Razor for hadronic top pair
    TVector3 metl;
    metl.SetPtEtaPhi(d.met.metPt[0], 0, d.met.metPhi[0]);
    d.evt.tt_MR = CalcMR(top1, top2);
    d.evt.tt_MTR = CalcMTR(top1, top2, metl);
    d.evt.tt_R = d.evt.tt_MTR / d.evt.tt_MR;
    d.evt.tt_R2 = pow(d.evt.tt_R, 2);
  }
}

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
  
  // Initialize TreeReader
  B2GTreeReader reader;
  
  // Class to Loop on files and read the Trees
  B2GTreeLooper looper(NTHFILE,1);
  
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
  sh.AddNewPostfix("ttbar,qcd",                   [&looper](){ return looper.it_sample; }, "ttbar;qcd", "t#bar{t};QCD", "2,6");
  sh.AddNewPostfix("ttbar,qcd,Susy3,Susy4",       [&looper](){ return looper.it_sample; }, "ttbar;qcd;susy3body;susy4body", "t#bar{t};QCD;T5tttt - 3body;T5tttt - 4body", "2,6,4,3");
  sh.AddNewPostfix("SideBand,Signal",             [&d](){ return d.evt.HTall > 1500; }, "SideBand;Signal", "H_{T,all} < 1500 GeV/c;H_{T,all} > 1500 GeV/c", "4,2");
  sh.AddNewPostfix("RBelow0p25,RAbove0p25",       [&d](){ return d.jetAK8.R > 0.25; },  "RBelow0p25;RAbove0p25", "R < 0.25;R > 0.25", "4,2");
  sh.AddNewPostfix("DPhiBelow0p28,DPhiAbove0p28", [&d](){ return fabs(d.evt.tt_dPhi) > 0.28; },  "DPhiBelow0p28;DPhiAbove0p28", "#Delta#phi_{t#bar{t}} < 0.28;#Delta#phi_{t#bar{t}} > 0.28", "4,2");
  
  //sh.AddNewPostfix("ttbar,Susy3,Susy4", &looper.it_sample, "ttbar;susy3body;susy4body", "SM t#bar{t};T5tttt - 3body;T5tttt - 4body", "2,4,3");
  sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("AK8JetsPtOrdered", [&d](){ return d.jetAK8.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:
  sh.AddNewFillParam("MuonEnergy",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.E;              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuonPt",            { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.mu.Pt;             }, .axis_title="Muon p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.E;             }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.ele.Pt;            }, .axis_title="Electron p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("MetPt",             { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.met.Pt;                      }, .axis_title="MET p_{T} (GeV/c)"});
  
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK4.E;                    }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK4.Pt;                   }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK4.Mass;                 }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  
  sh.AddNewFillParam("AK8JetEnergy",      { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK8.E;                    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("AK8JetPt",          { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK8.Pt;                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMass",        { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK8.Mass;                 }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetPrunedMass",  { .nbin= 400, .low=   0,   .high=2000, .fill=[&d](){ return d.jetAK8.prunedMass;           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("AK8JetTau1",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau1;                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau2",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau2;                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("AK8JetTau3",        { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3;                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("AK8JetTau21",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau2/d.jetAK8.tau1;   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau31",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3/d.jetAK8.tau1;   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("AK8JetTau32",       { .nbin= 100, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.tau3/d.jetAK8.tau2;   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("AK8JetMR",          { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetAK8.MR;                   }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("AK8JetMTR",         { .nbin=  50, .low=   0,   .high=5000, .fill=[&d](){ return d.jetAK8.MTR;                  }, .axis_title="M_{T}^{R} (GeV/c)"});

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
  sh.AddNewFillParam("AK8JetR",           { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.R;                    }, .axis_title="R"});
  sh.AddNewFillParam("AK8JetR2",          { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.R2;                   }, .axis_title="R^{2}"});
  sh.AddNewFillParam("H_TttFraction",     { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTttFraction;            }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("H_TexFraction",     { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTexFraction;            }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("H_T",               { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HT;                      }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("H_Tevt",            { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTall;                   }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("H_Ttt",             { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTtt;                    }, .axis_title="H_{T,tops} (GeV/c)"});
  sh.AddNewFillParam("H_Tex",             { .nbin=  25, .low=   0,   .high=5000, .fill=[&d](){ return d.evt.HTex;                    }, .axis_title="H_{T,extra} (GeV/c)"});
  
  // 2D binning of Signal cut variables
  sh.AddNewFillParam("DeltaPhi",          { .nbin=  16, .low=   0,   .high= 3.2, .fill=[&d](){ return fabs(d.evt.tt_dPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("Rtt",               { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.tt_R;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("R",                 { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.jetAK8.R;                    }, .axis_title="R"});
  sh.AddNewFillParam("HT",                { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HT;                      }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("HTall",             { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HTall;                   }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HTex",              { .nbin=  20, .low=   0,   .high=6000, .fill=[&d](){ return d.evt.HTex;                    }, .axis_title="H_{T,extra} (GeV/c)"});
  sh.AddNewFillParam("HTttFraction",      { .nbin=  20, .low=   0,   .high=   1, .fill=[&d](){ return d.evt.HTttFraction;            }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  
  // Define Cuts here:
  sh.AddNewCut("AK4Highest2Jet",   [&d](){ return d.jetAK4.jets_size>=2 && d.jetAK4.it<2; });
  sh.AddNewCut("AK4Highest3Jet",   [&d](){ return d.jetAK4.jets_size>=3 && d.jetAK4.it<3; });
  sh.AddNewCut("AK8Highest2Jet",   [&d](){ return d.jetAK8.jetsAK8_size>=2 && d.jetAK8.it<2; });
  sh.AddNewCut("AK8Highest3Jet",   [&d](){ return d.jetAK8.jetsAK8_size>=3 && d.jetAK8.it<3; });
  
  sh.AddNewCut("HadTop",           [&d](){ return d.jetAK8.tau1>0 && d.jetAK8.tau2>0 ? d.jetAK8.Pt > 400 && d.jetAK8.prunedMass > 140 && (d.jetAK8.tau2/d.jetAK8.tau1) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoPtCut",    [&d](){ return d.jetAK8.tau1>0 && d.jetAK8.tau2>0 ? d.jetAK8.prunedMass > 140 && (d.jetAK8.tau2/d.jetAK8.tau1) < 0.75 : 0; });
  sh.AddNewCut("HadTopNoTauCut",   [&d](){ return d.jetAK8.Pt > 400 && d.jetAK8.prunedMass > 140; });
  sh.AddNewCut("HadTopNoMassCut",  [&d](){ return d.jetAK8.tau1>0 && d.jetAK8.tau2>0 ? d.jetAK8.Pt > 400 && (d.jetAK8.tau2/d.jetAK8.tau1) < 0.75 : 0; });

  sh.AddNewCut("NHadTop==2",       [&d](){ return d.evt.nhadtops==2; });
  sh.AddNewCut("|DPhi|<2.8",       [&d](){ return fabs(d.evt.tt_dPhi)<2.8; });
  sh.AddNewCut("NTop==2",          [&d](){ return d.evt.nhadtops+d.evt.nleptops==2; });
  sh.AddNewCut("ttbar",            [&looper](){ return looper.it_sample==0; });
  sh.AddNewCut("ttbar,qcd",        [&looper](){ return looper.it_sample<2; });
  sh.AddNewCut("noqcd",            [&looper](){ return looper.it_sample!=1; });
  
  // Set Histogram weight (empty = 1)
  sh.SetHistoWeights({});
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  sh.AddHistos("mu",     { .fill="MuonEnergy",       .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("mu",     { .fill="MuonPt",           .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="EleEnergy",        .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  sh.AddHistos("ele",    { .fill="ElePt",            .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="", .ranges={0,1000, 0,0} });
  sh.AddHistos("met",    { .fill="MetPt",            .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0} });
  
  sh.AddHistos("jetAK4", { .fill="AK4JetEnergy",     .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK4", { .fill="AK4JetPt",         .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK4", { .fill="AK4JetMass",       .pfs={"AK4JetsPtOrdered"}, .cuts={"ttbar","AK4Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,500, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetEnergy",     .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });

  sh.AddHistos("jetAK8", { .fill="AK8JetPt",         .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetMass",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,250, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau1",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau2",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau3",       .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest2Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetAK8", { .fill="AK8JetTau21",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetAK8", { .fill="AK8JetTau31",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("jetAK8", { .fill="AK8JetTau32",      .pfs={"AK8JetsPtOrdered"}, .cuts={"ttbar","AK8Highest3Jet"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });

  // Top Tagging
  //  N-1 plots
  sh.AddHistos("jetAK8", { .fill="AK8JetPrunedMass", .pfs={"ttbar,qcd"}, .cuts={"ttbar,qcd","HadTopNoMassCut"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetPt",         .pfs={"ttbar,qcd"}, .cuts={"ttbar,qcd","HadTopNoPtCut"},   .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau21",      .pfs={"ttbar,qcd"}, .cuts={"ttbar,qcd","HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau31",      .pfs={"ttbar,qcd"}, .cuts={"ttbar,qcd","HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetAK8", { .fill="AK8JetTau32",      .pfs={"ttbar,qcd"}, .cuts={"ttbar,qcd","HadTopNoTauCut"},  .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  sh.AddHistos("evt",   { .fill="NHadTop",          .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="Norm", .ranges={0,4, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetMR",         .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetMTR",        .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="Norm", .ranges={0,3000, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_MR",            .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_MTR",           .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R",             .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,2500, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaPhi",      .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaEta",      .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,4, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_DeltaR",        .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Mass",          .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Pz",            .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_Hz",            .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_dPz",           .pfs={"ttbar,qcd,Susy3,Susy4"}, .cuts={"noqcd","NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // 01/16
  //// Signal cut variables
  //sh.AddHistos("evt",   { .fill="H_T",              .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Tevt",           .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Ttt",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_Tex" ,           .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_TttFraction",    .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="H_TexFraction",    .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_R",             .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_R2",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="AK8JetR2",         .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi",   .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  //
  //// 2D Correlation plots
  //sh.AddHistos("evt",   { .fill="HTall_vs_DeltaPhi",        .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTall_vs_R",               .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTall_vs_Rtt",             .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_DeltaPhi", .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_R",        .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="HTttFraction_vs_Rtt",      .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="R_vs_DeltaPhi",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("evt",   { .fill="Rtt_vs_DeltaPhi",          .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  // 01/17
  // Signal cut variables
  sh.AddHistos("evt",   { .fill="H_T",              .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Tevt",           .pfs={"ttbar,qcd,Susy3,Susy4"},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Ttt",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_Tex" ,           .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_TttFraction",    .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="H_TexFraction",    .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R",             .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_R2",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",          .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR2",         .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi",   .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="", .opt="Norm", .ranges={0,0, 0,0} });
  
  // 2D Correlation plots
  sh.AddHistos("evt",   { .fill="HTall_vs_DeltaPhi",        .pfs={"ttbar,qcd,Susy3,Susy4"},                   .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTall_vs_R",               .pfs={"ttbar,qcd,Susy3,Susy4"},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTall_vs_Rtt",             .pfs={"ttbar,qcd,Susy3,Susy4"},                   .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_DeltaPhi", .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_R",        .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HTttFraction_vs_Rtt",      .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2","|DPhi|<2.8"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="R_vs_DeltaPhi",            .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="Rtt_vs_DeltaPhi",          .pfs={"ttbar,qcd,Susy3,Susy4","SideBand,Signal"}, .cuts={"NHadTop==2"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  // 01/20
  // Signal cut variables
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi", .pfs={"RBelow0p25,RAbove0p25"},                         .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TT_AbsDeltaPhi", .pfs={"RBelow0p25,RAbove0p25","SideBand,Signal"},       .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"DPhiBelow0p28,DPhiAbove0p28"},                   .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"DPhiBelow0p28,DPhiAbove0p28","SideBand,Signal"}, .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="AK8JetR",        .pfs={"SideBand,Signal","DPhiBelow0p28,DPhiAbove0p28"}, .cuts={"ttbar","NHadTop==2"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 0,0} });
  
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
      looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/ttbar/*.root");
      //looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/qcd/*.root");
      //looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/susy3body/*.root");
      //looper.AddFile("/data/gridout/jkarancs/SusyAnalysis/B2G/ttree/susy4body/*.root");
    }
    
    while (looper.LoopOnSamples()) {
      while (looper.LoopOnFiles()) {
        TFile *curr_file = looper.CurrentFile();
        reader.Load_Tree(*curr_file,looper.TreeName());
        while(looper.LoopOnEntries()) {
          reader.GetEntry(looper.CurrentEntry());
          d = reader.data;
	  calcVariables(d);
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

