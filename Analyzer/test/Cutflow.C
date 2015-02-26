#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"
#include <iostream>

void get_yields(TFile* f, std::string name,
	     std::string cut1, std::string cut2, std::string cut3, 
	     std::string Cut1, std::string Cut2, std::string Cut3) {
  double signalxsec_pb = 0.0460525;
  double intlumi_invfb = 10;
  
  std::vector<std::string> samples;
  samples.push_back("T5ttttDeg_4bodydec");
  samples.push_back("T5ttttDeg_3bodydec");
  samples.push_back("TTBar");
  samples.push_back("QCD_HT");
  samples.push_back("WJets_HT");
  samples.push_back("SingleTop_tW");
  samples.push_back("DYJets");
  samples.push_back("GGJets_M");
  std::vector<std::string> cuts;
  cuts.push_back("");
  cuts.push_back(std::string("_")+cut3);
  cuts.push_back(std::string("_")+cut2);
  cuts.push_back(std::string("_")+cut1);
  cuts.push_back(std::string("_")+cut2+"_"+cut3);
  cuts.push_back(std::string("_")+cut1+"_"+cut3);
  cuts.push_back(std::string("_")+cut1+"_"+cut2);
  cuts.push_back(std::string("_")+cut1+"_"+cut2+"_"+cut3);
  std::vector<std::string> cutnames;
  cutnames.push_back("No Cut");
  cutnames.push_back(Cut3);
  cutnames.push_back(Cut2);
  cutnames.push_back(Cut1);
  cutnames.push_back(Cut2+" + "+Cut3);
  cutnames.push_back(Cut1+" + "+Cut3);
  cutnames.push_back(Cut1+" + "+Cut2);
  cutnames.push_back(Cut1+" + "+Cut2+" + "+Cut3);
  std::vector<double> nevt(samples.size()+1,0);
  printf("| *Cut* ");
  for (auto sample : samples) {
    printf("|  *%s*  ", sample.c_str());
  } printf("|  *All Background*  |  *Smin = Eff./(5/2+sqrt(B)) *  |\n");
  for (size_t i=0; i<cuts.size(); ++i) {
    printf("| %s ", cutnames[i].c_str());
    nevt[samples.size()]=0;
    for (size_t j=0; j<samples.size(); ++j) {
      TH1D* h = ((TH1D*)f->Get((name+samples[j]+cuts[i]).c_str()));
      nevt[j] = h->Integral();
      if (j>=2) nevt[samples.size()]+=nevt[j];
      printf("|  %.2f ", nevt[j]);
    } 
    double sig = std::min(nevt[0],nevt[1]);
    double bkg = nevt[samples.size()];
    //double sig_to_bkg = sig/sqrt(sig+bkg);
    double smin_punzi_5sig = sig/(signalxsec_pb * intlumi_invfb * 1000) / (5.0/2 + sqrt(bkg) );
    printf("|  %.2f |  %.5f | \n", bkg, smin_punzi_5sig);
  } printf("\n");
}

void Cutflow() {
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("CutTest_R0p35_DPhi2p8_HTall1500_fullstat.root");
  get_yields(f, "H_TexFraction/", "RAbove0p35", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  f->Close();

  f = TFile::Open("CutTest_R0p3_DPhi2p8_HTall1500_fullstat.root");
  get_yields(f, "H_TexFraction/", "RAbove0p3", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  f->Close();
  
}
