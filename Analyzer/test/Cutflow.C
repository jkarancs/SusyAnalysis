#include "TCanvas.h"
#include "TFile.h"
#include "TH3.h"
#include "TStyle.h"
#include <iostream>

void get_yields(std::string fname, std::vector<std::string>& samples, std::string name,
		std::string cut1, std::string cut2, std::string cut3, 
		std::string Cut1, std::string Cut2, std::string Cut3, bool full_eff=1, double signalxsec_pb = 0.0460525) {
  double intlumi_invfb = 10;
  
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
  TFile* f = TFile::Open(fname.c_str());
  printf("| *Cut* ");
  for (auto sample : samples) {
    printf("|  *%s*  ", sample.c_str());
  } printf("|  *All Background*  |  *Smin = Eff./(5/2+sqrt(B)) *  |\n");
  double allsig = 0;
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
    if (i==0) allsig = sig;
    double sig_eff = sig / (full_eff ? signalxsec_pb * intlumi_invfb * 1000 : allsig);
    double bkg = nevt[samples.size()];
    //double sig_to_bkg = sig/sqrt(sig+bkg);
    double smin_punzi_5sig = sig_eff / (5.0/2 + sqrt(bkg) );
    printf("|  %.2f |  %.5f | \n", bkg, smin_punzi_5sig);
  } printf("\n");
  f->Close();
}

void get_best_cuts_(std::string fname, std::vector<std::string>& samples,
		    std::string xname, std::string yname, std::string zname,
		    bool xlowcut, bool ylowcut, bool zlowcut, bool full_eff = 1, double signalxsec_pb = 0.0460525) {
  double intlumi_invfb = 10;
  
  std::string name = zname + "_vs_" + yname + "_vs_" + xname+"/";
  //std::cout<<name<<std::endl;
  TFile* f = TFile::Open(fname.c_str());
  TH3D *sig1 = 0, *sig2 = 0, *totbg = 0;
  int nx = 0, ny = 0, nz = 0;
  double allsig = 0;
  for (size_t i=0; i<samples.size(); ++i) {
    TH3D* h = (TH3D*)f->Get((name+samples[i]).c_str());
    if (i==0||(i==1&&allsig>h->Integral())) allsig = h->Integral();
    //std::cout<<samples[i]<<": "<<h->Integral()<<std::endl;
    if (i==0) sig1 = (TH3D*)h->Clone();
    if (i==1) sig2 = (TH3D*)h->Clone();
    if (i==2) totbg = (TH3D*)h->Clone();
    nx = h->GetNbinsX();
    ny = h->GetNbinsY();
    nz = h->GetNbinsZ();
    // Make cumulative Histos depending on the applied cut
    // lowcut --> bin low edge is the cut, everything else is summed up
    // highcut--> bin up edge is the cut
    for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
      for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
	double sumz = 0;
	for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
	  double bincont = h->GetBinContent(binx,biny,binz);
	  int binx_prev = xlowcut ? binx+1 : binx-1;
	  int biny_prev = ylowcut ? biny+1 : biny-1;
	  double sumx = (binx==(xlowcut ? nx : 1)) ? 0 : h->GetBinContent(binx_prev,biny,binz);
	  double sumy = ((binx==(xlowcut ? nx : 1))||(biny==(ylowcut ? ny : 1))) ? 0 : 
	    h->GetBinContent(binx,biny_prev,binz) - h->GetBinContent(binx_prev,biny_prev,binz);
	  sumz += bincont;
	  double sum_3d = sumx + sumy + sumz;
	  h->SetBinContent(binx,biny,binz,sum_3d);
	  if (i==0) sig1->SetBinContent(binx,biny,binz,sum_3d);
	  else if (i==1) sig2->SetBinContent(binx,biny,binz,sum_3d);
	  else if (i==2) totbg->SetBinContent(binx,biny,binz,sum_3d);
	  else totbg->SetBinContent(binx,biny,binz,sum_3d+totbg->GetBinContent(binx,biny,binz));
	}
      }
    }
  }
  double bestcutx = -9999;
  double bestcuty = -9999;
  double bestcutz = -9999;
  double best_smin = -9999;
  for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
    for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
      for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
	double cutx = xlowcut ? sig1->GetXaxis()->GetBinLowEdge(binx) : sig1->GetXaxis()->GetBinUpEdge(binx);
	double cuty = ylowcut ? sig1->GetYaxis()->GetBinLowEdge(biny) : sig1->GetYaxis()->GetBinUpEdge(biny);
	double cutz = zlowcut ? sig1->GetZaxis()->GetBinLowEdge(binz) : sig1->GetZaxis()->GetBinUpEdge(binz);
	double Sig1 = sig1->GetBinContent(binx,biny,binz);
	double Sig2 = sig2->GetBinContent(binx,biny,binz);
	double sig = std::min(Sig1, Sig2);
	double sig_eff = sig / (full_eff ? signalxsec_pb * intlumi_invfb * 1000 : allsig);
	double bkg = totbg->GetBinContent(binx,biny,binz);
	//double sig_to_bkg = sig/sqrt(sig+bkg);
	double smin_punzi_5sig = sig_eff / (5.0/2 + sqrt(bkg) );
	if (best_smin<smin_punzi_5sig) {
	  best_smin = smin_punzi_5sig;
	  bestcutx = cutx;
	  bestcuty = cuty;
	  bestcutz = cutz;
	  //std::cout<<xname<<(xlowcut ? ">" : "<")<<bestcutx<<" && "
          //         <<yname<<(ylowcut ? ">" : "<")<<bestcuty<<" && "
          //         <<zname<<(zlowcut ? ">" : "<")<<bestcutz<<" Smin = "<<best_smin<<std::endl;
	}
	//if ( (cutx<0.001&&cuty>3.19&&cutz<1) ||
	//     (cutx<0.001&&cuty>2.79&&cuty<2.81&&cutz<1) ||
	//     (cutx>0.349&&cutx<0.351&&cuty>2.79&&cuty<2.81&&cutz>1499&&cutz<1501) ) {
	//  std::cout<<xname<<" cut="<<cutx<<" "<<yname<<" cut="<<cuty<<" "<<zname<<" cut="<<cutz
	//	   <<" Sig1: "<<Sig1<<" Sig2: "<<Sig2<<" Tot BG: "<<bkg<<" Smin: "<<smin_punzi_5sig<<std::endl;
	//}
      }
    }
  }
  std::cout<<"Best cut is ( "
           <<xname<<(xlowcut ? ">" : "<")<<bestcutx<<" && "
	   <<yname<<(ylowcut ? ">" : "<")<<bestcuty<<" && "
	   <<zname<<(zlowcut ? ">" : "<")<<bestcutz<<" ) Smin = "<<best_smin<<std::endl;
  f->Close();
}

void Cutflow() {
  gStyle->SetOptStat(0);
  
  std::vector<std::string> samples;
  samples.push_back("T5ttttDeg_4bodydec");
  samples.push_back("T5ttttDeg_3bodydec");
  samples.push_back("TTBar");
  samples.push_back("QCD_HT");
  samples.push_back("WJets_HT");
  samples.push_back("SingleTop_tW");
  samples.push_back("DYJets");
  samples.push_back("GGJets_M");
  
  std::vector<std::string> newsamples;
  newsamples.push_back("T5ttttDeg_4bodydec_mGo1300");
  newsamples.push_back("T5ttttDeg_3bodydec_mGo1300");
  //newsamples.push_back("T5ttttDeg_4bodydec_mGo1000");
  //newsamples.push_back("T5ttttDeg_3bodydec_mGo1000");
  newsamples.push_back("TTBar");
  newsamples.push_back("T_tW");
  newsamples.push_back("TToLep_s_t");
  newsamples.push_back("QCD_HT");
  newsamples.push_back("QCD_Pt_bcToE");
  newsamples.push_back("WJets_HT");
  newsamples.push_back("ZJets_HT");
  newsamples.push_back("DYJets_HT");
  newsamples.push_back("GJets_HT");
  newsamples.push_back("GGJets_M");

  std::vector<std::string> newsamples_mGo1000;
  //newsamples.push_back("T5ttttDeg_4bodydec_mGo1300");
  //newsamples.push_back("T5ttttDeg_3bodydec_mGo1300");
  newsamples_mGo1000.push_back("T5ttttDeg_4bodydec_mGo1000");
  newsamples_mGo1000.push_back("T5ttttDeg_3bodydec_mGo1000");
  newsamples_mGo1000.push_back("TTBar");
  newsamples_mGo1000.push_back("T_tW");
  newsamples_mGo1000.push_back("TToLep_s_t");
  newsamples_mGo1000.push_back("QCD_HT");
  newsamples_mGo1000.push_back("QCD_Pt_bcToE");
  newsamples_mGo1000.push_back("WJets_HT");
  newsamples_mGo1000.push_back("ZJets_HT");
  newsamples_mGo1000.push_back("DYJets_HT");
  newsamples_mGo1000.push_back("GJets_HT");
  newsamples_mGo1000.push_back("GGJets_M");
  
  // Calculate best set of cuts
  //get_best_cuts_("ROOT_output/BestCuts_fullstat.root", samples, "TT_RFine",    "TT_AbsDeltaPhiFine", "TT_MRFine",  1, 0, 1, 1);
  //get_best_cuts_("ROOT_output/BestCuts_fullstat.root", samples, "TT_RFine",    "TT_AbsDeltaPhiFine", "H_TallFine", 1, 0, 1, 1);
  get_best_cuts_("ROOT_output/BestCuts_fullstat.root",       samples,            "AK8JetRFine", "TT_AbsDeltaPhiFine", "H_TallFine", 1, 0, 1, 0);
  get_best_cuts_("ROOT_output/CutTest_NewSamples_10th.root", newsamples,         "RFine",       "TTHadDPhiFine",      "HtAllFine",  1, 0, 1, 0);
  get_best_cuts_("ROOT_output/CutTest_NewSamples_10th.root", newsamples_mGo1000, "RFine",       "TTHadDPhiFine",      "HtAllFine",  1, 0, 1, 0);
  
  // Make tables
  //get_yields("ROOT_output/CutTest_R0p35_DPhi2p8_HTall1500_fullstat.root", samples, "H_TexFraction/", "RAbove0p35", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_R0p3_DPhi2p8_HTall1500_fullstat.root", samples, "H_TexFraction/", "RAbove0p3", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_R_DPhi_HTall_fullstat.root", samples, "H_TexFraction/", "RAbove0p35", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall", 1);
  get_yields("ROOT_output/BestCuts_fullstat.root",       samples,             "H_TexFraction/", "RAbove0p32", "DPhiBelow2p8", "HTallAbove1450", "R", "DPhi", "HTall", 0);
  get_yields("ROOT_output/CutTest_NewSamples_10th.root", newsamples,          "HtExFraction/",  "RAbove0p32", "DPhiBelow2p8", "HtAllAbove1450", "R", "DPhi", "HTall", 0);
  get_yields("ROOT_output/CutTest_NewSamples_10th.root", newsamples_mGo1000,  "HtExFraction/",  "RAbove0p32", "DPhiBelow2p8", "HtAllAbove1450", "R", "DPhi", "HTall", 0);
}
