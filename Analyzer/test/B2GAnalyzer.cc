#include "TCanvas.h"
#include "TFile.h"
#include "TH3.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include "../interface/Data.h"
#include "../interface/SmartHistos.h"

void get_yields(std::string fname, std::vector<std::string>& samples, size_t nsig, std::string name, 
		std::string cut1, std::string cut2, std::string cut3, 
		std::string Cut1, std::string Cut2, std::string Cut3, std::string pf = "") {
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
  std::vector<double> nevt(samples.size(),0);
  std::vector<double> nevt_sig(nsig,0);
  TFile* f = TFile::Open(fname.c_str());
  printf("| *Cut* ");
  for (auto sample : samples) {
    printf("|  *%s*  ", sample.c_str());
  } printf("|  *All Background*  |  *Smin = Eff./(5/2+sqrt(B)) *  |\n");
  for (size_t i=0; i<cuts.size(); ++i) {
    printf("| %s ", cutnames[i].c_str());
    double min_sig_eff = 9999;
    double tot_bg = 0;
    for (size_t j=0; j<samples.size(); ++j) {
      TH1D* h = ((TH1D*)f->Get((name+samples[j]+pf+cuts[i]).c_str()));
      nevt[j] = h->Integral();
      if (j<nsig) {
	if (i==0) nevt_sig[j]=nevt[j];
	double sig_eff = nevt[j] / nevt_sig[j];
	if (sig_eff<min_sig_eff) min_sig_eff = sig_eff;
      } else {
	tot_bg+=nevt[j];
      }
      printf("|  %.2f ", nevt[j]);
    }
    double smin_punzi_5sig = min_sig_eff / (5.0/2 + sqrt(tot_bg) );
    printf("|  %.2f |  %.5f | \n", tot_bg, smin_punzi_5sig);
  } printf("\n");
  f->Close();
}

std::vector<TH3D*> get_best_cuts_(std::string fname, std::vector<std::string>& samples, size_t nsig,
                                  std::string xname, std::string yname, std::string zname,
                                  std::string xvarname, std::string yvarname, std::string zvarname,
                                  bool xlowcut, bool ylowcut, bool zlowcut, std::string pf="") {
  std::vector<TH3D*> vh;
  std::string name = zname + "_vs_" + yname + "_vs_" + xname+"/";
  TFile* f = TFile::Open(fname.c_str());
  int nx=0, ny=0, nz=0;
  // Make cumulative plots for each samples
  for (size_t i=0; i<samples.size(); ++i) {
    TH3D* h = (TH3D*)f->Get((name+samples[i]+pf).c_str());
    vh.push_back((TH3D*)h->Clone());
    // Make cumulative Histos depending on the applied cut
    // lowcut --> bin low edge is the cut, everything else is summed up
    // highcut--> bin up edge is the cut
    nx = h->GetNbinsX(); ny = h->GetNbinsY(); nz = h->GetNbinsZ();
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
          vh[vh.size()-1]->SetBinContent(binx,biny,binz,sum_3d);
        }
      }
    }
  }
  std::vector<TH3D*> vh_sig_smin;
  for (size_t i=0; i<nsig; ++i) vh_sig_smin.push_back((TH3D*)vh[i]->Clone());
  // Calculate Smin for each sample
  for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
    for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
      for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
	double tot_bg = 0;
	for (int i=samples.size()-1; i>=0; --i) {
	  if (i>=((int)nsig)) {
	    tot_bg += vh[i]->GetBinContent(binx,biny,binz);
	  } else {
	    double sig = vh[i]->GetBinContent(binx,biny,binz);
	    double all_sig = vh[i]->GetBinContent(xlowcut? 1 : nx, ylowcut? 1 : ny, zlowcut? 1 : nz);
	    double eff_sig = sig / all_sig;
	    double smin_punzi_5sig = eff_sig / (5.0/2 + sqrt(tot_bg));
	    vh_sig_smin[i]->SetBinContent(binx,biny,binz, smin_punzi_5sig);
	  }
	}
      }
    }
  }
  // Get the best Smin for each bin but select the worst signal efficiency
  // Do for all combinations of 3 cuts
  std::vector<std::vector<bool> > calc_xyz;
  calc_xyz.push_back(std::vector<bool>({ 0, 0, 0 }));
  calc_xyz.push_back(std::vector<bool>({ 0, 0, 1 }));
  calc_xyz.push_back(std::vector<bool>({ 0, 1, 0 }));
  calc_xyz.push_back(std::vector<bool>({ 1, 0, 0 }));
  calc_xyz.push_back(std::vector<bool>({ 0, 1, 1 }));
  calc_xyz.push_back(std::vector<bool>({ 1, 0, 1 }));
  calc_xyz.push_back(std::vector<bool>({ 1, 1, 0 }));
  calc_xyz.push_back(std::vector<bool>({ 1, 1, 1 }));
  std::vector<double> best_smin;
  std::vector<std::vector<double> > bestcut_xyz;
  std::vector<std::vector<int> > bestbin_xyz;
  for (size_t i_cut=0; i_cut<calc_xyz.size(); ++i_cut) {
    best_smin.push_back(-9999);
    bestcut_xyz.push_back(std::vector<double>({ -9999, -9999, -9999 }));
    bestbin_xyz.push_back(std::vector<int>({ 0, 0, 0 }));
    for (int binx=(xlowcut==calc_xyz[i_cut][0]) ? nx : 1; xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
      for (int biny=(ylowcut==calc_xyz[i_cut][1]) ? ny : 1; ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
        for (int binz=(zlowcut==calc_xyz[i_cut][2]) ? nz : 1; zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
	  double smin_worst_signal = 9999;
	  for (size_t i=0; i<nsig; ++i) {
	    double smin_punzi_5sig = vh_sig_smin[i]->GetBinContent(binx,biny,binz);
	    if (smin_punzi_5sig<smin_worst_signal) smin_worst_signal = smin_punzi_5sig;
	  }
	  double cutx = xlowcut ? vh[0]->GetXaxis()->GetBinLowEdge(binx) : vh[0]->GetXaxis()->GetBinUpEdge(binx);
	  double cuty = ylowcut ? vh[0]->GetYaxis()->GetBinLowEdge(biny) : vh[0]->GetYaxis()->GetBinUpEdge(biny);
	  double cutz = zlowcut ? vh[0]->GetZaxis()->GetBinLowEdge(binz) : vh[0]->GetZaxis()->GetBinUpEdge(binz);
	  if (best_smin[i_cut]<smin_worst_signal) {
	    best_smin[i_cut] = smin_worst_signal;
	    bestcut_xyz[i_cut][0] = cutx; bestcut_xyz[i_cut][1] = cuty; bestcut_xyz[i_cut][2] = cutz;
	    bestbin_xyz[i_cut][0] = binx; bestbin_xyz[i_cut][1] = biny; bestbin_xyz[i_cut][2] = binz;
	  }
	}
      }
    }
  }
  //for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
  //  for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
  //    for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
  //      double smin_worst_signal = 9999;
  //      for (size_t i=0; i<nsig; ++i) {
  //        double smin_punzi_5sig = vh_sig_smin[i]->GetBinContent(binx,biny,binz);
  //        if (smin_punzi_5sig<smin_worst_signal) smin_worst_signal = smin_punzi_5sig;
  //      }
  //      double cutx = xlowcut ? vh[0]->GetXaxis()->GetBinLowEdge(binx) : vh[0]->GetXaxis()->GetBinUpEdge(binx);
  //      double cuty = ylowcut ? vh[0]->GetYaxis()->GetBinLowEdge(biny) : vh[0]->GetYaxis()->GetBinUpEdge(biny);
  //      double cutz = zlowcut ? vh[0]->GetZaxis()->GetBinLowEdge(binz) : vh[0]->GetZaxis()->GetBinUpEdge(binz);
  //      if (best_smin<smin_worst_signal) {
  //        best_smin = smin_worst_signal;
  //        best_cutx = cutx; best_cuty = cuty; best_cutz = cutz;
  //        best_binx = binx; best_biny = biny; best_binz = binz;
  //        //std::cout<<xname<<(xlowcut ? ">" : "<")<<best_cutx<<" && "
  //        //         <<yname<<(ylowcut ? ">" : "<")<<best_cuty<<" && "
  //        //         <<zname<<(zlowcut ? ">" : "<")<<best_cutz<<" Smin = "<<best_smin<<std::endl;
  //      }
  //      //if ( (cutx<0.001&&cuty>3.19&&cutz<1) ||
  //      //     (cutx<0.001&&cuty>2.79&&cuty<2.81&&cutz<1) ||
  //      //     (cutx>0.349&&cutx<0.351&&cuty>2.79&&cuty<2.81&&cutz>1499&&cutz<1501) ) {
  //      //  std::cout<<xname<<" cut="<<cutx<<" "<<yname<<" cut="<<cuty<<" "<<zname<<" cut="<<cutz
  //      //	   <<" Sig1: "<<Sig1<<" Sig2: "<<Sig2<<" Tot BG: "<<bkg<<" Smin: "<<smin_punzi_5sig<<std::endl;
  //      //}
  //    }
  //  }
  //}
  std::vector<std::string> cutnames;
  for (size_t i_cut=0; i_cut<calc_xyz.size(); ++i_cut) {
    std::stringstream ss1; ss1<<xvarname<<(xlowcut?">":"<")<<bestcut_xyz[i_cut][0]; std::string BestCut1 = ss1.str();
    std::stringstream ss2; ss2<<yvarname<<(ylowcut?">":"<")<<bestcut_xyz[i_cut][1]; std::string BestCut2 = ss2.str();
    std::stringstream ss3; ss3<<zvarname<<(zlowcut?">":"<")<<bestcut_xyz[i_cut][2]; std::string BestCut3 = ss3.str();
    //std::cout<<"Best cuts are ( "<<BestCut1<<" && "<<BestCut2<<" && "<<BestCut3<<" ) Smin = "<<best_smin[best_smin.size()-1]<<std::endl;
    std::string cutname = "";
    if (calc_xyz[i_cut][0]) 
      cutname += BestCut1;
    if (calc_xyz[i_cut][1]) {
      if (cutname.size()) cutname += " && ";
      cutname += BestCut2;
    }
    if (calc_xyz[i_cut][2]) {
      if (cutname.size()) cutname += " && ";
      cutname += BestCut3;
    }
    if (cutname.size()==0) cutname = "No cut";
    cutnames.push_back(cutname);
  }
  std::vector<double> nevt(samples.size(),0);
  std::vector<double> nevt_sig(nsig,0);
  printf("| *Best cut* ");
  for (auto sample : samples) {
    printf("|  *%s*  ", sample.c_str());
  } printf("|  *All bg*  |  *Smin = Eff./(5/2+sqrt(B)) *  |\n");
  for (size_t i_cut=0; i_cut<cutnames.size(); ++i_cut) {
    printf("| %s ", cutnames[i_cut].c_str());
    double min_sig_eff = 9999;
    double tot_bg = 0;
    for (size_t j=0; j<samples.size(); ++j) {
      nevt[j] = vh[j]->GetBinContent(bestbin_xyz[i_cut][0], bestbin_xyz[i_cut][1], bestbin_xyz[i_cut][2]);
      if (j<nsig) {
        if (i_cut==0) nevt_sig[j]=nevt[j];
        double sig_eff = nevt[j] / nevt_sig[j];
        if (sig_eff<min_sig_eff) min_sig_eff = sig_eff;
      } else {
        tot_bg+=nevt[j];
      }
      printf("|  %.2f ", nevt[j]);
    }
    double smin_punzi_5sig = min_sig_eff / (5.0/2 + sqrt(tot_bg) );
    printf("|  %.2f |  %.5f | \n", tot_bg, smin_punzi_5sig);
  } printf("\n");
  
  f->Close();
  return vh;
}

int main() {
  gStyle->SetOptStat(0);
  
  std::vector<std::string> samples;
  samples.push_back("T5ttttDeg_4bodydec");
  samples.push_back("T5ttttDeg_3bodydec");
  samples.push_back("QCD_HT");
  samples.push_back("TTBar");
  samples.push_back("WJets_HT");
  samples.push_back("SingleTop_tW");
  //samples.push_back("DYJets");
  //samples.push_back("GGJets_M");
  
  std::vector<std::string> newsamples;
  newsamples.push_back("T5ttttDeg_4bodydec_mGo1300");
  newsamples.push_back("T5ttttDeg_3bodydec_mGo1300");
  //newsamples.push_back("T5ttttDeg_4bodydec_mGo1000");
  //newsamples.push_back("T5ttttDeg_3bodydec_mGo1000");
  newsamples.push_back("TTJets");
  //newsamples.push_back("GJets_HT");
  //newsamples.push_back("GGJets_M");
  //newsamples.push_back("WJets");
  newsamples.push_back("WJets_HT");
  newsamples.push_back("ZJets_HT");
  newsamples.push_back("QCD_HT");
  //newsamples.push_back("QCD_Pt_bcToE");
  newsamples.push_back("T_tW");
  //newsamples.push_back("TToLep_s_t");
  //newsamples.push_back("DYJets");
  //newsamples.push_back("DYJets_HT");
  
  std::vector<std::string> newsamples2;
  //newsamples2.push_back("T5ttttDeg_4bodydec_mGo1300");
  //newsamples2.push_back("T5ttttDeg_3bodydec_mGo1300");
  newsamples2.push_back("T5ttttDeg_4bodydec_mGo1000");
  newsamples2.push_back("T5ttttDeg_3bodydec_mGo1000");
  newsamples2.push_back("TTJets");
  //newsamples2.push_back("GJets_HT");
  //newsamples2.push_back("GGJets_M");
  //newsamples2.push_back("WJets");
  newsamples2.push_back("WJets_HT");
  newsamples2.push_back("ZJets_HT");
  newsamples2.push_back("QCD_HT");
  //newsamples2.push_back("QCD_Pt_bcToE");
  newsamples2.push_back("T_tW");
  //newsamples2.push_back("TToLep_s_t");
  //newsamples2.push_back("DYJets");
  //newsamples2.push_back("DYJets_HT");

  std::string pfnames = "";
  for (auto sample : newsamples) pfnames += sample +";";
  pfnames = pfnames.substr(0, pfnames.size()-1);
  
  std::string latexnames = "";
  latexnames +="T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.3TeV);";
  latexnames +="T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.3TeV);";
  //latexnames +="T5tttt (#tilde{g}#rightarrowt#tilde{t}_{4body}, M_{#tilde{g}}=1.0TeV);";
  //latexnames +="T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.0TeV);";
  latexnames +="t#bar{t}+Jets;";
  //latexnames +="G+Jets (HT bins);";
  //latexnames +="GG+Jets (M bins);";
  //latexnames +="W+Jets #rightarrow l+#nu;";
  latexnames +="W+Jets #rightarrow l+#nu (HT bins);";
  latexnames +="Z+Jets #rightarrow #nu#nu (HT bins);";
  latexnames +="QCD (HT bins);";
  //latexnames +="QCD (Pt bins, b/c#rightarrowe);";
  latexnames +="single t/#bar{t} (tW channel);";
  //latexnames +="single t/#bar{t}#rightarrowl (s,t channel);";
  //latexnames +="DY+Jets #rightarrow l+l;";
  //latexnames +="DY+Jets #rightarrow l+l (HT bins);";
  latexnames = latexnames.substr(0, latexnames.size()-1);
  //std::string colors = "1,14,4,2,804,3,617,401";
  //std::string colors = "1,14,4,38,2,804,3,617,882,401,797";
  //std::string colors = "1,14,4,38,804,2,3,617,882,401,797";
  std::string colors = "1,14,2,3,401,4,617";
  
  // Remake Cut dependant distributions
#define CUTR 0.4
#define CUTDPHI 2.8
#define CUTHTALL 0
  
  std::stringstream ss1_pf; ss1_pf<<"RBelow"<<CUTR<<";RAbove"<<CUTR; std::string pf1 = ss1_pf.str(); std::replace(pf1.begin(), pf1.end(), '.', 'p');
  std::stringstream ss2_pf; ss2_pf<<"DPhiBelow"<<CUTDPHI<<";DPhiAbove"<<CUTDPHI; std::string pf2 = ss2_pf.str(); std::replace(pf2.begin(), pf2.end(), '.', 'p');
  std::stringstream ss3_pf; ss3_pf<<"HtAllBelow"<<CUTHTALL<<";HtAllAbove"<<CUTHTALL; std::string pf3 = ss3_pf.str();
  std::stringstream ss1_lat; ss1_lat<<"R < "<<CUTR<<";R > "<<CUTR; std::string lat1 = ss1_lat.str();
  std::stringstream ss2_lat; ss2_lat<<"#Delta#phi_{t#bar{t}} < "<<CUTDPHI<<";#Delta#phi_{t#bar{t}} > "<<CUTDPHI; std::string lat2 = ss2_lat.str();
  std::stringstream ss3_lat; ss3_lat<<"H_{T,all} < "<<CUTHTALL<<";H_{T,all} > "<<CUTHTALL; std::string lat3 = ss3_lat.str();
  
  Data d;
  SmartHistos sh;
  sh.AddHistoType("post");
  sh.AddHistoType("post2");
  
  double weight = 0;
  sh.SetHistoWeights({[&weight](){ return weight; }});
  
  size_t iSample = 0;
  sh.AddNewPostfix("AllSamples",     [&iSample](){ return iSample; }, pfnames, latexnames, colors);
  sh.AddNewPostfix("CutR",           [&d](){ return d.evt.R > CUTR;                  },  pf1, lat1, "4,2");
  sh.AddNewPostfix("CutDPhi",        [&d](){ return fabs(d.evt.TTHadDPhi) > CUTDPHI; },  pf2, lat2, "4,2");
  sh.AddNewPostfix("CutHtAll",       [&d](){ return d.evt.HTall > CUTHTALL;          },  pf3, lat3, "4,2");
  sh.AddNewPostfix("NTopHad",        [&d](){ return (d.evt.NTopHad>1)+(d.evt.NTopHad>2); }, "0To1HadTop;2HadTop;3HadTop", "N_{top,hadronic}<2;N_{top,hadronic}=2;N_{top,hadronic}=3", "4,2,3");
  
  sh.AddNewFillParam("R",           { .nbin=  20, .bins={   0,     1.0}, .fill=[&d](){ return d.evt.R;                }, .axis_title="R"});
  sh.AddNewFillParam("DPhi",        { .nbin=  16, .bins={   0,     3.2}, .fill=[&d](){ return fabs(d.evt.TTHadDPhi);  }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("HtAll",       { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HTall;            }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAllCoarse", { .nbin=  20, .bins={   0,    6000}, .fill=[&d](){ return d.evt.HTall;            }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("RCoarse",     { .nbin=  10, .bins={   0,    1.25}, .fill=[&d](){ return d.evt.R;                }, .axis_title="R"});
  sh.AddNewFillParam("MR",          { .nbin=  16, .bins={300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1200, 1600, 2000, 2800, 3500}, .fill=[&d](){ return d.evt.MR;               }, .axis_title="M_{R} (GeV)"});
  
  // N-(1,2,3) plots
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutR"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutHtAll"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutR","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutR","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutDPhi","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",               .pfs={"AllSamples","NTopHad","CutR","CutDPhi","CutHtAll"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutR"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutHtAll"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutR","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutR","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutDPhi","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",           .pfs={"AllSamples","NTopHad","CutR","CutDPhi","CutHtAll"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutR"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutHtAll"},                .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutR","CutDPhi"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutR","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutDPhi","CutHtAll"},           .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi",       .pfs={"AllSamples","NTopHad","CutR","CutDPhi","CutHtAll"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  // 2D Correlation plots
  sh.AddHistos("post",   { .fill="HtAllCoarse_vs_TTHadDPhi",   .pfs={"AllSamples","NTopHad"},      .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("post",   { .fill="HtAllCoarse_vs_R",           .pfs={"AllSamples","NTopHad"},      .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("post",   { .fill="R_vs_TTHadDPhi",             .pfs={"AllSamples","NTopHad"},      .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("post",   { .fill="HtAllCoarse_vs_TTHadDPhi",   .pfs={"AllSamples","NTopHad","CutR"}, .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("post",   { .fill="HtAllCoarse_vs_R",           .pfs={"AllSamples","NTopHad","CutDPhi"}, .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("post",   { .fill="R_vs_TTHadDPhi",             .pfs={"AllSamples","NTopHad","CutHtAll"}, .cuts={}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  // Plots for Background estimation
  sh.AddHistos("post",   { .fill="R",         .pfs={"CutDPhi",           "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"CutHtAll",          "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"CutDPhi","CutHtAll","NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"CutHtAll","CutDPhi","NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"CutR",              "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"CutHtAll",          "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"CutR","CutHtAll",   "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"CutHtAll","CutR",   "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"CutR",              "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"CutDPhi",           "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"CutR","CutDPhi",    "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"CutDPhi","CutR",    "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"NTopHad","CutDPhi",           "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"NTopHad","CutHtAll",          "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"NTopHad","CutDPhi","CutHtAll","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="R",         .pfs={"NTopHad","CutHtAll","CutDPhi","AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"NTopHad","CutR",              "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"NTopHad","CutHtAll",          "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"NTopHad","CutR","CutHtAll",   "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="TTHadDPhi", .pfs={"NTopHad","CutHtAll","CutR",   "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"NTopHad","CutR",              "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"NTopHad","CutDPhi",           "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"NTopHad","CutR","CutDPhi",    "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("post",   { .fill="HtAll",     .pfs={"NTopHad","CutDPhi","CutR",    "AllSamples"}, .cuts={}, .draw="", .opt="Log", .ranges={0,0, 1e-3,1e5} });
  
  // Razor fit 2D plots
  sh.AddHistos("post2",   { .fill="RCoarse_vs_MR", .pfs={"CutDPhi", "NTopHad","AllSamples"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  
  std::string R = "RFine";
  
  // Make tables
  //get_yields("ROOT_output/CutTest_R0p35_DPhi2p8_HTall1500_fullstat.root", samples, "H_Tall/", "RAbove0p35", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_R0p3_DPhi2p8_HTall1500_fullstat.root", samples, "H_Tall/", "RAbove0p3", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_R_DPhi_HTall_fullstat.root", samples, "H_Tall/", "RAbove0p35", "DPhiBelow2p8", "HTallAbove1500", "R", "DPhi", "HTall", 1);
  //get_yields("ROOT_output/BestCuts_fullstat.root",                samples,     2, "H_Tall/", "RAbove0p32", "DPhiBelow2p8", "HTallAbove1450", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_NewSamples_NoLepTop_10th.root", newsamples,  2, "HtAll/",  "RAbove0p32", "DPhiBelow2p8", "HtAllAbove1450", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_NewSamples_NoLepTop_fullstat.root", newsamples2, 2, "HtAll/",  "RAbove0p32", "DPhiBelow2p8", "HtAllAbove1450", "R", "DPhi", "HTall");
  //get_yields("ROOT_output/CutTest_NewSamples_fullstat.root",          newsamples2, 2, "HtAll/",  "RAbove0p32", "DPhiBelow2p8", "HtAllAbove1450", "R", "DPhi", "HTall");
  
  // Calculate best set of cuts
  //get_best_cuts_("ROOT_output/BestCuts_fullstat.root", samples, "TT_RFine",    "TT_AbsDeltaPhiFine", "TT_MRFine",  1, 0, 1);
  //get_best_cuts_("ROOT_output/BestCuts_fullstat.root", samples, "TT_RFine",    "TT_AbsDeltaPhiFine", "H_TallFine", 1, 0, 1);
  //get_best_cuts_("ROOT_output/BestCuts_fullstat.root",                    samples,     2, "AK8JetRFine", "TT_AbsDeltaPhiFine", "H_TallFine", "R", "DPhi", "HTall", 1, 0, 1);
  //get_best_cuts_("ROOT_output/CutTest_NewSamples_NoLepTop_fullstat.root", newsamples2, 2, "RFine",       "TTHadDPhiFine",      "HtAllFine",  "R", "DPhi", "HTall", 1, 0, 1); // Mgluino = 1000
  //get_best_cuts_("ROOT_output/CutTest_NewSamples_fullstat.root",          newsamples2, 2, "RFine",       "TTHadDPhiFine",      "HtAllFine",  "R", "DPhi", "HTall", 1, 0, 1);
  //++ get_best_cuts_("ROOT_output/CutTest_NewSamples_NoLepTop_fullstat.root",                        newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1); // Mgluino = 1300
  //++ get_best_cuts_("ROOT_output/CutTest_NewSamples_fullstat.root",                                 newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1);
  //++ get_best_cuts_("ROOT_output/BackGroundEstimate_NTopHadSideBand_TightLeptonVeto_fullstat.root", newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1);
  //++ get_best_cuts_("ROOT_output/BackGroundEstimate_NTopHadSideBand_LeptonVeto_fullstat.root",      newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1);
  //++ get_best_cuts_("ROOT_output/BackGroundEstimate_NTopHadSideBand_NoLeptonVeto_50MassCut_fullstat.root",      newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1, "_2HadTop");
  //get_best_cuts_("ROOT_output/BackGroundEstimate_Apr02_140MassCut_fullstat.root",      newsamples,  2, R, "TTHadDPhiFine", "HtAllFine", "R", "DPhi", "HTall", 1, 0, 1, "_2HadTop");
  get_best_cuts_("ROOT_output/BackGroundEstimate_Apr02_MRtt_fullstat.root",      newsamples,  2, R, "TTHadDPhiFine", "TTHadMRFine", "R", "DPhi", "MRtt", 1, 0, 1, "_2HadTop");
  
  //TFile* f = TFile::Open("ROOT_output/CutTest_NewSamples_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/CutTest_NewSamples_NoLepTop_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_NTopHadSideBand_TightLeptonVeto_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_NTopHadSideBand_LeptonVeto_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_NTopHadSideBand_TightLeptonVeto_50MassCut_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_NTopHadSideBand_NoLeptonVeto_50MassCut_fullstat.root");
  //TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_Apr02_140MassCut_fullstat.root");
  TFile* f = TFile::Open("ROOT_output/BackGroundEstimate_Apr02_MRtt_fullstat.root");
  int xlowcut = 1, ylowcut = 0, zlowcut=1;
  for (iSample=0; iSample<newsamples.size(); ++iSample) {
  //for (iSample=2; iSample==2; ++iSample) {
    std::cout<<"Start filling plots for "<<newsamples[iSample]<<"\n";
    //TH3D* h1 = (TH3D*)f->Get((std::string("HtAllFine_vs_TTHadDPhiFine_vs_")+R+"/"+newsamples[iSample]+"_0HadTop").c_str());
    //TH3D* h2 = (TH3D*)f->Get((std::string("HtAllFine_vs_TTHadDPhiFine_vs_")+R+"/"+newsamples[iSample]+"_1HadTop").c_str());
    //TH3D* h3 = (TH3D*)f->Get((std::string("HtAllFine_vs_TTHadDPhiFine_vs_")+R+"/"+newsamples[iSample]+"_2HadTop").c_str());
    TH3D* h1 = (TH3D*)f->Get((std::string("HtAllFine_vs_TTHadDPhiFine_vs_")+R+"/"+newsamples[iSample]+"_0To1HadTop").c_str());
    TH3D* h2 = (TH3D*)f->Get((std::string("HtAllFine_vs_TTHadDPhiFine_vs_")+R+"/"+newsamples[iSample]+"_2HadTop").c_str());
    int nx = h1->GetNbinsX(), ny = h1->GetNbinsY(), nz = h1->GetNbinsZ();
    for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
      for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
        for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
          weight = h1->GetBinContent(binx,biny,binz);
          d.evt.R = h1->GetXaxis()->GetBinCenter(binx);
          d.evt.TTHadDPhi = h1->GetYaxis()->GetBinCenter(biny);
          d.evt.HTall = h1->GetZaxis()->GetBinCenter(binz);
          d.evt.NTopHad = 0;
          sh.Fill("post");
          weight = h2->GetBinContent(binx,biny,binz);
          d.evt.R = h2->GetXaxis()->GetBinCenter(binx);
          d.evt.TTHadDPhi = h2->GetYaxis()->GetBinCenter(biny);
          d.evt.HTall = h2->GetZaxis()->GetBinCenter(binz);
          d.evt.NTopHad = 2;
          sh.Fill("post");
          //weight = h3->GetBinContent(binx,biny,binz);
          //d.evt.R = h3->GetXaxis()->GetBinCenter(binx);
          //d.evt.TTHadDPhi = h3->GetYaxis()->GetBinCenter(biny);
          //d.evt.HTall = h3->GetZaxis()->GetBinCenter(binz);
          //d.evt.NTopHad = 2;
          //sh.Fill("post");
        }
      }
    }
    //h1 = (TH3D*)f->Get((std::string("MRFine_vs_TTHadDPhiFine_vs_RFine/")+newsamples[iSample]+"_0HadTop").c_str());
    //h2 = (TH3D*)f->Get((std::string("MRFine_vs_TTHadDPhiFine_vs_RFine/")+newsamples[iSample]+"_1HadTop").c_str());
    //h3 = (TH3D*)f->Get((std::string("MRFine_vs_TTHadDPhiFine_vs_RFine/")+newsamples[iSample]+"_2HadTop").c_str());
    //nx = h1->GetNbinsX(), ny = h1->GetNbinsY(), nz = h1->GetNbinsZ();
    //for (int binx=(xlowcut ? nx : 1); xlowcut ? (binx>=1) : (binx<=nx); xlowcut ? --binx : ++binx) {
    //  for (int biny=(ylowcut ? ny : 1); ylowcut ? (biny>=1) : (biny<=ny); ylowcut ? --biny : ++biny) {
    //    for (int binz=(zlowcut ? nz : 1); zlowcut ? (binz>=1) : (binz<=nz); zlowcut ? --binz : ++binz) {
    //      // MR plots
    //      weight = h1->GetBinContent(binx,biny,binz);
    //      d.evt.R = h1->GetXaxis()->GetBinCenter(binx);
    //      d.evt.TTHadDPhi = h1->GetYaxis()->GetBinCenter(biny);
    //      d.evt.MR = h1->GetZaxis()->GetBinCenter(binz);
    //      d.evt.NTopHad = 0;
    //      sh.Fill("post2");
    //      weight = h2->GetBinContent(binx,biny,binz);
    //      d.evt.R = h2->GetXaxis()->GetBinCenter(binx);
    //      d.evt.TTHadDPhi = h2->GetYaxis()->GetBinCenter(biny);
    //      d.evt.MR = h2->GetZaxis()->GetBinCenter(binz);
    //      d.evt.NTopHad = 1;
    //      sh.Fill("post2");
    //      weight = h3->GetBinContent(binx,biny,binz);
    //      d.evt.R = h3->GetXaxis()->GetBinCenter(binx);
    //      d.evt.TTHadDPhi = h3->GetYaxis()->GetBinCenter(biny);
    //      d.evt.MR = h3->GetZaxis()->GetBinCenter(binz);
    //      d.evt.NTopHad = 2;
    //      sh.Fill("post2");
    //    }
    //  }
    //}
  }
  std::cout<<"Plots filled.\n\n";
  
  std::string cutx = pf1.substr((xlowcut ? pf1.find(";")+1 : 0), (xlowcut ? pf1.size()-pf1.find(";")-1 : pf1.find(";")));
  std::string cuty = pf2.substr((ylowcut ? pf2.find(";")+1 : 0), (ylowcut ? pf2.size()-pf2.find(";")-1 : pf2.find(";")));
  std::string cutz = pf3.substr((zlowcut ? pf3.find(";")+1 : 0), (zlowcut ? pf3.size()-pf3.find(";")-1 : pf3.find(";")));
  std::stringstream ss1_cut; ss1_cut<<"R "<<(xlowcut ? ">":"<")<<" "<<CUTR; std::string cut1 = ss1_cut.str();
  std::stringstream ss2_cut; ss2_cut<<"DPhi "<<(ylowcut ? ">":"<")<<" "<<CUTDPHI; std::string cut2 = ss2_cut.str();
  std::stringstream ss3_cut; ss3_cut<<"HTall "<<(zlowcut ? ">":"<")<<" "<<CUTHTALL; std::string cut3 = ss3_cut.str();
  
  TFile* file = new TFile((std::string("ROOT_output/SignalCutPlots_"+cutx+"_"+cuty+"_"+cutz+".root")).c_str(),"recreate");
  sh.DrawPlots();
  sh.Write();
  file->Close();
  
  std::cout<<"Manual cuts are ( "<<cut1<<" && "<<cut2<<" && "<<cut3<<" )"<<std::endl;
  //get_yields(std::string("ROOT_output/SignalCutPlots_"+cutx+"_"+cuty+"_"+cutz+".root"), newsamples, 2, "HtAll/", cutx, cuty, cutz, cut1, cut2, cut3, "");
  get_yields(std::string("ROOT_output/SignalCutPlots_"+cutx+"_"+cuty+"_"+cutz+".root"), newsamples, 2, "HtAll/", cutx, cuty, cutz, cut1, cut2, cut3, "_2HadTop");
  
  return 1;
}
