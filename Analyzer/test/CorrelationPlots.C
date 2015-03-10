#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

void do_plot(TCanvas* c, int n, TFile *f, const char* name, std::string title, float limit) {
  TVirtualPad* p = c->cd(n);
  p->SetLeftMargin(0.15);
  p->SetRightMargin(0.2);
  if (limit>0) p->SetLogz(1);
  TH2D* h = ((TH2D*)f->Get(name));
  h->SetTitle(title.c_str());
  h->GetYaxis()->SetTitleOffset(1.2);
  if (h->GetEntries()) h = (TH2D*)h->DrawNormalized("COLZ");
  else h->Draw("COLZ");
  if (limit>0) h->GetZaxis()->SetRangeUser(0.0001,limit);
}


void make_comparison_plot(TFile* f, std::string name, std::string canname, float limit=0, std::string pf="") {
  TCanvas* c = new TCanvas((canname+pf).c_str(),name.c_str(), 900,900);
  c->SetLeftMargin(0.2);
  c->Divide(2,2);
  std::string Pf = pf;
  Pf.erase(0,1);
  do_plot(c, 1, f, (name+"/ttbar"+pf).c_str(), std::string("ttbar ")+Pf, limit);
  do_plot(c, 2, f, (name+"/susy3body"+pf).c_str(), std::string("T5tttt 3-body decay ")+Pf, limit);
  do_plot(c, 3, f, (name+"/qcd"+pf).c_str(), std::string("QCD ")+Pf, limit);
  do_plot(c, 4, f, (name+"/susy4body"+pf).c_str(), std::string("T5tttt 4-body decay ")+Pf, limit);
  c->SaveAs((std::string("Plots/2015_01_17/Corr2D/")+canname+pf+".png").c_str());
}

void save_1d(TFile* f, float limit, std::string name, std::string canname) {
  TCanvas* c = ((TCanvas*)f->Get(name.c_str()));
  TH1D* h = (TH1D*)c->GetListOfPrimitives()->At(0);
  if (limit>0) h->GetYaxis()->SetRangeUser(0.00001,limit);
  c->SetLogy(1);
  c->Draw();
  //c->SaveAs((std::string("Plots/2015_01_17/")+canname+".png").c_str());
}

void CorrelationPlots() {
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("feb23_AllSamples_HTBinned_fullstat.root");
  
  // 01/16
  //save_1d(f, 0.1, "HT/ttbar,qcd,Susy3,Susy4",           "H_T");
  //save_1d(f, 0.1, "HTevt/ttbar,qcd,Susy3,Susy4",        "H_Tevt");
  //save_1d(f, 0.3, "HTttFraction/ttbar,qcd,Susy3,Susy4", "H_TttFraction");
  //save_1d(f, 0.3, "HTexFraction/ttbar,qcd,Susy3,Susy4", "H_TexFraction");
  //save_1d(f, 0.2, "R/ttbar,qcd,Susy3,Susy4",            "AK8JetR");
  //save_1d(f, 0.7, "R2/ttbar,qcd,Susy3,Susy4",           "AK8JetR2");
  //save_1d(f, 0.3, "Rtt/ttbar,qcd,Susy3,Susy4",          "TT_R");
  //save_1d(f, 0.8, "R2tt/ttbar,qcd,Susy3,Susy4",         "TT_R2");
  //save_1d(f, 0.4, "DeltaPhi/ttbar,qcd,Susy3,Susy4",     "TT_AbsDeltaPhi");
  //
  //make_comparison_plot(f, "HTevt_vs_DeltaPhi",        "HTevt_vs_DeltaPhi",        0.15);
  //make_comparison_plot(f, "HTevt_vs_R",               "HTevt_vs_R",               0.15);
  //make_comparison_plot(f, "HTevt_vs_Rtt",             "HTevt_vs_Rtt",             0.1);
  //make_comparison_plot(f, "HTttFraction_vs_DeltaPhi", "HTttFraction_vs_DeltaPhi", 0.3);
  //make_comparison_plot(f, "HTttFraction_vs_R",        "HTttFraction_vs_R",        0.25);
  //make_comparison_plot(f, "HTttFraction_vs_Rtt",      "HTttFraction_vs_Rtt",      0.25);
  //make_comparison_plot(f, "R_vs_DeltaPhi",            "R_vs_DeltaPhi",            0.25);
  //make_comparison_plot(f, "Rtt_vs_DeltaPhi",          "Rtt_vs_DeltaPhi",          0.25);

  // 01/17
  //save_1d(f, 0.6, "H_Tevt/ttbar,qcd,Susy3,Susy4",                 "HTevt");
  //save_1d(f, 0.7, "H_T/ttbar,qcd,Susy3,Susy4_SideBand",           "HT_SideBand");
  //save_1d(f, 0.7, "H_T/ttbar,qcd,Susy3,Susy4_Signal",             "HT_Signal");
  //save_1d(f, 0.4, "H_TttFraction/ttbar,qcd,Susy3,Susy4_SideBand", "HTttFraction_SideBand");
  //save_1d(f, 0.4, "H_TttFraction/ttbar,qcd,Susy3,Susy4_Signal",   "HTttFraction_Signal");
  //save_1d(f, 0.4, "H_TexFraction/ttbar,qcd,Susy3,Susy4_SideBand", "HTexFraction_SideBand");
  //save_1d(f, 0.4, "H_TexFraction/ttbar,qcd,Susy3,Susy4_Signal",   "HTexFraction_Signal");
  //save_1d(f, 0.6, "AK8JetR/ttbar,qcd,Susy3,Susy4_SideBand",       "R_SideBand");
  //save_1d(f, 0.6, "AK8JetR/ttbar,qcd,Susy3,Susy4_Signal",         "R_Signal");
  //save_1d(f, 1.0, "AK8JetR2/ttbar,qcd,Susy3,Susy4_SideBand",      "R2_SideBand");
  //save_1d(f, 1.0, "AK8JetR2/ttbar,qcd,Susy3,Susy4_Signal",        "R2_Signal");
  //save_1d(f, 0.6, "TT_R/ttbar,qcd,Susy3,Susy4_SideBand",          "Rtt_SideBand");
  //save_1d(f, 0.6, "TT_R/ttbar,qcd,Susy3,Susy4_Signal",            "Rtt_Signal");
  //save_1d(f, 1.0, "TT_R2/ttbar,qcd,Susy3,Susy4_SideBand",         "R2tt_SideBand");
  //save_1d(f, 1.0, "TT_R2/ttbar,qcd,Susy3,Susy4_Signal",           "R2tt_Signal");
  //save_1d(f, 0.8, "TT_AbsDeltaPhi/ttbar,qcd,Susy3,Susy4_SideBand","DeltaPhi_SideBand");
  //save_1d(f, 0.8, "TT_AbsDeltaPhi/ttbar,qcd,Susy3,Susy4_Signal",  "DeltaPhi_Signal");
  //
  //make_comparison_plot(f, "HTevt_vs_DeltaPhi",        "HTevt_vs_DeltaPhi",        0.2,  "");
  //make_comparison_plot(f, "HTevt_vs_R",               "HTevt_vs_R",               0.5,  "");
  //make_comparison_plot(f, "HTevt_vs_Rtt",             "HTevt_vs_Rtt",             0.25,  "");
  //make_comparison_plot(f, "HTttFraction_vs_DeltaPhi", "HTttFraction_vs_DeltaPhi_SideBand", 0.6,  "_SideBand");
  //make_comparison_plot(f, "HTttFraction_vs_DeltaPhi", "HTttFraction_vs_DeltaPhi_Signal",   0.6,  "_Signal");
  //make_comparison_plot(f, "HTttFraction_vs_R",        "HTttFraction_vs_R_SideBand",        0.25,  "_SideBand");
  //make_comparison_plot(f, "HTttFraction_vs_R",        "HTttFraction_vs_R_Signal",          0.25,  "_Signal");
  //make_comparison_plot(f, "HTttFraction_vs_Rtt",      "HTttFraction_vs_Rtt_SideBand",      0.25,  "_SideBand");
  //make_comparison_plot(f, "HTttFraction_vs_Rtt",      "HTttFraction_vs_Rtt_Signal",        0.25,  "_Signal");
  //make_comparison_plot(f, "R_vs_DeltaPhi",            "R_vs_DeltaPhi_SideBand",            0.4,  "_SideBand");
  //make_comparison_plot(f, "R_vs_DeltaPhi",            "R_vs_DeltaPhi_Signal",              0.4,  "_Signal");
  //make_comparison_plot(f, "Rtt_vs_DeltaPhi",          "Rtt_vs_DeltaPhi_SideBand",          0.4,  "_SideBand");
  //make_comparison_plot(f, "Rtt_vs_DeltaPhi",          "Rtt_vs_DeltaPhi_Signal",            0.4,  "_Signal");
  
  // 02/27
  save_1d(f, 1000000, "AK8JetR/AllSamples",       "R");
  save_1d(f, 1000000, "AK8JetR2/AllSamples",      "R2");
  save_1d(f, 1000000, "TT_R/AllSamples",          "Rtt");
  save_1d(f, 1000000, "TT_R2/AllSamples",         "R2tt");
  save_1d(f, 1000000, "TT_AbsDeltaPhi/AllSamples","DeltaPhi");
  save_1d(f, 1000000, "H_Tevt/AllSamples",        "HTevt");
  save_1d(f, 1000000, "H_T/AllSamples",           "HT");
  save_1d(f, 1000000, "H_TttFraction/AllSamples", "HTttFraction");
  save_1d(f, 1000000, "H_TexFraction/AllSamples", "HTexFraction");
}
