#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

void do_plot(TCanvas* c, int n, TFile *f, const char* name, const char* title, float limit) {
  TVirtualPad* p = c->cd(n);
  p->SetLeftMargin(0.15);
  p->SetRightMargin(0.2);
  p->SetLogz(1);
  TH2D* h = ((TH2D*)f->Get(name));
  h->SetTitle(title);
  h->GetYaxis()->SetTitleOffset(1.2);
  h = (TH2D*)h->DrawNormalized("COLZ");
  if (limit>0) h->GetZaxis()->SetRangeUser(0.0001,limit);
}


void make_comparison_plot(TFile* f, std::string canname, std::string name, float limit=0) {
  TCanvas* c = new TCanvas(canname.c_str(),name.c_str(), 900,900);
  c->SetLeftMargin(0.2);
  c->Divide(2,2);
  do_plot(c, 1, f, (name+"/ttbar").c_str(), "ttbar", limit);
  do_plot(c, 2, f, (name+"/susy3body").c_str(), "T5tttt 3-body decay", limit);
  do_plot(c, 3, f, (name+"/qcd").c_str(), "QCD", limit);
  do_plot(c, 4, f, (name+"/susy4body").c_str(), "T5tttt 4-body decay", limit);
  c->SaveAs((std::string("Plots/2015_01_14/Correction/Corr2D/")+canname+".png").c_str());
}

void save_1d(TFile* f, float limit, std::string canname, std::string name) {
  TCanvas* c = ((TCanvas*)f->Get((name+"/ttbar,qcd,Susy3,Susy4").c_str()));
  TH1D* h = (TH1D*)c->GetListOfPrimitives()->At(0);
  if (limit>0) h->GetYaxis()->SetRangeUser(0,limit);
  //c->SetLogy(1);
  c->Draw();
  c->SaveAs((std::string("Plots/2015_01_14/Correction/")+canname+".png").c_str());
}

void CorrelationPlots() {
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("plots.root");
  
  make_comparison_plot(f, "HTevt_vs_DeltaPhi",        "HTevt_vs_DeltaPhi",        0.15);
  make_comparison_plot(f, "HTevt_vs_R",               "HTevt_vs_R",               0.15);
  make_comparison_plot(f, "HTevt_vs_Rtt",             "HTevt_vs_Rtt",             0.1);
  make_comparison_plot(f, "HTttFraction_vs_DeltaPhi", "HTttFraction_vs_DeltaPhi", 0.3);
  make_comparison_plot(f, "HTttFraction_vs_R",        "HTttFraction_vs_R",        0.25);
  make_comparison_plot(f, "HTttFraction_vs_Rtt",      "HTttFraction_vs_Rtt",      0.25);
  make_comparison_plot(f, "R_vs_DeltaPhi",            "R_vs_DeltaPhi",            0.25);
  make_comparison_plot(f, "Rtt_vs_DeltaPhi",          "Rtt_vs_DeltaPhi",          0.25);
  
  //make_comparison_plot(f,"DefR_vs_HT",         "ak8jetR_vs_ht",                     0.0);
  //make_comparison_plot(f,"DefR_vs_HTextra",    "ak8jetR_vs_ht_extra_per_ht",        0.0);
  //make_comparison_plot(f,"TTR_vs_HT",          "tt_R_vs_ht",                        0.0);
  //make_comparison_plot(f,"TTR_vs_HTextra",     "tt_R_vs_ht_extra_per_ht",           0.0);
  //make_comparison_plot(f,"DeltaPhi_vs_HT",     "tt_AbsDeltaPhi_vs_ht",              0.0);
  //make_comparison_plot(f,"DeltaPhi_vs_HTextra","tt_AbsDeltaPhi_vs_ht_extra_per_ht", 0.0);
  //make_comparison_plot(f,"DeltaPhi_vs_DefR",   "tt_AbsDeltaPhi_vs_ak8jetR",         0.0);
  //make_comparison_plot(f,"DeltaPhi_vs_TTR",    "tt_AbsDeltaPhi_vs_tt_R",            0.0);
  
  save_1d(f, 0.1, "HT",           "H_T");
  save_1d(f, 0.1, "HTevt",        "H_Tevt");
  save_1d(f, 0.3, "HTttFraction", "H_TttFraction");
  save_1d(f, 0.3, "HTexFraction", "H_TexFraction");
  save_1d(f, 0.2, "R",            "AK8JetR");
  save_1d(f, 0.7, "R2",           "AK8JetR2");
  save_1d(f, 0.3, "Rtt",          "TT_R");
  save_1d(f, 0.8, "R2tt",         "TT_R2");
  save_1d(f, 0.4, "DeltaPhi",     "TT_AbsDeltaPhi");
  
}
