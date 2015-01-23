#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

void make_ratio_plot(TFile* f, std::string name, std::string canname, float limit) {
  int y1 = 350;
  int y2 = 150;
  int mid2 = 10;
  TCanvas* c = ((TCanvas*)f->Get(name.c_str()));
  if (c->GetListOfPrimitives()->GetEntries()>2) {
    TH1D* h1 = (TH1D*)c->GetListOfPrimitives()->At(0);
    TH1D* h2 = (TH1D*)c->GetListOfPrimitives()->At(1);
    TLegend* leg = (TLegend*)c->GetListOfPrimitives()->At(2);
    leg->SetX1(0.4);leg->SetX2(0.6);leg->SetY1(0.65);leg->SetY2(0.8);
    if (limit>0) h1->GetYaxis()->SetRangeUser(0.0002,limit);
    h1->SetLabelSize(20.0/(y1+40),"xyz");
    c = new TCanvas(canname.c_str(), c->GetTitle(), 604,626+(y1-500)+y2+mid2); // 600, 600
    c->Divide(1,2);
    // Pad 1 (80+500+20 x 40+500)
    TVirtualPad* p = c->cd(1);
    p->SetPad(0,(y2+60+mid2)/(y1+y2+100.0+mid2),1,1);
    p->SetTopMargin(40.0/(y1+40));
    p->SetBottomMargin(0);
    p->SetLogy(1);
    h1->Draw("SAMEHIST");
    h2->Draw("SAMEHIST");
    leg->Draw();
    // Pad 2 (80+500+20 x 200+60)
    p = c->cd(2);
    p->SetGrid(0,1);
    p->SetPad(0,0,1,(y2+60+mid2)/(y1+y2+100.0+mid2));
    p->SetTopMargin(((float)mid2)/(y2+60+mid2));
    p->SetBottomMargin(60.0/(y2+60+mid2));
    TH1D* ratio = (TH1D*)h1->Clone();
    TH1D* div = (TH1D*)h2->Clone();
    ratio->Scale(1/ratio->GetSumOfWeights());
    ratio->SetTitleSize(32.0/(y2+60+mid2),"xyz");
    ratio->SetLabelSize(20.0/(y2+60+mid2),"xyz");
    div->Scale(1/div->GetSumOfWeights());
    ratio->Divide(div);
    //ratio->GetYaxis()->SetRangeUser(0,2);
    ratio->GetYaxis()->SetNdivisions(305);
    ratio->GetYaxis()->SetTitle("Ratio");
    ratio->GetYaxis()->SetTitleOffset(0.45);
    ratio->SetTitleSize(24.0/(y2+60+mid2),"y");
    ratio->SetTitle("");
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(1);
    ratio->SetLineColor(1);
    ratio->Draw("PE1");
    TLine* l = new TLine(ratio->GetXaxis()->GetXmin(), 1, ratio->GetXaxis()->GetXmax(), 1);
    l->SetLineWidth(2);
    //l->SetLineColor(2);
    l->SetLineStyle(2);
    l->Draw();
  } else {
    c->Draw();
  }
  c->SaveAs((std::string("Plots/2015_01_20/")+canname+".png").c_str());
}

void RatioPlots() {
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("plots.root");
  
  // 01/20
  make_ratio_plot(f, "TT_AbsDeltaPhi/RBelow0p25,RAbove0p25",          "DeltaPhi_RBins", 1);
  make_ratio_plot(f, "TT_AbsDeltaPhi/RBelow0p25,RAbove0p25_SideBand", "DeltaPhi_RBins_SideBand", 1);
  make_ratio_plot(f, "TT_AbsDeltaPhi/RBelow0p25,RAbove0p25_Signal",   "DeltaPhi_RBins_Signal", 1);
  make_ratio_plot(f, "AK8JetR/DPhiBelow2p8,DPhiAbove2p8",             "R_DPhiBins", 1);
  make_ratio_plot(f, "AK8JetR/DPhiBelow2p8,DPhiAbove2p8_SideBand",    "R_DPhiBins_SideBand", 1);
  make_ratio_plot(f, "AK8JetR/DPhiBelow2p8,DPhiAbove2p8_Signal",      "R_DPhiBins_Signal", 1);
  make_ratio_plot(f, "AK8JetR/SideBand,Signal_DPhiBelow2p8",          "R_HTallBins_DPhiBelow2p8", 1);
  make_ratio_plot(f, "AK8JetR/SideBand,Signal_DPhiAbove2p8",          "R_HTallBins_DPhiAbove2p8", 1);
}

