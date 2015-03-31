#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

void EstimateBg() {
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  
  std::vector<std::string> samples;
  samples.push_back("TTBar");
  //samples.push_back("WJets");
  samples.push_back("WJets_HT");
  samples.push_back("ZJets_HT");
  samples.push_back("QCD_HT");
  //samples.push_back("QCD_Pt_bcToE");
  samples.push_back("T_tW");
  //samples.push_back("TToLep_s_t");
  samples.push_back("DYJets_HT");
  //samples.push_back("DYJets");
  //samples.push_back("GJets_HT");
  
  //int rebin[] = { 2, 4, 5, 4, 2, 4};
  //double sideband_fit_low_range[] = { 0.14, 0.16, 0.10, 0.28, 0.12, 0.12 };
  //int rebin[] = { 5, 5, 5, 5, 5, 5};
  //double sideband_fit_low_range[] = { 0.15, 0.20, 0.10, 0.30, 0.10, 0.15 };
  //int rebin[] = { 5, 5, 5, 5, 5, 5};
  int rebin[] = { 2, 2, 2, 2, 2, 2};
  double sideband_fit_low_range[] = { 0.14, 0.14, 0.30, 0.13, 0.13, 0.13 };
  bool dphi_sideband[] = { 1, 0, 0, 1, 0, 0 };
  bool baderror = false;
  double weight[] = { 2.69454, 0.0505037, 0.00921411, 6.80717, 0.354934, 0.00484915 };
  
  //TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p4_DPhiBelow2p8_HtAllAbove0_TightLeptonVeto.root");
  //TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p4_DPhiBelow2p8_HtAllAbove0.root");
  TFile *f = TFile::Open("ROOT_output/BackGroundEstimate_NTopHadSideBand_NoLeptonVeto_50MassCut_fullstat.root");
  double Rranges_ABCD[3] = { 0, 0.4, 1.0 };
  //TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p35_DPhiBelow2p8_HtAllAbove0.root");
  //double Rranges_ABCD[3] = { 0, 0.35, 1.0 };
  bool doFitting = false;
  double sum_a = 0, sum_b = 0, sum_c = 0, sum_d = 0, sum_d_abcd = 0, sum_d_nevt = 0;
  double sum_a_err = 0, sum_b_err = 0, sum_c_err = 0, sum_d_err = 0, sum_d_abcd_err = 0;
  double sum_b_fit = 0, sum_d_fit = 0, sum_d_fit_comb = 0;
  double sum_b_fit_err = 0, sum_d_fit_err = 0, sum_d_fit_comb_err = 0;
  //printf("| *Sample* | *A (R<0.4, NTop<2)* | *B (R>0.4, NTop<2)* | *C (R<0.4, NTop==2)* | *D (R>0.4, NTop==2) obs.* | *D = B*C/A pred.* | *A->B (R-shape fit) pred.* | *C->D (R-shape fit) pred.* | *D pred (From sideband R-shape fit)* | \n");
  printf("| *Sample* | *A (R<0.4, Sideband)* | *B (R fit in Sideband) pred.* | *B (R>0.4, Sideband)* | *C (R<0.4, Signal band)* | *D (R fit in C) pred.* | *D = B*C/A pred.* | *D = B (R fit in Sideband) * C/A* | *D (R>0.4, Signal band) obs.* | *Nevent in D* | \n");
  for (size_t iSample = 0; iSample<samples.size(); ++iSample) {
    std::string canname = std::string(dphi_sideband[iSample] ? "RFine/DPhi2p8_2HadTop_" : "RFine/NTopHad_DPhiBelow2p8_")+samples[iSample];
    //std::string canname = std::string(dphi_sideband[iSample] ? "R/CutDPhi_2HadTop_" : "R/NTopHad_DPhiBelow2p8_")+samples[iSample];
    TCanvas *can = (TCanvas*)f->Get(canname.c_str()); can->Draw();
    TH1D *h_side = (TH1D*)can->GetListOfPrimitives()->At(dphi_sideband[iSample]); h_side->Rebin(rebin[iSample]);
    TH1D *h_signal = (TH1D*)can->GetListOfPrimitives()->At(1-dphi_sideband[iSample]); h_signal->Rebin(rebin[iSample]);
    TLegend *leg = (TLegend*)can->GetListOfPrimitives()->At(2);
    leg->SetX1(0.35); leg->SetX2(0.65); leg->SetY1(0.6);

    // Add ratio plot
    int y1 = 350;
    int y2 = 150;
    int mid2 = 10;
    //TH1D* h1 = (TH1D*)c->GetListOfPrimitives()->At(0);
    //TH1D* h2 = (TH1D*)c->GetListOfPrimitives()->At(1);
    //TLegend* leg = (TLegend*)c->GetListOfPrimitives()->At(2);
    //leg->SetX1(0.4);leg->SetX2(0.6);leg->SetY1(0.65);leg->SetY2(0.8);
    //h1->SetLabelSize(20.0/(y1+40),"xyz");
    //c = new TCanvas(canname.c_str(), c->GetTitle(), 604,626+(y1-500)+y2+mid2); // 600, 600
    can->Divide(1,2);
    // Pad 1 (80+500+20 x 40+500)
    TVirtualPad* p = can->cd(1);
    p->SetPad(0,(y2+60+mid2)/(y1+y2+100.0+mid2),1,1);
    p->SetTopMargin(40.0/(y1+40));
    p->SetBottomMargin(0);
    p->SetLogy(1);
    h_side->GetYaxis()->SetRangeUser(1.00001e-4,1e4);
    h_side->Draw("HIST");
    h_signal->Draw("SAMEHIST");
    leg->Draw();
    // Pad 2 (80+500+20 x 200+60)
    p = can->cd(2);
    p->SetGrid(0,1);
    p->SetPad(0,0,1,(y2+60+mid2)/(y1+y2+100.0+mid2));
    p->SetTopMargin(((float)mid2)/(y2+60+mid2));
    p->SetBottomMargin(60.0/(y2+60+mid2));
    TH1D* ratio = (TH1D*)h_side->Clone();
    TH1D* div = (TH1D*)h_signal->Clone();
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
    p = can->cd(1);
    
    // calculate integrals
    double integral[2][2] = { { 0, 0 }, { 0, 0 } };
    double integral_error[2][2] = { { 0, 0 }, { 0, 0 } };
    double nevt[2][2] = { { 0, 0 }, { 0, 0 } };
    for (int i=0; i<2; ++i) {
      for (int bin=1; bin<=h_side->GetNbinsX(); ++bin) {
        if (h_signal->GetXaxis()->GetBinLowEdge(bin)>=Rranges_ABCD[i] && 
            h_signal->GetXaxis()->GetBinUpEdge(bin)<Rranges_ABCD[i+1]) {
	  double c0 = h_side->GetBinContent(bin), c1 = h_signal->GetBinContent(bin);
	  double e0 = h_side->GetBinError(bin),   e1 = h_signal->GetBinError(bin);
	  //std::cout<<h_signal->GetBinError(bin)<<" "<<sqrt(c1*weight[iSample])<<std::endl;
	  if (baderror) {
	    e0 = sqrt(c0*weight[iSample]);
	    e1 = sqrt(c1*weight[iSample]);
	  }
	  nevt[0][i] += (int)(c0*c0/(e0*e0) + 0.5);
	  nevt[1][i] += (int)(c1*c1/(e1*e1) + 0.5);
          integral[0][i] += c0;
          integral[1][i] += c1;
          integral_error[0][i] += e0*e0;
          integral_error[1][i] += e1*e1;
	  //if (e1>0) std::cout<<bin<<" "<<c1<<" +- "<<e1*e1<<" nevt = "<<((int)(c1*c1/(e1*e1) + 0.5))<<std::endl;
        }
      }
      //if (iSample==1&&i==1) std::cout<<integral[1][i]<<" +- "<<integral_error[1][i]<<std::endl;
      integral_error[0][i] = sqrt(integral_error[0][i]);
      integral_error[1][i] = sqrt(integral_error[1][i]);
    }
    //std::cout<<nevt[1][1]<<std::endl;
    
    // predict yields using 2 methods (ABCD and constrained R-shape fit method combined)
    // ABCD method
    double a = integral[0][0], b = integral[0][1], c = integral[1][0], d = integral[1][1];
    double a_err = integral_error[0][0], b_err = integral_error[0][1], c_err = integral_error[1][0], d_err = integral_error[1][1];
    // Calculate error
    // z = x / y -> z_err = sqrt( (x*x*y_err*y_err + y*y*x_err*x_err)/(y*y*y*y) )
    // z = x * y -> z_err = sqrt ( x*x*y_err*y_err + y*y*x_err*x_err )
    double c_per_a_err = sqrt((c*c*a_err*a_err + a*a*c_err*c_err)/(a*a*a*a));
    double d_abcd = b * (c/a), d_abcd_err = sqrt(b*b*c_per_a_err*c_per_a_err + (c/a)*(c/a)*b_err*b_err);
    double d_nevt = nevt[1][1];
    
    // Fit in the full range of NTop Sideband
    // Do fitting and calculate integrals
    TF1 *fit_side = new TF1("NTopSide_fit","exp([0]+[1]*x)", sideband_fit_low_range[iSample], Rranges_ABCD[2]);
    fit_side->SetLineColor(h_side->GetLineColor());
    h_side->Fit("NTopSide_fit","QRE");
    fit_side->Draw("SAME");
    double fit_integral[2][2], fit_integral_error[2][2], d_fit_comb = 0, d_fit_comb_err = 0;
    //double Rranges_ACfit[3] = { samples[iSample]=="ZJets_HT" ? 0.2 : 0.15, Rranges_ABCD[1], Rranges_ABCD[2] };
    double Rranges_ACfit[3] = { sideband_fit_low_range[iSample], Rranges_ABCD[1], Rranges_ABCD[2] };
    for (int i=0; i<2; ++i) {
      fit_integral[0][i] = fit_side->Integral(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
      fit_integral_error[0][i] = fit_side->IntegralError(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
    }
    double par0 = fit_side->GetParameter(0), par0_error = fit_side->GetParError(0);
    double par1 = fit_side->GetParameter(1), par1_error = fit_side->GetParError(1);
    double par1min, par1max; fit_side->GetParLimits(1, par1min, par1max);
    
    // Razor 2D fit
    //+ TCanvas *can2 = new TCanvas("can2","2D Razor fit",1200,600);
    //+ can2->Divide(2); can2->cd(1);
    //+ TH2D *h2_side = (TH2D*)f->Get((std::string("RCoarse_vs_MR/DPhiBelow2p8_0To1HadTop_")+samples[iSample]).c_str());
    //+ h2_side->Draw("COLZ");
    //+ 
    //+ //TF1 *Razor_fit = new TF1("Razor_fit","([0]*((x*x-[2])^(1/[1]))-1)*exp(-[0]*[1]*((x*x-[2])^(1/[1])))", Rranges_ABCD[0], Rranges_ABCD[2]);
    //+ TF2 *Razor_fit = new TF2("Razor_fit","([0]*((x-[2])^(1/[1]))*((y*y-[3])^(1/[1]))-1)*exp(-[0]*[1]*((x-[2])^(1/[1]))*((y*y-[3])^(1/[1])))",0,3500,0,1.25);
    //+ Razor_fit->SetParName(0, "b");
    //+ Razor_fit->SetParName(1, "n");
    //+ Razor_fit->SetParName(2, "MR0");
    //+ Razor_fit->SetParName(3, "R20");
    //+ Razor_fit->SetParameter(0,0.01);
    //+ Razor_fit->SetParameter(1,1.1);
    //+ Razor_fit->SetParameter(2,-1350);
    //+ Razor_fit->SetParameter(3,-0.01);
    //+ h2_side->Fit("Razor_fit","RM0");
    //+ Razor_fit->Draw("cont1 same");
    //+ can2->cd(2);
    //+ h2_side = new TH2D("hehe","",25,0,2500,25,0,1.25);
    //+ h2_side->FillRandom("Razor_fit",100000);
    //+ h2_side->Draw("COLZ");
    
    // Fit in the Signal region
    // Fitting in sideband, get B area under curve and scale by C/A
    TF1 *fit_signal = new TF1("NTopSignal_RSide_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ABCD[2]);
    fit_signal->SetLineColor(h_signal->GetLineColor());
    //fit_signal->SetParameter(1, par1); 
    //fit_signal->SetParLimits(1, par1min, par1max);
    h_signal->Fit("NTopSignal_RSide_fit","QREB");
    fit_signal->Draw("SAME");
    for (int i=0; i<2; ++i) {
      fit_integral[1][i] = fit_signal->Integral(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
      fit_integral_error[1][i] = fit_signal->IntegralError(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
    }
    d_fit_comb = fit_integral[0][1] * (c/a);
    d_fit_comb_err = sqrt(fit_integral[0][1]*fit_integral[0][1]*c_per_a_err*c_per_a_err + (c/a)*(c/a)*fit_integral_error[0][1]*fit_integral_error[0][1]);
    TF1 *fit_pred = new TF1("Predicted_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ACfit[2]);
    fit_pred->SetLineColor(1);
    fit_pred->SetLineStyle(2);
    fit_pred->FixParameter(0, par0+std::log(c/a)); 
    fit_pred->FixParameter(1, par1); 
    h_signal->Fit("Predicted_fit","QREB+");
    fit_pred->Draw("SAME");
    leg->AddEntry(fit_pred, "N_{top,hadronic}<2 fit scaled (ABCD)", "l");
    
    // Constrained R fit method, using fixed slope of fit in NTop sideband
    // Free parameter is the height - Do not really work, due to limited statistics
    // Full A area needs to be used
    //TF1 *fit_signal_constr = new TF1("NTopSignal_RSide_Constr_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ACfit[1]);
    //fit_signal_constr->SetLineStyle(2);
    //fit_signal_constr->FixParameter(1, par1); 
    //fit_signal_constr->SetParError(1, par1_error); 
    //h_signal->Fit("NTopSignal_RSide_Constr_fit","QREB+");
    //d_fit_comb = fit_signal_constr->Integral(Rranges_ACfit[1], Rranges_ACfit[2])/h_signal->GetXaxis()->GetBinWidth(1);
    //d_fit_comb_err = fit_signal_constr->IntegralError(Rranges_ACfit[1], Rranges_ACfit[2])/h_signal->GetXaxis()->GetBinWidth(1);
    
    printf("| %s |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |", samples[iSample].c_str(), a, a_err, fit_integral[0][1], fit_integral_error[0][1], b, b_err, c, c_err);
    printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %d |\n", fit_integral[1][1], fit_integral_error[1][1], d_abcd, d_abcd_err, d_fit_comb, d_fit_comb_err, d, d_err, d_nevt);
    sum_a += a; sum_b += b; sum_c += c; sum_d += d;
    sum_a_err += a_err*a_err; sum_b_err += b_err*b_err; sum_c_err += c_err*c_err; sum_d_err += d_err*d_err;
    sum_d_abcd += d_abcd; sum_d_abcd_err += d_abcd_err*d_abcd_err;
    sum_b_fit += fit_integral[0][1]; sum_b_fit_err += fit_integral_error[0][1]*fit_integral_error[0][1];
    sum_d_fit += fit_integral[1][1]; sum_d_fit_err += fit_integral_error[1][1]*fit_integral_error[1][1];
    sum_d_fit_comb += d_fit_comb; sum_d_fit_comb_err += d_fit_comb_err*d_fit_comb_err;
    sum_d_nevt += d_nevt;
    
    //std::cout<<samples[iSample]<<std::endl;
    //std::cout<<"Intagrals (NTop sideband,   R sideband)   - MC-sum: "<<integral[0][0]<<" +- "<<integral_error[0][0]<<" Fit: "<<fit_integral[0][0]<<" +- "<<fit_integral_error[0][0]<<std::endl;
    //std::cout<<"Intagrals (NTop sideband,   R signalband) - MC-sum: "<<integral[0][1]<<" +- "<<integral_error[0][1]<<" Fit: "<<fit_integral[0][1]<<" +- "<<fit_integral_error[0][1]<<std::endl;
    //std::cout<<"Intagrals (NTop signalband, R sideband)   - MC-sum: "<<integral[1][0]<<" +- "<<integral_error[1][0]<<" Fit: "<<fit_integral[1][0]<<" +- "<<fit_integral_error[1][0]<<std::endl;
    //std::cout<<"Intagrals (NTop signalband, R signalband) - MC-sum: "<<integral[1][1]<<" +- "<<integral_error[1][1]<<" Fit: "<<fit_integral[1][1]<<" +- "<<fit_integral_error[1][1]<<std::endl;
    //can->SaveAs((std::string("/afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_03_29/")+can->GetName()+".png").c_str());
  }
  sum_a_err = sqrt(sum_a_err); sum_b_err = sqrt(sum_b_err); sum_c_err = sqrt(sum_c_err); sum_d_err = sqrt(sum_d_err);
  sum_b_fit_err = sqrt(sum_b_fit_err); sum_d_fit_err = sqrt(sum_d_fit_err); sum_d_fit_comb_err = sqrt(sum_d_fit_comb_err);
  printf("| SUM |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |",
	 sum_a, sum_a_err, sum_b_fit, sum_b_fit_err, sum_b, sum_b_err, sum_c, sum_c_err);
  printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %d |\n", sum_d_fit, sum_d_fit_err, sum_d_abcd, sum_d_abcd_err, sum_d_fit_comb, sum_d_fit_comb_err, sum_d, sum_d_err, sum_d_nevt); // R-shape fit and ABCD combined result
}

