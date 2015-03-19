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
  samples.push_back("QCD_HT");
  //samples.push_back("QCD_Pt_bcToE");
  samples.push_back("TTBar");
  //samples.push_back("GJets_HT");
  //samples.push_back("WJets");
  samples.push_back("WJets_HT");
  samples.push_back("T_tW");
  //samples.push_back("TToLep_s_t");
  samples.push_back("ZJets_HT");
  //samples.push_back("DYJets");
  samples.push_back("DYJets_HT");
  
  //int rebin[] = { 2, 2, 4, 5, 4, 4};
  //double sideband_fit_low_range[] = { 0.12, 0.14, 0.16, 0.10, 0.28, 0.12 };
  //int rebin[] = { 5, 5, 5, 5, 5, 5};
  //double sideband_fit_low_range[] = { 0.10, 0.15, 0.20, 0.10, 0.30, 0.15 };
  int rebin[] = { 1, 1, 1, 1, 1, 1};
  //int rebin[] = { 5, 5, 5, 5, 5, 5};
  double sideband_fit_low_range[] = { 0.15, 0.15, 0.20, 0.15, 0.15, 0.15 };
  
  //TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p4_DPhiBelow2p8_HtAllAbove0_TightLeptonVeto.root");
  TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p4_DPhiBelow2p8_HtAllAbove0.root");
  double Rranges_ABCD[3] = { 0, 0.4, 1.0 };
  //TFile *f = TFile::Open("ROOT_output/SignalCutPlots_RAbove0p35_DPhiBelow2p8_HtAllAbove0.root");
  //double Rranges_ABCD[3] = { 0, 0.35, 1.0 };
  bool doFitting = false;
  double sum_a = 0, sum_b = 0, sum_c = 0, sum_d = 0, sum_d_abcd = 0;
  double sum_a_err = 0, sum_b_err = 0, sum_c_err = 0, sum_d_err = 0, sum_d_abcd_err = 0;
  double sum_b_fit = 0, sum_d_fit = 0, sum_d_fit_comb = 0;
  double sum_b_fit_err = 0, sum_d_fit_err = 0, sum_d_fit_comb_err = 0;
  //printf("| *Sample* | *A (R<0.4, NTop<2)* | *B (R>0.4, NTop<2)* | *C (R<0.4, NTop==2)* | *D (R>0.4, NTop==2) obs.* | *D = B*C/A pred.* | *A->B (R-shape fit) pred.* | *C->D (R-shape fit) pred.* | *D pred (From sideband R-shape fit)* | \n");
  printf("| *Sample* | *A (R<0.4, NTop<2)* | *B (R fit in AB) pred.* | *B (R>0.4, NTop<2)* | *C (R<0.4, NTop==2)* | *D (R fit in C) pred.* | *D = B*C/A pred.* | *D = B (R fit in AB) * C/A* | *D (R>0.4, NTop==2) obs.* | \n");
  for (size_t iSample = 0; iSample<samples.size(); ++iSample) {
    TCanvas* can = ((TCanvas*)f->Get((std::string("R/NTopHad_DPhiBelow2p8_")+samples[iSample]).c_str())); can->Draw();
    TH1D* h_side = (TH1D*)can->GetListOfPrimitives()->At(0); h_side->Rebin(rebin[iSample]);
    TH1D* h_signal = (TH1D*)can->GetListOfPrimitives()->At(1); h_signal->Rebin(rebin[iSample]);
    TLegend* leg = (TLegend*)can->GetListOfPrimitives()->At(2);
    leg->SetX1(0.4); leg->SetX2(0.6); leg->SetY1(0.7);
    
    // calculate integrals
    double integral[2][2] = { { 0, 0 }, { 0, 0 } };
    double integral_error[2][2] = { { 0, 0 }, { 0, 0 } };
    for (int i=0; i<2; ++i) {
      for (int bin=1; bin<=h_side->GetNbinsX(); ++bin) {
        if (h_signal->GetXaxis()->GetBinLowEdge(bin)>=Rranges_ABCD[i] && 
            h_signal->GetXaxis()->GetBinUpEdge(bin)<Rranges_ABCD[i+1]) {
	  double c0 = h_side->GetBinContent(bin), c1 = h_signal->GetBinContent(bin);
	  double e0 = h_side->GetBinError(bin),   e1 = h_signal->GetBinError(bin);
          integral[0][i] += c0;
          integral[1][i] += c1;
          integral_error[0][i] += e0*e0;
          integral_error[1][i] += e1*e1;
        }
      }
      integral_error[0][i] = sqrt(integral_error[0][i]);
      integral_error[1][i] = sqrt(integral_error[1][i]);
    }
    
    // predict yields using 2 methods (ABCD and constrained R-shape fit method combined)
    // ABCD method
    double a = integral[0][0], b = integral[0][1], c = integral[1][0], d = integral[1][1];
    double a_err = integral_error[0][0], b_err = integral_error[0][1], c_err = integral_error[1][0], d_err = integral_error[1][1];
    // Calculate error
    // z = x / y -> z_err = sqrt( (x*x*y_err*y_err + y*y*x_err*x_err)/(y*y*y*y) )
    // z = x * y -> z_err = sqrt ( x*x*y_err*y_err + y*y*x_err*x_err )
    double c_per_a_err = sqrt((c*c*a_err*a_err + a*a*c_err*c_err)/(a*a*a*a));
    double d_abcd = b * (c/a), d_abcd_err = sqrt(b*b*c_per_a_err*c_per_a_err + (c/a)*(c/a)*b_err*b_err);
    
    // Fit in the full range of NTop Sideband
    // Do fitting and calculate integrals
    TF1* fit_side = new TF1("NTopSide_fit","exp([0]+[1]*x)", sideband_fit_low_range[iSample], Rranges_ABCD[2]);
    fit_side->SetLineColor(4);
    h_side->Fit("NTopSide_fit","QRE");
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
    
    // Fit in the Signal region
    // Fitting in sideband, get B area under curve and scale by C/A
    TF1* fit_signal = new TF1("NTopSignal_RSide_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ABCD[2]);
    //fit_signal->SetParameter(1, par1); 
    //fit_signal->SetParLimits(1, par1min, par1max);
    h_signal->Fit("NTopSignal_RSide_fit","QREB");
    for (int i=0; i<2; ++i) {
      fit_integral[1][i] = fit_signal->Integral(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
      fit_integral_error[1][i] = fit_signal->IntegralError(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
    }
    d_fit_comb = fit_integral[0][1] * (c/a);
    d_fit_comb_err = sqrt(fit_integral[0][1]*fit_integral[0][1]*c_per_a_err*c_per_a_err + (c/a)*(c/a)*fit_integral_error[0][1]*fit_integral_error[0][1]);
    TF1* fit_pred = new TF1("Predicted_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ACfit[2]);
    fit_pred->SetLineColor(1);
    fit_pred->SetLineStyle(2);
    fit_pred->FixParameter(0, par0+std::log(c/a)); 
    fit_pred->FixParameter(1, par1); 
    h_signal->Fit("Predicted_fit","QREB+");
    leg->AddEntry(fit_pred, "N_{top,hadronic}<2 fit scaled (ABCD)", "l");
    
    // Constrained R fit method, using fixed slope of fit in NTop sideband
    // Free parameter is the height - Do not really work, due to limited statistics
    // Full A area needs to be used
    //TF1* fit_signal_constr = new TF1("NTopSignal_RSide_Constr_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ACfit[1]);
    //fit_signal_constr->SetLineStyle(2);
    //fit_signal_constr->FixParameter(1, par1); 
    //fit_signal_constr->SetParError(1, par1_error); 
    //h_signal->Fit("NTopSignal_RSide_Constr_fit","QREB+");
    //d_fit_comb = fit_signal_constr->Integral(Rranges_ACfit[1], Rranges_ACfit[2])/h_signal->GetXaxis()->GetBinWidth(1);
    //d_fit_comb_err = fit_signal_constr->IntegralError(Rranges_ACfit[1], Rranges_ACfit[2])/h_signal->GetXaxis()->GetBinWidth(1);
    
    printf("| %s |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |", samples[iSample].c_str(), a, a_err, fit_integral[0][1], fit_integral_error[0][1], b, b_err, c, c_err);
    printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |\n", fit_integral[1][1], fit_integral_error[1][1], d_abcd, d_abcd_err, d_fit_comb, d_fit_comb_err, d, d_err);
    sum_a += a; sum_b += b; sum_c += c; sum_d += d;
    sum_a_err += a_err*a_err; sum_b_err += b_err*b_err; sum_c_err += c_err*c_err; sum_d_err += d_err*d_err;
    sum_d_abcd += d_abcd; sum_d_abcd_err += d_abcd_err*d_abcd_err;
    sum_b_fit += fit_integral[0][1]; sum_b_fit_err += fit_integral_error[0][1]*fit_integral_error[0][1];
    sum_d_fit += fit_integral[1][1]; sum_d_fit_err += fit_integral_error[1][1]*fit_integral_error[1][1];
    sum_d_fit_comb += d_fit_comb; sum_d_fit_comb_err += d_fit_comb_err*d_fit_comb_err;
    
    //std::cout<<samples[iSample]<<std::endl;
    //std::cout<<"Intagrals (NTop sideband,   R sideband)   - MC-sum: "<<integral[0][0]<<" +- "<<integral_error[0][0]<<" Fit: "<<fit_integral[0][0]<<" +- "<<fit_integral_error[0][0]<<std::endl;
    //std::cout<<"Intagrals (NTop sideband,   R signalband) - MC-sum: "<<integral[0][1]<<" +- "<<integral_error[0][1]<<" Fit: "<<fit_integral[0][1]<<" +- "<<fit_integral_error[0][1]<<std::endl;
    //std::cout<<"Intagrals (NTop signalband, R sideband)   - MC-sum: "<<integral[1][0]<<" +- "<<integral_error[1][0]<<" Fit: "<<fit_integral[1][0]<<" +- "<<fit_integral_error[1][0]<<std::endl;
    //std::cout<<"Intagrals (NTop signalband, R signalband) - MC-sum: "<<integral[1][1]<<" +- "<<integral_error[1][1]<<" Fit: "<<fit_integral[1][1]<<" +- "<<fit_integral_error[1][1]<<std::endl;
    //can->SaveAs((std::string("/afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_03_14/")+can->GetName()+".png").c_str());
  }
  sum_a_err = sqrt(sum_a_err); sum_b_err = sqrt(sum_b_err); sum_c_err = sqrt(sum_c_err); sum_d_err = sqrt(sum_d_err);
  sum_b_fit_err = sqrt(sum_b_fit_err); sum_d_fit_err = sqrt(sum_d_fit_err); sum_d_fit_comb_err = sqrt(sum_d_fit_comb_err);
  printf("| SUM |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |",
	 sum_a, sum_a_err, sum_b_fit, sum_b_fit_err, sum_b, sum_b_err, sum_c, sum_c_err);
  printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |\n", sum_d_fit, sum_d_fit_err, sum_d_abcd, sum_d_abcd_err, sum_d_fit_comb, sum_d_fit_comb_err, sum_d, sum_d_err); // R-shape fit and ABCD combined result
}

