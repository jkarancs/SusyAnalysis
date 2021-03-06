#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

void EstimateBg() {
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  bool latex = false;
  bool RunIIMC = true;
  
  //const char* filename = RunIIMC ? "ROOT_output/RunIIMC_Bands_OldTopDef_fullstat.root" : "ROOT_output/NotePlots_AllSamples_fullstat.root";
  const char* filename = RunIIMC ? "ROOT_output/RunIIMC_Bands_NewTopDef_fullstat.root" : "ROOT_output/NotePlots_AllSamples_fullstat.root";
  
  std::vector<std::string> samples[4];
  if (RunIIMC) {
    //samples[0].push_back("TTJets_NLO");
    samples[0].push_back("TTJets_LO_madgraph");
    //samples[0].push_back("TTJets_LO_pythia8");
    
    //samples[1].push_back("TTJets_NLO");
    samples[1].push_back("TTJets_LO_madgraph");
    //samples[1].push_back("TTJets_LO_pythia8");
    samples[1].push_back("WJets");
    samples[1].push_back("ST");
    //samples[1].push_back("Diboson");
    
    //samples[2].push_back("All Bkg.");
  } else {
    // PHYS14 Samples
    // DPhi sideband - Better estimate for TTbar
    samples[0].push_back("TTJets");
    //samples[0].push_back("TT");
    
    // Main method - NTop Sideband - separately
    samples[1].push_back("TTJets");
    //samples[1].push_back("TT");
    //samples[1].push_back("WJets");
    samples[1].push_back("WJets");
    samples[1].push_back("ZJets");
    samples[1].push_back("QCD");
    //samples[1].push_back("QCD_Pt_bcToE");
    samples[1].push_back("T_tW");
    //samples[1].push_back("TToLep_s_t");
    //samples[1].push_back("DYJets_HT");
    //samples[1].push_back("DYJets");
    //samples[1].push_back("GJets_HT");
    
    // NTop Sideband All background summed
    samples[2].push_back("All Bkg.");
    
    // Signal in NTop bins
    samples[3].push_back("T5ttttDeg_4bodydec_mGo1300");
    samples[3].push_back("T5ttttDeg_3bodydec_mGo1300");
    samples[3].push_back("T5ttttDeg_4bodydec_mGo1000");
    samples[3].push_back("T5ttttDeg_3bodydec_mGo1000");
    samples[3].push_back("T1tttt_mGl1500_mLSP100");
  }
  
  bool baderror = false;
  double weight[] = { 0.32686, 0.0505037, 0.00921411, 6.80717, 0.354934, 0.00484915 };
  int rebin[] = { 10, 10, 10, 10, 10 };
  double sideband_fit_low_range[] = { 0.15, 0.15, 0.15, 0.15, 0.10, 0.15 };
  bool dphi_sideband[] = { 0, 0, 0, 0, 0, 0 };
  
  int i_h_side[4]   = { 2, 0, 0, 0 };
  int i_h_signal[4] = { 3, 1, 1, 1 };
  
  Double_t Rranges_ABCD[4][4] = 
    { { 2.8, 3.2, 0.0, 2.8  },
      { 0.0, 0.4, 0.4 - 1e-10, 1.20 },
      { 0.0, 0.4, 0.4 - 1e-10, 1.20 },
      { 0.0, 0.4, 0.4 - 1e-10, 1.20 } };
  bool doFitting = false;
  double sum_a = 0, sum_b = 0, sum_c = 0, sum_d = 0, sum_d_abcd = 0, sum_d_nevt = 0;
  double sum_a_err = 0, sum_b_err = 0, sum_c_err = 0, sum_d_err = 0, sum_d_abcd_err = 0;
  double sum_b_fit = 0, sum_d_fit = 0, sum_d_fit_comb = 0;
  double sum_b_fit_err = 0, sum_d_fit_err = 0, sum_d_fit_comb_err = 0;
  double comb_d = 0, comb_d_err = 0, comb_d_abcd = 0, comb_d_abcd_err = 0;
  if (latex) {
    printf("\\begin{table*}[htbH]\n");
    printf("\\small\n");
    printf("\\begin{center}\n");
    printf("\\topcaption{Estimated Standard Model background yields in ABCD regions. A, B is in the sideband, C and D is the signal band, D is the signal region.\\label{tab:SMBkgEstimate}}\n");
    printf("\\begin{tabular}{lrrrrrrrr}\n");
  }
  for (size_t iMethod = 0; iMethod<4; ++iMethod) {
    TFile *f = TFile::Open(filename);
    if (latex) {
      if (iMethod==0) {
	printf("\\hline\n");
	printf("Method 2 & A & B & C & D = B*C/A & D obs. & Ratio pred./obs. & Pull & KS test\\\\\n");
      }	else if (iMethod==1){
	printf("\\hline\n");
	printf("Method 1 & A & B & C & D = B*C/A & D obs. & Ratio pred./obs. & Pull & KS test\\\\\n");
      }
      printf("\\hline\n");
    } else {
      //if (iMethod==0) printf("| *Sample* | *A (DPhi>2.8, SB)* | *B (DPhi<2.8, SB)* | *C (DPhi>2.8, Sig.B.)* | - | *D = B*C/A pred.* | *D (DPhi<2.8, Sig.B.) obs.* | *Ratio pred./obs.* |\n");
      //else if (iMethod==1) printf("| *Sample* | *A (R<0.4, SB)* | *B (R>0.4, SB)* | *C (R<0.4, Sig.B.)* | *D = B (R fit, SB) * C/A pred.* | *D = B*C/A pred.* | *D (R>0.4, Sig.B.) obs.* | *Ratio pred./obs.* | \n");
      if (iMethod==0) printf("| *Sample* | *A (DPhi>2.8, SB)* | *B (DPhi<2.8, SB)* | *C (DPhi>2.8, Sig.B.)* | *D = B*C/A pred.* | *D (DPhi<2.8, Sig.B.) obs.* | *Ratio pred./obs.* | *Pull (pred-obs)/error* | *KS test* |\n");
      else if (iMethod==1) printf("| *Sample* | *A (R<0.4, SB)* | *B (R>0.4, SB)* | *C (R<0.4, Sig.B.)* | *D = B*C/A pred.* | *D (R>0.4, Sig.B.) obs.* | *Ratio pred./obs.* | *Pull (pred-obs)/error* | *KS test* |\n");
    }
    for (size_t iSample = 0; iSample<samples[iMethod].size(); ++iSample) {
      std::string canname = 
	iMethod==0 ? std::string("DPhiBins/RBins_2HadTop_")+samples[iMethod][iSample] :
	iMethod==1 ? (std::string(dphi_sideband[iSample] ? "RFine/DPhi2p8_2HadTop_" : "RFine/NTopHad_DPhiBelow2p8_")+samples[iMethod][iSample]) :
	iMethod==2 ? (dphi_sideband[iSample] ? "RBins/CutDPhi_2HadTop" : "RBins/NTopHad_DPhiBelow2p8") : 
	iMethod==3 ? (std::string("RBins/NTopHad_DPhiBelow2p8_")+samples[iMethod][iSample]) : "";
      if (RunIIMC) {
        canname = 
          iMethod==0 ? std::string("DPhiBinsNew/RBins_2HadTop_")+samples[iMethod][iSample] :
          iMethod==1 ? (std::string(dphi_sideband[iSample] ? "RFine/DPhi2p8_2HadTop_" : "RFine/NTopBandsNew_DPhiBelow2p8_")+samples[iMethod][iSample]) :
          iMethod==2 ? (dphi_sideband[iSample] ? "RBins/DPhiBandsNew_2HadTop" : "RBins/NTopBandsNew_DPhiBelow2p8") : 
          iMethod==3 ? (std::string("RBins/NTopBandsNew_DPhiBelow2p8_")+samples[iMethod][iSample]) : "";
      }
      TCanvas *can = (TCanvas*)(f->Get(canname.c_str())); 
      can = (TCanvas*)can->Clone();
      can->Draw();
      TH1D *h_side = (TH1D*)can->GetListOfPrimitives()->At(i_h_side[iMethod]);
      TH1D *h_signal = (TH1D*)can->GetListOfPrimitives()->At(i_h_signal[iMethod]);
      TH1D *h_pred =(TH1D*)h_side->Clone();
      if (iMethod==1&&rebin[iSample]>1) {
        h_side->Rebin(rebin[iSample]);
        h_signal->Rebin(rebin[iSample]);
        h_pred->Rebin(rebin[iSample]);
      }
      TLegend *leg = (TLegend*)can->GetListOfPrimitives()->At(can->GetListOfPrimitives()->GetEntries()-1);
      leg->SetX1(0.35); leg->SetX2(0.65); leg->SetY1(0.6);
      
      // Add ratio plot
      int y1 = 350;
      int y2 = 150;
      int mid2 = 10;
      can->Divide(1,2);
      // Pad 1 (80+500+20 x 40+500)
      TVirtualPad* p = can->cd(1);
      p->SetPad(0,(y2+60+mid2)/(y1+y2+100.0+mid2),1,1);
      p->SetTopMargin(40.0/(y1+40));
      p->SetBottomMargin(0);
      p->SetRightMargin(0.05);
      p->SetLogy(1);
      h_side->GetYaxis()->SetRangeUser(1.00001e-4,1e4);
      h_side->Draw("HIST");
      h_signal->Draw("SAMEHISTE1");
      leg->Draw();
      // Pad 2 (80+500+20 x 200+60)
      p = can->cd(2);
      p->SetGrid(0,1);
      p->SetPad(0,0,1,(y2+60+mid2)/(y1+y2+100.0+mid2));
      p->SetTopMargin(((float)mid2)/(y2+60+mid2));
      p->SetBottomMargin(60.0/(y2+60+mid2));
      p->SetRightMargin(0.05);
      TH1D* ratio = (TH1D*)h_signal->Clone();
      TH1D* div = (TH1D*)h_side->Clone();
      ratio->Scale(1/ratio->GetSumOfWeights());
      ratio->SetTitleSize(32.0/(y2+60+mid2),"xyz");
      ratio->SetLabelSize(20.0/(y2+60+mid2),"xyz");
      div->Scale(1/div->GetSumOfWeights());
      ratio->Divide(div);
      //ratio->GetYaxis()->SetRangeUser(0,2);
      ratio->GetXaxis()->SetTitleOffset(0.7);
      ratio->GetYaxis()->SetNdivisions(305);
      ratio->GetYaxis()->SetTitle("Ratio (Norm.)");
      ratio->GetYaxis()->SetTitleOffset(0.4);
      ratio->SetTitleSize(24.0/(y2+60+mid2),"y");
      ratio->SetTitle("");
      ratio->SetMarkerStyle(20);
      ratio->SetMarkerColor(1);
      ratio->SetLineColor(1);
      ratio->GetYaxis()->SetRangeUser(0,4);
      ratio->Draw("PE1");
      TLine* l = new TLine(ratio->GetXaxis()->GetXmin(), 1, ratio->GetXaxis()->GetXmax(), 1);
      l->SetLineWidth(2);
      //l->SetLineColor(2);
      l->SetLineStyle(2);
      l->Draw();
      gPad->Update();
      // Fit ratio
      //TF1 *fit_ratio;
      //if (iMethod==1) {
      //  fit_ratio = new TF1("fit_ratio","pol0", Rranges_ABCD[iMethod][0], Rranges_ABCD[iMethod][1]);
      //  fit_ratio->SetLineColor(1);
      //  ratio->Fit("fit_ratio","RE");
      //  fit_ratio->SetRange(Rranges_ABCD[iMethod][0], Rranges_ABCD[iMethod][2]);
      //  fit_ratio->Draw("SAME");
      //}
      p = can->cd(1);
      
      
      // calculate integrals
      double integral[2][2] = { { 0, 0 }, { 0, 0 } };
      double integral_error[2][2] = { { 0, 0 }, { 0, 0 } };
      double nevt[2][2] = { { 0, 0 }, { 0, 0 } };
      for (int i=0; i<2; ++i) {
        for (int bin=1; bin<=h_side->GetNbinsX(); ++bin) {
          if (h_signal->GetXaxis()->GetBinLowEdge(bin)>=Rranges_ABCD[iMethod][i*2] && 
              h_signal->GetXaxis()->GetBinUpEdge(bin)<=Rranges_ABCD[iMethod][i*2+1]) {
            double c0 = h_side->GetBinContent(bin), c1 = h_signal->GetBinContent(bin);
            double e0 = h_side->GetBinError(bin),   e1 = h_signal->GetBinError(bin);
            //std::cout<<h_signal->GetBinError(bin)<<" "<<sqrt(c1*weight[iSample])<<std::endl;
            if (baderror&&iMethod==1) {
              e0 = sqrt(c0*weight[iSample]);
              e1 = sqrt(c1*weight[iSample]);
            }
            nevt[0][i] += (int)(c0*c0/(e0*e0) + 0.5);
            nevt[1][i] += (int)(c1*c1/(e1*e1) + 0.5);
            integral[0][i] += c0;
            integral[1][i] += c1;
            //if (iMethod==1) { // weight bin by projected ratio (correction)
            //  double bincent = h_signal->GetXaxis()->GetBinLowEdge(bin);
            //  integral[0][1] *= fit_ratio->Eval(bincent);
            //}
            integral_error[0][i] += e0*e0;
            integral_error[1][i] += e1*e1;
            //if (iSample==0&&e1>0) std::cout<<bin<<" "<<c1<<" +- "<<e1*e1<<" nevt = "<<((int)(c1*c1/(e1*e1) + 0.5))<<std::endl;
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
      double d_ratio = d_abcd / d;
      double d_ratio_err = sqrt((d_err/d)*(d_err/d) + (d_abcd_err/d_abcd)*(d_abcd_err/d_abcd))*d_ratio;
      double d_pull = (d_abcd-d)/sqrt(d_abcd_err*d_abcd_err + d_err*d_err);
      // Scaled plot
      h_pred->Scale(c/a);
      h_pred->SetLineColor(1);
      h_pred->SetLineStyle(2);
      h_pred->Draw("SAMEHISTE1");
      leg->AddEntry(h_pred, "Prediction (ABCD)", "l");
      
      double fit_integral[2][2], fit_integral_error[2][2], d_fit_comb = 0, d_fit_comb_err = 0;
      if (iMethod==1) {
	// Fit in the full range of NTop Sideband
	// Do fitting and calculate integrals
	TF1 *fit_side = new TF1("NTopSide_fit","exp([0]+[1]*x)", iMethod==2 ? 0.2 : sideband_fit_low_range[iSample], Rranges_ABCD[iMethod][3]);
        fit_side->SetLineColor(h_side->GetLineColor());
        h_side->Fit("NTopSide_fit","QRE");
        //fit_side->Draw("SAME");
        double Rranges_ACfit[3] = { sideband_fit_low_range[iSample], Rranges_ABCD[iMethod][2], Rranges_ABCD[iMethod][3] };
        for (int i=0; i<2; ++i) {
          fit_integral[0][i] = fit_side->Integral(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
          fit_integral_error[0][i] = fit_side->IntegralError(Rranges_ACfit[i], Rranges_ACfit[i+1])/h_signal->GetXaxis()->GetBinWidth(1);
        }
        double par0 = fit_side->GetParameter(0), par0_error = fit_side->GetParError(0);
        double par1 = fit_side->GetParameter(1), par1_error = fit_side->GetParError(1);
        double par1min, par1max; fit_side->GetParLimits(1, par1min, par1max);
	
        // Fit in the Signal region
        // Fitting in sideband, get B area under curve and scale by C/A
        TF1 *fit_signal = new TF1("NTopSignal_RSide_fit","exp([0]+[1]*x)", Rranges_ACfit[0], Rranges_ABCD[iMethod][3]);
        fit_signal->SetLineColor(h_signal->GetLineColor());
        //fit_signal->SetParameter(1, par1); 
        //fit_signal->SetParLimits(1, par1min, par1max);
        h_signal->Fit("NTopSignal_RSide_fit","QREB");
        //fit_signal->Draw("SAME");
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
        //fit_pred->Draw("SAME");
      }
      
      // Check compatibility of prediction to observed distribution
      double ks_test = h_pred->KolmogorovTest(h_signal);
      
      if (iMethod==1) {
	sum_a += a; sum_b += b; sum_c += c; sum_d += d;
	sum_a_err += a_err*a_err; sum_b_err += b_err*b_err; sum_c_err += c_err*c_err; sum_d_err += d_err*d_err;
	sum_d_abcd += d_abcd; sum_d_abcd_err += d_abcd_err*d_abcd_err;
	sum_b_fit += fit_integral[0][1]; sum_b_fit_err += fit_integral_error[0][1]*fit_integral_error[0][1];
	sum_d_fit += fit_integral[1][1]; sum_d_fit_err += fit_integral_error[1][1]*fit_integral_error[1][1];
	sum_d_fit_comb += d_fit_comb; sum_d_fit_comb_err += d_fit_comb_err*d_fit_comb_err;
	sum_d_nevt += d_nevt;
	//printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |\n", d_fit_comb, d_fit_comb_err, d_abcd, d_abcd_err, d, d_err, d_ratio, d_ratio_err);
      }
      if (latex) {
	printf("%s &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &", samples[iMethod][iSample].c_str(), a, a_err, b, b_err, c, c_err);
	printf("  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  %.2f &  %.2f \\\\\n", d_abcd, d_abcd_err, d, d_err, d_ratio, d_ratio_err, d_pull, ks_test);
      } else {
	printf("| %s |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |", samples[iMethod][iSample].c_str(), a, a_err, b, b_err, c, c_err);
	printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f |  %.2f |\n", d_abcd, d_abcd_err, d, d_err, d_ratio, d_ratio_err, d_pull, ks_test);
      }
      // Combining best methods
      if ((iMethod==0&&iSample==0)||(iMethod==1&&iSample!=0)) {
	comb_d_abcd += d_abcd; comb_d += d; comb_d_abcd_err += d_abcd_err*d_abcd_err; comb_d_err += d_err*d_err;
      }
    }
    if (iMethod==1) {
      sum_a_err = sqrt(sum_a_err); sum_b_err = sqrt(sum_b_err); sum_c_err = sqrt(sum_c_err); sum_d_err = sqrt(sum_d_err);
      sum_b_fit_err = sqrt(sum_b_fit_err); sum_d_fit_err = sqrt(sum_d_fit_err); sum_d_fit_comb_err = sqrt(sum_d_fit_comb_err);
      double sum_d_ratio = sum_d_abcd / sum_d;
      double sum_d_ratio_err = sqrt((sum_d_err/sum_d)*(sum_d_err/sum_d) + (sum_d_abcd_err/sum_d_abcd)*(sum_d_abcd_err/sum_d_abcd))*sum_d_ratio;
      double sum_d_pull = (sum_d_abcd-sum_d)/sqrt(sum_d_abcd_err*sum_d_abcd_err + sum_d_err*sum_d_err);
      if (latex) {
	printf("\\hline\n");
	printf("Sum Bkg. &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &", sum_a, sum_a_err, sum_b, sum_b_err, sum_c, sum_c_err);
	printf("  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  %.2f & \\\\\n", sum_d_abcd, sum_d_abcd_err, sum_d, sum_d_err, sum_d_ratio, sum_d_ratio_err, sum_d_pull);
      } else {
	//printf("| Sum Bkg.|  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |", sum_a, sum_a_err, sum_b_fit, sum_b_fit_err, sum_b, sum_b_err, sum_c, sum_c_err);
	printf("| Sum Bkg.|  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |", sum_a, sum_a_err, sum_b, sum_b_err, sum_c, sum_c_err);
	//printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |\n", sum_d_fit_comb, sum_d_fit_comb_err, sum_d_abcd, sum_d_abcd_err, sum_d, sum_d_err, sum_d_ratio, sum_d_ratio_err);
	printf("  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f |  - |\n", sum_d_abcd, sum_d_abcd_err, sum_d, sum_d_err, sum_d_ratio, sum_d_ratio_err, sum_d_pull);
      }
    } else if (iMethod==2) {
      double comb_d_ratio = comb_d_abcd / comb_d;
      double comb_d_ratio_err = sqrt((comb_d_err/comb_d)*(comb_d_err/comb_d) + (comb_d_abcd_err/comb_d_abcd)*(comb_d_abcd_err/comb_d_abcd))*comb_d_ratio;
      double comb_d_pull = (comb_d_abcd-comb_d)/sqrt(comb_d_abcd_err*comb_d_abcd_err + comb_d_err*comb_d_err);
      if (latex) {
	printf("\\hline\n");
	printf("\\hline\n");
	printf("Combined Bkg. & & & &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  $%.2f \\pm %.2f$ &  %.2f & \\\\\n", comb_d_abcd, comb_d_abcd_err, comb_d, comb_d_err, comb_d_ratio, comb_d_ratio_err, comb_d_pull);
	printf("\\hline\n");
      } else {
	printf("| Combined Bkg.| | | |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f +- %.2f |  %.2f | - |\n", comb_d_abcd, comb_d_abcd_err, comb_d, comb_d_err, comb_d_ratio, comb_d_ratio_err, comb_d_pull);
      }
    }
  }
  if (latex) {
    printf("\\hline\n");
    printf("\\end{tabular}\n");
    printf("\\end{center}\n");
    printf("\\end{table*}\n");
  }
}
/*
mv RBins_NTopHad_DPhiBelow2p8_T1tttt_mGl1500_mLSP100.png     /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T1tttt_mGl1500_mLSP100.png	 
mv RBins_NTopHad_DPhiBelow2p8_T5ttttDeg_3bodydec_mGo1000.png /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T5ttttDeg_3bodydec_mGl1000.png
mv RBins_NTopHad_DPhiBelow2p8_T5ttttDeg_4bodydec_mGo1000.png /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T5ttttDeg_4bodydec_mGl1000.png
mv RBins_NTopHad_DPhiBelow2p8_T5ttttDeg_3bodydec_mGo1300.png /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T5ttttDeg_3bodydec_mGl1300.png
mv RBins_NTopHad_DPhiBelow2p8_T5ttttDeg_4bodydec_mGo1300.png /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T5ttttDeg_4bodydec_mGl1300.png
mv RBins_NTopHad_DPhiBelow2p8.png        /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut.png	     
mv RFine_NTopHad_DPhiBelow2p8_T_tW.png   /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_T_tW.png  
mv RFine_NTopHad_DPhiBelow2p8_QCD.png    /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_QCD.png   
mv RFine_NTopHad_DPhiBelow2p8_ZJets.png  /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_ZJets.png 
mv RFine_NTopHad_DPhiBelow2p8_WJets.png  /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_WJets.png 
mv RFine_NTopHad_DPhiBelow2p8_TTJets.png /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/RBins_NTopSB_DPhiCut_TTJets.png
mv DPhiBins_RBins_2HadTop_TTJets.png     /afs/cern.ch/user/j/jkarancs/public/SUSY/Analysis/Plots/2015_05_28/DPhiBins_RSB_2TopCut_TTJets.png
*/
