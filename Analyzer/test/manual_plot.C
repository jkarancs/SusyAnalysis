void manual_plot() {
  TChain *c = new TChain("c");
  //c->Add("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb28_edm_Feb20/TT/*.root/B2GTTreeMaker/B2GTree");
  c->Add("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb28_edm_Feb20/WJetsToLNu_HT-600toInf/*.root/B2GTTreeMaker/B2GTree");
  //c->Add("/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Feb28_edm_Feb20/T5ttttDeg_mGo1300_4bodydec/*.root/B2GTTreeMaker/B2GTree");
  //c->Draw("evt_TTHadR","evt_NTopHad==2");
  TH1D* h[4];
  h[0] = new TH1D("h1","",120,0,1.2);
  h[1] = new TH1D("h2","",120,0,1.2);
  h[2] = new TH1D("h3","",120,0,1.2);
  h[3] = new TH1D("h4","",120,0,1.2);
  TCanvas *can = new TCanvas("c","",1000,500);
  can->Divide(2);
  can->cd(1)->SetLogy(1);
  c->Draw("evt_R>>h1","jetAK8_Pt[0]>400&&jetAK8_Pt[1]>400&&evt_NTopHad<2");
  c->Draw("evt_R>>h2","jetAK8_Pt[0]>400&&jetAK8_Pt[1]>400&&jetAK8_prunedMass[0]>50&&jetAK8_prunedMass[1]>50&&evt_NTopHad<2","SAME");
  c->Draw("evt_R>>h3","jetAK8_Pt[0]>400&&jetAK8_Pt[1]>400&&jetAK8_prunedMass[0]>100&&jetAK8_prunedMass[1]>100&&evt_NTopHad<2","SAME");
  c->Draw("evt_R>>h4","evt_NTopHad==2","SAME");
  h[0]->GetYaxis()->SetRangeUser(0.0001,1000000);

  TH1D* slope = new TH1D("slope","",4,0,4);
  slope->GetXaxis()->SetBinLabel(1,"#color[1]{NTop<2 p_{T}>400}");
  slope->GetXaxis()->SetBinLabel(2,"#color[2]{NTop<2 p_{T}>400 Mass>50}");
  slope->GetXaxis()->SetBinLabel(3,"#color[3]{NTop<2 p_{T}>400 Mass>100}");
  slope->GetXaxis()->SetBinLabel(4,"#color[4]{NTop==2}");
  for (int i=0; i<4; ++i) {
    h[i]->SetLineColor(i+1);
    TF1* fit = new TF1((std::string(h[i]->GetName())+"_fit").c_str(),"exp([0]+[1]*x)", 0.2, 1.2);
    h[i]->Fit((std::string(h[i]->GetName())+"_fit").c_str(),"RE0");
    fit->SetLineColor(i+1);
    fit->Draw("SAME");
    slope->SetBinContent(i+1,fit->GetParameter(1));
    slope->SetBinError(i+1,fit->GetParError(1));
  }
  can->cd(2);
  slope->Draw("");
}
