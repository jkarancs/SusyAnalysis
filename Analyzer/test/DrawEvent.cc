#include "TApplication.h"
#include "TLatex.h"
#include "TArc.h"

#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "../interface/Samples.h"
#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"


int main(int argc, char **argv) {
  TApplication theApp("App", &argc, argv);
  
  std::map<int, std::string> particle;
  particle[1] = "d";
  particle[2] = "u";
  particle[3] = "s";
  particle[4] = "c";
  particle[5] = "b";
  particle[6] = "t";
  
  particle[-1] = "_d";
  particle[-2] = "_u";
  particle[-3] = "_s";
  particle[-4] = "_c";
  particle[-5] = "_b";
  particle[-6] = "_t";
  
  particle[11] = "e-";
  particle[12] = "ve";
  particle[13] = "mu-";
  particle[14] = "vmu";
  particle[15] = "tau-";
  particle[16] = "vtau";
  
  particle[-11] = "e+";
  particle[-12] = "_ve";
  particle[-13] = "mu+";
  particle[-14] = "_vmu";
  particle[-15] = "tau+";
  particle[-16] = "_vtau";
  
  particle[-24] = "W-";
  
  particle[2212] = "p";
  
  particle[1000001] = "~d";
  particle[1000002] = "~u";
  particle[1000003] = "~s";
  particle[1000004] = "~c";
  particle[1000005] = "~b";
  particle[1000006] = "~d";
  particle[1000011] = "~e-";
  particle[1000012] = "~ve";
  particle[1000013] = "~mu-";
  particle[1000014] = "~vmu";
  particle[1000015] = "~tau-";
  particle[1000016] = "~vtau";
  
  Samples samples;
  std::string sample_dir = "/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Apr02_edm_Apr01/";
  samples.AddSample("T5ttttDeg_3bodydec_mGo1300", "T5tttt (#tilde{g}#rightarrowt#tilde{t}_{2,3body}, M_{#tilde{g}}=1.3TeV)", "12",
		    { { .dir=sample_dir+"T5ttttDeg_mGo1300_23bodydec/*.root", .xsec_pb=0.0460525 } });
  samples.AddSample("TTJets", "t#bar{t}+jets", "633", { { .dir=sample_dir+"TTJets/*.root", .xsec_pb=831.76 } }); //
  // Initialize TreeReader
  B2GTreeReader reader;
  
  // Class to Loop on files and read the Trees
  B2GTreeLooper looper(1,1);
  
  // Data variable
  Data d;
  
  // Histogram storage class
  SmartHistos sh;
  sh.AddHistoType("evt");
  
  gStyle->SetOptStat(0);

  std::vector<std::string> dirs = samples.GetListOfDirectories();
  for ( std::string dir : dirs ) looper.AddFile(dir);
  looper.LoopOnSamples();
  looper.LoopOnSamples();
  looper.LoopOnFiles();
  TFile *curr_file = looper.CurrentFile();
  reader.Load_Tree(*curr_file,looper.TreeName());
  while (looper.LoopOnEntries()) {
    reader.GetEntry(looper.CurrentEntry());
    d = reader.data;
    TH2D *h = new TH2D("h","",100,-5,5, 64,-3.2,3.2);
    std::vector<double> top_pt;
    for (int i=0; i<d.gen.size; ++i) {
      if ((abs(d.gen.ID[i])==5||abs(d.gen.ID[i])==6||abs(d.gen.ID[i])==11||abs(d.gen.ID[i])==13||abs(d.gen.ID[i])==15||abs(d.gen.ID[i])>1e6||abs(d.gen.ID[i])==24)&&d.gen.Pt[i]>10) {
        //std::cout<<i<<" Particle: "<<particle[d.gen.ID[i]]<<" ID: "<<d.gen.ID[i]<<" MomID: "<<d.gen.MomID[i]<<" Status: "<<d.gen.Status[i]<<" Pt: "<<d.gen.Pt[i]<<" Eta: "<<d.gen.Eta[i]<<" Phi: "<<d.gen.Phi[i]<<std::endl;
        int fillid = abs(d.gen.ID[i]) > 1e6 ? abs(d.gen.ID[i])-999970 : abs(d.gen.ID[i]);
        if (d.gen.ID[i]!=d.gen.MomID[i]) {
          if (fillid==6) {
	    top_pt.push_back(d.gen.Pt[i]);
	    //h->Fill(d.gen.Eta[i],d.gen.Phi[i],4);
	  }
          //if (fillid==5) h->Fill(d.gen.Eta[i],d.gen.Phi[i],1);
          //if (fillid==24) h->Fill(d.gen.Eta[i],d.gen.Phi[i],7);
        }
      }
    }
    if (top_pt.size()>1) if (top_pt[0]>800&&top_pt[1]>800) {
      h->Draw("COLZ");
      for (int i=0; i<d.jetsAK8.size; ++i) {
	TLatex *jet = new TLatex(d.jetsAK8.Eta[i], d.jetsAK8.Phi[i], "#times"); jet->SetTextAlign(22); jet->Draw();
	//h->Fill(d.jetsAK8.Eta[i],d.jetsAK8.Phi[i],12);
	std::stringstream ss; 
	ss.precision(2); ss<<"#tau_{32} = "<<(d.jetsAK8.tau3[i]/d.jetsAK8.tau2[i]);
	ss.precision(3); ss<<" M = "<<d.jetsAK8.prunedMass[i];
	TLatex *info = new TLatex(d.jetsAK8.Eta[i]+0.8, d.jetsAK8.Phi[i], ss.str().c_str());
	info->SetTextAlign(12); info->SetTextSizePixels(16); info->Draw();
	TArc *a = new TArc(d.jetsAK8.Eta[i], d.jetsAK8.Phi[i], 0.8); 
	a->SetFillStyle(0); a->SetLineColor((d.jetsAK8.tau3[i]/d.jetsAK8.tau2[i])<0.75 && d.jetsAK8.prunedMass[i]>140 ? 3 : 1 ); 
	a->Draw();
      }
      for (int i=0; i<d.gen.size; ++i) if (d.gen.ID[i]!=d.gen.MomID[i]) {
	if (abs(d.gen.ID[i])==6) {
	  TLatex *lat = new TLatex(d.gen.Eta[i], d.gen.Phi[i], "#color[3]{t}"); lat->SetTextAlign(22); lat->Draw();
	}
	if (abs(d.gen.ID[i])==5&&abs(d.gen.MomID[i])==6) {
	  TLatex *lat = new TLatex(d.gen.Eta[i], d.gen.Phi[i], "#color[4]{b}"); lat->SetTextAlign(22); lat->Draw();
	}
	if (abs(d.gen.ID[i])==24&&abs(d.gen.MomID[i])==6) {
	  TLatex *lat = new TLatex(d.gen.Eta[i], d.gen.Phi[i], "#color[2]{W}"); lat->SetTextAlign(22); lat->Draw();
	}
	if ((abs(d.gen.ID[i])==11||abs(d.gen.ID[i])==13||abs(d.gen.ID[i])==15)&&abs(d.gen.MomID[i])==24) {
	  TLatex *lat = new TLatex(d.gen.Eta[i], d.gen.Phi[i], "l"); lat->SetTextAlign(22); lat->Draw();
	}
      }
      gPad->Update();
      sleep(10);
    }
    delete h;
  }
  
  return 1;
}
