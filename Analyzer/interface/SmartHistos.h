#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TF1.h"
#include "TFile.h"
#include "TH3.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TProfile2D.h"

class Postfixes {
public:
  Postfixes() {}
  ~Postfixes() {}
  typedef struct Postfix { std::function<size_t()> sel; std::vector<std::string> vec; std::vector<std::string> leg; std::string colz; } Postfix;
  
private:
  std::map<const char*, Postfix> pf_map_;
  
  double get_dbl_(std::string& str) {
    std::stringstream ss(str); 
    double d; ss>>d; 
    std::stringstream size; size<<d;
    str.erase(0,size.str().size());
    return d; 
  }
  
  std::vector<double> str_to_vec_dbl_(std::string s) {
    std::vector<double> v_dbl;
    std::istringstream ss(s);
    std::string sub;
    while(std::getline(ss, sub, ',')) {
      if (sub.find("to")!=std::string::npos) {
	double i,j,k=1;
	i = get_dbl_(sub);
	sub.erase(0,2);
	j = get_dbl_(sub);
	if (sub.size()) {
	  sub.erase(0,1);
	  k = get_dbl_(sub);
	} else k=(i<j)?1:-1;
	if (k>=0) for (double ii=i; ii<=j; ii+=k) v_dbl.push_back(ii);
	else for (double ii=i; ii>=j; ii+=k) v_dbl.push_back(ii);
      } else v_dbl.push_back(get_dbl_(sub));
    }
    return v_dbl;
  }
  
  std::vector<std::string> interpret_substring_(std::string s, bool format=0) {
    std::vector<std::string> v_str;
    std::istringstream ss(s);
    std::string sub;
    std::string elem;
    while(std::getline(ss, sub, '[')) {
      size_t f = sub.find("]");
      if (f==std::string::npos) {
	elem += sub;
	for (size_t i=0; i<v_str.size(); ++i) v_str[i]+=sub;
      } else {
	std::vector<double> v_dbl = str_to_vec_dbl_(sub.substr(0,f));
	if (!v_str.size()) for (size_t i=0; i<v_dbl.size(); ++i) {
	  std::stringstream ss2; ss2<<elem<<v_dbl[i]<<sub.substr(f+1,sub.size()-f-1);
	  std::string str = ss2.str();
	  if (format) if (str.find(".")!=std::string::npos) str.replace(str.find("."),1,"p");
	  v_str.push_back(str);
	} else for (size_t i=0; i<std::min(v_dbl.size(),v_str.size()); ++i) {
	  std::stringstream ss2; ss2<<v_str[i]<<v_dbl[i];
	  std::string str = ss2.str();
	  if (format) if (str.find(".")!=std::string::npos) str.replace(str.find("."),1,"p");
	  v_str[i]=str;
	}
      }
    }
    if (!v_str.size()) v_str.push_back(elem);
    return v_str;
  }
  
  std::vector<std::string> interpret_string_(std::string s, bool format=0) {
    std::vector<std::string> v_str;
    std::istringstream ss(s);
    std::string sub;
    while(std::getline(ss, sub, ';')) {
      std::vector<std::string> subs = interpret_substring_(sub,format);
      for (size_t i=0; i<subs.size(); ++i) v_str.push_back(subs[i]);
    }
    return v_str;
  }
public:
  void AddNew(const char* name, std::function<size_t()> sel, std::string pf, std::string leg, std::string colz) { 
    pf_map_.insert(std::pair<const char*, Postfix>(name, { .sel=sel, .vec=interpret_string_(pf,1), .leg=interpret_string_(leg), .colz=colz })); 
  }
  
  Postfix GetPostfix(const char* name) { 
    size_t count = pf_map_.count(name);
    if (!count) std::cout<<"!!! ERROR: Postfixes::GetPostfix: Postfix with name = "<<name<<" was not found."<<std::endl;
    return (count) ? pf_map_[name] : Postfix({0,std::vector<std::string>(),std::vector<std::string>(),""});
  }
};


class Cuts {
public:
  Cuts() {}
  ~Cuts() {}
  
private:
  std::map<const char*, std::function<bool()>> cut_map_;
  
public:
  void AddNew(const char* name, std::function<bool()> cut) { cut_map_.insert(std::pair<const char*, std::function<bool()>>(name, cut )); }
  
  std::function<bool()> GetCut(const char* name) {
    size_t count = cut_map_.count(name);
    if (!count) std::cout<<"!!! ERROR: Cuts::GetCut: Cut with name = "<<name<<" was not found."<<std::endl;
    return (count) ? cut_map_[name] : [](){ return 0; };
  }
};


class SmartHisto {
  
public:
  // constructors, destructor
  SmartHisto(std::string name, std::vector<const char*>& pf_names, std::vector<Postfixes::Postfix> pfs,
	     std::vector<std::function<double()> >& ffs, std::vector<std::function<double()> >& weights, std::vector<std::function<bool()> >& cuts, 
	     std::string draw, std::string opt, std::vector<double> axis_ranges) { 
    name_=name;
    pf_names_=pf_names;
    pfs_ = pfs;
    npf_=pfs_.size();
    ndim_=ffs.size();
    if (ndim_>0) fill_1d_ = ffs[0];
    if (ndim_>1) fill_2d_ = ffs[1];
    if (ndim_>2) fill_3d_ = ffs[2];
    nweight_=weights.size();
    if (nweight_>0) weight1_ = weights[0];
    if (nweight_>1) weight2_ = weights[1];
    if (nweight_>2) weight3_ = weights[2];
    ncut_=cuts.size();
    if (ncut_>0) cut1_ = cuts[0];
    if (ncut_>1) cut2_ = cuts[1];
    if (ncut_>2) cut3_ = cuts[2];
    if (ncut_>3) cut4_ = cuts[3];
    if (ncut_>4) cut5_ = cuts[4];
    draw_=draw;
    norm_  = opt.find("Norm")!=std::string::npos;
    keep_  = opt.find("Keep")!=std::string::npos;
    stat_  = opt.find("Stat")!=std::string::npos ? 1110 : 0;
    sumw2_ = opt.find("Sumw2")!=std::string::npos;
    axis_ranges_=axis_ranges;
    if (npf_>4) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 4 postfixes, only use max 4, or redefine functions!\n";
    if (ndim_>3) std::cout<<"!!! ERROR: SmartHisto::constructor: More than 4 dimension, define a maximum of 3!\n";
    if (ncut_>5) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 5 cuts specified, please add new variables to store them!\n";
    if (nweight_>3) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 3 weights specified, please add new variables to store them!\n";
    init_();
  }
  ~SmartHisto() {}
  
private:
  std::string name_;
  std::vector<const char*> pf_names_;
  std::vector<Postfixes::Postfix> pfs_;
  size_t npf_;
  
  // Fill functions specify how to fill histos
  size_t ndim_;
  std::function<double()> fill_1d_;
  std::function<double()> fill_2d_;
  std::function<double()> fill_3d_;
  
  // Weights: pointers to floats
  size_t nweight_;
  std::function<double()> weight1_;
  std::function<double()> weight2_;
  std::function<double()> weight3_;
  
  // Cuts: pointers to bools
  size_t ncut_;
  std::function<bool()> cut1_;
  std::function<bool()> cut2_;
  std::function<bool()> cut3_;
  std::function<bool()> cut4_;
  std::function<bool()> cut5_;
  
  // Draw arguments and plot options
  // options:
  // NORM - normalize
  // STAT - draw stat boxes)
  std::string draw_;
  bool norm_; // Normalize histo
  bool keep_; // Keep colors/markers in order
  Int_t stat_; // Draw Stat boxes
  bool sumw2_;
  
  // axis ranges: xlow, xhigh, ylow, yhigh, zlow, zhigh
  // if low==high -> do not set
  std::vector<double> axis_ranges_;
  
  // histo containers
  TH1D* h1d_0p_;
  TH2D* h2d_0p_;
  TH3D* h3d_0p_;
  std::vector<TH1D*> h1d_1p_;
  std::vector<TH2D*> h2d_1p_;
  std::vector<TH3D*> h3d_1p_;
  std::vector<std::vector<TH1D*> > h1d_2p_;
  std::vector<std::vector<TH2D*> > h2d_2p_;
  std::vector<std::vector<TH3D*> > h3d_2p_;
  std::vector<std::vector<std::vector<TH1D*> > > h1d_3p_;
  std::vector<std::vector<std::vector<TH2D*> > > h2d_3p_;
  std::vector<std::vector<std::vector<TH3D*> > > h3d_3p_;
  std::vector<std::vector<std::vector<std::vector<TH1D*> > > > h1d_4p_;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > h2d_4p_;
  std::vector<std::vector<std::vector<std::vector<TH3D*> > > > h3d_4p_;
  
  // --------------------------------------------------------------------------
  //                      Special Histogram Calculations:
  
  // Define plots that are calculated from +1 dimensional objects
  // Eg.: MPV from cluster charge distribution
  std::vector<std::vector<std::string> > spec_;
  std::vector<std::vector<std::string> > spec2_;
  TF1* ring_fit_[3][4];
  void init_() {
    //                HistoParam name, +1D name, Axis Title, +1D Axis Title
    spec_.push_back({"HitEfficiency","ValidHit","Hit Efficiency", "Valid Hit"});
    spec_.push_back({"DColEfficiency","Validhit","DCol Efficiency", "Valid Hit"});
    //                Name Prefix, Title Prefix
    spec2_.push_back({"Avg","Avg. "});
    
    // Get Hit Eff vs DCol Eff functions (previously measured)
    std::string fname[3] = {"hiteff_vs_dcol_l1", "hiteff_vs_dcol_l2", "hiteff_vs_dcol_l3"};
    std::string ringname[4] = {"_ring1", "_ring2", "_ring3", "_ring4"};
    TFile *f_func = TFile::Open("/data/jkarancs/CMSSW/TimingStudy/CMSSW_7_1_0_pre1/src/DPGAnalysis/PixelTimingStudy/test/DynIneff_scale_factors/HitEffvsDColFunctions.root");
    if (f_func) {
      for (int lay=0; lay<3; ++lay) for (int ring=0; ring<4; ring++) {
	TF1* f = (TF1*)f_func->Get((fname[lay]+ringname[ring]).c_str());
	ring_fit_[lay][ring] = (TF1*)f->Clone();
      }
      f_func->Close();
    } else std::cout<<"SmartHisto::init_() - File not found: "<<f_func->GetName()<<std::endl;
    
    // Initialize some variables
    h1d_0p_ = 0;
    h2d_0p_ = 0;
    h3d_0p_ = 0;
  }
  
  size_t find_spec_(std::string name) {
    size_t find = -1;
    for (size_t s=0; s<spec_.size(); ++s) if (name.find(spec_[s][0])!=std::string::npos) find = s;
    return find;
  }
  size_t find_spec2_(std::string name) {
    size_t find = -1;
    for (size_t s=0; s<spec2_.size(); ++s) if (name.find(spec2_[s][0])==0) find = s;
    return find;
  }
  
  // Define which functions to use for each special histo
  void calc_spec_1d_(TH1D* h1d, TH2D* h2d) {
    if (find_spec_(name_)==0) calc_eff_1d_(h1d, h2d); // HitEfficiency
    if (find_spec_(name_)==1) calc_dcol_1d_(h1d, h2d); // DColEfficiency
    if (find_spec2_(name_)==0) calc_eff_1d_(h1d, h2d, 0); // Average (Use Profile)
  }
  
  void calc_spec_2d_(TH2D* h2d, TH3D* h3d) {
    if (find_spec_(name_)==0) calc_eff_2d_(h2d, h3d); // HitEfficiency
    //if (find_spec_(name_)==1) calc_dcol_2d_(h2d, h3d); // DColEfficiency
    if (find_spec2_(name_)==0) calc_eff_2d_(h2d, h3d); // Average (Use 3DProfile)
  }
  
  // Define functions below: (eg: Efficiency, MPV etc...)
  
  // ********************  HitEfficiency *********************
  void calc_eff_1d_(TH1D* h1d, TH2D* h2d, int err_type = 2) {
    if (err_type==0) { // Use TProfile
      TProfile* p = h2d->ProfileX();
      for (int i=1; i<=p->GetNbinsX(); ++i) {
	h1d->SetBinContent(i,p->GetBinContent(i));
	h1d->SetBinError(i,p->GetBinError(i));
      }
      delete p;
      //if (h2d->GetNbinsY()==2) {
      //  TProfile* p = h2d->ProfileX();
      //  for (int i=1; i<=p->GetNbinsX(); ++i) {
      //    h1d->SetBinContent(i,p->GetBinContent(i));
      //    h1d->SetBinError(i,p->GetBinError(i));
      //  }
      //  delete p;
      //} else {
      //  for (int i=1; i<=h2d->GetNbinsX(); ++i) {
      //    double mean=0, n=0, err=0;
      //    for (int j=1; j<=h2d->GetNbinsY(); ++j) {
      //      double cont = h2d->GetBinContent(i,j);
      //      double bincent = h2d->GetYaxis()->GetBinCenter(j);
      //      mean += bincent * cont;
      //      n += cont;
      //    }
      //    mean = (n>0) ? mean/n : 0;
      //    for (int j=1; j<=h2d->GetNbinsY(); ++j) {
      //      double cont = h2d->GetBinContent(i,j);
      //      double bincent = h2d->GetYaxis()->GetBinCenter(j);
      //      err += (bincent-mean)*(bincent-mean)*cont;
      //    }
      //    err = (n>0) ? sqrt(err/n) : 0;
      //    if (n>0) {
      //      std::cout<<h1d->GetName()<<" bin"<<i<<" "<<mean<<" +- "<<err<<std::endl;
      //      h1d->SetBinContent(i,mean);
      //      h1d->SetBinError(i,err);
      //    }
      //  }
      //}
    } else if (err_type==1) { // Normal Approximation interval (Same result as TProfile)
      double z = 1; // N Sigma confidence
      for (int i=1; i<=h2d->GetNbinsX(); ++i) {
	double val = h2d->GetBinContent(i,2), mis = h2d->GetBinContent(i,1);
	if (val+mis>0) {
	  double eff = val / (val + mis);
	  double err = z*sqrt(eff*(1.0-eff)/(val+mis));
	  h1d->SetBinContent(i,eff);
	  h1d->SetBinError(i,err);
	}
      }
    } else if (err_type==2) { // Wilson Score Interval
      double z = 1; // N Sigma confidence
      for (int i=1; i<=h2d->GetNbinsX(); ++i) {
	double val = h2d->GetBinContent(i,2), mis = h2d->GetBinContent(i,1);
	if (val+mis>0) {
	  double eff = val / (val + mis);
	  double n = val + mis;
	  double cen = (eff+(z*z/(2*n))) / (1.0 + (z*z/n));
	  double halfwidth = z*sqrt( eff*(1.0-eff)/n + (z*z/(4*n*n)) ) / (1.0 + (z*z/n));
	  double err = halfwidth + fabs(cen-eff);
	  // Assymmetric error -> Choose larger for a conservative error estimate
	  h1d->SetBinContent(i,eff);
	  h1d->SetBinError(i,err);
	}
      }
    }
    h1d->SetEntries(h2d->GetEntries());
  }
  
  void calc_eff_2d_(TH2D* h2d, TH3D* h3d) {
    TProfile2D* p = h3d->Project3DProfile("yx");
    for (int i=1; i<=p->GetNbinsX(); ++i) for (int j=1; j<=p->GetNbinsY(); ++j) {
      h2d->SetBinContent(i,j,p->GetBinContent(i,j));
      h2d->SetBinError(i,j,p->GetBinError(i,j));
    }
    delete p;
    h2d->SetEntries(h3d->GetEntries());
  }
  
  // ******************** DColEfficiency *********************
  void calc_dcol_1d_(TH1D* h1d, TH2D* h2d) {
    std::string name =h1d->GetName();
    // Check which Layer/Ring is plotted (or if it is a module plot)
    bool mod_plot = (std::string(h1d->GetName()).find("Modules")!=std::string::npos);
    int Lay = 0; // Assume layer 1 if no such postfix is given
    for (int lay=0; lay<3; ++lay) {
      std::stringstream ss; ss<<"Lay"<<lay+1;
      if (name.find(ss.str())!=std::string::npos) Lay = lay;
    }
    int Ring = -1;
    for (int ring=0; ring<4; ring++) {
      std::stringstream ss; ss<<"Mod"<<ring+1;
      if (name.find(ss.str())!=std::string::npos) Ring = ring;
    }
    // Get HitEfficiency and use hiteff_vs_dcol functions to get DCol Efficiency
    calc_eff_1d_(h1d, h2d);
    for (int bin=1; bin<=h1d->GetNbinsX(); ++bin) {
      if (mod_plot) Ring = abs(bin-5)-1;
      if (Lay==0&&Ring==3) Ring = 2;
      double hiteff = h1d->GetBinContent(bin);
      double dcoleff = (Ring!=-1) ? ( hiteff>ring_fit_[Lay][Ring]->Eval(0.75) ? ring_fit_[Lay][Ring]->GetX(hiteff) : 0 ): 0;
      h1d->SetBinContent(bin, dcoleff);
      h1d->SetBinError(bin, 0);
    }
    h1d->SetEntries(h2d->GetEntries());
  }
  
  //TF1* get_eff_vs_dcol_func_(TH1* h) {
  //  // Check which Layer, Ring is plotted
  //  std::string name =h->GetName();
  //  int Lay = 0; // Assume layer 1 if no such postfix is given
  //  for (int lay=0; lay<3; ++lay) {
  //    std::stringstream ss; ss<<"Lay"<<lay+1;
  //    if (name.find(ss.str())!=std::string::npos) Lay = lay;
  //  }
  //  int Ring = -1;
  //  for (int ring=0; ring<4; ring++) {
  //    std::stringstream ss; ss<<"Ring"<<ring+1;
  //    if (name.find(ss.str())!=std::string::npos) Ring = ring;
  //  }
  //  if (Lay==0&&Ring==3) Ring = 2;
  //  if (Ring!=-1) {
  //    // Get Hit Eff vs DCol Eff function (previously measured)
  //    std::string name[3] = {"hiteff_vs_dcol_l1", "hiteff_vs_dcol_l2", "hiteff_vs_dcol_l3"};
  //    std::string ringname[4] = {"_mod1", "_mod2", "_mod3", "_mod4"};
  //    TFile *f_func = TFile::Open("/data/jkarancs/CMSSW/TimingStudy/CMSSW_7_1_0_pre1/src/DPGAnalysis/PixelTimingStudy/test/DynIneff_scale_factors/HitEffvsDColFunctions.root");
  //    TF1* ring_fit = (TF1*)f_func->Get((name[Lay]+ringname[Ring]).c_str());
  //    f_func->Close();
  //    std::cout<<name<<" "<<name[Lay]+ringname[Ring]<<std::endl;
  //    return ring_fit;
  //  } else return NULL;
  //}
  //
  //void calc_dcol_2d_(TH2D* h2d, TH3D* h3d) {
  //  TF1* ring_fit = get_eff_vs_dcol_func_(h2d);
  //  if (ring_fit) {
  //    // Get HitEfficiency and use hiteff(dcoleff) to get DCol Efficiency
  //    calc_eff_2d_(h2d, h3d);
  //    //for (int binx=1; binx<=h2d->GetNbinsX(); ++binx) for (int biny=1; biny<=h2d->GetNbinsY(); ++biny) {
  //    //  double hiteff = h2d->GetBinContent(binx, biny);
  //    //  double dcoleff = ring_fit->GetX(hiteff);
  //    //  h2d->SetBinContent(binx, biny, dcoleff);
  //    //}
  //  }
  //}
  
  void calc_specials_() {
    if (find_spec_(name_)!=(size_t)-1||find_spec2_(name_)!=(size_t)-1) {
      switch (npf_) {
      case 0:
	switch (ndim_) {
	case 2: calc_spec_1d_(h1d_0p_, h2d_0p_); break;
	case 3: calc_spec_2d_(h2d_0p_, h3d_0p_); break;
	} break;
      case 1:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_1p_.size(); ++i)
	    calc_spec_1d_(h1d_1p_[i], h2d_1p_[i]); break;
	case 3: 
	  for (size_t i=0; i<h3d_1p_.size(); ++i)
	    calc_spec_2d_(h2d_1p_[i], h3d_1p_[i]); break;
	} break;
      case 2:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j)
	    calc_spec_1d_(h1d_2p_[i][j], h2d_2p_[i][j]); break;
	case 3:
	  for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j)
	    calc_spec_2d_(h2d_2p_[i][j], h3d_2p_[i][j]); break;
	} break;
      case 3:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
	    for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) calc_spec_1d_(h1d_3p_[i][j][k], h2d_3p_[i][j][k]); break;
	case 3:
	  for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j)
	    for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) calc_spec_2d_(h2d_3p_[i][j][k], h3d_3p_[i][j][k]); break;
	} break;
      case 4:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) for (size_t k=0; k<h2d_4p_[i][j].size(); ++k)
	    for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) calc_spec_1d_(h1d_4p_[i][j][k][l], h2d_4p_[i][j][k][l]); break;
	case 3:
	  for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) for (size_t k=0; k<h3d_4p_[i][j].size(); ++k)
	    for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) calc_spec_2d_(h2d_4p_[i][j][k][l], h3d_4p_[i][j][k][l]); break;
	} break;
      }
    }
  }
  
  std::string spec_dirname_(TObject* obj) {
    bool is_spec = false;
    std::string dirname = name_;
    std::string hname = obj->GetName();
    if (hname.find(name_)==std::string::npos) for (size_t s=0; s<spec_.size(); ++s) if (hname.find(spec_[s][1])!=std::string::npos) {
      is_spec = true;
      dirname.replace(dirname.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
    }
    if (!is_spec) dirname = "";
    else dirname += "/";
    return dirname;
  }
  
  std::string name_only_pf_(TObject* obj) {
    std::string name = obj->GetName();
    if (npf_>0) {
      std::string del = spec_dirname_(obj);
      if (del.size()) del.replace(del.size()-1,1,"_");
      else del = name_ + "_";
      if (name.find(del)!=std::string::npos)
	name.erase(name.find(del),del.size());
    }
    return name;
  }
  
  void write_(TObject* obj) {
    if (obj) {
      std::string main_dir = name_;
      std::string spec_dir = spec_dirname_(obj);
      bool is_spec = (spec_dir.size()!=0);
      if (!gDirectory->GetKey(main_dir.c_str())) gDirectory->mkdir(main_dir.c_str());
      gDirectory->cd(main_dir.c_str());
      if (is_spec) {
        if (!gDirectory->GetKey(spec_dir.c_str())) gDirectory->mkdir(spec_dir.c_str()); 
        gDirectory->cd(spec_dir.c_str());
      }
      obj->Write(name_only_pf_(obj).c_str()); 
      gDirectory->cd(is_spec?"../..":"..");
    }
  }
  
  void load_(TFile* f, TH1D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH1D* h_temp = (TH1D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) { h->Add(h_temp); h->SetEntries(h->GetEntries()+h_temp->GetEntries()); }
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  void load_(TFile* f, TH2D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH2D* h_temp = (TH2D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) { h->Add(h_temp); h->SetEntries(h->GetEntries()+h_temp->GetEntries()); }
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  void load_(TFile* f, TH3D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH3D* h_temp = (TH3D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) { h->Add(h_temp); h->SetEntries(h->GetEntries()+h_temp->GetEntries()); }
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  // Load/Add all histos from a file
  void load_all_(TFile* f, bool add = 0) {
    if (npf_==0) {
      load_(f,h1d_0p_,add);
      load_(f,h2d_0p_,add);
      load_(f,h3d_0p_,add);
    } else if (npf_==1) {
      for (size_t i=0; i<h1d_1p_.size(); ++i) load_(f,h1d_1p_[i],add);
      for (size_t i=0; i<h2d_1p_.size(); ++i) load_(f,h2d_1p_[i],add);
      for (size_t i=0; i<h3d_1p_.size(); ++i) load_(f,h3d_1p_[i],add);
    } else if (npf_==2) {
      for (size_t i=0; i<h1d_2p_.size(); ++i) for (size_t j=0; j<h1d_2p_[i].size(); ++j) load_(f,h1d_2p_[i][j],add);
      for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j) load_(f,h2d_2p_[i][j],add);
      for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j) load_(f,h3d_2p_[i][j],add);
    } else if (npf_==3) {
      for (size_t i=0; i<h1d_3p_.size(); ++i) for (size_t j=0; j<h1d_3p_[i].size(); ++j)
        for (size_t k=0; k<h1d_3p_[i][j].size(); ++k) load_(f,h1d_3p_[i][j][k],add);
      for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) load_(f,h2d_3p_[i][j][k],add);
      for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) load_(f,h3d_3p_[i][j][k],add);
    } else if (npf_==4) {
      for (size_t i=0; i<h1d_4p_.size(); ++i) for (size_t j=0; j<h1d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_4p_[i][j].size(); ++k) for (size_t l=0; l<h1d_4p_[i][j][k].size(); ++l) load_(f,h1d_4p_[i][j][k][l],add);
      for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_4p_[i][j].size(); ++k) for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) load_(f,h2d_4p_[i][j][k][l],add);
      for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_4p_[i][j].size(); ++k) for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) load_(f,h3d_4p_[i][j][k][l],add);
    }
  }
  
  bool pass_cuts_() {
    switch (ncut_) {
    case 0: return 1; break;
    case 1: return (cut1_()); break;
    case 2: return (cut1_())&&(cut2_()); break;
    case 3: return (cut1_())&&(cut2_())&&(cut3_()); break;
    case 4: return (cut1_())&&(cut2_())&&(cut3_())&&(cut4_()); break;
    case 5: return (cut1_())&&(cut2_())&&(cut3_())&&(cut4_())&&(cut5_()); break;
    default: return 1; // Warning in constructor
    }
  }
  
public:
  const std::string& GetName() { return name_; }
  const std::vector<const char*>& GetPFNames() { return pf_names_; }
  
  // Add New SmartHisto to the container vector
  // and set Filling properties, Postfixes, titles etc...
  // 2D/3D: If one of the Variable is special, also create
  // -1D object for the result
  // Eg: 2D: ClusterCharge vs Variable
  // --> 1D: MPV vs ClusteCharge
  //
  // 1D
  void AddNew(std::string name, std::string title,
	      int nbin1, double low1, double high1) {
    if (npf_==0) {
      h1d_0p_ = new TH1D(name.c_str(), title.c_str(), nbin1, low1, high1);
      h1d_0p_->Sumw2();
    } else if (npf_==1) {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h1d_1p_.push_back(new TH1D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(), nbin1, low1, high1));
	if (sumw2_) h1d_1p_[i]->Sumw2();
      }
    } else if (npf_==2) {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h1d_2p_.push_back(std::vector<TH1D*>());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h1d_2p_[i].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
					title.c_str(), nbin1, low1, high1));
	  if (sumw2_) h1d_2p_[i][j]->Sumw2();
	}
      }
    } else if (npf_==3) {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h1d_3p_.push_back(std::vector<std::vector<TH1D*> > ());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h1d_3p_[i].push_back(std::vector<TH1D*>());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	    h1d_3p_[i][j].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
					     title.c_str(), nbin1, low1, high1));
	    if (sumw2_) h1d_3p_[i][j][k]->Sumw2();
	  }
	}
      }
    } else if (npf_==4) {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h1d_4p_.push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h1d_4p_[i].push_back(std::vector<std::vector<TH1D*> >());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	    h1d_4p_[i][j].push_back(std::vector<TH1D*>());
	    for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
	      h1d_4p_[i][j][k].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
						  title.c_str(), nbin1, low1, high1));
	      if (sumw2_) h1d_4p_[i][j][k][l]->Sumw2();
	    }
	  }
	}
      }
    }
  }
  // 2D
  void AddNew(std::string name, std::string title,
	      int nbin1, double low1, double high1,
	      int nbin2, double low2, double high2) { 
    size_t s =find_spec_(name); 
    size_t s2 =find_spec2_(name); 
    if (npf_==0) {
      if (s!=(size_t)-1) {
	h1d_0p_ = new TH1D(name.c_str(), title.c_str(), nbin1, low1, high1);
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      } else if (s2!=(size_t)-1) {
	h1d_0p_ = new TH1D(name.c_str(), std::string(title).insert(1,spec2_[s2][1]).c_str(), nbin1, low1, high1);
	name.erase(0,spec2_[s2][0].size());
      }
      h2d_0p_ = new TH2D(name.c_str(), title.c_str(), nbin1, low1, high1, nbin2, low2, high2);
    } else if (npf_==1) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	  h1d_1p_.push_back(new TH1D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(), nbin1, low1, high1));
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      } else if (s2!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	  h1d_1p_.push_back(new TH1D((name+"_"+pfs_[0].vec[i]).c_str(), std::string(title).insert(1,spec2_[s2][1]).c_str(), nbin1, low1, high1));
	name.erase(0,spec2_[s2][0].size());
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	h2d_1p_.push_back(new TH2D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(), 
				   nbin1, low1, high1, nbin2, low2, high2));
    } else if (npf_==2) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_2p_.push_back(std::vector<TH1D*>());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	    h1d_2p_[i].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), title.c_str(), nbin1, low1, high1));
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      } else if (s2!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_2p_.push_back(std::vector<TH1D*>());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	    h1d_2p_[i].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), std::string(title).insert(1,spec2_[s2][1]).c_str(), nbin1, low1, high1));
	}
	name.erase(0,spec2_[s2][0].size());
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h2d_2p_.push_back(std::vector<TH2D*>());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	  h2d_2p_[i].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), title.c_str(),
					nbin1, low1, high1, nbin2, low2, high2));
      }
    } else if (npf_==3) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_3p_.push_back(std::vector<std::vector<TH1D*> > ());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h1d_3p_[i].push_back(std::vector<TH1D*>());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	      h1d_3p_[i][j].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(), title.c_str(), nbin1, low1, high1));
	  }
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      } else if (s2!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_3p_.push_back(std::vector<std::vector<TH1D*> > ());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h1d_3p_[i].push_back(std::vector<TH1D*>());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	      h1d_3p_[i][j].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
					       std::string(title).insert(1,spec2_[s2][1]).c_str(), nbin1, low1, high1));
	  }
	}
	name.erase(0,spec2_[s2][0].size());
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h2d_3p_.push_back(std::vector<std::vector<TH2D*> > ());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h2d_3p_[i].push_back(std::vector<TH2D*>());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	    h2d_3p_[i][j].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(), title.c_str(),
					     nbin1, low1, high1, nbin2, low2, high2));
	}
      }
    } else if (npf_==4) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_4p_.push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h1d_4p_[i].push_back(std::vector<std::vector<TH1D*> >());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	      h1d_4p_[i][j].push_back(std::vector<TH1D*>());
	      for (size_t l=0; l<pfs_[3].vec.size(); ++l)
		h1d_4p_[i][j][k].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(), title.c_str(), nbin1, low1, high1));
	    }
	  }
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      } else if (s2!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h1d_4p_.push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h1d_4p_[i].push_back(std::vector<std::vector<TH1D*> >());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	      h1d_4p_[i][j].push_back(std::vector<TH1D*>());
	      for (size_t l=0; l<pfs_[3].vec.size(); ++l)
		h1d_4p_[i][j][k].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
						    std::string(title).insert(1,spec2_[s2][1]).c_str(), nbin1, low1, high1));
	    }
	  }
	}
	name.erase(0,spec2_[s2][0].size());
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h2d_4p_.push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h2d_4p_[i].push_back(std::vector<std::vector<TH2D*> >());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	    h2d_4p_[i][j].push_back(std::vector<TH2D*>());
	    for (size_t l=0; l<pfs_[3].vec.size(); ++l)
	      h2d_4p_[i][j][k].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(), title.c_str(),
						  nbin1, low1, high1, nbin2, low2, high2));
	  }
	}
      }
    }
  }
  // 3D
  void AddNew(std::string name, std::string title,
	      int nbin1, double low1, double high1,
	      int nbin2, double low2, double high2,
	      int nbin3, double low3, double high3) { 
    size_t s =find_spec_(name); 
    if (npf_==0) {
      if (s!=(size_t)-1) {
	h2d_0p_ = new TH2D(name.c_str(), title.c_str(), nbin1, low1, high1, nbin2, low2, high2);
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      }
      h3d_0p_ = new TH3D(name.c_str(), title.c_str(), nbin1, low1, high1, nbin2, low2, high2, nbin3, low3, high3);
    } else if (npf_==1) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) 
	  h2d_1p_.push_back(new TH2D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(),
				     nbin1, low1, high1, nbin2, low2, high2));
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) 
	h3d_1p_.push_back(new TH3D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(),
				   nbin1, low1, high1, nbin2, low2, high2, nbin3, low3, high3));
    } else if (npf_==2) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h2d_2p_.push_back(std::vector<TH2D*>());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	    h2d_2p_[i].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), title.c_str(),
					  nbin1, low1, high1, nbin2, low2, high2));
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h3d_2p_.push_back(std::vector<TH3D*>());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	  h3d_2p_[i].push_back(new TH3D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), title.c_str(),
					nbin1, low1, high1, nbin2, low2, high2, nbin3, low3, high3));
      }
    } else if (npf_==3) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h2d_3p_.push_back(std::vector<std::vector<TH2D*> > ());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h2d_3p_[i].push_back(std::vector<TH2D*>());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	      h2d_3p_[i][j].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(), title.c_str(),
					       nbin1, low1, high1, nbin2, low2, high2));
	  }
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h3d_3p_.push_back(std::vector<std::vector<TH3D*> > ());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h3d_3p_[i].push_back(std::vector<TH3D*>());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	    h3d_3p_[i][j].push_back(new TH3D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(), title.c_str(),
					     nbin1, low1, high1, nbin2, low2, high2, nbin3, low3, high3));
	}
      }
    } else if (npf_==4) {
      if (s!=(size_t)-1) {
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	  h2d_4p_.push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    h2d_4p_[i].push_back(std::vector<std::vector<TH2D*> >());
	    for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	      h2d_4p_[i][j].push_back(std::vector<TH2D*>());
	      for (size_t l=0; l<pfs_[3].vec.size(); ++l)
		h2d_4p_[i][j][k].push_back(new TH2D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(), title.c_str(),
						    nbin1, low1, high1, nbin2, low2, high2));
	    }
	  }
	}
	name.replace(name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	title.replace(title.find(spec_[s][2]),spec_[s][2].size(),spec_[s][3]);
      }
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	h3d_4p_.push_back(std::vector<std::vector<std::vector<TH3D*> > >());
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	  h3d_4p_[i].push_back(std::vector<std::vector<TH3D*> >());
	  for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	    h3d_4p_[i][j].push_back(std::vector<TH3D*>());
	    for (size_t l=0; l<pfs_[3].vec.size(); ++l)
	      h3d_4p_[i][j][k].push_back(new TH3D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(), title.c_str(),
						  nbin1, low1, high1, nbin2, low2, high2, nbin3, low3, high3));
	  }
	}
      }
    }
  }
  
  // Fill Histograms using the std::function<double()>
  void Fill() {
    //std::cout<<name_<<std::endl;
    //if (name_=="OnCluChargeNorm"&&npf_==2&&pfs_[1].sel()==0&&pass_cuts_()==1) std::cout<<pfs_[0].sel()<<" "<<pfs_[1].sel()<<" "<<fill_1d_()<<" "<<*weight1_<<std::endl;
    //if (name_=="HitEfficiency_vs_LayersDisks"&&npf_==1&&pass_cuts_()==1) std::cout<<pfs_[0].sel()<<" "<<fill_1d_()<<" "<<*weight1_<<std::endl;
    if (pass_cuts_()) {
      double weight = 0;
      switch (nweight_) {
      case 0: weight = 1; break;
      case 1: weight = (weight1_()); break;
      case 2: weight = (weight1_())*(weight2_()); break;
      case 3: weight = (weight1_())*(weight2_())*(weight3_()); break;
      default: break;
      }
      switch (npf_) {
      case 0:
        switch (ndim_) {
        case 1: h1d_0p_->Fill(fill_1d_(),weight); break;
        case 2: h2d_0p_->Fill(fill_1d_(),fill_2d_(),weight); break;
        case 3: h3d_0p_->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
        } break;
      case 1:
        if (pfs_[0].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 2:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 3:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1&&pfs_[2].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 4:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1&&pfs_[2].sel()!=(size_t)-1&&pfs_[3].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      }
    }
  }
  
  void Load(TFile* f) { load_all_(f); calc_specials_(); }
  
  void Add(TFile* f) { load_all_(f, 1); calc_specials_(); }
  
  void CalcSpecials() { calc_specials_(); }

  // Calculate special histos and write them in a file
  void Write() {
    //calc_specials_();
    if (npf_==0) {
      write_(h1d_0p_);
      write_(h2d_0p_);
      write_(h3d_0p_);
    } else if (npf_==1) {
      for (size_t i=0; i<h1d_1p_.size(); ++i) write_(h1d_1p_[i]);
      for (size_t i=0; i<h2d_1p_.size(); ++i) write_(h2d_1p_[i]);
      for (size_t i=0; i<h3d_1p_.size(); ++i) write_(h3d_1p_[i]);
    } else if (npf_==2) {
      for (size_t i=0; i<h1d_2p_.size(); ++i) for (size_t j=0; j<h1d_2p_[i].size(); ++j) write_(h1d_2p_[i][j]);
      for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j) write_(h2d_2p_[i][j]);
      for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j) write_(h3d_2p_[i][j]);
    } else if (npf_==3) {
      for (size_t i=0; i<h1d_3p_.size(); ++i) for (size_t j=0; j<h1d_3p_[i].size(); ++j)
        for (size_t k=0; k<h1d_3p_[i][j].size(); ++k) write_(h1d_3p_[i][j][k]);
      for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) write_(h2d_3p_[i][j][k]);
      for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) write_(h3d_3p_[i][j][k]);
    } else if (npf_==4) {
      for (size_t i=0; i<h1d_4p_.size(); ++i) for (size_t j=0; j<h1d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_4p_[i][j].size(); ++k) for (size_t l=0; l<h1d_4p_[i][j][k].size(); ++l) write_(h1d_4p_[i][j][k][l]);
      for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_4p_[i][j].size(); ++k) for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) write_(h2d_4p_[i][j][k][l]);
      for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_4p_[i][j].size(); ++k) for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) write_(h3d_4p_[i][j][k][l]);
    }
  }
  
  //______________________________________________________________________
  //                       Multidraw functions
  
  TCanvas* custom_can_(TH1D* h, std::string canname, int gx = 0, int gy = 0,
        	       int histosize_x = 500, int histosize_y = 500,
        	       int mar_left = 80, int mar_right = 20, int mar_top = 20, int mar_bottom = 60) {
    mar_top += (std::string(h->GetTitle()).size()>0)*20;
    int titlefontsize = 32;
    int labelfontsize = 20;
    int yoffset_x = mar_left - titlefontsize - 4;
    int xoffset_y = mar_bottom - titlefontsize - 4;
    int zoffset_x = mar_right - titlefontsize - 4;
    int padsize_x = histosize_x + mar_left + mar_right;
    int padsize_y = histosize_y + mar_top + mar_bottom;
    int padsize = ((padsize_x<=padsize_y) ? padsize_x : padsize_y);
    float padratio_yx = (float)padsize_y/padsize_x > 1 ? 1 : (float)padsize_y/padsize_x;
    float padratio_xy = (float)padsize_x/padsize_y > 1 ? 1 : (float)padsize_x/padsize_y;
    Float_t xoffset = ((Float_t)xoffset_y/titlefontsize+0.5) * padratio_xy /1.6;
    Float_t yoffset = ((Float_t)yoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    Float_t zoffset = ((Float_t)zoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    Float_t titlesize = (Float_t)titlefontsize/padsize;
    Float_t labelsize = (Float_t)labelfontsize/padsize;
    h->SetTitleFont(42,"xyz");
    h->SetLabelFont(42,"xyz");
    h->SetTitleSize(titlesize,"xyz");
    h->SetLabelSize(labelsize,"xyz");
    h->GetXaxis()->SetTitleOffset(xoffset);
    h->GetYaxis()->SetTitleOffset(yoffset);
    h->GetZaxis()->SetTitleOffset(zoffset);
    h->GetYaxis()->SetDecimals(1);
    h->GetZaxis()->SetDecimals(1);
    //gStyle->SetPadLeftMargin((Float_t)mar_left/padsize_x);
    //gStyle->SetPadRightMargin((Float_t)mar_right/padsize_x);
    //gStyle->SetPadTopMargin((Float_t)mar_top/padsize_y);
    //gStyle->SetPadBottomMargin((Float_t)mar_bottom/padsize_y);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleH(titlefontsize/padsize);
    gStyle->SetTitleFontSize(titlesize);
    TCanvas* canvas = new TCanvas(canname.c_str(), h->GetTitle(), padsize_x + 4, padsize_y + 26);
    TVirtualPad* pad = canvas->cd(1);
    pad->SetLeftMargin((Float_t)mar_left/padsize_x);
    pad->SetRightMargin((Float_t)mar_right/padsize_x);
    pad->SetTopMargin((Float_t)mar_top/padsize_y);
    pad->SetBottomMargin((Float_t)mar_bottom/padsize_y);
    canvas->SetGrid(gx,gy);
    return canvas;
  }
  
  TCanvas* custom_can_(TH2D* h, std::string canname, int gx = 0, int gy = 0, 
		       int histosize_x = 500, int histosize_y = 500, int mar_left = 80, int mar_right = 120) {
    int mar_top = 20 + (std::string(h->GetTitle()).size()>0)*20;
    int mar_bottom = 60;
    int titlefontsize = 32;
    int labelfontsize = 20;
    int pal_offset_x = 5;
    int pal_width_x = 25;
    int xoffset_y = mar_bottom - titlefontsize - 4;
    int yoffset_x = mar_left - titlefontsize - 4;
    int zoffset_x = mar_right - pal_offset_x - pal_width_x - titlefontsize;
    int padsize_x = histosize_x + mar_left + mar_right;
    int padsize_y = histosize_y + mar_top + mar_bottom;
    int padsize = ((padsize_x<=padsize_y) ? padsize_x : padsize_y);
    float padratio_yx = (Float_t)padsize_y/padsize_x > 1 ? 1 : (Float_t)padsize_y/padsize_x;
    float padratio_xy = (Float_t)padsize_x/padsize_y > 1 ? 1 : (Float_t)padsize_x/padsize_y;
    Float_t xoffset = ((Float_t)xoffset_y/titlefontsize+0.5) * padratio_xy /1.6;
    Float_t yoffset = ((Float_t)yoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    Float_t zoffset = ((Float_t)zoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    Float_t titlesize = (Float_t)titlefontsize/padsize;
    Float_t labelsize = (Float_t)labelfontsize/padsize;
    h->SetTitleFont(42,"xyz");
    h->SetLabelFont(42,"xyz");
    h->SetTitleSize(titlesize,"xyz");
    h->SetLabelSize(labelsize,"xyz");
    h->GetXaxis()->SetTitleOffset(xoffset);
    h->GetYaxis()->SetTitleOffset(yoffset);
    h->GetZaxis()->SetTitleOffset(zoffset);
    h->GetZaxis()->RotateTitle(1);
    h->GetYaxis()->SetDecimals(1);
    h->GetZaxis()->SetDecimals(1);
    if (histosize_y<250) h->GetZaxis()->SetNdivisions(505);
    //gStyle->SetPadLeftMargin((Float_t)mar_left/padsize_x);
    //gStyle->SetPadRightMargin((Float_t)mar_right/padsize_x);
    //gStyle->SetPadTopMargin((Float_t)mar_top/padsize_y);
    //gStyle->SetPadBottomMargin((Float_t)mar_bottom/padsize_y);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleH(titlefontsize/padsize);
    gStyle->SetTitleFontSize(titlesize);
    TCanvas* canvas = new TCanvas(canname.c_str(), h->GetTitle(), padsize_x + 4, padsize_y + 26);
    TVirtualPad* pad = canvas->cd(1);
    pad->SetLeftMargin((Float_t)mar_left/padsize_x);
    pad->SetRightMargin((Float_t)mar_right/padsize_x);
    pad->SetTopMargin((Float_t)mar_top/padsize_y);
    pad->SetBottomMargin((Float_t)mar_bottom/padsize_y);
    canvas->SetGrid(gx,gy);
    if (norm_&&h->Integral()>0) h = (TH2D*)h->DrawNormalized(draw_.c_str());
    else h->Draw(draw_.c_str());
    if (h->Integral()>0&&draw_=="COLZ") {
      gPad->Update();
      TPaletteAxis* palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
      if (palette) {
	palette->SetX1NDC(1 - (Float_t)(mar_right - pal_offset_x)/padsize_x);
	palette->SetX2NDC(1 - (Float_t)(mar_right - pal_offset_x - pal_width_x)/padsize_x);
	palette->SetY1NDC((Float_t)mar_bottom/padsize_y);
	palette->SetY2NDC(1 - (Float_t)mar_top/padsize_y);
      }
    }
    gStyle->SetOptTitle(0);
    return canvas;
  }
  
  std::vector<int> string_to_vector_(std::string val) {
    std::vector<int> vec;
    // Reading values from string separated by commas, (ranges x-y can be specified too)
    // eg. "3-6,20,22,24-27,11," - add comma at the end of the string
    std::stringstream ss(val);
    while (ss) {
      while (ss && !isdigit(ss.peek ())) ss.get ();
      if (! ss)	break;
      int i = 0;
      int j = 0;
      ss >> i;
      while (ss && !ispunct(ss.peek ())) ss.get ();
      char a; ss>>a;
      if (a=='-') {
        while (ss && !isdigit(ss.peek ())) ss.get ();
        ss >> j;
      } else j=i;
      for (int n=i; n<=j; n++) vec.push_back(n);
    }
    return vec;
  }
  
  void set_stat_(TH1D* h, Color_t col, int i) {
    if (i<7) {
      gPad->Update();
      TPaveStats* stats = (TPaveStats*)h->FindObject("stats");
      stats->SetLineColor(col);
      stats->SetTextColor(col);
      stats->SetX1NDC(0.74);
      stats->SetX2NDC(0.94);
      stats->SetY1NDC(0.8-i*0.11);
      stats->SetY2NDC(0.9-i*0.11);
    }
  }
  
  void multidraw_with_legend_(size_t skip, std::vector<TH1D*>& hvec, std::vector<std::string> pf, std::string colz,
			      std::string legtitle="", float x1=0.15, float y2=0.9) {
    // Draw multiple histograms, set their marker/line color/style
    // Then Draw legend for all histo with titles from a postfix
    std::vector<int> col = string_to_vector_(colz);
    std::vector<Int_t> marker = { 20, 21, 22, 23, 29, //full circle, square, tri-up, tri-down, star
				  20, 21, 22, 23, 29,
				  20, 21, 22, 23, 29 };
    int nleg = 0;
    //size_t i_highest = -1;
    //double highest = 0;
    for (size_t i=skip; i<hvec.size(); i++) if (hvec[i]->GetEntries()>0) {
      // Determine the highest histo
      // we draw this histo first, so default range is ok for all histo
      // then redraw all histo in predetedmined order
      //double max = hvec[i]->GetMaximum() / (norm_ ? hvec[i]->GetSumOfWeights() : 1 );
      //if (max > highest) {
      //  highest = max;
      //  i_highest = i;
      //}
      if ((draw_.find("P")!=std::string::npos)) { 
	hvec[i]->SetMarkerColor((Color_t)col[i-(keep_?skip:0)]); 
	hvec[i]->SetMarkerStyle(marker[i-(keep_?skip:0)]); 
      } else { 
	hvec[i]->SetLineColor((Color_t)col[i-(keep_?skip:0)]);
	hvec[i]->SetLineWidth(2);
      }
      if (stat_) hvec[i]->SetStats(1);
      ++nleg;
    }
    if (legtitle.size()>0) ++nleg;
    float x2 = x1 + 0.2;
    float y1 = y2 - nleg * 0.05;
    TLegend *leg = new TLegend(x1,y1,x2,y2,legtitle.c_str());
    std::string same = (stat_ ? "SAMES" : "SAME") + draw_;
    for (size_t i=skip; i<hvec.size(); i++) if (hvec[i]->GetEntries()>0) {
      // Draw (highest plot first, then redraw all incorrect order)
      //if (i==skip) {
      //  if (norm_) hvec[i_highest]->DrawNormalized(draw_.c_str());
      //  else hvec[i_highest]->Draw(draw_.c_str());	  
      //}
      if (norm_&&hvec[i]->Integral()>0) hvec[i]->DrawNormalized((i==skip) ? draw_.c_str() : same.c_str());
      else hvec[i]->Draw((i==skip) ? draw_.c_str() : same.c_str());
      if (stat_) set_stat_(hvec[i], (Color_t)col[i-(keep_?skip:0)], i-skip);
      std::stringstream colored_text;
      colored_text<<"#color["<<(Color_t)col[i-(keep_?skip:0)]<<"]{"<<pf[i]<<"}";
      leg->AddEntry(hvec[i], colored_text.str().c_str(), draw_.find("P")!=std::string::npos ? "P" : "L");
    }
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw("SAME");
  }
  
  void draw_one_(TH1D* h) {
    if ((draw_.find("P")!=std::string::npos)) { 
      h->SetMarkerColor(1); 
      h->SetMarkerStyle(20); 
    } else { 
      h->SetLineColor(1);
      h->SetLineWidth(2);
    }
    if (stat_) h->SetStats(1);
    if (norm_&&h->Integral()>0) h->DrawNormalized(draw_.c_str());
    else h->Draw(draw_.c_str());
    if (stat_) set_stat_(h,1,0);
  }
  
  typedef struct DrawParams1D { std::vector<TH1D*> hvec; std::string canname; std::string legtitle;  } DrawParams1D;
  std::vector<DrawParams1D> get_dps_1stpf_1d_() {
    std::vector<DrawParams1D> dps;
    bool is_spec = (find_spec_(name_)!=(size_t)-1)||(find_spec2_(name_)!=(size_t)-1);
    if (ndim_==1||(ndim_==2&&is_spec)) {
      if (npf_==0) {
	dps.push_back({ .hvec={h1d_0p_}, .canname=name_, .legtitle="" });
      } else if (npf_==1) {
	dps.push_back({ .hvec={}, .canname=name_+"_"+pf_names_[0], .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_1p_[i]);
      } else if (npf_==2) for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_2p_[i][j]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j];
	dps[dps.size()-1].legtitle=pfs_[1].leg[j];
      } else if (npf_==3) for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_3p_[i][j][k]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k];
	dps[dps.size()-1].legtitle=pfs_[1].leg[j]+" "+pfs_[2].leg[k];
      } else if (npf_==4) for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_4p_[i][j][k][l]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l];
	dps[dps.size()-1].legtitle=pfs_[1].leg[j]+" "+pfs_[2].leg[k]+" "+pfs_[3].leg[l];
      }
    }
    return dps;
  }
  
  typedef struct DrawParams2D { TH2D* h; std::string canname; } DrawParams2D;
  std::vector<DrawParams2D> get_dps_2d_() {
    std::vector<DrawParams2D> dps;
    bool is_spec = (find_spec_(name_)!=(size_t)-1)||(find_spec2_(name_)!=(size_t)-1);
    if ((ndim_==2&&!is_spec)||(ndim_==3&&is_spec)) {
      if (npf_==0)
	dps.push_back({ .h=h2d_0p_, .canname=name_});
      else if (npf_==1) for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	dps.push_back({ .h=h2d_1p_[i], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i] });
      else if (npf_==2) for (size_t i=0; i<pfs_[0].vec.size(); ++i) for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	dps.push_back({ .h=h2d_2p_[i][j], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j] });
      else if (npf_==3) for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k)
	  dps.push_back({ .h=h2d_3p_[i][j][k], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k] });
      else if (npf_==4) for (size_t i=0; i<pfs_[0].vec.size(); ++i) for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	for (size_t k=0; k<pfs_[2].vec.size(); ++k) for (size_t l=0; l<pfs_[3].vec.size(); ++l)
	  dps.push_back({ .h=h2d_4p_[i][j][k][l], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l] });
    }
    return dps;
  }
  
  void DrawPlots() {
    gStyle->SetOptStat(stat_);
    // 1D plots
    std::vector<DrawParams1D> dps1d = get_dps_1stpf_1d_();
    for (size_t i=0; i<dps1d.size(); ++i) {
      size_t skip = 0;
      while (dps1d[i].hvec[skip]->GetEntries()==0) { 
	++skip; 
	if (skip==dps1d[i].hvec.size()) break; 
      }
      if (skip<dps1d[i].hvec.size()) {
	if (axis_ranges_.size()>=2) if (axis_ranges_[0]!=axis_ranges_[1]) 
	  dps1d[i].hvec[skip]->GetXaxis()->SetRangeUser(axis_ranges_[0],axis_ranges_[1]);
	if (axis_ranges_.size()>=4) if (axis_ranges_[2]!=axis_ranges_[3]) 
	  dps1d[i].hvec[skip]->GetYaxis()->SetRangeUser(axis_ranges_[2],axis_ranges_[3]);
	TCanvas *c = custom_can_(dps1d[i].hvec[skip], dps1d[i].canname);
	if (npf_) multidraw_with_legend_(skip, dps1d[i].hvec, pfs_[0].leg, pfs_[0].colz, dps1d[i].legtitle);
	else draw_one_(dps1d[i].hvec[0]);
	write_(c);
      }
    }
    // 2D Plots
    std::vector<DrawParams2D> dps2d = get_dps_2d_();
    for (size_t i=0; i<dps2d.size(); ++i) {
      if (dps2d[i].h->GetEntries()) {
	if (axis_ranges_.size()>=2) if (axis_ranges_[0]!=axis_ranges_[1]) 
	  dps2d[i].h->GetXaxis()->SetRangeUser(axis_ranges_[0],axis_ranges_[1]);
	if (axis_ranges_.size()>=4) if (axis_ranges_[2]!=axis_ranges_[3]) 
	  dps2d[i].h->GetYaxis()->SetRangeUser(axis_ranges_[2],axis_ranges_[3]);
	if (axis_ranges_.size()==6) if (axis_ranges_[4]!=axis_ranges_[5]) 
	  dps2d[i].h->GetZaxis()->SetRangeUser(axis_ranges_[4],axis_ranges_[5]);
	TCanvas *c = custom_can_(dps2d[i].h, dps2d[i].canname);
	write_(c);
      }
    }
  }
  
  // histo containers
  TH1D*& GetH1D() { return h1d_0p_; }
  TH2D*& GetH2D() { return h2d_0p_; }
  TH3D*& GetH3D() { return h3d_0p_; }
  std::vector<TH1D*>& GetH1D1P() { return h1d_1p_; }
  std::vector<TH2D*>& GetH2D1P() { return h2d_1p_; }
  std::vector<TH3D*>& GetH3D1P() { return h3d_1p_; }
  std::vector<std::vector<TH1D*> >& GetH1D2P() { return h1d_2p_; }
  std::vector<std::vector<TH2D*> >& GetH2D2P() { return h2d_2p_; }
  std::vector<std::vector<TH3D*> >& GetH3D2P() { return h3d_2p_; }
  std::vector<std::vector<std::vector<TH1D*> > >& GetH1D3P() { return h1d_3p_; }
  std::vector<std::vector<std::vector<TH2D*> > >& GetH2D3P() { return h2d_3p_; }
  std::vector<std::vector<std::vector<TH3D*> > >& GetH3D3P() { return h3d_3p_; }
  std::vector<std::vector<std::vector<std::vector<TH1D*> > > >& GetH1D4P() { return h1d_4p_; }
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > >& GetH2D4P() { return h2d_4p_; }
  std::vector<std::vector<std::vector<std::vector<TH3D*> > > >& GetH3D4P() { return h3d_4p_; }
  
};

class SmartHistos {
  
public:
  // constructor, destructor
  SmartHistos() { pf_ = new Postfixes(); cut_ = new Cuts(); weights_={}; }
  ~SmartHistos() {}
  typedef struct FillParams { int nbin; double low; double high; std::function<double()> fill; std::string axis_title; } FillParams;
  
private:
  // Postfix container
  Postfixes* pf_;
  
  // Cut container
  Cuts* cut_;
  
  // SmartHisto containers
  std::map<std::string, std::vector<SmartHisto*> > sh_;
  
  // FillParams container
  std::map<std::string, FillParams> hp_map_;
  
  std::vector<std::function<double()> > weights_;
  
  FillParams get_hp_(std::string name) {
    size_t count = hp_map_.count(name);
    if (count) {
      return hp_map_[name];
    } else {
      // Check if name has a special prefix (Average etc)
      size_t find = std::string::npos;
      std::string match = "";
      for ( auto el : hp_map_) {
	size_t f = name.find(el.first);
	if (f<find) {
	  find = f;
	  match = el.first;
	}
      }
      if (find!=std::string::npos) {
	return hp_map_[match];
      } else {
	std::cout<<"!!! ERROR: SmartHistos::get_hp_: FillParams with name = "<<name<<" was not found."<<std::endl;
	return FillParams({.nbin=0,.low=0,.high=0,.fill={ [](){return -9999.9;} },.axis_title=""});
      }
    }
  }
  
  std::vector<FillParams> get_hp_vec_(std::string fill) {
    std::vector<FillParams> vec;
    size_t found1 = 0;
    size_t found2 = fill.find("_vs_",found1);
    while (found2 != std::string::npos) {
      fill.erase(found2,4);
      vec.push_back(get_hp_(fill.substr(found1,found2-found1)));
      found1=found2;
      found2 = fill.find("_vs_",found1);
    }
    vec.push_back(get_hp_(fill.substr(found1,fill.length()-found1)));
    return vec;
  }
  
  void set_default_style_() {
    gStyle->SetPaperSize(20.,20.);
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetErrorX(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineWidth(2);
    //gStyle->SetLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(1);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleX(1);
    gStyle->SetTitleY(1);
    gStyle->SetTitleAlign(33);
  }

public:
  void AddNewFillParam(std::string name, FillParams hp) { hp_map_.insert(std::pair<std::string, FillParams>(name, hp )); }
  
  void AddNewPostfix(const char* name, std::function<size_t()> sel, std::string pf, std::string leg, std::string colz) { pf_->AddNew(name, sel, pf, leg, colz); }
  
  void AddNewCut(const char* name, std::function<bool()> cut) { cut_->AddNew(name, cut); }
  
  void SetHistoWeights(std::vector<std::function<double()> > weights) { weights_ = weights; }
  
  void AddHistoType(std::string type) { sh_[type] = std::vector<SmartHisto*>(); }

  typedef struct HistoParams { std::string fill; std::vector<const char*> pfs; std::vector<const char*> cuts; std::string draw; std::string opt; std::vector<double> ranges; } HistoParams;
  void AddHistos(std::string name, HistoParams hp, bool AddCutsToTitle = true) {
    if (sh_.count(name)) {
      std::vector<FillParams> hp_vec = get_hp_vec_(hp.fill);
      bool valid = 1;
      for (size_t i=0; i<hp_vec.size(); ++i) valid = (valid&&(hp_vec[i].nbin!=0));
      if (valid) {
        std::vector<Postfixes::Postfix> pfs;
        for (size_t i=0; i<hp.pfs.size(); ++i) pfs.push_back(pf_->GetPostfix(hp.pfs[i]));
        std::string axis_titles = "";
        std::vector<std::function<double()> > fillfuncs;
	if (AddCutsToTitle) for (size_t icut=0; icut<hp.cuts.size(); ++icut) {
	  if (icut>0) axis_titles += ", ";
	  axis_titles += hp.cuts[icut];
	}
        for (size_t i=hp_vec.size(); i>0; --i) {
          axis_titles += ";" + hp_vec[i-1].axis_title;
          fillfuncs.push_back(hp_vec[i-1].fill);
        }
	std::vector<std::function<bool()>> cuts;
        for (size_t i=0; i<hp.cuts.size(); ++i) cuts.push_back(cut_->GetCut(hp.cuts[i]));
        sh_[name].push_back(new SmartHisto(hp.fill.c_str(), hp.pfs, pfs, fillfuncs, weights_, cuts, hp.draw, hp.opt, hp.ranges));
        if (hp_vec.size()==1) sh_[name][sh_[name].size()-1]->AddNew(hp.fill, axis_titles, 
								    hp_vec[0].nbin, hp_vec[0].low, hp_vec[0].high);
        else if (hp_vec.size()==2) sh_[name][sh_[name].size()-1]->AddNew(hp.fill, axis_titles,
									 hp_vec[1].nbin, hp_vec[1].low, hp_vec[1].high,
									 hp_vec[0].nbin, hp_vec[0].low, hp_vec[0].high);
        else if (hp_vec.size()==3) sh_[name][sh_[name].size()-1]->AddNew(hp.fill, axis_titles,
									 hp_vec[2].nbin, hp_vec[2].low, hp_vec[2].high,
									 hp_vec[1].nbin, hp_vec[1].low, hp_vec[1].high,
									 hp_vec[0].nbin, hp_vec[0].low, hp_vec[0].high);
      }
    }
  }
  
  void PrintNames() { 
    for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it) {
      for (size_t i=0; i<it->second.size(); ++i)  {
	std::cout<<it->second[i]->GetName();
	size_t npf = it->second[i]->GetPFNames().size();
	if (npf>0) {
	  std::cout<<" - ";
	  for (size_t j=0; j<npf; ++j) {
	    std::cout<<it->second[i]->GetPFNames()[j];
	    if (j+1!=npf) std::cout<<", ";
	  }
	}
	std::cout<<"\n";
      }
    }
  }
  
  void Load(const char* filename) {
    TFile* f = TFile::Open(filename); 
    for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
      for (size_t i=0; i<it->second.size(); ++i) it->second[i]->Load(f); 
    f->Close();
  }
  
  void Add(const char* filenames) {
    TChain* fc = new TChain("fc");
    fc->Add(filenames);
    TObjArray* fileElements=fc->GetListOfFiles();
    TIter next(fileElements); TChainElement* chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      TFile* f = TFile::Open(chEl->GetTitle()); 
      for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
	for (size_t i=0; i<it->second.size(); ++i) it->second[i]->Add(f); 
      f->Close();
    }
  }

  void CalcSpecials() { for (auto vs : sh_) for (auto s : vs.second) s->CalcSpecials(); }
  
  void Fill(std::string name) { 
    for (size_t i=0; i<sh_[name].size(); ++i) sh_[name][i]->Fill(); 
  }
  
  void DrawPlots() {
    set_default_style_();
    for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
      for (size_t i=0; i<it->second.size(); ++i) it->second[i]->DrawPlots(); 
  }
  
  void Write(std::string name = "") { 
    if (name.size()) for (size_t i=0; i<sh_[name].size(); ++i) sh_[name][i]->Write();
    else for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
      for (size_t i=0; i<it->second.size(); ++i) it->second[i]->Write();
  }
  
  SmartHisto* GetHistos(std::string name, size_t i) { return sh_[name][i]; }
  
};
