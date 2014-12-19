#ifndef B2GTreeLooper_h
#define B2GTreeLooper_h

#include "TFile.h"
#include "TH3.h"
#include "TStopwatch.h"
#include "TChain.h"
#include <string>
#include <iomanip>
#include <iostream>

class B2GTreeLooper {
public:
  B2GTreeLooper(int nthfile = 1, bool show_progress = 0) :
    nthfile_(nthfile),show_progress_(show_progress) {
    init_();
  }
  ~B2GTreeLooper() {
    delete sw_;
    delete chain_;
  }
  
private:
  const int nthfile_;
  const bool show_progress_;
  
  void init_() {
    chain_ = new TChain("chain");
    if (show_progress_) {
      sw_ = new TStopwatch;
      processed_entries_=0;
      progress_=0;
      step_size_ = 0.0001;
    }
  }
  
  // For B2GTreeLooper
  TChain* chain_;
  TObjArray* files_;
  TFile* file_;
  int fileindex_;
  Long64_t entry_;
  Long64_t all_entries_;
  const char* treename_;
  
  // For ProgressEstimator
  TStopwatch* sw_;
  Long64_t total_entries_;
  Long64_t processed_entries_;
  double progress_;
  double step_size_;
  

  void ProgressEstimator_() {
    if (double(processed_entries_)/double(total_entries_)>progress_+step_size_) {
      progress_ += step_size_;
      sw_->Stop();
      double totaltime = sw_->RealTime();
      sw_->Continue();
      double entriespersec = double(processed_entries_)/totaltime;
      double timeleft = double(total_entries_-processed_entries_)/ entriespersec;
      int time_h = int(timeleft/3600);
      int time_m = (timeleft-(time_h*3600))/60;
      int time_s = timeleft-(time_h*3600)-(time_m)*60;
      step_size_ = ((time_h)?30:1)*(entriespersec/double(total_entries_));
      std::cout<<"\r"<<std::flush;
      printf("----------- %03.1f%% ( %02d:%02d:%02d Left ) -----------", progress_*100, time_h, time_m, time_s);
    }
    processed_entries_++;
  }
  
public:
  //Modifiers
  void AddFiles(std::string address) {
    chain_->Add(address.c_str());
    
    files_=chain_->GetListOfFiles();
    fileindex_=0;
    
    // Calculate total number of entries - For ProgressEstimator
    total_entries_ = (nthfile_==1) ? chain_->GetEntries() : 0;
    if (show_progress_&&nthfile_!=1) {
      for (int nf=0; nf<files_->GetEntries(); ++nf) {
        if (nf%nthfile_==nthfile_/2) {
          TFile f(files_->At(nf)->GetTitle());
          total_entries_ += ((TTree*)f.Get(files_->At(nf)->GetName()))->GetEntries();
        }
      }
      step_size_ = 500000/total_entries_;
    }
  }
  
  bool LoopOnFiles() {
    if (fileindex_==0&&show_progress_) {
      std::cout<<"--------- Started Looping on Trees ---------"<<std::flush;
      sw_->Start(0);
    }
    while (fileindex_%nthfile_!=(nthfile_-1)/2) ++fileindex_; // Skip files if specified
    if (fileindex_<files_->GetEntries()) {
      file_ = TFile::Open(files_->At(fileindex_)->GetTitle());
      treename_ = files_->At(fileindex_)->GetName();
      all_entries_ = ((TTree*)file_->Get(treename_))->GetEntries();
      if (!show_progress_) std::cout<<"Opened file: "<<files_->At(fileindex_)->GetTitle()<<std::endl;
      ++fileindex_;
      entry_=-1;
      return 1;
    } else {
      return 0;
    }
  }
  
  bool LoopOnEntries() { 
    ++entry_; 
    if (show_progress_) ProgressEstimator_();
    return entry_<all_entries_;
  }
  
  // Accessors
  const char* TreeName() { return treename_; }
  
  TFile* CurrentFile() { return file_; }
  
  Long64_t CurrentEntry() { return entry_; }


};

#endif