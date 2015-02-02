#ifndef B2GTreeLooper_h
#define B2GTreeLooper_h

#include "TFile.h"
#include "TH3.h"
#include "TStopwatch.h"
#include "TChain.h"
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>

class B2GTreeLooper {
public:
  B2GTreeLooper(int nthfile = 1, bool show_progress = 0) :
    nthfile_(nthfile),show_progress_(show_progress) {
    init_();
  }
  ~B2GTreeLooper() {
    delete sw_;
  }
  
private:
  const int nthfile_;
  const bool show_progress_;
  
  void init_() {
    if (show_progress_) {
      sw_ = new TStopwatch;
      processed_entries_=0;
      progress_=0;
      step_size_ = 0.0001;
      total_entries_=0;
    }
    it_sample=-1;
  }
  
  // For B2GTreeLooper
  std::vector<TChain*> samples_;
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
  void AddFile(std::string fileaddress, int new_loop=1) {
    std::string treename = "/DMTreesDumper/TreeBase";
    size_t size = samples_.size();
    if (new_loop!=0) {
      std::stringstream ss;
      ss<<"_"<<size+1;
      samples_.push_back(new TChain((std::string("filechain")+ss.str()).c_str()));
      size = samples_.size();
    }
    if (new_loop!=-1) samples_[size-1]->Add((fileaddress+treename).c_str());

    // Calculate total number of entries - For ProgressEstimator
    total_entries_ += (nthfile_==1) ? samples_[size-1]->GetEntries() : 0;
    if (show_progress_&&nthfile_!=1) {
      int nadded = samples_[size-1]->GetListOfFiles()->GetEntries();
      int nth = std::min(nadded, nthfile_);
      for (int nf=0; nf<nadded; ++nf) {
        if (nf%nth==nth/2) {
          TFile f(samples_[size-1]->GetListOfFiles()->At(nf)->GetTitle());
          total_entries_ += ((TTree*)f.Get(samples_[size-1]->GetListOfFiles()->At(nf)->GetName()))->GetEntries();
        }
      }
      step_size_ = 50000/total_entries_;
    }
  }
  
  size_t it_sample;
  bool LoopOnSamples() {
    ++it_sample;
    if (it_sample<samples_.size()) {
      files_ = samples_[it_sample]->GetListOfFiles();
      fileindex_=0;
      return 1;
    } else {
      return 0;
    }
  }
  
  bool LoopOnFiles() {
    if (fileindex_==0&&show_progress_) {
      //std::cout<<"--------- Started Looping on Trees ---------"<<std::flush;
      sw_->Start(0);
    }
    while (fileindex_%nthfile_!=(nthfile_-1)/2) ++fileindex_; // Skip files if specified (/2 instead 0 --> distribute files evenly)
    if (fileindex_<files_->GetEntries()) {
      //std::cout<<"Opening file: "<<files_->At(fileindex_)->GetTitle()<<std::endl;
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
