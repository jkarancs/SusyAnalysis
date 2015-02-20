#ifndef Samples_h
#define Samples_h

#include <vector>
#include <string>

class Samples {
  
public:
  
  Samples() { NSample_=0; }
  ~Samples() {}
  
private:
  
  size_t NSample_;
  std::vector<std::string> PFNames_;
  std::vector<std::string> LatexNames_;
  std::vector<std::vector<std::string> > DirNames_;
  std::vector<std::vector<double> > CrossSections_;
  std::vector<size_t > DirToIndex_;
  
public:
  
  typedef struct SubSample { std::string dir; double xsec_pb; } SubSample;
  void AddSample(std::string name, std::string latex, std::vector<SubSample> sample) {
    PFNames_.push_back(name);
    LatexNames_.push_back(latex);
    DirNames_.push_back(std::vector<std::string>());
    CrossSections_.push_back(std::vector<double>());
    for ( SubSample subsample : sample ) {
      DirToIndex_.push_back(NSample_);
      DirNames_[DirNames_.size()-1].push_back(subsample.dir);
      CrossSections_[CrossSections_.size()-1].push_back(subsample.xsec_pb);
    }
    ++NSample_;
  }
  
  std::vector<std::string> GetListOfDirectories() {
    std::vector<std::string> dirs;
    for ( auto subdirs : DirNames_ ) for ( auto subdir : subdirs ) dirs.push_back(subdir);
    return dirs;
  }
  
  std::vector<double> GetCrossSections() {
    std::vector<double> v_xsec;
    for ( auto subsample_xsecs : CrossSections_ ) for ( auto xsec : subsample_xsecs) v_xsec.push_back(xsec);
    return v_xsec;
  }
  
  std::vector<size_t > GetDirToIndex() { return DirToIndex_; }
  
  std::string GetPFNames() {
    std::string names;
    for ( std::string name : PFNames_ ) names += name+";";
    return names.substr(0, names.size()-1);
  }
  
  std::string GetLatexNames() {
    std::string names;
    for ( std::string name : LatexNames_ ) names += name+";";
    return names.substr(0, names.size()-1);
  }
};

#endif
