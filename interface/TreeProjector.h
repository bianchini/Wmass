#ifndef TreeProjector_h
#define TreeProjector_h

#include <stdio.h>
#include <TTree.h>
#include <TFile.h>
#include <ROOT/TDataFrame.hxx>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <string>

using Farray_t = ROOT::Experimental::TDF::TArrayBranch<float>;
using Iarray_t = ROOT::Experimental::TDF::TArrayBranch<int>;

using namespace std;
using namespace ROOT::Experimental;

class TreeProjector {

 public:
  TreeProjector(const std::string&, const std::string&, const std::string&, const int&);
  vector<TDF::TResultProxy<TH1D>> add_histo_withcut(const vector<int>&, const string&, const vector<float>&, const string&, const vector<float>&);
  float get_std(const vector<float>&);
  float get_max(const vector<float>&);
  int run();
  ~TreeProjector();
  
 private:
  int status;
  ROOT::Experimental::TDataFrame* tdf;
  TFile* fout;
  vector<int> all_weights;
};

#endif
