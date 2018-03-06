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
#include<chrono>

using Farray_t = ROOT::Experimental::TDF::TArrayBranch<float>;
using Iarray_t = ROOT::Experimental::TDF::TArrayBranch<int>;

using namespace std;
using namespace ROOT::Experimental;
using namespace std::chrono;

class TreeProjector {

 public:
  TreeProjector(const std::string&, const std::string&, const std::string&, const int&);
  vector<TDF::TResultProxy<TH1D>> plot_pt_with_qt_cut(const vector<int>&, const string&, const vector<float>&, const string&, const vector<float>&, const string&, const vector<float>&);
  float get_std(const vector<float>&);
  float get_max(const vector<float>&);
  int run_pt_bias_vs_qt();
  int run_test();
  ~TreeProjector();
  
 private:
  bool verbose;
  int status;
  int count_filters;
  ROOT::Experimental::TDataFrame* tdf;
  TFile* fout;
  TTree* tout;
};

#endif
