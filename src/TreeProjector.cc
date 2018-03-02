#include "../interface/TreeProjector.h"

using namespace std;
using namespace ROOT::Experimental;

TreeProjector::TreeProjector(const string& file_name="", const string& tree_name="", const string& out_name="", const int& num_of_cores=64){

  cout << "TreeProjector::TreeProjector()" << endl;

  // output file
  fout = TFile::Open(out_name.c_str(), "RECREATE");

  // run parallel
  ROOT::EnableImplicitMT(num_of_cores); 

  // TDataFrame
  tdf = new ROOT::Experimental::TDataFrame(tree_name.c_str(), file_name.c_str()); 

  all_weights = {};
  // flag
  status = 0;
}

TreeProjector::~TreeProjector(){
  cout << "TreeProjector::~TreeProjector()" << endl;
  delete tdf;
}


float TreeProjector::get_std(const vector<float>& v){
  float sum = std::accumulate(v.begin(), v.end(), 0.0);
  float mean = sum / v.size();
  std::vector<float> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(), [&mean](float x) { return x - mean; });
  float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  float stdev = std::sqrt(sq_sum / v.size());
  return stdev;
}

float TreeProjector::get_max(const vector<float>& v){
  float min = *std::min_element(v.begin(), v.end());
  float max = *std::max_element(v.begin(), v.end());
  return (max-min)*0.5;
}


vector<TDF::TResultProxy<TH1D>> TreeProjector::add_histo_withcut(const vector<int>& weights, const string& cut_var, const vector<float>& qt_max, const string& plot_var, const vector<float>& pt_max){

  vector<TDF::TResultProxy<TH1D>> histos = {};

  // for string formatting
  char buffer[10];
  
  for(auto w : weights){
    string weight_name = "weight_"+to_string(w);
    for(auto qt : qt_max){
      sprintf(buffer, "qt%.0f", qt);
      string qt_name(buffer);
      for(auto pt : pt_max){
	sprintf(buffer, "pt%.0f", pt);
	string pt_name(buffer);
	string histo_name = qt_name+"_"+pt_name+"_"+weight_name;
	histos.push_back( tdf->Filter(cut_var+">="+to_string(qt)).Filter(plot_var+"<="+to_string(pt)).Define(histo_name, "weights["+to_string(w)+"]").Histo1D( TH1D(("h_"+histo_name).c_str(),"", 200, 25.0, pt), plot_var, histo_name));
      }
    }  
  }

  return histos;
}

int TreeProjector::run(){

  //vector<float> qt_max = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
  vector<float> qt_max = {0.0, 10.0, 20.0, 30.0};
  //vector<float> pt_max = {50.0, 55.0, 60.0, 65.0, 70.0};
  vector<float> pt_max = {55.0};

  cout << "Filling nominal histo" << endl;
  vector<int> weights_nominal= {0};
  vector<TDF::TResultProxy<TH1D>> histos_nominal = add_histo_withcut( weights_nominal, "Wdress_qt", qt_max, "Wbare_mu_pt", pt_max);

  cout << "Filling PDF" << endl;
  vector<int> weights_pdf = {}; for(int w = 9 ; w < 20; ++w) weights_pdf.push_back(w);
  std::vector<TDF::TResultProxy<TH1D>> histos_pdf = add_histo_withcut( weights_pdf, "Wdress_qt", qt_max, "Wbare_mu_pt", pt_max);

  cout << "Filling Scale" << endl;
  vector<int> weights_scale = {1,2,3,4,6,8};
  std::vector<TDF::TResultProxy<TH1D>> histos_scale  = add_histo_withcut( weights_scale, "Wdress_qt", qt_max, "Wbare_mu_pt", pt_max);

  fout->cd();

  cout << "Writing..." << endl;
  for(auto h : histos_nominal) h->Write();
  for(auto h : histos_pdf) h->Write();
  for(auto h : histos_scale) h->Write();

  /////////////////////////////////////

  char buffer[30];   
  for(auto pt : pt_max){
    sprintf(buffer, "h_qt%.0f_pt%0.f_weight_%d", 0.0, pt, 0 );
    float norm = ((TH1D*)fout->Get(buffer))->Integral();

    for(auto qt : qt_max){

      sprintf(buffer, "h_qt%.0f_pt%0.f_weight_%d", qt, pt, 0 );
      TH1D* h_nominal = (TH1D*)fout->Get(buffer);
      float mean_nominal = h_nominal->GetMean();
      float mean_err_nominal = h_nominal->GetMeanError();
      float integ = h_nominal->Integral();

      std::vector<float> mu_pdf = { mean_nominal };      
      for(auto w : weights_pdf){
	sprintf(buffer, "h_qt%.0f_pt%0.f_weight_%d", qt, pt, w );
	TH1D* h = (TH1D*)fout->Get(buffer);
	mu_pdf.push_back( h->GetMean() );
      }
      float rms_pdf = get_std(mu_pdf);

      std::vector<float> mu_scale = { mean_nominal };      
      for(auto w : weights_pdf){
	sprintf(buffer, "h_qt%.0f_pt%0.f_weight_%d", qt, pt, w );
	TH1D* h = (TH1D*)fout->Get(buffer);
	mu_scale.push_back( h->GetMean() );
      }
      float rms_scale = get_max(mu_scale);
      float rms = TMath::Sqrt(rms_pdf*rms_pdf + rms_scale*rms_scale);      
      cout << "pt<=" << pt << ", qt>=" << qt << ": " << " [frac=" << (integ/norm) << "] " << rms_pdf/mean_nominal*(integ/norm) << ", " << rms/mean_nominal << " (" << mean_err_nominal/mean_nominal << ")" << endl;      

    }    
  }


  fout->Close();
  
  
  return 0;
}
