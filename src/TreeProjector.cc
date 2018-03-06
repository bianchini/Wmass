#include "../interface/TreeProjector.h"

using namespace std;
using namespace ROOT::Experimental;

TreeProjector::TreeProjector(const string& file_name="", const string& tree_name="", const string& out_name="", const int& num_of_cores=64){

  cout << "TreeProjector::TreeProjector()" << endl;

  // output file
  fout = TFile::Open(out_name.c_str(), "RECREATE");
  tout = new TTree("tree", "tree");
  
  // run parallel
  ROOT::EnableImplicitMT(num_of_cores); 

  // TDataFrame
  tdf = new ROOT::Experimental::TDataFrame(tree_name.c_str(), file_name.c_str()); 

  verbose = true;

  // flag
  status = 0;
}

TreeProjector::~TreeProjector(){
  cout << "TreeProjector::~TreeProjector()" << endl;
  //fout->Close();
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


vector<TDF::TResultProxy<TH1D>> TreeProjector::add_histo_with_cut(const vector<int>& weights, 
								  const string& cut_var_qt, const vector<float>& qt_max, 
								  const string& cut_var_y, const vector<float>& y_max, 
								  const string& plot_var, const vector<float>& pt_max){

  vector<TDF::TResultProxy<TH1D>> histos = {};

  // for string formatting
  char buffer[10];
  
  for(auto w : weights){
    auto get_weight = [w](Farray_t weights){return weights[w];};
    string weight_name = "weight_"+to_string(w);
    for(auto qt : qt_max){
      sprintf(buffer, "qt%.0f", qt);
      auto cut_qt = [qt](float x){ return x>qt; };
      string qt_name(buffer);
      for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
	float y_down = y_max[iy];
	float y_up = y_max[iy+1];
	auto cut_y = [y_down, y_up](float x){ return (x>=y_down && x<y_up); };
	sprintf(buffer, "y%.0f_%.0f", y_down, y_up);
	string y_name(buffer);
	for(auto pt : pt_max){
	  sprintf(buffer, "pt%.0f", pt);
	  auto cut_pt = [pt](float x){ return x<pt; };
	  string pt_name(buffer);
	  string histo_name = qt_name+"_"+y_name+"_"+pt_name+"_"+weight_name;
	  //histos.push_back( tdf->Filter(cut_var+">="+to_string(qt)).Filter(plot_var+"<="+to_string(pt)).Define(histo_name, "weights["+to_string(w)+"]").Histo1D( TH1D(("h_"+histo_name).c_str(),"", 200, 25.0, pt), plot_var, histo_name));
	  histos.push_back( tdf->Filter( cut_qt, {cut_var_qt} ).Filter( cut_y, {cut_var_y} ).Filter(cut_pt, {plot_var} ).Define(histo_name, get_weight, {"weights"} ).Histo1D( TH1D(("h_"+histo_name).c_str(),"", 200, 25.0, pt), plot_var, histo_name));
	}
      }
    }  
  }

  return histos;
}

int TreeProjector::run_test(){ 

  cout << "TreeProjector::run_test()" << endl; 

  cout << "Filling histos" << endl;
  TDF::TResultProxy<TH2D> h2_plus = tdf->Filter("mu_charge==-13").Histo2D( TH2D("h_Wplus_test","", 50, -2.5, 2.5, 50, 25.0, 55.0), "Wbare_mu_eta", "Wbare_mu_pt");
  TDF::TResultProxy<TH2D> h2_minus = tdf->Filter("mu_charge==+13").Histo2D( TH2D("h_Wminus_test","", 50, -2.5, 2.5, 50, 25.0, 55.0), "Wbare_mu_eta", "Wbare_mu_pt");

  cout << "Write histos" << endl;
  fout->mkdir("test")->cd();
  h2_plus->Write();
  h2_minus->Write();

  //////////////////
  vector<string> histos = {"Wplus", "Wminus"};
  for(auto hname : histos){
    TH2D* h2 = (TH2D*)fout->Get(("test/h_"+hname+"_test").c_str());
    cout << "Histo2D " << hname << " filled with integral: " << h2->Integral() << endl;
  }

  fout->Close();
    
  return 0;
}

int TreeProjector::run_pt_bias_vs_qt(){

  // start the clock....
  auto t1 = high_resolution_clock::now();

  out.time = static_cast<int>(duration_cast<milliseconds>(t2-t1).count());
  cout << "TreeProjector::run_pt_bias_vs_qt()" << endl; 

  char buffer[50];   

  vector<float> qt_max = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
  vector<float> y_max  = {0.0, 2.0, 4.0};
  vector<float> pt_max = {50.0, 60.0, 70.0};

  //vector<float> qt_max = {0.0, 10.0, 20.0};
  //vector<float> y_max  = {0.0, 1.0, 2.0};
  //vector<float> pt_max = {55.0};

  map<string, float> variables;
  vector<string> res = {"frac", "stat_err", "scale", "pdf", "tot"};
  for(auto qt : qt_max){
    for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
      float y_down = y_max[iy];
      float y_up = y_max[iy+1];
      for(auto pt : pt_max){
	for(auto v : res){
	  sprintf(buffer, "qt%.0f_y%.0f_%.0f_pt%0.f_%s", qt, y_down, y_up, pt, v.c_str());      
	  variables[string(buffer)] = {0.0};
	  //cout << string(buffer) << endl;
	  tout->Branch(buffer, &(variables[buffer]), (string(buffer)+"/F").c_str());
	}
      }
    }
  }

  cout << "Filling nominal histo" << endl;
  vector<int> weights_nominal= {0};
  vector<TDF::TResultProxy<TH1D>> histos_nominal = add_histo_with_cut( weights_nominal, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  cout << "Filling PDF" << endl;
  vector<int> weights_pdf = {}; for(int w = 9 ; w < 109; ++w) weights_pdf.push_back(w);
  //vector<int> weights_pdf = {9};
  std::vector<TDF::TResultProxy<TH1D>> histos_pdf = add_histo_with_cut( weights_pdf, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  cout << "Filling Scale" << endl;
  vector<int> weights_scale = {1,2,3,4,6,8};
  //vector<int> weights_scale = {1};
  std::vector<TDF::TResultProxy<TH1D>> histos_scale  = add_histo_with_cut( weights_scale, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  fout->mkdir("pt_bias_vs_qt")->cd();

  cout << "Writing..." << endl;
  for(auto h : histos_nominal) h->Write();
  for(auto h : histos_pdf) h->Write();
  for(auto h : histos_scale) h->Write();

  /////////////////////////////////////

  for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
    float y_down = y_max[iy];
    float y_up = y_max[iy+1];

    for(auto pt : pt_max){
      sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.0f_%.0f_pt%0.f_weight_%d", 0.0, y_down, y_up, pt, 0 );
      float norm = ((TH1D*)fout->Get(buffer))->Integral();
      
      for(auto qt : qt_max){

	sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.0f_%.0f_pt%0.f_weight_%d", qt, y_down, y_up, pt, 0 );
	TH1D* h_nominal = (TH1D*)fout->Get(buffer);
	float mean_nominal = h_nominal->GetMean();
	float mean_err_nominal = h_nominal->GetMeanError();
	float integ = h_nominal->Integral();
	
	std::vector<float> mu_pdf = { mean_nominal };      
	for(auto w : weights_pdf){
	  sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.0f_%.0f_pt%0.f_weight_%d", qt, y_down, y_up, pt, w );
	  TH1D* h = (TH1D*)fout->Get(buffer);
	  mu_pdf.push_back( h->GetMean() );
	}
	float rms_pdf = get_std(mu_pdf);
	
	std::vector<float> mu_scale = { mean_nominal };      
	for(auto w : weights_pdf){
	  sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.0f_%.0f_pt%0.f_weight_%d", qt, y_down, y_up, pt, w );
	  TH1D* h = (TH1D*)fout->Get(buffer);
	  mu_scale.push_back( h->GetMean() );
	}
	float rms_scale = get_max(mu_scale);
	float rms = TMath::Sqrt(rms_pdf*rms_pdf + rms_scale*rms_scale);      
	
	sprintf(buffer, "qt%.0f_y%.0f_%.0f_pt%0.f", qt, y_down, y_up, pt);
	variables[string(buffer)+"_frac"] = integ/norm;
	variables[string(buffer)+"_stat_err"] = mean_err_nominal/mean_nominal;
	variables[string(buffer)+"_scale"] = rms_scale/mean_nominal; 
	variables[string(buffer)+"_pdf"] = rms_pdf/mean_nominal; 
	variables[string(buffer)+"_tot"] = rms/mean_nominal; 
	
	if(verbose){
	  cout << "pt<=" << pt << ", qt>=" << qt << ", " << y_down << "< y <" << y_up << ": "
	       << "[frac=" << (integ/norm) << "], bias = " << rms/mean_nominal*(integ/norm) <<  " (stat = " << mean_err_nominal/mean_nominal << ")" << endl;      
	}
      }
    }    
  }

  tout->Fill();
  tout->Write();
  fout->Close();

  // stop the clock!
  auto t2 = high_resolution_clock::now();
  cout << "Done in " << static_cast<int>(duration_cast<milliseconds>(t2-t1).count()) << endl;
   
  return 0;
}
