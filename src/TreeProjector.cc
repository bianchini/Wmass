#include "../interface/TreeProjector.h"

using namespace std;
using namespace ROOT::Experimental;

TreeProjector::TreeProjector(const string& file_name="", const string& tree_name="", const string& out_name="", const int& num_of_cores=64){

  cout << "TreeProjector::TreeProjector()" << endl;

  // output file
  fin = nullptr;
  fout = TFile::Open(out_name.c_str(), "RECREATE");
  tout = new TTree("tree", "tree");
  
  // run parallel
  ROOT::EnableImplicitMT(num_of_cores); 

  // TDataFrame
  tdf = new ROOT::Experimental::TDataFrame(tree_name.c_str(), file_name.c_str()); 

  // flags
  status = {-1};
  count_filters = {0};
  verbose = {true};
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

int TreeProjector::run_test(){ 

  cout << "TreeProjector::run_test()" << endl; 

  // start the clock....
  auto t1 = high_resolution_clock::now();

  cout << "Filling histos" << endl;

  TDF::TResultProxy<TH2D> h2_plus = (*tdf)
    .Filter("mu_charge==-13")
    .Histo2D( TH2D("h_Wplus_test","", 50, -2.5, 2.5, 50, 25.0, 55.0), "Wbare_mu_eta", "Wbare_mu_pt"); 
  count_filters++;
  TDF::TResultProxy<TH2D> h2_minus = (*tdf)
    .Filter("mu_charge==+13")
    .Histo2D( TH2D("h_Wminus_test","", 50, -2.5, 2.5, 50, 25.0, 55.0), "Wbare_mu_eta", "Wbare_mu_pt");  
  count_filters++;

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

  // stop the clock!
  auto t2 = high_resolution_clock::now();
  cout << "Total number of filters = " << count_filters << endl;
  cout << "Done in " << static_cast<int>(duration_cast<minutes>(t2-t1).count()) << " min." << endl;
       
  return 0;
}

vector<TDF::TResultProxy<TH1D>> TreeProjector::plot_pt_with_qt_cut(const vector<int>& weights, 
								   const string& cut_var_qt, const vector<float>& qt_max, 
								   const string& cut_var_y, const vector<float>& y_max, 
								   const string& plot_var, const vector<float>& pt_max){
  
  // to be returned
  vector<TDF::TResultProxy<TH1D>> histos = {};

  // for string formatting
  char buffer[10];
  
  // We define here a temporary variable holding the weights collections to workaround
  // a performance degradation being fixed in ROOT proper.
  auto toVector = [](Farray_t weights){std::vector<float> v(weights.begin(), weights.end()); return v;};
  static auto tdf_cachedWeights = tdf->Define("weights_cached", toVector, {"weights"} );

  for(auto w : weights){

    // define weight[w]
    auto get_weight = [w](Farray_t &weights){return weights[w];};

    string weight_name = "weight_"+to_string(w);

    auto tdf_tmp = tdf_cachedWeights.Define(weight_name, get_weight, {"weights_cached"} );

    for(auto qt : qt_max){
      sprintf(buffer, "qt%.0f", qt);
      string qt_name(buffer);
      for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
	float y_down = y_max[iy];
	float y_up = y_max[iy+1];
	auto cut = [qt, y_down, y_up](float x, float y){ return (x >= qt && TMath::Abs(y)>=y_down && TMath::Abs(y)<y_up); };
	sprintf(buffer, "y%.1f_%.1f", y_down, y_up);
	string y_name(buffer);
	for(auto pt : pt_max){
	  sprintf(buffer, "pt%.0f", pt);
	  auto cut_pt = [pt](float x){ return x<pt; };
	  string pt_name(buffer);
	  string histo_name = qt_name+"_"+y_name+"_"+pt_name+"_"+weight_name;
	  histos.push_back( //(*tdf)
			   tdf_tmp
			   .Filter( cut_pt, {plot_var} )
			   .Filter( cut, {cut_var_qt, cut_var_y}, qt_name+"_"+y_name )
			   .Histo1D( TH1D(("h_"+histo_name).c_str(),"", 200, 25.0, pt), plot_var, weight_name));
	  count_filters++;
	}
      }
    }  
  }

  return histos;
}


int TreeProjector::run_pt_bias_vs_qt(){

  cout << "TreeProjector::run_pt_bias_vs_qt()" << endl; 

  // start the clock....
  auto t1 = high_resolution_clock::now();

  //vector<float> qt_max = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
  //vector<float> y_max  = {0.0, 5.0};
  //vector<float> pt_max = {55.0, 60.0, 65.0, 70.0};

  vector<float> qt_max = {0.0, 30.0};
  vector<float> y_max  = {0.0, 5.0};
  vector<float> pt_max = {55.0};

  map<string, float> variables;
  vector<string> res = {"frac", "stat_err", "scale", "pdf", "tot"};
  char buffer[50];   
  for(auto qt : qt_max){
    for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
      float y_down = y_max[iy];
      float y_up   = y_max[iy+1];
      for(auto pt : pt_max){
	for(auto v : res){
	  sprintf(buffer, "qt%.0f_y%.1f_%.1f_pt%0.f_%s", qt, y_down, y_up, pt, v.c_str());      
	  variables[string(buffer)] = {0.0};
	  tout->Branch(buffer, &(variables[buffer]), (string(buffer)+"/F").c_str());
	}
      }
    }
  }


  cout << " > Filling nominal" << endl;
  vector<int> weights_nominal= {0};
  vector<TDF::TResultProxy<TH1D>> histos_nominal = plot_pt_with_qt_cut( weights_nominal, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  cout << " > Filling PDF" << endl;
  //vector<int> weights_pdf = {}; for(int w = 9 ; w < 109; ++w) weights_pdf.push_back(w);
  vector<int> weights_pdf = {9};
  std::vector<TDF::TResultProxy<TH1D>> histos_pdf = plot_pt_with_qt_cut( weights_pdf, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  cout << " > Filling Scale" << endl;
  //vector<int> weights_scale = {1,2,3,4,6,8};
  vector<int> weights_scale = {1};
  std::vector<TDF::TResultProxy<TH1D>> histos_scale =plot_pt_with_qt_cut( weights_scale, "Wdress_qt", qt_max, "Wdress_y", y_max, "Wbare_mu_pt", pt_max);

  cout << " > Total number of filters = " << count_filters << endl;
  fout->mkdir("pt_bias_vs_qt")->cd();

  cout << " > Writing to disk..." << endl;
  for(auto h : histos_nominal) h->Write();
  for(auto h : histos_pdf)     h->Write();
  for(auto h : histos_scale)   h->Write();
  cout << " > Done. Fill output tree" << endl;

  /////////////////////////////////////

  for(unsigned int iy=0; iy< y_max.size()-1 ; ++iy){
    float y_down = y_max[iy];
    float y_up = y_max[iy+1];

    for(auto pt : pt_max){
      sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.1f_%.1f_pt%0.f_weight_%d", 0.0, y_down, y_up, pt, 0 );
      float norm = ((TH1D*)fout->Get(buffer))->Integral();
      
      for(auto qt : qt_max){

	sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.1f_%.1f_pt%0.f_weight_%d", qt, y_down, y_up, pt, 0 );
	TH1D* h_nominal = (TH1D*)fout->Get(buffer);
	float mean_nominal = h_nominal->GetMean();
	float mean_err_nominal = h_nominal->GetMeanError();
	float integ = h_nominal->Integral();
	
	std::vector<float> mu_pdf = { mean_nominal };      
	for(auto w : weights_pdf){
	  sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.1f_%.1f_pt%0.f_weight_%d", qt, y_down, y_up, pt, w );
	  TH1D* h = (TH1D*)fout->Get(buffer);
	  mu_pdf.push_back( h->GetMean() );
	}
	float rms_pdf = get_std(mu_pdf);
	
	std::vector<float> mu_scale = { mean_nominal };      
	for(auto w : weights_pdf){
	  sprintf(buffer, "pt_bias_vs_qt/h_qt%.0f_y%.1f_%.1f_pt%0.f_weight_%d", qt, y_down, y_up, pt, w );
	  TH1D* h = (TH1D*)fout->Get(buffer);
	  mu_scale.push_back( h->GetMean() );
	}
	float rms_scale = get_max(mu_scale);
	float rms = TMath::Sqrt(rms_pdf*rms_pdf + rms_scale*rms_scale);      
	
	sprintf(buffer, "qt%.0f_y%.1f_%.1f_pt%0.f", qt, y_down, y_up, pt);
	variables[string(buffer)+"_frac"] = integ/norm;
	variables[string(buffer)+"_stat_err"] = mean_err_nominal/mean_nominal;
	variables[string(buffer)+"_scale"] = rms_scale/mean_nominal; 
	variables[string(buffer)+"_pdf"] = rms_pdf/mean_nominal; 
	variables[string(buffer)+"_tot"] = rms/mean_nominal; 
	
	if(verbose){
	  cout << "Selection: [pt<=" << pt << ", qt>=" << qt << ", " << y_down << "<y<" << y_up << "] : "
	       << "[frac=" << (integ/norm) << "], dpT/pT = " << rms/mean_nominal*(integ/norm) <<  " (stat = " << mean_err_nominal/mean_nominal << ")" << endl;      
	}
      }
    }    
  }

  tdf->Report();

  tout->Fill();
  tout->Write();
  fout->Close();

  // stop the clock!
  auto t2 = high_resolution_clock::now();
  cout << "Done in " << static_cast<int>(duration_cast<seconds>(t2-t1).count()) << " sec." << endl;
   
  return 0;
}
