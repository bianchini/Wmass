{
  ROOT::EnableImplicitMT(64);
  ROOT::RDataFrame d("tree", "/scratch/bertacch/wmass/pdf_peak_shift/MC_data/MC_Wjet2lnu.root");

  //"/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180220_102334/0000/tree*root",
  //"/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180220_102024/0000/tree*root"}

  using floats = ROOT::VecOps::RVec<float>;

  auto Wplus   = [](float x) {return x==-13;};
  auto Wminus  = [](float x) {return x==+13;};
  auto acc_rec = [](float x, float y) {return (TMath::Abs(y)<2.5 && x>25.0);};
  auto acc_gen = [](float x, float y) {return (x>20.0 && TMath::Abs(y)<10.0);};
  auto energy  = [](float x, float y) { return x*TMath::CosH(y); };
  auto abspz  = [](float x, float y) { return TMath::Abs(x*TMath::SinH(y)); };
  auto abspy  = [](float x, float y) { return TMath::Abs(x*TMath::Sin(y)); };
  auto abspx  = [](float x, float y) { return TMath::Abs(x*TMath::Cos(y)); };

  int nbins   = 100;
  float xlow  = 35.;
  float xhigh = 45.;
  //float xlow  = 0.;
  //float xhigh = 200.;

  auto d2 = d
    //.Range(0,10000)
    .Filter(acc_rec, {"WpreFSR_mu_pt","WpreFSR_mu_eta"})    
    .Define("WpreFSR_mu_e",    energy, {"WpreFSR_mu_pt","WpreFSR_mu_eta"})
    //.Define("WpreFSR_mu_abspz", abspz, {"WpreFSR_mu_pt","WpreFSR_mu_eta"})
    //.Define("WpreFSR_mu_abspx", abspx, {"WpreFSR_mu_pt","WpreFSR_mu_phi"})
    //.Define("WpreFSR_mu_abspy", abspy, {"WpreFSR_mu_pt","WpreFSR_mu_phi"});
    .Define("w",  [](floats &weights) {return weights[0];}, {"weights"})
    .Define("count",  []() {return 0.5;});

  //std::vector<int> weights = {0,1,2,3,4,5,6,7,8};
  std::vector<double> weights = {80.119, 80.169, 80.219, 80.269, 80.319, 80.369, 80.419, 80.469, 80.519, 80.569, 80.619, 80.669, 80.719};
  //std::vector<double> weights = {};
  //for(int i = 0; i < 9; ++i)   weights.push_back(i);
  //for(int i = 9; i < 109; ++i) weights.push_back(i);
  //for(int i = 9; i < 109; ++i) weights.push_back(i);

  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos;

  TFile* fout = TFile::Open("out_pt_weights_mass_acc.root","RECREATE");

  histos.push_back( d2.Histo1D<double>( {"norm", "", 1,0,1}, "count", "w"  ) );

  for(unsigned int iw=0; iw<weights.size(); ++iw){

    float w = weights[iw];
    //auto dw = d2.Define(std::string(Form("w%i",w)),  [w](floats &weights) {return weights[w];}, {"weights"});
    auto dw = d2.Define(std::string(Form("w%i",iw)),  [w](float &mass, float &weight) {return (weight*(TMath::Power(mass*mass-80.419*80.419,2) + 80.419*80.419*2.0476*2.0476)/(TMath::Power(mass*mass-w*w, 2) + w*w*2.0476*2.0476));}, {"WpreFSR_mass", "w"});

    histos.push_back( dw.Filter(Wplus, {"mu_charge"}).Histo1D<double>({("WpreFSR_mu_e_plus_all_"+std::string(Form("w%i",iw))).c_str(),"Energy", nbins ,xlow, xhigh }, "WpreFSR_mu_e", std::string(Form("w%i",iw))));
    histos.push_back( dw.Filter(Wplus, {"mu_charge"}).Histo1D<double>({("WpreFSR_mu_pt_plus_all_"+std::string(Form("w%i",iw))).c_str(),"Pt", nbins ,xlow, xhigh },    "WpreFSR_mu_pt", std::string(Form("w%i",iw))));
    
    histos.push_back( dw.Filter(Wminus, {"mu_charge"}).Histo1D<double>({("WpreFSR_mu_e_minus_all_"+std::string(Form("w%i",iw))).c_str(),"Energy", nbins ,xlow, xhigh }, "WpreFSR_mu_e", std::string(Form("w%i",iw))));
    histos.push_back( dw.Filter(Wminus, {"mu_charge"}).Histo1D<double>({("WpreFSR_mu_pt_minus_all_"+std::string(Form("w%i",iw))).c_str(),"Pt", nbins ,xlow, xhigh },    "WpreFSR_mu_pt", std::string(Form("w%i",iw))));

  }
  
  fout->cd();
  for(auto h : histos) h->Write("", TObject::kOverwrite);
  fout->Close();

}
