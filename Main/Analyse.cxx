#ifndef __MAIN_ANALYSE_CXX__
#define __MAIN_ANALYSE_CXX__

#include "Analyse.h"


namespace Main {


  void Analyse::SetBNBCosmicFile(std::string f) {
   
    mc_bnbcosmic_file_name = f;

  }

  void Analyse::SetInTimeCosmicFile(std::string f) {
   
    mc_intimecosmic_file_name = f;

  }

  void Analyse::SetDirtFile(std::string f) {
   
    mc_dirt_file_name = f;

  }

  void Analyse::SetBNBONFile(std::string f) {
   
    bnbon_file_name = f;

  }

  void Analyse::SetEXTBNBFile(std::string f) {
   
    extbnb_file_name = f;

  }

  void Analyse::SetBNBPOT(double v) {
   
    bnbon_pot_meas = v;

  }

  void Analyse::SetBNBONTriggers(double v) {
   
    bnbon_triggers = v;

  }

  void Analyse::SetEXTBNBTriggers(double v) {
   
    extbnb_triggers = v;

  }

  void Analyse::SetTargetFluxSystematic(std::string s) {
    _target_flux_syst = s;
  }

  void Analyse::SetPrefix(std::string p) {
    _prefix = p;
  }


  // void Analyse::AddExtraDiagonalUncertainty(TH2D & matrix, double frac_unc)
  // {

  //   for (int i = 0; i < matrix.GetNbinsX(); i++) {
  //     double current_content = matrix.GetBinContent(i+1, i+1);
  //     matrix.SetBinContent(i+1, i+1, current_content + )
  //   }

  // }




  void Analyse::DoAnalise() 
  {

  clock_t begin = clock();


  std::string analyser_outdir = std::getenv("MYSW_OUTDIR");
  analyser_outdir += "output_data_mc/";
  std::string mkdir_command = "mkdir -p ";
  mkdir_command += analyser_outdir;
  system(mkdir_command.c_str());

  LOG_NORMAL() << "Created output folder with name " << analyser_outdir << "." << std::endl;
  
  //gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  // *************************************
  // Opening files
  // *************************************
  TFile* mc_bnbcosmic_file = TFile::Open(mc_bnbcosmic_file_name.c_str(), "READ");
  TFile* mc_intimecosmic_file = TFile::Open(mc_intimecosmic_file_name.c_str(), "READ");
  TFile* mc_dirt_file = TFile::Open(mc_dirt_file_name.c_str(), "READ");
  TFile* bnbon_file = TFile::Open(bnbon_file_name.c_str(), "READ");
  TFile* extbnb_file = TFile::Open(extbnb_file_name.c_str(), "READ");

  if (!mc_dirt_file) 
    LOG_NORMAL() << "MC Dirt File not available. Will run without." << std::endl;
  

  
  // *************************************
  // Getting number of events for intimecosmic
  // *************************************
  //TH1D* h_nevts_intimecosmic = (TH1D*)mc_intimecosmic_file->Get("h_nevts");
  //intimecosmic_total_events = h_nevts_intimecosmic->GetBinContent(1);
  //std::cout << "Number of events (InTimeCosmic): " << intimecosmic_total_events << std::endl;
  
  LOG_NORMAL()<<"Got the histograms need to make plots " << std::endl;
  
  // *************************************
  // Getting the relevant histograms from MC file BNBCosmic
  // *************************************
  //std::map<std::string,TH1D*>* temp_map;
  //mc_bnbcosmic_file->GetObject("hmap_trkmom", temp_map);
  //std::map<std::string,TH1D*> hmap_trkmom_mc = *temp_map;

  


  // Create placeholders to get stuff from file
  std::map<std::string,std::map<std::string,TH1D*>>* temp_map_bs;
  BootstrapTH1D * temp_bs;

  mc_bnbcosmic_file->GetObject("hmap_trkmom_geant_pm1_bs", temp_map_bs);
  std::map<std::string,std::map<std::string,TH1D*>> map_bs = *temp_map_bs;

  mc_bnbcosmic_file->GetObject("hmap_trkcostheta_geant_pm1_bs", temp_map_bs);
  std::map<std::string,std::map<std::string,TH1D*>> map_bs_angle = *temp_map_bs;


  // Bootstrap efficiency - GEANT pm1sigma
  mc_bnbcosmic_file->GetObject("bs_geant_pm1_eff_mumom_num", temp_bs);
  BootstrapTH1D bs_genie_pm1_eff_mumom_num = *temp_bs;
  mc_bnbcosmic_file->GetObject("bs_geant_pm1_eff_mumom_den", temp_bs);
  BootstrapTH1D bs_genie_pm1_eff_mumom_den = *temp_bs;

  mc_bnbcosmic_file->GetObject("bs_geant_pm1_eff_mucostheta_num", temp_bs);
  BootstrapTH1D bs_genie_pm1_eff_mucostheta_num = *temp_bs;
  mc_bnbcosmic_file->GetObject("bs_geant_pm1_eff_mucostheta_den", temp_bs);
  BootstrapTH1D bs_genie_pm1_eff_mucostheta_den = *temp_bs;


  // Boostrap reco-true
  //std::map<std::string,TH2D*>* temp_map_bs2; 
  //mc_bnbcosmic_file->GetObject("bs_geant_pm1_true_reco_mom", temp_map_bs2); 
  //std::map<std::string,TH2D*> bs_true_reco_mom_mc = *temp_map_bs2;

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Create a list of the backgrounds that will be subtracted
  std::vector<std::string> bkg_names = {"chxbac", "reabac"};


  LOG_NORMAL()<<"Start to making the plots"<<std::endl;
  // Instantiate the GENIE reweighting plotter
  ReweightingPlotter genie_rw_plotter;

  if (_do_pm1sigma_plots) {

    // Bootstrap number of events per type
    std::map<std::string, BootstrapTH1D> bs;
    for (auto it : map_bs) {
      BootstrapTH1D temp;
      temp.SetAllHistograms(it.second);
      bs[it.first] = temp;
    }
    // Make +-1 sigma plots from GENIE
    genie_rw_plotter.SetEventBootstrapMap(bs);
    genie_rw_plotter.SetEfficiencyBootstraps(bs_genie_pm1_eff_mumom_num, bs_genie_pm1_eff_mumom_den);
    genie_rw_plotter.MakePlots(0, false, true);
    //genie_rw_plotter.MakePlots(2, false, true);
    genie_rw_plotter.MakeBackgroundPlots(0, false, true);  
  }

  if (_do_pm1sigma_plots) {

    // Bootstrap number of events per type
    std::map<std::string, BootstrapTH1D> bs;
    for (auto it : map_bs_angle) {
      BootstrapTH1D temp;
      temp.SetAllHistograms(it.second);
      bs[it.first] = temp;
    }
    // Make +-1 sigma plots from GENIE
    genie_rw_plotter.SetEventBootstrapMap(bs);
    genie_rw_plotter.SetEfficiencyBootstraps(bs_genie_pm1_eff_mucostheta_num, bs_genie_pm1_eff_mucostheta_den);
    genie_rw_plotter.MakePlots(1, false, true);
    //genie_rw_plotter.MakePlots(2, false, true);
    genie_rw_plotter.MakeBackgroundPlots(1, false, true);  
  }










  LOG_NORMAL() << "Preparing to load Event Histo from BNBCosmic file." << std::endl;

  UBXSecEventHisto1D * _event_histo_1d_mc = 0;
  //mc_bnbcosmic_file->GetObject("UBXSecEventHisto1D", _event_histo_1d_mc);


  Bool_t Object_exists;
  Object_exists=mc_bnbcosmic_file->GetListOfKeys()->Contains("UBXSecEventHisto1D");
  LOG_NORMAL()<<"Is this object exist ???? "<<Object_exists<<std::endl;
  mc_bnbcosmic_file->GetObject("UBXSecEventHisto1D", _event_histo_1d_mc);


  LOG_NORMAL() << "Event Histo correclty loaded from BNBCosmic file." << std::endl;




  gROOT->SetBatch(kTRUE);

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Create a list of the backgrounds that will be subtracted


  // Computing time<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
            
  //rootapp->Run();
  //rootapp->Terminate(0);
                 
  return;




  }//end of do analyze

  void Analyse::PlotMCTHStack(THStack *hs_trklen, std::map<std::string,TH1D*> themap, double scale_factor_mc_bnbcosmic)
{

  for (auto iter : themap) {
    iter.second->Scale(scale_factor_mc_bnbcosmic);
  }

  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["qe"]->SetLineColor(kGreen+2);
  themap["qe"]->SetFillColor(kGreen+2);
  themap["res"]->SetLineColor(kRed+1);
  themap["res"]->SetFillColor(kRed+1);
  themap["dis"]->SetLineColor(kBlue+2);
  themap["dis"]->SetFillColor(kBlue+2);
  themap["coh"]->SetLineColor(kMagenta+1);
  themap["coh"]->SetFillColor(kMagenta+1);
  themap["mec"]->SetLineColor(kOrange-3);
  themap["mec"]->SetFillColor(kOrange-3);

  hs_trklen->Add(themap["qe"]);
  hs_trklen->Add(themap["mec"]);
  hs_trklen->Add(themap["res"]);
  hs_trklen->Add(themap["coh"]);
  hs_trklen->Add(themap["dis"]);
  
  

  hs_trklen->Draw("hist");

  themap["total"]->Draw("E2 same");

  TLegend* leg2;

  TString xaxis_title = hs_trklen->GetXaxis()->GetTitle();
  if (xaxis_title.Contains("theta")) {
    leg2 = new TLegend(0.1747851,0.5431579,0.4355301,0.8231579,NULL,"brNDC");
    hs_trklen->SetMaximum(6200);
  } else if (xaxis_title.Contains("phi")) {
    leg2 = new TLegend(0.6203438,0.5905263,0.8810888,0.8705263,NULL,"brNDC");
    hs_trklen->SetMaximum(2500);
  } else if (xaxis_title.Contains("Energy")){
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
    hs_trklen->SetMaximum(3600);
  } else if (xaxis_title.Contains("Momentum")){
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
    hs_trklen->SetMaximum(4200);
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }

  std::stringstream sstm;

  sstm << "QE, " << std::setprecision(2)  << themap["qe"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["qe"],sstm.str().c_str(),"f");
  sstm.str("");

  sstm << "MEC, " << std::setprecision(2)  << themap["mec"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["mec"],sstm.str().c_str(),"f");
  sstm.str("");

  sstm << "RES, " << std::setprecision(2)  << themap["res"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["res"],sstm.str().c_str(),"f");
  sstm.str("");

  sstm << "COH, " << std::setprecision(2)  << themap["coh"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["coh"],sstm.str().c_str(),"f");
  sstm.str("");

  sstm << "DIS, " << std::setprecision(2)  << themap["dis"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["dis"],sstm.str().c_str(),"f");
  sstm.str("");

  

  

  leg2->AddEntry(themap["total"],"Stat. Unc. ","f");
  sstm.str("");

  leg2->Draw();

  return;

}

void Analyse::PrintFakeDataMessage() {
  for (int i = 0; i < 100; i++) {   
    std::cout << "****************************** RUNNING WITH FAKE DATA ******************************" << std::endl;
  }
}




}

#endif
