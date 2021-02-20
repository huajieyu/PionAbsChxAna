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
  std::vector<std::string> bkg_names = {"chxbac", "reabac", "other"};

  TH1D *h_Evttot_mc = (TH1D*)mc_bnbcosmic_file->Get("h_Evttot");
  TH1D *h_Evttot_data = (TH1D*)bnbon_file->Get("h_Evttot");

  Double_t Nevtmc= h_Evttot_mc->GetBinContent(1);
  Double_t Nevtdata = h_Evttot_data->GetBinContent(1);
 
  double scale_fac = Nevtdata/(Nevtmc*1.0);

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
    //calculate smeared efficiecy
    //get the histograms of efficiency

    int n_bins_mumom = 50;
    int n_bins_mucostheta = 50;
    int n_bins_muphi=50;

    double bins_mumom[51]; 
    double bins_mucostheta[51]; 
    double bins_muphi[51];

    for(int i=0; i<n_bins_mumom+1; i++){
         bins_mumom[i]=0.0+1.2/n_bins_mumom*i;
    }
    for(int i=0; i<n_bins_mucostheta+1; i++){
         bins_mucostheta[i]=-1.0+2.0/n_bins_mucostheta*i;
    }
    for(int i=0; i<n_bins_muphi+1; i++){
         bins_muphi[i]=-TMath::Pi()+2.0*TMath::Pi()/n_bins_muphi*i;
    }
    //Set the map with the signal and background histograms
    std::map<std::string, TH1D*> hmap_trkmom;
    hmap_trkmom["total"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sel_pmom");
    hmap_trkmom["signal"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sig_pmom");
    hmap_trkmom["chxbac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_chxbac_pmom");
    hmap_trkmom["reabac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_reabac_pmom");
    hmap_trkmom["other"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_other_pmom");

    std::map<std::string, TH1D*> hmap_trkmom_data;
    hmap_trkmom_data["total"]= (TH1D*)bnbon_file->Get("h_PiAbs_sel_pmom");


    std::map<std::string, TH1D*> hmap_trkcostheta;
    hmap_trkcostheta["total"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sel_pcostheta");
    hmap_trkcostheta["signal"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sig_pcostheta");
    hmap_trkcostheta["chxbac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_chxbac_pcostheta");
    hmap_trkcostheta["reabac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_reabac_pcostheta");
    hmap_trkcostheta["other"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_other_pcostheta");

    std::map<std::string, TH1D*> hmap_trkcostheta_data;
    hmap_trkcostheta_data["total"]= (TH1D*)bnbon_file->Get("h_PiAbs_sel_pcostheta");
 

    std::map<std::string, TH1D*> hmap_trkphi;
    hmap_trkphi["total"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sel_pphi");
    hmap_trkphi["signal"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_sig_pphi");
    hmap_trkphi["chxbac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_chxbac_pphi");
    hmap_trkphi["reabac"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_reabac_pphi");
    hmap_trkphi["other"]= (TH1D*)mc_bnbcosmic_file->Get("h_PiAbs_other_pphi");

    std::map<std::string, TH1D*> hmap_trkphi_data;
    hmap_trkphi_data["total"]= (TH1D*)bnbon_file->Get("h_PiAbs_sel_pphi");
 


    TH1D* h_eff_pmom_num; TH1D* h_eff_pmom_den;
    TH1D* h_eff_pcostheta_num; TH1D* h_eff_pcostheta_den;
    TH1D* h_eff_pphi_num; TH1D* h_eff_pphi_den;

    TH2D* h_true_reco_mom;
    TH2D* h_true_reco_costheta;
    TH2D* h_true_reco_phi;
    h_eff_pmom_den = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pmom_den");
    h_eff_pmom_num = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pmom_num");

    h_eff_pcostheta_den = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pcostheta_den");
    h_eff_pcostheta_num = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pcostheta_num");

    h_eff_pphi_den = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pphi_den");
    h_eff_pphi_num = (TH1D*)mc_bnbcosmic_file->Get("h_eff_pphi_num");

    h_true_reco_mom = (TH2D*)mc_bnbcosmic_file->Get("h_true_reco_mom");
    h_true_reco_costheta = (TH2D*)mc_bnbcosmic_file->Get("h_true_reco_costheta");
    h_true_reco_phi = (TH2D*)mc_bnbcosmic_file->Get("h_true_reco_phi");

    std::string varcalc[3] = {"mom", "costheta", "phi"}; 
    //Set Migration Matrix and print out
    TLatex *prelim_Right = new TLatex(0.9, 0.93, "protoDUNE Simulation");
 
    prelim_Right->SetTextFont(62);
    prelim_Right->SetTextColor(kGray+1); 
    prelim_Right->SetNDC(); 
    prelim_Right->SetTextSize(0.05); 
    prelim_Right->SetTextAlign(32); 
    //prelim->Draw();

    TLatex* prelim_Left = new TLatex(0.4,0.8, "protoDUNE Simulation");
    prelim_Left->SetTextFont(62);
    prelim_Left->SetTextColor(kGray+1); 
    prelim_Left->SetNDC(); 
    prelim_Left->SetTextSize(0.05); 
    prelim_Left->SetTextAlign(32); 
 
 
  if(_do_reweighting_plots){ 

    //
    // Proton Momentum  Momentum Cross Section
    //
    TMatrix _S_mom; _S_mom.Clear(); _S_mom.ResizeTo(n_bins_mumom + 1, n_bins_mumom + 1);
    MigrationMatrix2D migrationmatrix2d;
    migrationmatrix2d.SetOutDir("migration_matrix_2d_trkmom");
    migrationmatrix2d.SetNBins(n_bins_mumom, n_bins_mumom);
    migrationmatrix2d.SetTrueRecoHistogram(h_true_reco_mom);
    _S_mom = migrationmatrix2d.CalculateMigrationMatrix();
    migrationmatrix2d.PlotMatrix();
    migrationmatrix2d.SetOutputFileName("migration_matrix_2d_trkmom.tex");
    migrationmatrix2d.PrintSmearingMatrixLatex();

 
    TCanvas *c_mom=new TCanvas("c_mom", "c_mom", 900, 900);
    h_true_reco_mom->Draw("colz");


    int bins_x = h_true_reco_mom->GetNbinsX()+1;
    int bins_y = h_true_reco_mom->GetNbinsY()+1;
    int _m = bins_x - 1;
    int _n = bins_y - 1;
    TH2D * smearing_matrix_mom= new TH2D("smearing_matrix_mom", "", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
    smearing_matrix_mom->GetXaxis()->SetTitle("p_{#mu} (Truth) [GeV/c]");
    smearing_matrix_mom->GetYaxis()->SetTitle("p_{#mu} (Reco) [GeV/c]");
    smearing_matrix_mom->GetYaxis()->SetTitleOffset(1.4);
    for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_mom->SetBinContent(j+1, i+1, _S_mom[i][j]);

      }
    }
    //smearing_matrix_mom->GetXaxis()->SetRangeUser(0.0,2.5);
    //smearing_matrix_mom->GetYaxis()->SetRangeUser(0.0,2.5);
    smearing_matrix_mom->Draw("colz");
    prelim_Right->Draw("same");
    c_mom->SaveAs("Fractional_MigrationMatrix_trkmom.png");
    //------------------------------------------------------------------------------ 
    //Get the smeared efficiency of proton momentum
    CrossSectionCalculator1D _xsec_calc_mom;
    _xsec_calc_mom.Reset();
    _xsec_calc_mom.SetVerbose(true);
    _xsec_calc_mom.SetOutDir("output_data_mc");
    _xsec_calc_mom.SetNameAndLabel("mom", ";P_{proton}[GeV]; Selected Events");
    _xsec_calc_mom.SetMigrationMatrix(_S_mom);
    _xsec_calc_mom.SetTruthHistograms(h_eff_pmom_num, h_eff_pmom_den);
    _xsec_calc_mom.SetHistograms(hmap_trkmom, hmap_trkmom["total"], hmap_trkmom["total"]);
     _xsec_calc_mom.SetHistograms(hmap_trkmom, hmap_trkmom["total"],  hmap_trkmom["total"], hmap_trkmom, hmap_trkmom["total"]);
    _xsec_calc_mom.Smear(n_bins_mumom, n_bins_mumom);
 
    TH1D * xsec_mumom = _xsec_calc_mom.ExtractCrossSection_PiAbs(bkg_names, "p_{proton}^{reco} [GeV]", "mb", scale_fac);

 
    //-------------------------------------------------------------------------------
    TMatrix _S_costheta; _S_costheta.Clear(); _S_costheta.ResizeTo(n_bins_mucostheta + 1, n_bins_mucostheta + 1);
    //MigrationMatrix2D migrationmatrix2d;
    migrationmatrix2d.SetOutDir("migration_matrix_2d_trkcostheta");
    migrationmatrix2d.SetNBins(n_bins_mucostheta, n_bins_mucostheta);
    migrationmatrix2d.SetTrueRecoHistogram(h_true_reco_costheta);
    _S_costheta = migrationmatrix2d.CalculateMigrationMatrix();
    migrationmatrix2d.PlotMatrix();
    migrationmatrix2d.SetOutputFileName("migration_matrix_2d_trkcostheta.tex");
    migrationmatrix2d.PrintSmearingMatrixLatex();
    
 

    TCanvas *c_costheta=new TCanvas("c_costheta", "c_costheta", 900, 900);
    h_true_reco_costheta->Draw("colz");


    bins_x = h_true_reco_costheta->GetNbinsX()+1;
    bins_y = h_true_reco_costheta->GetNbinsY()+1;
    _m = bins_x - 1;
    _n = bins_y - 1;
    TH2D * smearing_matrix_costheta= new TH2D("smearing_matrix_costheta", "", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
    smearing_matrix_costheta->GetXaxis()->SetTitle("cos#theta_{#mu} (Truth)");
    smearing_matrix_costheta->GetYaxis()->SetTitle("cos#theta_{#mu} (Reco)");
    smearing_matrix_costheta->GetYaxis()->SetTitleOffset(1.4);
    for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_costheta->SetBinContent(j+1, i+1, _S_costheta[i][j]);

      }
    }
    smearing_matrix_costheta->Draw("colz");
    prelim_Right->Draw("same");
    c_costheta->SaveAs("Fractional_MigrationMatrix_trkcostheta.png");
    //------------------------------------------------------------------------------ 
    //Get the smeared efficiency of proton momentum
    CrossSectionCalculator1D _xsec_calc_costheta;
    _xsec_calc_costheta.Reset();
    _xsec_calc_costheta.SetVerbose(true);
    _xsec_calc_costheta.SetOutDir("output_data_mc");
    _xsec_calc_costheta.SetNameAndLabel("costheta", ";cos#theta_{proton}; Selected Events");
    _xsec_calc_costheta.SetMigrationMatrix(_S_costheta);
    _xsec_calc_costheta.SetTruthHistograms(h_eff_pcostheta_num, h_eff_pcostheta_den);
    _xsec_calc_costheta.Smear(n_bins_mucostheta, n_bins_mucostheta);
 
    //_xsec_calc_costheta.SetBkgToSubtract(bkg_names);

    //------------------------------------------------------------

    TMatrix _S_phi; _S_phi.Clear(); _S_phi.ResizeTo(n_bins_muphi + 1, n_bins_muphi + 1);
    //MigrationMatrix2D migrationmatrix2d;
    migrationmatrix2d.SetOutDir("migration_matrix_2d_trkphi");
    migrationmatrix2d.SetNBins(n_bins_muphi, n_bins_muphi);
    migrationmatrix2d.SetTrueRecoHistogram(h_true_reco_phi);
    _S_phi = migrationmatrix2d.CalculateMigrationMatrix();
    migrationmatrix2d.PlotMatrix();
    migrationmatrix2d.SetOutputFileName("migration_matrix_2d_trkphi.tex");
    migrationmatrix2d.PrintSmearingMatrixLatex();


    TCanvas *c_phi=new TCanvas("c_phi", "c_phi", 900, 900);
    h_true_reco_phi->Draw("colz");


    bins_x = h_true_reco_phi->GetNbinsX()+1;
    bins_y = h_true_reco_phi->GetNbinsY()+1;
    _m = bins_x - 1;
    _n = bins_y - 1;
    TH2D * smearing_matrix_phi= new TH2D("smearing_matrix_phi", "", n_bins_muphi, bins_muphi, n_bins_muphi, bins_muphi);
    smearing_matrix_phi->GetXaxis()->SetTitle("#phi_{proton} (Truth)");
    smearing_matrix_phi->GetYaxis()->SetTitle("#phi_{proton} (Reco)");
    smearing_matrix_phi->GetYaxis()->SetTitleOffset(1.4);
    for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_phi->SetBinContent(j+1, i+1, _S_phi[i][j]);

      }
    }
    smearing_matrix_phi->Draw("colz");
    prelim_Right->Draw("same");
    c_phi->SaveAs("Fractional_MigrationMatrix_trkphi.png");
    //------------------------------------------------------------------------------ 
    //Get the smeared efficiency of proton momentum
    CrossSectionCalculator1D _xsec_calc_phi;
    _xsec_calc_phi.Reset();
    _xsec_calc_phi.SetVerbose(true);
    _xsec_calc_phi.SetOutDir("output_data_mc");
    _xsec_calc_phi.SetNameAndLabel("phi", ";#phi_{proton}; Selected Events");
    _xsec_calc_phi.SetMigrationMatrix(_S_phi);
    _xsec_calc_phi.SetTruthHistograms(h_eff_pphi_num, h_eff_pphi_den);
    _xsec_calc_phi.Smear(n_bins_muphi, n_bins_muphi);
 
    //_xsec_calc_phi.SetBkgToSubtract(bkg_names);
  }//end of if _do_reweighting_plots
    //------------------------------------------------------------------------------------------

   LOG_NORMAL()<<"Start to get the basic graphs and do calculation "<<std::endl;

   TGraph* gr_selected_tot_slc = (TGraph*)mc_dirt_file->Get("gr_selected_tot_slc");
   TGraph* gr_selected_pibkg_slc = (TGraph*)mc_dirt_file->Get("gr_selected_pibkg_slc");
   TGraph* gr_selected_mubkg_slc = (TGraph*)mc_dirt_file->Get("gr_selected_mubkg_slc");
   TGraph* gr_selected_pibkg_elastic_slc = (TGraph*)mc_dirt_file->Get("gr_selected_pibkg_elastic_slc");


   TGraph* gr_intabs_slc = (TGraph*)mc_dirt_file->Get("gr_intabs_slc");
   TGraph* gr_recoabs_slc = (TGraph*)mc_dirt_file->Get("gr_recoabs_slc");

   TGraph* gr_intabs_sel_slc = (TGraph*)mc_dirt_file->Get("gr_intabs_sel_slc");
   TGraph* gr_recoabs_sel_slc = (TGraph*)mc_dirt_file->Get("gr_recoabs_sel_slc");
   LOG_NORMAL()<<"Successfully got the graphs for further calculation"<<std::endl;
   Double_t intabs[nslices+1]; Double_t intabs_sel[nslices+1]; 
   Double_t recoabs[nslices+1]; Double_t recoabs_sel[nslices+1]; 

   Int_t recoabs_nbins = gr_recoabs_slc->GetN();

   Double_t slcid[nslices+1];
   Double_t rslcid[nslices+1];
   Double_t sslcid[nslices+1];
   Double_t tslcid[nslices+1];
   Double_t uslcid[nslices+1];


   Double_t Purity_reco[nslices+1];

   Double_t Purity_reco_err[nslices+1];


   Double_t selected_tot[nslices+1];
   Double_t selected_pibkg[nslices+1];
   Double_t selected_mubkg[nslices+1];
   Double_t selected_pibkg_elastic[nslices+1];
   Double_t Err_x[nslices+1];

   TH1D *h_Purity=new TH1D("h_Purity", "h_Purity", nslices+1, -0.5, 24.5);
   TH1D *h_Selected=new TH1D("h_Selected", "h_Selected", nslices+1, -0.5, 24.5);
   TH1D *h_Selected_BKG=new TH1D("h_Selected_BKG", "h_Selected_BKG", nslices+1, -0.5, 24.5);

   Double_t total_mubkg=0.0;
   Double_t total_pibkg=0.0;   
   Double_t total_pibkg_elastic = 0.0;
   Double_t total_signal_true=0.0;
   Double_t total_signal_reco=0.0;
   Double_t total_sel=0.0;
   for(int ind=0; ind<recoabs_nbins; ind++){ 
            Err_x[ind]=0;
            gr_recoabs_slc->GetPoint(ind, slcid[ind], recoabs[ind]);   
            gr_intabs_slc->GetPoint(ind, rslcid[ind], intabs[ind]);   

            gr_recoabs_sel_slc->GetPoint(ind, uslcid[ind], recoabs_sel[ind]);
            gr_intabs_sel_slc->GetPoint(ind, rslcid[ind], intabs_sel[ind]);   

            gr_selected_tot_slc->GetPoint(ind, sslcid[ind], selected_tot[ind]);
            h_Selected->SetBinContent(ind+1, selected_tot[ind]);
            h_Selected->SetBinError(ind+1, TMath::Sqrt(selected_tot[ind]));
 
            gr_selected_pibkg_slc->GetPoint(ind, tslcid[ind], selected_pibkg[ind]);
            gr_selected_mubkg_slc->GetPoint(ind, tslcid[ind], selected_mubkg[ind]);
            gr_selected_pibkg_elastic_slc->GetPoint(ind, tslcid[ind], selected_pibkg_elastic[ind]);

            h_Selected_BKG->SetBinContent(ind+1, selected_mubkg[ind]+selected_pibkg[ind]);
            h_Selected_BKG->SetBinError(ind+1, TMath::Sqrt(selected_mubkg[ind])+selected_pibkg[ind]);

            total_sel +=selected_tot[ind];
            total_pibkg +=selected_pibkg[ind];
            total_pibkg_elastic +=selected_pibkg_elastic[ind];
            total_mubkg +=selected_mubkg[ind];
            total_signal_true +=intabs_sel[ind];
            total_signal_reco +=recoabs_sel[ind];
            Purity_reco[ind]=recoabs_sel[ind]/selected_tot[ind];

            Purity_reco_err[ind] = Purity_reco[ind]*TMath::Sqrt(1/selected_tot[ind]+1/recoabs_sel[ind]);
            h_Purity->SetBinContent(ind+1, Purity_reco[ind]);
            h_Purity->SetBinError(ind+1, Purity_reco_err[ind]);
   }

   LOG_NORMAL()<<"Total Selected "<<total_sel<<std::endl;
   LOG_NORMAL()<<"Total Pion Background is "<<total_pibkg<<" Total Muon Background is "<<total_mubkg<<"   Total Pion Elastic Background is "<<total_pibkg_elastic<<std::endl;
   LOG_NORMAL()<<"Total True Signal is "<<total_signal_true<<" Total Reco Signal is "<<total_signal_reco<<std::endl;


   TH1D *h_efficiency_reco_num = new TH1D("h_efficiency_reco_num", "h_efficiency_reco_num", nslices+1, -0.5, 24.5);
   TH1D *h_efficiency_reco_den = new TH1D("h_efficiency_reco_den", "h_efficiency_reco_den", nslices+1, -0.5, 24.5);
    
   TH1D *h_efficiency_true_num = new TH1D("h_efficiency_true_num", "h_efficiency_true_num", nslices+1, -0.5, 24.5);
   TH1D *h_efficiency_true_den = new TH1D("h_efficiency_true_den", "h_efficiency_true_den", nslices+1, -0.5, 24.5);
 
   for(int i=0; i<nslices+1; i++){
         h_efficiency_reco_num->SetBinContent(i+1, recoabs_sel[i]);
         h_efficiency_reco_num->SetBinError(i+1, TMath::Sqrt(recoabs_sel[i]));
         h_efficiency_reco_den->SetBinContent(i+1, recoabs[i]);
         h_efficiency_reco_den->SetBinError(i+1, TMath::Sqrt(recoabs[i]));
 
         h_efficiency_true_num->SetBinContent(i+1, intabs_sel[i]);
         h_efficiency_true_num->SetBinError(i+1, TMath::Sqrt(intabs_sel[i]));
         h_efficiency_true_den->SetBinContent(i+1, intabs[i]);
         h_efficiency_true_den->SetBinError(i+1, TMath::Sqrt(intabs[i]));
   }

   h_efficiency_reco_num->Sumw2();
   h_efficiency_reco_den->Sumw2();

   h_efficiency_true_num->Sumw2();
   h_efficiency_true_den->Sumw2();

   std::string outFileName = "output_calcualted.root";
   TFile *outputFile;
   outputFile = new TFile(outFileName.c_str(),  "RECREATE");

   h_efficiency_reco_num->Write();
   h_efficiency_reco_den->Write();

   h_efficiency_true_num->Write();
   h_efficiency_true_den->Write();

   h_Purity->Write();
   h_Selected->Write();
   h_Selected_BKG->Write();


   TEfficiency *pEff_reco=new TEfficiency(*h_efficiency_reco_num, *h_efficiency_reco_den);
   pEff_reco->SetStatisticOption(TEfficiency::kFNormal);
   pEff_reco->SetMarkerStyle(8);
   pEff_reco->Write("pEff_reco");

   TEfficiency *pEff_true=new TEfficiency(*h_efficiency_true_num, *h_efficiency_true_den);
   pEff_true->SetStatisticOption(TEfficiency::kFNormal);
   pEff_true->SetMarkerStyle(8);
   pEff_true->Write("pEff_true");

    
   auto gr_test_Purity = new TGraphAsymmErrors(25, slcid, Purity_reco, Err_x, Err_x, Purity_reco_err, Purity_reco_err);
   gr_test_Purity->GetXaxis()->SetTitle("sliceID(reco)");
   gr_test_Purity->GetYaxis()->SetTitle("Purity");
   gr_test_Purity->Write("gr_test_Purity");
   LOG_NORMAL()<<"Calculated Purity and Efficiency and Saved to Histograms"<<std::endl;

   TH1D *h_CalcSignal_Selected = new TH1D("h_CalcSignal_Selected", "h_CalcSignal_Selected", nslices+1, -0.5, 24.5);
   *h_CalcSignal_Selected = (*h_Selected)*(*h_Purity);
   h_CalcSignal_Selected->SetName("h_CalcSignal_Selected");
   LOG_NORMAL()<<"closure test for efficiency "<<std::endl; 
   TH1D *h_CalcSignal_Generated = new TH1D("h_CalcSignal_Generated", "h_CalcSignal_Generated", nslices+1, -0.5, 24.5);
   Double_t temp_siggen[nslices+1];

   for(int i=0; i<nslices+1; i++){
          temp_siggen[i]=selected_tot[i]*Purity_reco[i]/(recoabs_sel[i]/recoabs[i]);
          std::cout<<"i= "<<i<<"  temp_siggen is "<<temp_siggen[i]<<"  true abs generated is "<<recoabs[i]<<std::endl;
          h_CalcSignal_Generated->SetBinContent(i+1, temp_siggen[i]);       

   }
   

  

   h_CalcSignal_Generated->Write();




   outputFile->Close();


   /*
   TH2D *sliceIDmat_abs_den = (TH2D*)mc_intimecosmic_file->Get("sliceIDmat_abs_den");

 
    TMatrix _S_sliceIDmat_den; _S_sliceIDmat_den.Clear(); _S_sliceIDmat_den.ResizeTo(nslices + 2, nslices + 2);
    MigrationMatrix2D migrationmatrix2d;
    migrationmatrix2d.SetOutDir("migration_matrix_2d_sliceIDmat");
    migrationmatrix2d.SetNBins(nslices+1, nslices+1);
    migrationmatrix2d.SetTrueRecoHistogram(sliceIDmat_abs_den);
    _S_sliceIDmat_den = migrationmatrix2d.CalculateMigrationMatrix_uf();
    migrationmatrix2d.PlotMatrix();
    migrationmatrix2d.SetOutputFileName("migration_matrix_2d_sliceIDmat_den.tex");
    migrationmatrix2d.PrintSmearingMatrixLatex();


    TCanvas *c_sliceIDmat_den=new TCanvas("c_sliceIDmat_den", "c_sliceIDmat_den", 1800, 1200);
    sliceIDmat_abs_den->Draw("colz");


    int bins_x = sliceIDmat_abs_den->GetNbinsX()+1;
    int bins_y = sliceIDmat_abs_den->GetNbinsY()+1;
    int _m = bins_x - 1;
    int _n = bins_y - 1;
    TH2D * smearing_matrix_sliceIDmat_den= new TH2D("smearing_matrix_sliceIDmat_den", "", nslices+1, -0.5, 24.5, nslices+1, -0.5, 24.5);
    smearing_matrix_sliceIDmat_den->GetXaxis()->SetTitle("#sliceID (Truth)");
    smearing_matrix_sliceIDmat_den->GetYaxis()->SetTitle("#sliceID (Reco)");
    smearing_matrix_sliceIDmat_den->GetYaxis()->SetTitleOffset(1.4);
    for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_sliceIDmat_den->SetBinContent(j+1, i+1, _S_sliceIDmat_den[i][j]);

      }
    }
    smearing_matrix_sliceIDmat_den->Draw("colz, TEXT");
    prelim_Right->Draw("same");
    c_sliceIDmat_den->SaveAs("Fractional_MigrationMatrix_sliceIDmat_den.png");
    Double_t trueabs_test[nslices+1];
    for(int idx=0; idx<nslices+1; idx++){
        trueabs_test[idx]=0;
        for(int idy=0; idy<nslices+1; idy++){
             trueabs_test[idx] +=recoabs[idy]*smearing_matrix_sliceIDmat_den->GetBinContent(idx+1, idy+1);   
        }
        std::cout<<"test value is "<<trueabs_test[idx] <<"       true value is  "<<intabs[idx]<<std::endl;

    }



    TH2D *sliceIDmat_abs_num = (TH2D*)mc_intimecosmic_file->Get("sliceIDmat_abs_num");

 
    TMatrix _S_sliceIDmat_num; _S_sliceIDmat_num.Clear(); _S_sliceIDmat_num.ResizeTo(nslices + 2, nslices + 2);
    //MigrationMatrix2D migrationmatrix2d;
    migrationmatrix2d.SetOutDir("migration_matrix_2d_sliceIDmat");
    migrationmatrix2d.SetNBins(nslices+1, nslices+1);
    migrationmatrix2d.SetTrueRecoHistogram(sliceIDmat_abs_num);
    _S_sliceIDmat_num = migrationmatrix2d.CalculateMigrationMatrix_uf();
    migrationmatrix2d.PlotMatrix();
    migrationmatrix2d.SetOutputFileName("migration_matrix_2d_sliceIDmat_num.tex");
    migrationmatrix2d.PrintSmearingMatrixLatex();


    TCanvas *c_sliceIDmat_num=new TCanvas("c_sliceIDmat_num", "c_sliceIDmat_num", 1800, 1200);
    sliceIDmat_abs_num->Draw("colz");


    bins_x = sliceIDmat_abs_num->GetNbinsX()+1;
    bins_y = sliceIDmat_abs_num->GetNbinsY()+1;
    _m = bins_x - 1;
    _n = bins_y - 1;
    TH2D * smearing_matrix_sliceIDmat_num= new TH2D("smearing_matrix_sliceIDmat_num", "", nslices+1, -0.5, 24.5, nslices+1, -0.5, 24.5);
    smearing_matrix_sliceIDmat_num->GetXaxis()->SetTitle("#sliceID (Truth)");
    smearing_matrix_sliceIDmat_num->GetYaxis()->SetTitle("#sliceID (Reco)");
    smearing_matrix_sliceIDmat_num->GetYaxis()->SetTitleOffset(1.4);
    for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_sliceIDmat_num->SetBinContent(j+1, i+1, _S_sliceIDmat_num[i][j]);

      }
    }
    smearing_matrix_sliceIDmat_num->Draw("colz,TEXT");
    prelim_Right->Draw("same");
    c_sliceIDmat_num->SaveAs("Fractional_MigrationMatrix_sliceIDmat_num.png");
  */ 




  //=======================================================================================
  LOG_NORMAL() << "Preparing to load Event Histo from BNBCosmic file." << std::endl;
  /*
  UBXSecEventHisto1D * _event_histo_1d_mc = 0;
  mc_bnbcosmic_file->GetObject("UBXSecEventHisto1D", _event_histo_1d_mc);
  LOG_NORMAL() << "Get the object of UBXSecEventHisto1D"<<std::endl;

  Bool_t Object_exists;
  Object_exists=mc_bnbcosmic_file->GetListOfKeys()->Contains("UBXSecEventHisto1D");
  LOG_NORMAL()<<"Is this object exist ???? "<<Object_exists<<std::endl;
  mc_bnbcosmic_file->GetObject("UBXSecEventHisto1D", _event_histo_1d_mc);


  LOG_NORMAL() << "Event Histo correclty loaded from BNBCosmic file." << std::endl;
  */



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
