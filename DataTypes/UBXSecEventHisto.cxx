#ifndef __DATATYPES_UBXSECEVENTHISTO_CXX__
#define __DATATYPES_UBXSECEVENTHISTO_CXX__

#include "UBXSecEventHisto.h"


namespace DataTypes {


    void UBXSecEventHisto::InitializeBootstraps()
    {

      /*double n_bins_double_mucostheta = _n_bins_double_mucostheta;
      double * bins_double_mucostheta = _bins_double_mucostheta;
      double n_bins_double_mumom = _n_bins_double_mumom;
      double * bins_double_mumom = _bins_double_mumom;
      */
      int n_bins_pmom = 50; 
      int n_bins_pcostheta = 50;
      int n_bins_pphi = 50;
      double bins_pmom[51];
      double bins_pcostheta[51];
      double bins_pphi[51];
      for(int i=0; i<51; i++){
          bins_pmom[i]=0.0 + 1.2/50.0*i;
          bins_pcostheta[i]=-1.0 + 2.0/50.0*i;
          bins_pphi[i]=-TMath::Pi() + 2*TMath::Pi()*i/50.0;
      }



      h_eff_pmom_num = new TH1D("h_eff_pmom_num", "h_eff_pmom_num", n_bins_pmom, bins_pmom);
      h_eff_pmom_den = new TH1D("h_eff_pmom_den", "h_eff_pmom_den", n_bins_pmom, bins_pmom);

      h_eff_pcostheta_num = new TH1D("h_eff_pcostheta_num", "h_eff_pcostheta_num", n_bins_pcostheta, bins_pcostheta);
      h_eff_pcostheta_den = new TH1D("h_eff_pcostheta_den", "h_eff_pcostheta_den", n_bins_pcostheta, bins_pcostheta);
     
      h_eff_pphi_num = new TH1D("h_eff_pphi_num", "h_eff_pphi_num", n_bins_pphi, bins_pphi);
      h_eff_pphi_den = new TH1D("h_eff_pphi_den", "h_eff_pphi_den", n_bins_pphi, bins_pphi);





      h_true_pphi_test = new TH1D("h_true_pphi_test", "h_true_pphi_test", n_bins_pphi, -TMath::Pi(), TMath::Pi());

      htotinc_reco_beamwire = new TH1D("htotinc_reco_beamwire", "htotinc_reco_beamwire", 601, -0.5, 600.5);
 

      htotinc_pion_reco_beamz = new TH1D("htotinc_pion_reco_beamz", "htotinc_pion_reco_beamz", 300, 0.0, 300);
      htotinc_pion_reco_beamwire = new TH1D("htotinc_pion_reco_beamwire", "htotinc_pion_reco_beamwire", 601, -0.5, 600.5);
         
      htotinc_muon_reco_beamz = new TH1D("htotinc_muon_reco_beamz", "htotinc_muon_reco_beamz", 300, 0.0, 300);
      htotinc_muon_reco_beamwire = new TH1D("htotinc_muon_reco_beamwire", "htotinc_muon_reco_beamwire", 601, -0.5, 600.5);

      htotinc_true_beamwire = new TH1D("htotinc_true_beamwire", "htotinc_true_beamwire", 601, -0.5, 600.5);


      hsel_reco_beamz = new TH1D("hsel_reco_beamz", "hsel_reco_beamz", 300, 0.0, 300.0);
      hsel_reco_beamwire = new TH1D("hsel_reco_beamwire", "hsel_reco_beamwire", 601, -0.5, 600.5);

      h_gen_recovstrue_beamz = new TH2D("h_gen_recovstrue_beamz","h_gen_recovstrue_beamz", 300, 0.0, 300.0, 300, 0.0, 300.0);
      h_gen_recovstrue_beamwire = new TH2D("h_gen_recovstrue_beamwire","h_gen_recovstrue_beamwire", 601, -0.5, 600.5, 601, -0.5, 600.5);

      htotinc_reco_beamz_reso = new TH1D("htotinc_reco_beamz_reso", "htotinc_reco_beamz_reso", 100, -30., 50.); 
      htotinc_reco_beamz_reso->GetXaxis()->SetTitle("End Z resolution");      
      htotinc_reco_beamz_reso->GetYaxis()->SetTitle("incident particles");      

      htotinc_reco_beamwire_reso = new TH1D("htotinc_reco_beamwire_reso", "htotinc_reco_beamwire_reso", 151, -50.5, 100.5); 
      htotinc_reco_beamwire_reso->GetXaxis()->SetTitle("Wire # resolution");
      htotinc_reco_beamwire_reso->GetYaxis()->SetTitle("incident particles");      

      htotinc_reco_beam_intE_reso = new TH1D("htotinc_reco_beam_intE_reso", "htotinc_reco_beam_intE_reso", 100, -300., 100.); 
      htotinc_reco_beam_intE_reso->GetXaxis()->SetTitle("Interacting Energy resolution");
      htotinc_reco_beam_intE_reso->GetYaxis()->SetTitle("interacting #pi");      

      htotdecay_reco_beamz_reso = new TH1D("htotdecay_reco_beamz_reso", "htotdecay_reco_beamz_reso", 100, -100., 100.); 
      htotdecay_reco_beamz_reso->GetXaxis()->SetTitle("End Z resolution");      
      htotdecay_reco_beamz_reso->GetYaxis()->SetTitle("decayident particles");      

      htotdecay_reco_beamwire_reso = new TH1D("htotdecay_reco_beamwire_reso", "htotdecay_reco_beamwire_reso", 151, -50.5, 100.5); 
      htotdecay_reco_beamwire_reso->GetXaxis()->SetTitle("Wire # resolution");
      htotdecay_reco_beamwire_reso->GetYaxis()->SetTitle("decayident particles");      

      htotdecay_reco_beam_intE_reso = new TH1D("htotdecay_reco_beam_intE_reso", "htotdecay_reco_beam_intE_reso", 100, -300., 100.); 
      htotdecay_reco_beam_intE_reso->GetXaxis()->SetTitle("Interacting Energy resolution");
      htotdecay_reco_beam_intE_reso->GetYaxis()->SetTitle("interacting #pi");      

      htotmuon_reco_beamz_reso = new TH1D("htotmuon_reco_beamz_reso", "htotmuon_reco_beamz_reso", 100, -100., 100.); 
      htotmuon_reco_beamz_reso->GetXaxis()->SetTitle("End Z resolution");      
      htotmuon_reco_beamz_reso->GetYaxis()->SetTitle("muonident particles");      

      htotmuon_reco_beamwire_reso = new TH1D("htotmuon_reco_beamwire_reso", "htotmuon_reco_beamwire_reso", 151, -50.5, 100.5); 
      htotmuon_reco_beamwire_reso->GetXaxis()->SetTitle("Wire # resolution");
      htotmuon_reco_beamwire_reso->GetYaxis()->SetTitle("muonident particles");      

      htotmuon_reco_beam_intE_reso = new TH1D("htotmuon_reco_beam_intE_reso", "htotmuon_reco_beam_intE_reso", 100, -300., 100.); 
      htotmuon_reco_beam_intE_reso->GetXaxis()->SetTitle("Interacting Energy resolution");
      htotmuon_reco_beam_intE_reso->GetYaxis()->SetTitle("interacting #pi");      



      Int_t nbinse=12; 
      Int_t nbinsthickness = 100;

      for(long unsigned int i=0; i<sizeof(slicebins)/sizeof(slicebins[0])-1; i++){
          //hgen_true_beam_trajz[i]=new TH1D(Form("hgen_true_beam_trajz_%d", i), Form("hgen_true_beam_trajz_%d",i), 100, slicebins[i]*2-10, slicebins[i+1]*2+10); 
          //hgen_reco_beam_endz[i]=new TH1D(Form("hgen_reco_beam_endz_%d", i), Form("hgen_reco_beam_endz_%d",i), 100, slicebins[i]*2-10, slicebins[i+1]*2+10); 
          htotinc_true_beam_trajz[i]=new TH1D(Form("htotinc_true_beam_trajz_%d", i), Form("htotinc_true_beam_trajz_%d",i), 100, -30., 250.); 
          htotinc_reco_beam_endz[i]=new TH1D(Form("htotinc_reco_beam_endz_%d", i), Form("htotinc_reco_beam_endz_%d",i), 100, -30., 250.); 

          htotinc_reco_beam_incidentEnergy[i]=new TH1D(Form("htotinc_reco_beam_incidentEnergy_%d", i), Form("htotinc_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 

          htotinc_reco_beam_sliceThickness[i]=new TH1D(Form("htotinc_reco_beam_sliceThickness_%d", i), Form("htotinc_reco_beam_sliceThickness_%d",i), nbinsthickness, 0., 20.); 



 
          hgen_reco_beam_incidentEnergy[i]=new TH1D(Form("hgen_reco_beam_incidentEnergy_%d", i), Form("hgen_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
          hsel_reco_beam_incidentEnergy[i]=new TH1D(Form("hsel_reco_beam_incidentEnergy_%d", i), Form("hsel_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
          hsig_reco_beam_incidentEnergy[i]=new TH1D(Form("hsig_reco_beam_incidentEnergy_%d", i), Form("hsig_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
          hchxbac_reco_beam_incidentEnergy[i]=new TH1D(Form("hchxbac_reco_beam_incidentEnergy_%d", i), Form("hchxbac_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
          hreabac_reco_beam_incidentEnergy[i]=new TH1D(Form("hreabac_reco_beam_incidentEnergy_%d", i), Form("hreabac_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
          hotherbac_reco_beam_incidentEnergy[i]=new TH1D(Form("hotherbac_reco_beam_incidentEnergy_%d", i), Form("hotherbac_reco_beam_incidentEnergy_%d",i), nbinse, 0.0, 1200.); 
               
      }//end of looping over all the slices

      h_muon_beamendzvsp=new TH2D("h_muon_beamendzvsp", "h_muon_beamendzvsp", 100, 0.0, 500, 100, 0, 1.5);
      h_muon_beamendz_true=new TH1D("h_muon_beamendz_true", "h_muon_beamendz_true", 100, 0.0, 500.);
      h_muon_beamendz_reco=new TH1D("h_muon_beamendz_reco", "h_muon_beamendz_reco", 100, 0.0, 500.);
      h_pion_beamendz_true=new TH1D("h_pion_beamendz_true", "h_pion_beamendz_true", 100, 0.0, 500.);
      h_pion_beamendz_reco=new TH1D("h_pion_beamendz_reco", "h_pion_beamendz_reco", 100, 0.0, 500.);
      h_pion_decay_beamendz_true=new TH1D("h_pion_decay_beamendz_true", "h_pion_decay_beamendz_true", 100, 0.0, 500.);
      h_pion_decay_beamendz_reco=new TH1D("h_pion_decay_beamendz_reco", "h_pion_decay_beamendz_reco", 100, 0.0, 500.);
      h_upstream_beamendz_true=new TH1D("h_upstream_beamendz_true", "h_upstream_beamendz_true", 100, 0.0, 500.);
      h_upstream_beamendz_reco=new TH1D("h_upstream_beamendz_reco", "h_upstream_beamendz_reco", 100, 0.0, 500.);





      h_truepion_beam_deltax = new TH1D("h_truepion_beam_deltax", "h_truepion_beam_deltax", 100, -20, 40);
      h_truepion_beam_deltay = new TH1D("h_truepion_beam_deltay", "h_truepion_beam_deltay", 100, -20, 40);
      h_truepion_beam_deltaz = new TH1D("h_truepion_beam_deltaz", "h_truepion_beam_deltaz", 100, 0, 40);
      h_truepion_beam_cos = new TH1D("h_truepion_beam_cos", "h_truepion_beam_cos", 100, -1., 1.);

      h_truemuon_beam_deltax = new TH1D("h_truemuon_beam_deltax", "h_truemuon_beam_deltax", 100, -20, 40);
      h_truemuon_beam_deltay = new TH1D("h_truemuon_beam_deltay", "h_truemuon_beam_deltay", 100, -20, 40);
      h_truemuon_beam_deltaz = new TH1D("h_truemuon_beam_deltaz", "h_truemuon_beam_deltaz", 100, 0, 40);
      h_truemuon_beam_cos = new TH1D("h_truemuon_beam_cos", "h_truemuon_beam_cos", 100, -1., 1.);

      h_trueproton_beam_deltax = new TH1D("h_trueproton_beam_deltax", "h_trueproton_beam_deltax", 100, -20, 40);
      h_trueproton_beam_deltay = new TH1D("h_trueproton_beam_deltay", "h_trueproton_beam_deltay", 100, -20, 40);
      h_trueproton_beam_deltaz = new TH1D("h_trueproton_beam_deltaz", "h_trueproton_beam_deltaz", 100, 0, 40);
      h_trueproton_beam_cos = new TH1D("h_trueproton_beam_cos", "h_trueproton_beam_cos", 100, -1., 1.);

      h_trueelectron_beam_deltax = new TH1D("h_trueelectron_beam_deltax", "h_trueelectron_beam_deltax", 100, -20, 40);
      h_trueelectron_beam_deltay = new TH1D("h_trueelectron_beam_deltay", "h_trueelectron_beam_deltay", 100, -20, 40);
      h_trueelectron_beam_deltaz = new TH1D("h_trueelectron_beam_deltaz", "h_trueelectron_beam_deltaz", 100, 0, 40);
      h_trueelectron_beam_cos = new TH1D("h_trueelectron_beam_cos", "h_trueelectron_beam_cos", 100, -1., 1.);


      h_truecosmic_beam_deltax = new TH1D("h_truecosmic_beam_deltax", "h_truecosmic_beam_deltax", 100, -20, 40);
      h_truecosmic_beam_deltay = new TH1D("h_truecosmic_beam_deltay", "h_truecosmic_beam_deltay", 100, -20, 40);
      h_truecosmic_beam_deltaz = new TH1D("h_truecosmic_beam_deltaz", "h_truecosmic_beam_deltaz", 100, 0, 40);
      h_truecosmic_beam_cos = new TH1D("h_truecosmic_beam_cos", "h_truecosmic_beam_cos", 100, -1., 1.);

      h_truenottrig_beam_deltax = new TH1D("h_truenottrig_beam_deltax", "h_truenottrig_beam_deltax", 100, -20, 40);
      h_truenottrig_beam_deltay = new TH1D("h_truenottrig_beam_deltay", "h_truenottrig_beam_deltay", 100, -20, 40);
      h_truenottrig_beam_deltaz = new TH1D("h_truenottrig_beam_deltaz", "h_truenottrig_beam_deltaz", 100, 0, 40);
      h_truenottrig_beam_cos = new TH1D("h_truenottrig_beam_cos", "h_truenottrig_beam_cos", 100, -1., 1.);

      h_trueother_beam_deltax = new TH1D("h_trueother_beam_deltax", "h_trueother_beam_deltax", 100, -20, 40);
      h_trueother_beam_deltay = new TH1D("h_trueother_beam_deltay", "h_trueother_beam_deltay", 100, -20, 40);
      h_trueother_beam_deltaz = new TH1D("h_trueother_beam_deltaz", "h_trueother_beam_deltaz", 100, 0, 40);
      h_trueother_beam_cos = new TH1D("h_trueother_beam_cos", "h_trueother_beam_cos", 100, -1., 1.);
   
      h_data_beam_deltax = new TH1D("h_data_beam_deltax", "h_data_beam_deltax", 100, -20, 40);
      h_data_beam_deltay = new TH1D("h_data_beam_deltay", "h_data_beam_deltay", 100, -20, 40);
      h_data_beam_deltaz = new TH1D("h_data_beam_deltaz", "h_data_beam_deltaz", 100, 0, 40);
      h_data_beam_cos = new TH1D("h_data_beam_cos", "h_data_beam_cos", 100, -1., 1.);

      h_ispiinelastic_beam_deltax = new TH1D("h_ispiinelastic_beam_deltax", "h_ispiinelastic_beam_deltax", 100, -10., 10.);
      h_ispiinelastic_beam_deltay = new TH1D("h_ispiinelastic_beam_deltay", "h_ispiinelastic_beam_deltay", 100, -10., 10.);
      h_ispiinelastic_beam_deltaz = new TH1D("h_ispiinelastic_beam_deltaz", "h_ispiinelastic_beam_deltaz", 100, 20., 35.);
      h_ispiinelastic_beam_cos = new TH1D("h_ispiinelastic_beam_cos", "h_ispiinelastic_beam_cos", 100, 0.5, 1.);

      h_ispidecay_beam_deltax = new TH1D("h_ispidecay_beam_deltax", "h_ispidecay_beam_deltax", 100, -10., 10.);
      h_ispidecay_beam_deltay = new TH1D("h_ispidecay_beam_deltay", "h_ispidecay_beam_deltay", 100, -10., 10.);
      h_ispidecay_beam_deltaz = new TH1D("h_ispidecay_beam_deltaz", "h_ispidecay_beam_deltaz", 100, 20., 35.);
      h_ispidecay_beam_cos = new TH1D("h_ispidecay_beam_cos", "h_ispidecay_beam_cos", 100, 0.5, 1.);

      h_ismuon_beam_deltax = new TH1D("h_ismuon_beam_deltax", "h_ismuon_beam_deltax", 100, -10., 10.);
      h_ismuon_beam_deltay = new TH1D("h_ismuon_beam_deltay", "h_ismuon_beam_deltay", 100, -10., 10.);
      h_ismuon_beam_deltaz = new TH1D("h_ismuon_beam_deltaz", "h_ismuon_beam_deltaz", 100, 20., 35.);
      h_ismuon_beam_cos = new TH1D("h_ismuon_beam_cos", "h_ismuon_beam_cos", 100, 0.5, 1.);

      h_isupstream_beam_deltax = new TH1D("h_isupstream_beam_deltax", "h_isupstream_beam_deltax", 100, -10., 10.);
      h_isupstream_beam_deltay = new TH1D("h_isupstream_beam_deltay", "h_isupstream_beam_deltay", 100, -10., 10.);
      h_isupstream_beam_deltaz = new TH1D("h_isupstream_beam_deltaz", "h_isupstream_beam_deltaz", 100, 20., 35.);
      h_isupstream_beam_cos = new TH1D("h_isupstream_beam_cos", "h_isupstream_beam_cos", 100, 0.5, 1.);



      h_reco_beam_tmdqdx = new TH1D("h_reco_beam_tmdqdx", "h_reco_beam_tmdqdx", 100, 0.0, 10.0);
      h_reco_beam_pion_tmdqdx = new TH1D("h_reco_beam_pion_tmdqdx", "h_reco_beam_pion_tmdqdx", 100, 0.0, 10.0);
      h_reco_beam_muon_tmdqdx = new TH1D("h_reco_beam_muon_tmdqdx", "h_reco_beam_muon_tmdqdx", 100, 0.0, 10.0);
      h_reco_beam_proton_tmdqdx = new TH1D("h_reco_beam_proton_tmdqdx", "h_reco_beam_proton_tmdqdx", 100, 0.0, 10.0);
      h_reco_beam_gamma_tmdqdx = new TH1D("h_reco_beam_gamma_tmdqdx", "h_reco_beam_gamma_tmdqdx", 100, 0.0, 10.0);
      h_reco_beam_other_tmdqdx = new TH1D("h_reco_beam_other_tmdqdx", "h_reco_beam_other_tmdqdx", 100, 0.0, 10.0);

      h_reco_beam_costheta_pion = new TH1D("h_reco_beam_costheta_pion", "h_reco_beam_costheta_pion", 100, -1., 1.0);
      h_reco_beam_costheta_muon = new TH1D("h_reco_beam_costheta_muon", "h_reco_beam_costheta_muon", 100, -1., 1.0);
      h_reco_beam_costheta_proton = new TH1D("h_reco_beam_costheta_proton", "h_reco_beam_costheta_proton", 100, -1., 1.0);
      h_reco_beam_costheta_gamma = new TH1D("h_reco_beam_costheta_gamma", "h_reco_beam_costheta_gamma", 100, -1., 1.0);
      h_reco_beam_costheta_other = new TH1D("h_reco_beam_costheta_other", "h_reco_beam_costheta_other", 100, -1., 1.0);




      h_reco_bdangle_other = new TH1D("h_reco_bdangle_other", "h_reco_bdangle_other", 50, 0., TMath::Pi());
      h_reco_bdangle_bkmuon = new TH1D("h_reco_bdangle_bkmuon", "h_reco_bdangle_bkmuon", 50, 0., TMath::Pi());
      h_geteta_same = new TH1D("h_geteta_same", "h_geteta_same", 50, -1.0, 1.0);
      h_geteta_diff = new TH1D("h_geteta_diff", "h_geteta_diff", 50, -1.0, 1.0);

      h_reco_true_sliceID_ori = new TH2D("h_reco_true_sliceID_ori", "h_reco_true_sliceID_ori", 28, -1.5, 26.5, 28, -1.5, 26.5);
      h_reco_true_sliceID_ori->GetXaxis()->SetTitle("True sliceID");
      h_reco_true_sliceID_ori->GetYaxis()->SetTitle("sliceID");

      h_reco_true_sliceID_corr = new TH2D("h_reco_true_sliceID_corr", "h_reco_true_sliceID_corr", 28, -1.5, 26.5, 28, -1.5, 26.5);
      h_reco_true_sliceID_corr->GetXaxis()->SetTitle("True sliceID");
      h_reco_true_sliceID_corr->GetYaxis()->SetTitle("sliceID");

      h_reco_true_incident_piinelastic = new TH2D("h_reco_true_incident_piinelastic", "h_reco_true_incident_piinelastic", 28, -1.5, 26.5, 28, -1.5, 26.5);
      h_reco_true_incident_piinelastic->GetXaxis()->SetTitle("True sliceID");
      h_reco_true_incident_piinelastic->GetYaxis()->SetTitle("sliceID");


      
     h_upstream_cosmic = new TH1D("h_upstream_cosmic", "h_upstream_cosmic", 100, -70.0, 230);
     h_upstream_pion = new TH1D("h_upstream_pion", "h_upstream_pion", 100, -70.0, 230);
     h_upstream_muon = new TH1D("h_upstream_muon", "h_upstream_muon", 100, -70.0, 230);
     h_upstream_proton = new TH1D("h_upstream_proton", "h_upstream_proton", 100, -70.0, 230);
     h_upstream_gamma = new TH1D("h_upstream_gamma", "h_upstream_gamma", 100, -70.0, 230);
     h_upstream_other = new TH1D("h_upstream_other", "h_upstream_other", 100, -70.0, 230);



    }




    void UBXSecEventHisto::AddPolyBins() {

      /*_AddPolyBins_(h_eff_muangle_mumom_poly_num);
      _AddPolyBins_(h_eff_muangle_mumom_poly_den);
      _AddPolyBins_(hmap_trktheta_trkmom_poly);
      _AddPolyBins_(h_poly_reco_per_true);

      _AddPolyBins_(bs_genie_multisim_eff_poly_muangle_mumom_num);
      _AddPolyBins_(bs_genie_multisim_eff_poly_muangle_mumom_den);
      _AddPolyBins_(bs_flux_multisim_eff_poly_muangle_mumom_num);
      _AddPolyBins_(bs_flux_multisim_eff_poly_muangle_mumom_den);
      _AddPolyBins_(bs_extra_syst_multisim_eff_poly_muangle_mumom_num);
      _AddPolyBins_(bs_extra_syst_multisim_eff_poly_muangle_mumom_den);
      _AddPolyBins_(bs_mc_stat_multisim_eff_poly_muangle_mumom_num);
      _AddPolyBins_(bs_mc_stat_multisim_eff_poly_muangle_mumom_den);

      _AddPolyBins_(hmap_trktheta_trkmom_poly_genie_multisim_bs);
      _AddPolyBins_(hmap_trktheta_trkmom_poly_flux_multisim_bs);
      _AddPolyBins_(hmap_trktheta_trkmom_poly_extra_syst_multisim_bs);
      _AddPolyBins_(hmap_trktheta_trkmom_poly_mc_stat_multisim_bs);

      UBTH2Poly* h_poly_binnumber = h_eff_muangle_mumom_poly_num->GetCopyWithBinNumbers("bs");
      TCanvas* canvas_binnumber_poly = new TCanvas("canvas_binnumber_poly", "canvas", 800, 700);
      // gStyle->SetPaintTextFormat("4.0f");
      h_poly_binnumber->Draw("text");
      canvas_binnumber_poly->SaveAs("bin_numbers.pdf");
      */
    }





    void UBXSecEventHisto::_AddPolyBins_(std::map<std::string,UBTH2Poly*> h) {
      /*for (auto it : h) {
        _AddPolyBins_(it.second);
      }*/
    }

    void UBXSecEventHisto::_AddPolyBins_(std::vector<UBTH2Poly*> h) {
      /*for (auto h1 : h) {
        _AddPolyBins_(h1);
      }*/
    }

    void UBXSecEventHisto::_AddPolyBins_(std::vector<std::vector<UBTH2Poly*>> h) {
      /*for (int i = 0; i < (int)h.size(); i++) {
        for (int m = 0; m < _n_poly_bins; m++) {
          _AddPolyBins_(h.at(i).at(m));
        }
      }*/
    }

    void UBXSecEventHisto::_AddPolyBins_(std::map<std::string,std::vector<UBTH2Poly*>> h) {
      /*for (auto it1 : h) {
        for (int m = 0; m < _n_poly_bins; m++) {
          _AddPolyBins_(it1.second.at(m));
        }
      }*/
    }


    void UBXSecEventHisto::_AddPolyBins_(std::map<std::string,std::map<std::string,UBTH2Poly*>> h) {
      /*for (auto it1 : h) { 
        for (auto it2 : it1.second) { 
          _AddPolyBins_(it2.second); 
        } 
      }*/
    }


    void UBXSecEventHisto::_AddPolyBins_(BootstrapTH2DPoly * h) {
/*
      std::map<int, std::pair<int, int>> _exclusion_map2 = { {3, std::make_pair(0, 1)},
                                                             {4, std::make_pair(0, 1)},
                                                             {5, std::make_pair(0, 1)}, };

      h->SetNBinsX(_n_bins_double_mucostheta);

      int separator_counter = 0;
      std::vector<int> separators;
      separators.clear();

      for (int x = 0; x < _n_bins_double_mucostheta; x++) {
        for (int y = 0; y < _n_bins_double_mumom; y++) {

          auto it = _exclusion_map1.find(x);
          if (it != _exclusion_map1.end()) {
            if (y == it->second.first) {
              separator_counter++;
              h->AddBin(_bins_double_mucostheta[it->first], _bins_double_mumom[it->second.first], _bins_double_mucostheta[it->first+1], _bins_double_mumom[it->second.second+1]);
              continue;
            } else if (y == it->second.second) {
              continue;
            }
          }

          it = _exclusion_map2.find(x);
          if (it != _exclusion_map2.end()) {
            if (y == it->second.first) {
              separator_counter++;
              h->AddBin(_bins_double_mucostheta[it->first], _bins_double_mumom[it->second.first], _bins_double_mucostheta[it->first+1], _bins_double_mumom[it->second.second+1]);
              continue;
            } else if (y == it->second.second) {
              continue;
            }
          }

          separator_counter++;
          h->AddBin(_bins_double_mucostheta[x], _bins_double_mumom[y], _bins_double_mucostheta[x+1], _bins_double_mumom[y+1]);
        }
        separators.emplace_back(separator_counter);
        separator_counter = 0;
      }

      h->SetSeparators(separators);
*/    }


    void UBXSecEventHisto::_AddPolyBins_(UBTH2Poly * h) {
/*
      std::map<int, std::pair<int, int>> _exclusion_map2 = { {3, std::make_pair(0, 1)},
                                                             {4, std::make_pair(0, 1)},
                                                             {5, std::make_pair(0, 1)}, };

      h->SetNBinsX(_n_bins_double_mucostheta);

      int separator_counter = 0;
      std::vector<int> separators;
      separators.clear();

      for (int x = 0; x < _n_bins_double_mucostheta; x++) {
        for (int y = 0; y < _n_bins_double_mumom; y++) {

          auto it = _exclusion_map1.find(x);
          if (it != _exclusion_map1.end()) {
            if (y == it->second.first) {
              separator_counter++;
              h->AddBin(_bins_double_mucostheta[it->first], _bins_double_mumom[it->second.first], _bins_double_mucostheta[it->first+1], _bins_double_mumom[it->second.second+1]);
              // std::cout << "h->AddBin(" << _bins_double_mucostheta[it->first] << ", " << _bins_double_mumom[it->second.first] << ", " <<  _bins_double_mucostheta[it->first+1] << ", " << _bins_double_mumom[it->second.second+1] << ")" << std::endl;
              continue;
            } else if (y == it->second.second) {
              continue;
            }
          }

          it = _exclusion_map2.find(x);
          if (it != _exclusion_map2.end()) {
            if (y == it->second.first) {
              separator_counter++;
              h->AddBin(_bins_double_mucostheta[it->first], _bins_double_mumom[it->second.first], _bins_double_mucostheta[it->first+1], _bins_double_mumom[it->second.second+1]);
              // std::cout << "h->AddBin(" << _bins_double_mucostheta[it->first] << ", " << _bins_double_mumom[it->second.first] << ", " <<  _bins_double_mucostheta[it->first+1] << ", " << _bins_double_mumom[it->second.second+1] << ")" << std::endl;
              continue;
            } else if (y == it->second.second) {
              continue;
            }
          }

          separator_counter++;
          h->AddBin(_bins_double_mucostheta[x], _bins_double_mumom[y], _bins_double_mucostheta[x+1], _bins_double_mumom[y+1]);
          // std::cout << "h->AddBin(" << _bins_double_mucostheta[x] << ", " << _bins_double_mumom[y] << ", " <<  _bins_double_mucostheta[x+1] << ", " << _bins_double_mumom[y+1] << ")" << std::endl;
        }
        separators.emplace_back(separator_counter);
        separator_counter = 0;
      }

      h->SetSeparators(separators);
*/    }

  }

#endif
