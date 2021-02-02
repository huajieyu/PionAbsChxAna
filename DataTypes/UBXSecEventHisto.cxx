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



        
      htotinc_reco_beamz = new TH1D("htotinc_reco_beamz", "htotinc_reco_beamz", 300, 0.0, 300);
      htotinc_reco_beamwire = new TH1D("htotinc_reco_beamwire", "htotinc_reco_beamwire", 601, -0.5, 600.5);

      hsel_reco_beamz = new TH1D("hsel_reco_beamz", "hsel_reco_beamz", 300, 0.0, 300.0);
      hsel_reco_beamwire = new TH1D("hsel_reco_beamwire", "hsel_reco_beamwire", 601, -0.5, 600.5);

      h_gen_recovstrue_beamz = new TH2D("h_gen_recovstrue_beamz","h_gen_recovstrue_beamz", 300, 0.0, 300.0, 300, 0.0, 300.0);
      h_gen_recovstrue_beamwire = new TH2D("h_gen_recovstrue_beamwire","h_gen_recovstrue_beamwire", 601, -0.5, 600.5, 601, -0.5, 600.5);

      htotinc_reco_beamz_reso = new TH1D("htotinc_reco_beamz_reso", "htotinc_reco_beamz_reso", 100, -20., 100.); 
      htotinc_reco_beamz_reso->GetXaxis()->SetTitle("End Z resolution");      
      htotinc_reco_beamz_reso->GetYaxis()->SetTitle("incident particles");      

      htotinc_reco_beamwire_reso = new TH1D("htotinc_reco_beamwire_reso", "htotinc_reco_beamwire_reso", 151, -50.5, 100.5); 
      htotinc_reco_beamwire_reso->GetXaxis()->SetTitle("Wire # resolution");
      htotinc_reco_beamz_reso->GetXaxis()->SetTitle("incident particles");      


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
