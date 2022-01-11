/**
 * \file UBXSecEventHisto1D.h
 *
 * \ingroup DataTypes
 * 
 * \brief Class def header for a class UBXSecEventHisto1D
 *
 * @author deltutto
 */

/** \addtogroup DataTypes

    @{*/
#ifndef __DATATYPES_UBXSECEVENTHISTO1D_H__
#define __DATATYPES_UBXSECEVENTHISTO1D_H__

#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <TChain.h>
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TFile.h>
#include "TH2Poly.h"

#include "UBTH2Poly.h"
#include "BootstrapTH1D.h"
#include "BootstrapTH2D.h"
#include "BootstrapTH2DPoly.h"

namespace DataTypes {

  /**
     \class UBXSecEventHisto1D
     User defined class UBXSecEventHisto1D ... these comments are used to generate
     doxygen documentation!
  */
  class UBXSecEventHisto1D{
    
  public:
    
    /// Default constructor
    UBXSecEventHisto1D() {}
    
    /// Default destructor
    ~UBXSecEventHisto1D(){}

    ///
    void InitializeBootstraps();


    Int_t nbins_beamE = 6;
    Double_t bins_beamE[7] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    
    double bin_size_int = 50.;
    double bin_size_inc = 50.; //2 
    double eStart = 1000.; //1200
    double eEnd = 0.; 
    int nBin_int = (eStart - eEnd) / bin_size_int;
    int nBin_inc = (eStart - eEnd) / bin_size_inc;

    TH1D *h_orimom_proton;
    TH1D *h_orimom_neutron;
    TH1D *h_orimom_pionpm;
    TH1D *h_orimom_pion0;
    TH1D *h_orimom_kaon;
    TH1D *h_orimom_electron;
    TH1D *h_orimom_muon;
    TH1D *h_orimom_photon;
    TH1D *h_orimom_other;


    TH1D *h_chi2_phypo_proton;
    TH1D *h_chi2_phypo_pionpm;
    TH1D *h_chi2_phypo_muon;
    TH1D *h_chi2_phypo_electron;
    TH1D *h_chi2_phypo_kaon;
    TH1D *h_chi2_phypo_other;
    TH1D *h_chi2_phypo_mc;

    TH1D *htmdqdx_chi2_phypo_proton;
    TH1D *htmdqdx_chi2_phypo_pionpm;
    TH1D *htmdqdx_chi2_phypo_muon;
    TH1D *htmdqdx_chi2_phypo_electron;
    TH1D *htmdqdx_chi2_phypo_kaon;
    TH1D *htmdqdx_chi2_phypo_other;
    TH1D *htmdqdx_chi2_phypo_mc;


    TH1D *h_nhits_proton;
    TH1D *h_nhits_pionpm;
    TH1D *h_nhits_muon;
    TH1D *h_nhits_electron;
    TH1D *h_nhits_kaon;
    TH1D *h_nhits_other;
    TH1D *h_reconhits_proton;
    TH1D *h_tmdqdxnhits_proton;
 
    TH1D *h_recomom_selected_proton;
    TH1D *h_recomom_selected_pionpm;
    TH1D *h_recomom_selected_muon;
    TH1D *h_recomom_selected_electron;
    TH1D *h_recomom_selected_kaon;
    TH1D *h_recomom_selected_other;

    TH1D *h_recopimom_selected_proton;
    TH1D *h_recopimom_selected_pionpm;
    TH1D *h_recopimom_selected_muon;
    TH1D *h_recopimom_selected_electron;
    TH1D *h_recopimom_selected_kaon;
    TH1D *h_recopimom_selected_other;




    TH1D *h_mom_gentruep;
    TH1D *h_thetax_gentruep;
    TH1D *h_mom_recotruep;
    TH1D *h_mom_recotruep_test;
    TH1D *h_thetax_recotruep_test;
    TH1D *h_mom_nonrecop;

    TH1D *h_mom_selectedtruep;
    TH1D *h_mom_trkscoretruep;
    TH1D *h_mom_tmdqdxtruep;
    TH1D *h_mom_chi2truep;

    TH1D *h_mom_selectedtruepipm;
    TH1D *h_mom_trkscoretruepipm;
    TH1D *h_mom_tmdqdxtruepipm;
    TH1D *h_mom_chi2truepipm;

 
    TH1D *h_tmdqdx_new;
    TH1D *h_tmdqdx_proton_new;
    TH1D *h_tmdqdx_pionpm_new;
    TH1D *h_tmdqdx_photon_new;
    TH1D *h_tmdqdx_other_new;

    TH1D *h_tmdqdx;
    TH1D *h_tmdqdx_proton;
    TH1D *h_tmdqdx_pionpm;
    TH1D *h_tmdqdx_photon;
    TH1D *h_tmdqdx_other;

    TH2D *h_tmdqdxvsrange_proton;
    TH2D *h_tmdqdxvsrange_pionpm;
    TH2D *h_tmdqdxvsrange_photon;
    TH2D *h_tmdqdxvsrange_other;

    TH1D  *h_trackscore_all;
    TH1D  *h_trackscore_proton;
    TH1D  *h_trackscore_electron;
    TH1D  *h_trackscore_pionpm;
    TH1D  *h_trackscore_photon;
    TH1D  *h_trackscore_muon;
    TH1D  *h_trackscore_pion0;
    TH1D  *h_trackscore_other;

    TH1D  *h_trackscore_shwlike_photon;
    TH1D  *h_trackscore_shwlike_electron;
    TH1D  *h_trackscore_shwlike_nonphoton;
    TH1D  *h_trackscore_shwlike_all;

    TH1D  *h_nhits_shwlike_photon;
    TH1D  *h_nhits_shwlike_electron;
    TH1D  *h_nhits_shwlike_nonphoton;
    TH1D  *h_nhits_shwlike_all;

    TH1D  *h_shweng_photon;
    TH1D  *h_shweng_electron;
    TH1D  *h_shweng_nonphoton;
    TH1D  *h_shweng_all;

    TH1D  *h_shwang_photon;
    TH1D  *h_shwang_electron;
    TH1D  *h_shwang_nonphoton;
    TH1D  *h_shwang_all;

    TH1D  *h_shwdis_photon;
    TH1D  *h_shwdis_electron;
    TH1D  *h_shwdis_nonphoton;
    TH1D  *h_shwdis_all;




    TH2D  *h_nonrecop_momvsnhits;

    TH1D  *h_lownhitsp_thetax;
    TH1D  *h_lownhitsp_thetay;
    TH1D  *h_lownhitsp_thetaz;

    TH1D  *h_lownhitsp_startX;
    TH1D  *h_lownhitsp_startY;;
    TH1D  *h_lownhitsp_startZ;

    TH1D  *h_lownhitsp_endX;
    TH1D  *h_lownhitsp_endY;;
    TH1D  *h_lownhitsp_endZ;

    TH1D  *h_chi2cutp_thetax;
    TH1D  *h_chi2cutp_thetay;
    TH1D  *h_chi2cutp_thetaz;



    TH1D  *h_mom_gentruepion0;
    TH1D  *h_mom_recotruepion0;
    TH1D  *h_mom_recotruenonpi0;

    TH1D  *h_mom_gentruepionpm;
    TH1D  *h_mom_recotruepionpm;
    TH1D  *h_mom_gentrueother;
    TH1D  *h_mom_recotrueother;


    TH1D  *h_Tk_selectedp;

    TH1D  *h_sel_Emissing;
    TH1D  *h_sel_Pmissing;
    TH1D  *h_sel_Ptmissing;
    TH1D  *h_sel_Plongit;   

    TH2D  *h_sig_Ptmissing_vs_pcand;
 
    TH1D  *h_sig_Emissing;
    TH1D  *h_sig_Pmissing;
    TH1D  *h_sig_Ptmissing;
    TH1D  *h_sig_Plongit;

    TH2D  *h_sig_pvsnmult;
    TH2D  *h_PiAbs_sig_true_totalKE_protonvsneutron;

    TH1D  *h_PiAbs_sig_proton_mom;
    TH1D  *h_PiAbs_sig_proton_costheta;


    TH1D  *h_PiAbs_sig_neutron_mom;
    TH1D  *h_PiAbs_sig_neutron_costheta;
  

    TH1D  *h_PiAbs_sig_true_Ptmissing;
    TH1D  *h_PiAbs_sig_true_Ptmissing_onlyproton;

    TH1D  *h_chxbac_Emissing;
    TH1D  *h_chxbac_Pmissing;
    TH1D  *h_chxbac_Ptmissing;
    TH1D  *h_chxbac_Plongit;

    TH1D  *h_reabac_Emissing;
    TH1D  *h_reabac_Pmissing;
    TH1D  *h_reabac_Ptmissing;
    TH1D  *h_reabac_Plongit;
 
    TH1D  *h_other_Emissing;
    TH1D  *h_other_Pmissing;
    TH1D  *h_other_Ptmissing;
    TH1D  *h_other_Plongit;

    TH1D  *h_PiAbs_sel_pmom;
    TH1D  *h_PiAbs_sig_pmom;
    TH1D  *h_PiAbs_chxbac_pmom;
    TH1D  *h_PiAbs_reabac_pmom;
    TH1D  *h_PiAbs_other_pmom;

    TH1D  *h_PiAbs_sel_energetic_pmom;
    TH1D  *h_PiAbs_sig_energetic_pmom;
    TH1D  *h_PiAbs_chxbac_energetic_pmom;
    TH1D  *h_PiAbs_reabac_energetic_pmom;
    TH1D  *h_PiAbs_other_energetic_pmom;
     
    TH1D  *h_PiAbs_sel_pcostheta;
    TH1D  *h_PiAbs_sig_pcostheta;
    TH1D  *h_PiAbs_chxbac_pcostheta;
    TH1D  *h_PiAbs_reabac_pcostheta;
    TH1D  *h_PiAbs_other_pcostheta;

    TH1D  *h_PiAbs_sel_ptheta;
    TH1D  *h_PiAbs_sig_ptheta;
    TH1D  *h_PiAbs_chxbac_ptheta;
    TH1D  *h_PiAbs_reabac_ptheta;
    TH1D  *h_PiAbs_other_ptheta;


    TH1D  *h_PiAbs_sel_pphi;
    TH1D  *h_PiAbs_sig_pphi;
    TH1D  *h_PiAbs_chxbac_pphi;
    TH1D  *h_PiAbs_reabac_pphi;
    TH1D  *h_PiAbs_other_pphi;

    TH1D  *h_PiAbs_sel_pcosthetax;


    TH1D  *h_PiAbs_sel_pmult;
    TH1D  *h_PiAbs_sig_pmult;
    TH1D  *h_PiAbs_chxbac_pmult;
    TH1D  *h_PiAbs_reabac_pmult;

    TH1D  *h_PiAbs_sig_energeticproton_reco_mom;
    TH1D  *h_PiAbs_sig_energeticproton_reco_costheta;
    TH1D  *h_PiAbs_sig_energeticproton_reco_phi;

    TH2D  *h_PiAbs_sig_energeticproton_truevsreco_mom;
    TH2D  *h_PiAbs_sig_energeticproton_truevsreco_costheta;
    TH2D  *h_PiAbs_sig_energeticproton_truevsreco_phi;




    TH1D  *h_PiAbs_gen_sig_energeticproton_mom;
    TH1D  *h_PiAbs_gen_sig_energeticproton_costheta;
    TH1D  *h_PiAbs_gen_sig_energeticproton_phi;


   







 
    TH2D  *h_PiAbs_sig_truevsreco_pmult;

 
    TH1D  *h_PiAbs_sig_nmult;   

    TH1D  *h_reco_pionpm_rea;
    TH1D  *h_reco_photon_rea;
    TH1D  *h_reco_photon_chx;

    TH1D  *h_daughter_beam_dist;

    TH1D  *h_reco_grand_daughter_beam_dist;
    TH1D  *h_reco_daughter_beam_dist;
    TH1D  *h_reco_grand_daughter_angle;
    TH1D  *h_reco_daughter_angle;

    TH1D  *h_reco_daughter_deltax;
    TH1D  *h_reco_daughter_deltay;
    TH1D  *h_reco_daughter_deltaz;


    TH2D  *h_reco_daughter_distvsangle;

    TH1D  *h_gdfromproton;
    TH1D  *h_gdfromother;

    TH2D  *h_selproton_momreso;
    TH1D  *h_selproton_momreso_new;
    TH1D  *h_selproton_momPxreso_new;
    TH1D  *h_selproton_momPyreso_new;
    TH1D  *h_selproton_momPzreso_new;
    TH1D  *h_selproton_costhetareso_new;

    TH2D  *h_selproton_costhetareso;
 
    TH1D  *h_ngamma_frompi0;
    TH1D  *h_pgamma_frompi0;
    //TH1D  *selected_signal_percut;
    //TH1D  *selected_percut;
    TH1D  *h_true_sig_trkmult;
    TH1D  *h_true_chxbac_trkmult;
    TH1D  *h_true_reabac_trkmult;
    TH2D  *h_selproton_dedxRR;

    TH1D  *h_shwcut_mult;
    TH1D  *h_tmdqdxcut_mult;
    TH1D  *h_chi2cut_mult;

    TH1D  *hsig_shwcut_mult;
    TH1D  *hsig_tmdqdxcut_mult;
    TH1D  *hsig_chi2cut_mult;

    TH1D  *hchxbac_shwcut_mult;
    TH1D  *hchxbac_tmdqdxcut_mult;
    TH1D  *hchxbac_chi2cut_mult;

    TH1D  *hreabac_shwcut_mult;
    TH1D  *hreabac_tmdqdxcut_mult;
    TH1D  *hreabac_chi2cut_mult;

    TH1D  *hother_shwcut_mult;
    TH1D  *hother_tmdqdxcut_mult;
    TH1D  *hother_chi2cut_mult;


    TH1D  *hsig_nrecogamma;
    TH1D  *hchxbac_nrecogamma;
    TH1D  *hreabac_nrecogamma;
    TH1D  *hother_nrecogamma;



    TH1D  *h_totalp_mom;
    TH1D  *h_reintp_mom;    

    TH1D  *h_Emissing_true;
    TH1D  *h_Emissing_reco;
    TH1D  *h_Pmissing_true;
    TH1D  *h_Pmissing_reco;

    TH1D  *h_trklen_reso;
    
    TH1D  *h_seltrk_ptheta_cosmic;
    TH1D  *h_seltrk_ptheta_beam_daughter;
    TH1D  *h_seltrk_ptheta_beam_granddaughter;
    TH1D  *h_seltrk_pphi_cosmic;
    TH1D  *h_seltrk_pphi_beam_daughter;
    TH1D  *h_seltrk_pphi_beam_granddaughter;
    TH1D  *h_seltrk_distvtx_daughter;
    TH1D  *h_seltrk_distvtx_granddaughter;

    TH1D  *h_seltrk_angle_daughter;
    TH1D  *h_seltrk_angle_granddaughter;


    TH1D  *h_true_beam_endE_den;
    TH1D  *h_true_beam_endE_num;
    TH2D  *h_sel_gdvsd;

    TH1D  *h_michele_trackscore;
    TH1D  *h_gammae_trackscore;
    TH1D  *h_pione_trackscore;
    TH1D  *h_othere_trackscore;

    TH1D  *h_energetic_shower_eng_all;
    TH1D  *h_energetic_shower_eng_abs;
    TH1D  *h_energetic_shower_eng_chx;
    TH1D  *h_energetic_shower_eng_rea;
    TH1D  *h_energetic_shower_eng_other;

    TH1D  *h_energetic_shower_ang_all;
    TH1D  *h_energetic_shower_ang_abs;
    TH1D  *h_energetic_shower_ang_chx;
    TH1D  *h_energetic_shower_ang_rea;
    TH1D  *h_energetic_shower_ang_other;

    TH1D  *h_energetic_shower_dist_all;
    TH1D  *h_energetic_shower_dist_abs;
    TH1D  *h_energetic_shower_dist_chx;
    TH1D  *h_energetic_shower_dist_rea;
    TH1D  *h_energetic_shower_dist_other;

    TH2D  *h_pcostheta1st2nd_ptmissing;

    TH1D* h_trueE_truePion_inc_initE;
    TH1D* h_trueE_truePion_inc_interE;
    TH1D* h_trueE_truePionInel;
    TH1D* h_trueE_truePion_incident;
    TH1D* h_trueE_trueAbs_interacting;
    TH1D* h_trueE_trueChx_interacting;
    TH1D* h_betheMean_pion;
    TH1D* h_xs_trueE_trueInel;
    TH1D* h_xs_trueE_trueAbs;
    TH1D* h_xs_trueE_trueChx;

  protected:


    
  };
}

#endif
/** @} */ // end of doxygen group 

