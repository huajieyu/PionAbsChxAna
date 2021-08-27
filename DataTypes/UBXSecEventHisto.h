/**
 * \file UBXSecEventHisto.h
 *
 * \ingroup DataTypes
 * 
 * \brief Class def header for a class UBXSecEventHisto
 *
 * @author deltutto
 */

/** \addtogroup DataTypes

    @{*/
#ifndef __DATATYPES_UBXSECEVENTHISTO_H__
#define __DATATYPES_UBXSECEVENTHISTO_H__

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
     \class UBXSecEventHisto
     User defined class UBXSecEventHisto ... these comments are used to generate
     doxygen documentation!
  */
  class UBXSecEventHisto{
    
  public:
    
    /// Default constructor
    UBXSecEventHisto() {}
    
    /// Default destructor
    ~UBXSecEventHisto(){}

    ///
    void InitializeBootstraps();

    ///
    void AddPolyBins();

    TH1D*  h_eff_pmom_num;
    TH1D*  h_eff_pmom_den;

    TH1D*  h_eff_pcostheta_num;
    TH1D*  h_eff_pcostheta_den;

    TH1D*  h_eff_pphi_num;
    TH1D*  h_eff_pphi_den;
   
    TH1D*  h_true_pphi_test;

    TH1D*  htotinc_reco_beamwire;

    TH1D*  htotinc_pion_reco_beamz;
    TH1D*  htotinc_pion_reco_beamwire;
    TH1D*  htotinc_muon_reco_beamz;
    TH1D*  htotinc_muon_reco_beamwire;

    TH1D*  htotinc_true_beamwire;


    TH1D*  hsel_reco_beamz;
    TH1D*  hsel_reco_beamwire;
    TH2D*  h_gen_recovstrue_beamz; 
    TH2D*  h_gen_recovstrue_beamwire; 



    Double_t slicebins[25];

 
    TH1D*  htotinc_true_beam_trajz[24];
    TH1D*  htotinc_reco_beam_endz[24];

    TH1D*  hgen_reco_beam_incidentEnergy[24];

    TH1D*  htotinc_reco_beam_incidentEnergy[24];
    TH1D*  htotinc_reco_beam_sliceThickness[24];


    TH1D*  hsel_reco_beam_incidentEnergy[24];
    TH1D*  hsig_reco_beam_incidentEnergy[24];
    TH1D*  hchxbac_reco_beam_incidentEnergy[24];
    TH1D*  hreabac_reco_beam_incidentEnergy[24];
    TH1D*  hotherbac_reco_beam_incidentEnergy[24];


    TH2D* h_muon_beamendzvsp;

    TH1D* h_muon_beamendz_true;
    TH1D* h_muon_beamendz_reco;
    TH1D* h_pion_beamendz_true;
    TH1D* h_pion_beamendz_reco;

    TH1D* h_pion_decay_beamendz_true;
    TH1D* h_pion_decay_beamendz_reco;

    TH1D* h_upstream_beamendz_true;
    TH1D* h_upstream_beamendz_reco;

    TH1D* h_beamendz_data;



    TH1D* hcuts_muon_beamendz_true;
    TH1D* hcuts_muon_beamendz_reco;
    TH1D* hcuts_pion_beamendz_true;
    TH1D* hcuts_pion_beamendz_reco;

    TH1D* hcuts_pion_decay_beamendz_true;
    TH1D* hcuts_pion_decay_beamendz_reco;

    TH1D* hcuts_upstream_beamendz_true;
    TH1D* hcuts_upstream_beamendz_reco;



    TH1D* h_truepion_beam_deltax;
    TH1D* h_truepion_beam_deltay;
    TH1D* h_truepion_beam_deltaz;
    TH1D* h_truepion_beam_cos;

    TH1D* h_truemuon_beam_deltax;
    TH1D* h_truemuon_beam_deltay;
    TH1D* h_truemuon_beam_deltaz;
    TH1D* h_truemuon_beam_cos;

    TH1D* h_trueproton_beam_deltax;
    TH1D* h_trueproton_beam_deltay;
    TH1D* h_trueproton_beam_deltaz;
    TH1D* h_trueproton_beam_cos;

    TH1D* h_trueelectron_beam_deltax;
    TH1D* h_trueelectron_beam_deltay;
    TH1D* h_trueelectron_beam_deltaz;
    TH1D* h_trueelectron_beam_cos;


    TH1D* h_truecosmic_beam_deltax;
    TH1D* h_truecosmic_beam_deltay;
    TH1D* h_truecosmic_beam_deltaz;
    TH1D* h_truecosmic_beam_cos;

    TH1D* h_truenottrig_beam_deltax;
    TH1D* h_truenottrig_beam_deltay;
    TH1D* h_truenottrig_beam_deltaz;
    TH1D* h_truenottrig_beam_cos;

    TH1D* h_trueother_beam_deltax;
    TH1D* h_trueother_beam_deltay;
    TH1D* h_trueother_beam_deltaz;
    TH1D* h_trueother_beam_cos;

    TH1D* h_data_beam_deltax;
    TH1D* h_data_beam_deltay;
    TH1D* h_data_beam_deltaz;
    TH1D* h_data_beam_cos;

    TH1D* h_mc_beam_deltax;
    TH1D* h_mc_beam_deltay;
    TH1D* h_mc_beam_deltaz;
    TH1D* h_mc_beam_cos;

    TH1D* h_data_beamtype;
    TH1D* h_mc_beamtype_muon;
    TH1D* h_mc_beamtype_pion;

    TH1D *h_beam_trackscore_collection;
    TH1D *h_beam_trackscore_collection_pion;
    TH1D *h_beam_trackscore_collection_muon;
    TH1D *h_beam_trackscore_collection_proton;
    TH1D *h_beam_trackscore_collection_other;

    TH1D* h_ispiinelastic_beam_deltax;
    TH1D* h_ispiinelastic_beam_deltay;
    TH1D* h_ispiinelastic_beam_deltaz;
    TH1D* h_ispiinelastic_beam_cos;
    TH1D* h_ispidecay_beam_deltax;
    TH1D* h_ispidecay_beam_deltay;
    TH1D* h_ispidecay_beam_deltaz;
    TH1D* h_ispidecay_beam_cos;
    TH1D* h_ismuon_beam_deltax;
    TH1D* h_ismuon_beam_deltay;
    TH1D* h_ismuon_beam_deltaz;
    TH1D* h_ismuon_beam_cos;
    TH1D* h_isupstream_beam_deltax;
    TH1D* h_isupstream_beam_deltay;
    TH1D* h_isupstream_beam_deltaz;
    TH1D* h_isupstream_beam_cos;

    TH1D* h_test_deltax_piinc;
    TH1D* h_test_deltay_piinc;
    TH1D* h_test_deltaz_piinc;
    TH1D* h_test_cos_piinc;

    TH1D* h_test_deltax_pidec;
    TH1D* h_test_deltay_pidec;
    TH1D* h_test_deltaz_pidec;
    TH1D* h_test_cos_pidec;
	    
    TH1D* h_test_deltax_muon;
    TH1D* h_test_deltay_muon;
    TH1D* h_test_deltaz_muon;
    TH1D* h_test_cos_muon;
	    
    TH1D* h_test_deltax_ups;
    TH1D* h_test_deltay_ups;
    TH1D* h_test_deltaz_ups;
    TH1D* h_test_cos_ups;
	
    TH1D* h_test_deltax_data;
    TH1D* h_test_deltay_data;
    TH1D* h_test_deltaz_data;
    TH1D* h_test_cos_data;
	


    TH1D* h_reco_beam_tmdqdx;
    TH1D* h_reco_beam_pion_tmdqdx;
    TH1D* h_reco_beam_muon_tmdqdx;
    TH1D* h_reco_beam_proton_tmdqdx;
    TH1D* h_reco_beam_gamma_tmdqdx;
    TH1D* h_reco_beam_other_tmdqdx;

    //TH1D* h_reco_beam;
    TH1D* h_reco_beam_costheta_pion;
    TH1D* h_reco_beam_costheta_muon;
    TH1D* h_reco_beam_costheta_proton;
    TH1D* h_reco_beam_costheta_gamma;
    TH1D* h_reco_beam_costheta_other;


    TH1D* h_reco_bdangle_bkmuon;
    TH1D* h_reco_bdangle_other;
    TH1D* h_geteta_same;
    TH1D* h_geteta_diff;

    TH2D* h_reco_true_sliceID_ori;
    TH2D* h_reco_true_sliceID_corr;

    TH2D* h_reco_true_incident_piinelastic;

    TH1D* h_upstream_cosmic;
    TH1D* h_upstream_pion;
    TH1D* h_upstream_muon;
    TH1D* h_upstream_proton;
    TH1D* h_upstream_gamma;
    TH1D* h_upstream_other;
    TH1D* h_bendangle_total;
    TH1D* h_bendangle_cosmic;
    TH1D* h_bendangle_pion;
    TH1D* h_bendangle_muon;
    TH1D* h_bendangle_proton;
    TH1D* h_bendangle_gamma;
    TH1D* h_bendangle_other;

    TH1D* h_bendangle_piinelastic;
    TH1D* h_bendangle_pidecay;
    TH1D* h_bendangle_tmuon;

    TH1D* h_beam_trklen_total;
    TH1D* h_beam_trklen_muon;
    TH1D* h_beam_trklen_pion;

    TH1D* h_beamdeltaz_primarypion_afterbqxycut;
    TH1D* h_beamdeltaz_primarymuon_afterbqxycut;
    TH1D* h_beamdeltaz_primaryproton_afterbqxycut;
    TH1D* h_beamdeltaz_primaryelectron_afterbqxycut;
    TH1D* h_beamdeltaz_cosmic_afterbqxycut;
    TH1D* h_beamdeltaz_primarynottrig_afterbqxycut;
    TH1D* h_beamdeltaz_other_afterbqxycut;


    TH1D* h_dminms;
    TH1D* h_muon_dminms;
    TH1D* h_pion_dminms;
    TH1D* h_pion_decay_dminms;
    TH1D* h_upstream_dminms;


    TH1D* h_reco_beam_startx;
    TH1D* h_reco_beam_starty;
    TH1D* h_reco_beam_startz;

    TH1D* h_beam_intKE_true;
    TH1D* h_beam_intKE_true_piinelastic; 
    TH1D* h_beam_intKE_true_pidecay;
    TH1D* h_beam_intKE_true_muon;
    TH1D* h_beam_intKE_true_upstream;	  

    TH1D* h_beam_intKE_true_aftercut;
    TH1D* h_beam_intKE_true_piinelastic_aftercut; 
    TH1D* h_beam_intKE_true_pidecay_aftercut;
    TH1D* h_beam_intKE_true_muon_aftercut;
    TH1D* h_beam_intKE_true_upstream_aftercut; 


    TH1D* h_beam_intKE_reco;
    TH1D* h_beam_intKE_reco_piinelastic; 
    TH1D* h_beam_intKE_reco_pidecay;
    TH1D* h_beam_intKE_reco_muon;
    TH1D* h_beam_intKE_reco_upstream;

    TH1D* h_beam_intKE_reco_aftercut;
    TH1D* h_beam_intKE_reco_piinelastic_aftercut; 
    TH1D* h_beam_intKE_reco_pidecay_aftercut;
    TH1D* h_beam_intKE_reco_muon_aftercut;
    TH1D* h_beam_intKE_reco_upstream_aftercut; 


  protected:

    


    // TFile *_file = 0;

    ///
    void _AddPolyBins_(UBTH2Poly * h);

    ///
    void _AddPolyBins_(BootstrapTH2DPoly * h);

    ///
    void _AddPolyBins_(std::map<std::string,std::vector<UBTH2Poly*>> h);

    ///
    void _AddPolyBins_(std::vector<std::vector<UBTH2Poly*>> h);

    ///
    void _AddPolyBins_(std::map<std::string,std::map<std::string,UBTH2Poly*>> h);

    ///
    void _AddPolyBins_(std::map<std::string,UBTH2Poly*> h);

    ///
    void _AddPolyBins_(std::vector<UBTH2Poly*> h);

    
  };
}

#endif
/** @} */ // end of doxygen group 

