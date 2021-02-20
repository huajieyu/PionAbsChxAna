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
    TH1D*  htotinc_reco_beamz_reso;
    TH1D*  htotinc_reco_beamwire_reso;
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

