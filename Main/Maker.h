/**
 * \file Maker.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class Maker
 *
 * @author deltutto
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_MAKER_H__
#define __MAIN_MAKER_H__

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
#include <numeric>
#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TLatex.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TH2Poly.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "ubana/DataTypes/UBTH2Poly.h"
#include "ubana/DataTypes/BootstrapTH2DPoly.h"
#include "ubana/DataTypes/UBXSecEventHisto.h"
#include "ubana/DataTypes/UBXSecEventHisto1D.h"

#include "UBXSecEvent.h"
#include "ubana/DataTypes/BootstrapTH1D.h"
#include "ubana/DataTypes/BootstrapTH2D.h"
#include "ubana/Base/PlottingTools.h"

#include "ubana/Base/LoggerFeature.h"

using namespace DataTypes;
using namespace Base;

namespace Main{

  /**
     \class Maker
     User defined class Maker ... these comments are used to generate
     doxygen documentation!
  */
  class Maker : public LoggerFeature {
    
  public:
    
    /// Default constructor
    Maker(std::string name = "Maker") 
    : LoggerFeature(name) {}
    
    /// Default destructor
    ~Maker(){}

    /// Produces an output file with the relevant histograms
    void MakeFile();

    void SetMomThreshCut(double);
    void SetTrkScoreCut(double);
    void SetFVCut(bool); 
 


    /// Sets the name of the UBXSec input file (the file needs to contain a UBXSec/tree TTree)
    void SetInputFile(std::string);

    void SetInputFileSCE(std::string);

    void SetInputFile_add(std::string);

    /// Sets the name of the ouput file
    void SetOutputFile(std::string);

    /// Sets the number of entries to loop over (-1: all entries)
    void SetEntries(int);

    /// Sets the first entry that will be used in the tree loop (default is 0)
    void SetInitialEntry(int);

    /// Sets the start of beam spill (only used if selection is re-run)
    void SetBeamSpillStart(double);

    /// Sets the end of beam spill (only used if selection is re-run)
    void SetBeamSpillEnd(double);

    /// Sets the shift between beam-on and beam-off flash timing (just for plotting flash time)
    void SetFlashShift(double);

    /// Setd the gain calibration (for dQ/dx) (only used if selection is re-run)
    void SetGainCalibration(double);

    /// If set to True, calculated the POT for this file
    void SetCalculatePOT(bool);

    /// Sets if the file is a data one or not
    void SetIsData(bool);

    void SetSignalTypeAbs(bool);

    void SetSignalTypeChx(bool);
 

    /// If true, MEC events are turned off, and MA is scaled up
    void SetMaUpMECOff(bool option) {_maup_mecoff = option;}

    /// Call this to scale the kaon flux by 1.5
    void ReweighKaons(bool option = true, double factor = 1.5) {_reweigh_kaons = option; _kaon_reweigh_factor = factor;}

    ///
    void PrintConfig();

    /// Sets an extra weight to be applied to every MC event
    void SetExtraWeight(double w) {_extra_weight = w;};

    /// If called, will scaled all cosmics by "w"
    void ScaleCosmics(double w) {_scale_cosmics = true; _scale_factor_cosmic = w;};

    /// If True is passed, filles all the universes histogram for Flux
    void FillBootstrapFlux(bool option) {_fill_bootstrap_flux = option;}

    /// If True is passed, filles all the universes histogram for GENIE
    void FillBootstrapGenie(bool option) {_fill_bootstrap_genie = option;}

    void FillBootstrapGeant(bool option) {_fill_bootstrap_geant = option;}

    /// If True is passed, filles all the universes histogram for Extra Syst
    void FillBootstrapExtraSyst(bool option) {_fill_bootstrap_extra_syst = option;}

    /// If True is passed, filles all the universes histogram for MC Stat
    void FillBootstrapMCStat(bool option) {_fill_bootstrap_mc_stat = option;}

    /// Sets the number of universes to create for MC stat Poisson reweighing
    void SetNUniversesMCStat(int n) {_mc_stat_n_events = n;}

    /// Sets the target used for Flux systematics
    void SetTargetFluxSystematic(std::string s) { _target_flux_syst = s; }

    /// Sets the target used for extra systematics
    void SetTargetExtraSystematic(std::string s) { _extra_syst_target_syst = s; }

    //UBXSecEventHisto1D * _event_histo_1d;
    //UBXSecEventHisto   * _event_histo;


  private:

    /// Prints a warning message if running with MA up and MEC off
    void PrintMaUpMECOff();

    /// Prints a warning message if running with scaled kaon flux
    void PrintReweighKaons();

    void DrawProgressBar(double progress, double barWidth);

    void DrawPOT2(double pot, double target = 6.6e20);

    double eff_uncertainty(int _n, int _N);

    void FillBootstrap(double fill_value,
                       double evt_wgt,
                       std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_genie_pm1_bs, 
                       std::string channel_namel, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts_genie);

    void FillBootstrap(double fill_value1,
                       double fill_value2,
                       double evt_wgt,
                       std::map<std::string,std::map<std::string,TH2D*>> hmap_trkmom_genie_pm1_bs, 
                       std::string channel_namel, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts_genie);

    // void FillBootstrap(double fill_value,
    //                std::map<std::string,TH1D*> hmap_trkmom_genie_pm1_bs, 
    //                std::vector<std::string> fname, 
    //                std::vector<double> wgts_genie);

    void FillBootstrap(double fill_value1,
                       double fill_value2,
                       double evt_wgt,
                       std::map<std::string,TH2D*> hmap_trkmom_genie_pm1_bs, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts_genie);

    void FillBootstrap(double fill_value1, // reco value x (costheta)
                       double fill_value2, // reco value y (momentum)
                       int m, // true bin m (costheta)
                       int n, // true bin n (momentum)
                       double evt_wgt,
                       std::map<std::string,std::vector<std::vector<TH2D*>>> bs_reco_per_true, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts);

    void FillBootstrap(double fill_value1, // reco value x (costheta)
                       double fill_value2, // reco value y (momentum)
                       int n, // true bin n (1 number, unrolled)
                       double evt_wgt,
                       std::map<std::string,std::vector<UBTH2Poly*>> bs_poly_reco_per_true, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts);
        
    // void FillBootstrap(double fill_value1, // reco value x (costheta)
    //                    double fill_value2, // reco value y (momentum)
    //                    int n, // true bin n (1 number, unrolled)
    //                    double evt_wgt,
    //                    std::vector<std::vector<UBTH2Poly*>> bs_poly_reco_per_true, 
    //                    std::vector<std::string> fname, 
    //                    std::vector<double> wgts);

    void FillBootstrap(int m, // true bin m (1 number, unrolled)
                                int j, // reco bin i (1 number, unrolled)
                                double evt_wgt,
                                std::map<std::string,std::vector<std::vector<double>>> &bs_poly_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts);

    void FillBootstrap(double fill_value1,
                       double fill_value2,
                       double evt_wgt,
                       std::map<std::string,std::map<std::string,UBTH2Poly*>> hmap, 
                       std::string channel_namel, 
                       std::vector<std::string> fname, 
                       std::vector<double> wgts_genie);


    double Sce_Corrected_endZ(double a, double b, double c);

    void AddPolyBins(UBTH2Poly * h);

    void AddPolyBins(BootstrapTH2DPoly h);

    UBXSecEventHisto1D * _event_histo_1d;
    UBXSecEventHisto   * _event_histo;



    bool _maup_mecoff = false;
    bool _reweigh_kaons = false;
    double _kaon_reweigh_factor = 1.5;

    const bool _breakdownPlots = true;
    const bool _makePlots = false;

    //const bool _fill_bootstrap = true;
    bool _fill_bootstrap_flux = false;
    bool _fill_bootstrap_genie = false;
    bool _fill_bootstrap_geant = false;
    bool _fill_bootstrap_extra_syst = false;

    std::string _target_flux_syst = "";
    std::string _extra_syst_target_syst = "";

    const bool _check_duplicate_events = false;

    double _beamSpillStarts = 3.2;  // us
    double _beamSpillEnds   = 4.8;  // us
    double _flashShift      = 0.;   //4.06; //us
    double _gainCalib       = 198;  // e-/ADC

    std::string filen     = "ubxsec_output.root";


    std::string filen_add = "";
    std::string fileoutn  = "ubxsecana_output.root";

    bool evalPOT          = false;
    double Evttot         = 0.0;

    int maxEntries        = -1;
    int _initial_entry    = 0; ///< Entry in Tree to begin with
    bool isdata           = false;

    bool sel_abs          = false;
 
    bool sel_chx          = false;

    double _extra_weight = 1.; ///Extra weight to be applied to the events

    const double _pe_cut = 50;

    const double targetPOT = 4.95e19;

    double bins_mumom[7] = {0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50};
    double bins_mucostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};

    int n_bins_mumom = 6;
    int n_bins_mucostheta = 9;

    int n_bins_double_mumom = 6; ///< Number of momentum bins for double differential
    double bins_double_mumom[7] = {0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50}; ///< Momentum bins for double differential
    int n_bins_double_mucostheta = 9; ///< Number of costheta bins for double differential
    double bins_double_mucostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00}; ///< costheta bins for double differential

    int _n_poly_bins = 43;
    std::map<int, std::pair<int, int>> _exclusion_map = { {0, std::make_pair(4, 5)},
                                                          {1, std::make_pair(4, 5)}, };

    bool _scale_cosmics = false; ///< If true scales the cosmic background by _scale_factor_cosmic
    double _scale_factor_cosmic = 1.; ///< Factor used to scale the cosmic background (used only if _scale_cosmics is true)


    // These variables are filled in the reco-true TTree in the code
    std::vector<std::string> _wgtsnames_genie_multisim;
    std::vector<double> _wgts_genie_multisim;
    std::vector<std::string> _wgtsnames_extra_syst;
    std::vector<double> _wgts_extra_syst;
    std::vector<std::string> _wgtsnames_flux_multisim;
    std::vector<double> _wgts_flux_multisim;


    std::map<int, int> sliceID_incidentN;
    std::map<int, double> sliceID_incidentthickness;
    std::map<int, int> sliceID_interactN;
    std::map<int, double> sliceID_interactthickness;
    std::map<int, int> sliceID_backgroundN;
    std::map<int, double> sliceID_backgroundthickness;
          
    const static int nwires_in_slice = 20;
    const static int nslices = 480/nwires_in_slice;

    Int_t nbinse=12; 
    Int_t nbinsthickness = 100;

    double NA=6.02214076e23;
    double MAr=35.95; //gmol
    double Density = 1.39; // g/cm^3


    TH1D *dslcID = new TH1D("dslcID","reco slice ID - true slice ID",20,-10,10);
 
    double interaction[nslices];
    double signal[nslices];
    double incident[nslices];

    double true_incident[nslices];
    double true_interaction[nslices];
    double true_abs[nslices];
 
    double interaction_sel[nslices];
    double incident_sel[nslices];

    double selected_tot[nslices];
    double selected_bkg[nslices];


    double interaction_gen[nslices];
    double incident_gen[nslices];



    //TH1D *incE[nslices];
    //TH1D *pitch[nslices];

    double slicebins[nslices+1];


    bool _fill_bootstrap_mc_stat = false; ///< If true, fills bootstrap with poisson weights (mc stat)
    int _mc_stat_n_events = 100; ///< Number of universes uses for poisson weights (mc stat)

    TRandom _random_engine; ///< The engine to generate random numbers
    //=========================================================================
    bool inFV(double x, double y, double z);
    double thetax(double theta, double phi);
    bool FVcuton = false; 

    double Calc_reco_beam_energy(TVector3 momentum, double mass);
    TVector3 Calc_reco_beam_momentum(vector<double>*dedx, TVector3 *dir3);

    //std::string filen = "";
    //std::string filen_add = "";
    //std::string fileoutn = "";
    //int _initial_entry = 0;
    //int maxEntries = -1;
    double momthreshcut = 0.0;
    double trkscorecut=0.0;
    double cut_dEdX_high = 2.8;
    double cut_dEdX_low = 0.6;
    double cut_for_proton = 3.4;
    double libo_low = 0.16;
    double libo_high = 0.16;
    //PDEventHisto1D *_event_histo_1d;

    bool isProton = false;
    double Ecalcmiss(double Esum, double PTmiss, int np); 
    const double NeutronMass = 0.93956542; 
    const double ProtonMass = 0.938272;
    const double PionMass = 0.13957; 

    double cutAPA3_Z = 226.;
    double cut_trackScore = 0.4;
    //daughter Distance cut

    int temp_index_shw= -999;

    double cut_daughter_track_distance = 10.;

    double cut_daughter_shower_distance_low = 5.;
    double cut_daughter_shower_distance_high = 1000.;
    double cut_primary_chi2 = 140.;
    double cut_secondary_chi2 = 50.;

    int cut_nHits_shower_low = 40;
    int cut_nHits_shower_high = 1000;

    int cut_energy_shower_low = 80;
    int cut_energy_shower_high = 1000; 

    double cut_angle_shower_low = 0.1;
    double cut_angle_shower_high = 0.2;

    //
    //For MC from Owen Goodwins studies
    double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
    double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;
    //
    //For Data from Owen Goodwin
    double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
    double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

    bool isBeamType(int i);
    bool data_beam_PID(const std::vector<int> *pidCandidates);

    bool manual_beamPos_mc(double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ);
    
    bool manual_beamPos_data (int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks); 
    bool endAPA3(double reco_beam_endZ); 
    //
    //Tag PrimaryPion without elastic Scattering
    bool has_pi0shower(const std::vector<double> &track_score);
    const std::vector<double> truncatedMean(double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX); 
    const std::vector<double> truncatedMean_libo(double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX); 
    const std::vector<double> truncatedMean_xglu(double tmdqdx_test); 

    double GetTruncatedMean(const vector<double> tmparr, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac);
    bool has_shower_nHits_distance(const std::vector<double> &track_score,
                                    const std::vector<int> &nHits,
                                    const std::vector<double> &distance); 

    bool has_shower_energy_distance(const std::vector<double> &track_score,
                                    const std::vector<double> &energy,
                                    const std::vector<double> &distance); 

     bool has_shower_nHits(const std::vector<double> &track_score,
                                  const std::vector<int> &nHits);

    bool has_shower_Eng(const std::vector<double> &track_score,
                                  const std::vector<double> &shweng);

    bool has_shower_Ang(const std::vector<double> &track_score,
                                  const std::vector<double> &shwang);

    bool has_shower_Dist(const std::vector<double> &track_score,
                                  const std::vector<double> &distance);

    bool secondary_noPion( const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dEdX);

    bool secondary_noPion_libo( const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dEdX);





 
    const std::vector<double> compute_distanceVertex( double beam_endX,
                                 double beam_endY,
                                 double beam_endZ, 
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ);


    const std::vector<double> compute_angleVertex   ( double beam_endX,
                                 double beam_endY,
                                 double beam_endZ,
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ,
                                 const std::vector<double> &d_startTheta,
                                 const std::vector<double> &d_startPhi,
                                 const std::vector<double> &d_endTheta,
                                 const std::vector<double> &d_endPhi);

    const std::vector<double> compute_angleVertex_shower   ( double beam_endX,
                                 double beam_endY,
                                 double beam_endZ,
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ,
                                 const std::vector<double> &d_dirX,
                                 const std::vector<double> &d_dirY,
                                 const std::vector<double> &d_dirZ);

     //=========================================================================
    
  };
}

#endif
/** @} */ // end of doxygen group 

