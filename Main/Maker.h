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
#include "TSpline.h"
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

    void LoadHist();
    TSpline3* MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const;
    double InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const;
    bool IsInsideBoundaries(TVector3 const& point) const;
    bool IsTooFarFromBoundaries(TVector3 const& point) const;
    TVector3 PretendAtBoundary(TVector3 const& point) const;

    TH3F *hDz_sim_pos_orig_new;
    TH3F *hDz_sim_neg_orig_new;

    TSpline3 *spline_dx_fwd_neg[31][37];
    TSpline3 *spline_dy_fwd_neg[19][37];
    TSpline3 *spline_dz_fwd_neg[19][31];
   
    TSpline3 *spline_dx_bkwd_neg[31][37];
    TSpline3 *spline_dy_bkwd_neg[19][37];
    TSpline3 *spline_dz_bkwd_neg[19][31];
   
    TSpline3 *spline_dEx_neg[31][37];
    TSpline3 *spline_dEy_neg[19][37];
    TSpline3 *spline_dEz_neg[19][31];
   
    TSpline3 *spline_dx_fwd_pos[31][37];
    TSpline3 *spline_dy_fwd_pos[19][37];
    TSpline3 *spline_dz_fwd_pos[19][31];
   
    TSpline3 *spline_dx_bkwd_pos[31][37];
    TSpline3 *spline_dy_bkwd_pos[19][37];
    TSpline3 *spline_dz_bkwd_pos[19][31];
   
    TSpline3 *spline_dEx_pos[31][37];
    TSpline3 *spline_dEy_pos[19][37];
    TSpline3 *spline_dEz_pos[19][31];
    

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
    double Sce_Corrected_endZ_2nd(double a, double b, double c);

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
    
    const static int thinslicewidth = 10; //cm      
    const static int nwires_in_slice = 20;
    const static int nslices = 480/nwires_in_slice;
    //const static int nslices = 24;
    double intabs_array_den[nslices+3][nslices+3];  
    double intchx_array_den[nslices+3][nslices+3];

    double intabs_array_num[nslices+3][nslices+3];  
    double intchx_array_num[nslices+3][nslices+3];

    Int_t nbinse=12; 
    Int_t nbinsthickness = 100;
    Int_t nbinspmom=10;
    Int_t nbinspcostheta=12;
    Int_t nbinspphi=12;



    double NA=6.02214076e23;
    double MAr=39.95; //gmol
    double Density = 1.39; // g/cm^3


    TH1D *dslcID = new TH1D("dslcID","reco slice ID - true slice ID",20,-10,10);

    TH2D *sliceIDmat_abs_den = new TH2D("sliceIDmat_abs_den", "sliceIDmat_abs_den", 25, -0.5, 24.5, 25, -0.5, 24.5);
    TH2D *sliceIDmat_abs_num = new TH2D("sliceIDmat_abs_num", "sliceIDmat_abs_num", 25, -0.5, 24.5, 25, -0.5, 24.5);
    
    TH2D *sliceIDmat_chx_den = new TH2D("sliceIDmat_chx_den", "sliceIDmat_chx_den", 25, -0.5, 24.5, 25, -0.5, 24.5);
    TH2D *sliceIDmat_chx_num = new TH2D("sliceIDmat_chx_num", "sliceIDmat_chx_num", 25, -0.5, 24.5, 25, -0.5, 24.5);

    TH2D *sliceIDmat_piinelas = new TH2D("sliceIDmat_piinelas", "sliceIDmat_piinelas", 27, -1.5, 25.5, 27, -1.5, 25.5);   
    TH2D *sliceIDmat_pidecay = new TH2D("sliceIDmat_pidecay", "sliceIDmat_pidecay", 27, -1.5, 25.5, 27, -1.5, 25.5);   
    TH2D *sliceIDmat_mudecay = new TH2D("sliceIDmat_mudecay", "sliceIDmat_mudecay", 27, -1.5, 25.5, 27, -1.5, 25.5);   
    TH2D *sliceIDmat_upstream = new TH2D("sliceIDmat_upstream", "sliceIDmat_upstream", 27, -1.5, 25.5, 27, -1.5, 25.5);   
    
    double interaction[nslices+3];
    double signal[nslices+3];

    double incident[nslices+3];

    double incident_pion[nslices+3];
    double incident_pion_decay[nslices+3];
    double incident_pion_elastic[nslices+3];
    double incident_muon[nslices+3];
    double incident_upstream[nslices+3];
    double incident_upstream_pion[nslices+3];


    double true_incident_pandora_identified[nslices+3]; 

    double true_abs_pandora_identified[nslices+3];
    double true_abs_beam_qualified[nslices+3];

    double true_incident[nslices+3];
    double true_incident_pion_decay[nslices+3];
    double true_incident_pion_decay_reco_nnn[nslices+3];
    double true_incident_pion_decay_broken[nslices+3];
    double true_incident_pion_decay_normal[nslices+3];

    double true_incident_upstream[nslices+3];
    double true_incident_upstream_pion[nslices+3];
    double true_incident_pion_elastic[nslices+3];
    double true_incident_pion[nslices+3];
    double true_incident_muon[nslices+3];
    double true_incident_muon_decay_reco_nnn[nslices+3];
    double true_incident_muon_decay_broken[nslices+3];
    double true_incident_muon_decay_normal[nslices+3];

    double true_interaction[nslices+3];

    double true_abs[nslices+3];
    double true_abs_test[nslices+3];
    double reco_abs[nslices+3];
    double true_abs_neg = 0;
    double reco_abs_neg = 0;

    double true_chx[nslices+3];
    double reco_chx[nslices+3];

    double true_chx_neg = 0;
    double reco_chx_neg = 0;

    double true_rea_neg = 0;
    double reco_rea_neg = 0;

    double true_piinelastic_neg = 0;
    double reco_piinelastic_neg = 0;
    double true_pidecay_neg = 0;
    double reco_pidecay_neg = 0;
    double true_piupstream_neg = 0;
    double reco_piupstream_neg = 0;
    double true_muon_neg = 0;
    double reco_muon_neg = 0;
    double true_upstream_neg = 0;
    double reco_upstream_neg = 0;




    double true_abs_sel[nslices+3];
    double reco_abs_sel[nslices+3];
  
    double true_chx_sel[nslices+3];
    double reco_chx_sel[nslices+3];
 

    double selected_tot[nslices+3];
    double selected_mubkg[nslices+3];
    double selected_pibkg[nslices+3];
    double selected_pibkg_elastic[nslices+3];

    double slicebins[nslices+3];


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

    double cutAPA3_Z = 215.;
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
    //double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
    //double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;
    double xlow = -5., xhigh = 14, ylow= -8., yhigh = 10;
    double zlow = 27, zhigh = 37, coslow = 0.93;
    //For Data from Owen Goodwin
    //double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
    //double data_yhigh= 10.;
    //double data_zlow=30., data_zhigh=35., data_coslow=.93;
 
    double data_xlow = -5., data_xhigh = 14., data_ylow= -8., data_yhigh=10.;
    double data_zlow=27, data_zhigh=37, data_coslow = 0.93;
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

    void manual_beamPos_mc_vector(double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ, std::vector<double> *temp_vmc);
    
    void manual_beamPos_data_vector (int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks, std::vector<double> *temp_vdata); 


    std::string beam_particle_Identification(std::string &reco_beam_true_byHits_process, Bool_t &reco_beam_true_byHits_matched, Int_t &reco_beam_true_byHits_origin, Int_t &reco_beam_true_byHits_PDG);

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


   double getEta_broken(vector<vector<double>> canddQdx, vector<vector<double>> trkRR, vector<double> trklen,  int muind, vector<double> beamdQdx, vector<double> beamRR, double beamlen);
 

     //=========================================================================
    
  };
}

#endif
/** @} */ // end of doxygen group 

