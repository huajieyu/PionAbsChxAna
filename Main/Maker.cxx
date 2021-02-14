#ifndef __MAIN_MAKER_CXX__
#define __MAIN_MAKER_CXX__

#include "Maker.h"
using namespace Base;

void Main::Maker::SetInputFile(std::string in)
{
  filen = in;
}

void Main::Maker::SetInputFileSCE(std::string in)
{
 filen_add = in;
}

void Main::Maker::SetOutputFile(std::string in)
{
  fileoutn = in;
}

void Main::Maker::SetEntries(int e)
{
  maxEntries = e;
}

void Main::Maker::SetInitialEntry(int e)
{
  _initial_entry = e;
}

void Main::Maker::SetBeamSpillStart(double v)
{
  _beamSpillStarts = v;
}

void Main::Maker::SetBeamSpillEnd(double v)
{
  _beamSpillEnds = v;
}

void Main::Maker::SetFlashShift(double v)
{
  _flashShift = v;
}

void Main::Maker::SetGainCalibration(double v)
{
  _gainCalib = v;
}

void Main::Maker::SetCalculatePOT(bool v)
{
  evalPOT = v;
}


void Main::Maker::SetIsData(bool v)
{
  isdata = v;
}

void Main::Maker::SetSignalTypeAbs(bool v)
{
  sel_abs = v;
}

void Main::Maker::SetSignalTypeChx(bool v)
{
  sel_chx = v;
}
double Main::Maker::Sce_Corrected_endZ_2nd(double true_beam_endX, double true_beam_endY, double true_beam_endZ){
     double true_endz=true_beam_endZ;
     double offset_z = 0;
     TVector3 point(true_beam_endX, true_beam_endY, true_beam_endZ);
        if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);

        if (point.X()>0){
          //true_endz += RecoFwd_Displacement_Z_Pos->GetBinContent(RecoFwd_Displacement_Z_Pos->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
          //offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Pos, point.X(), point.Y(), point.Z(), 3, 1, 2);
          offset_z = InterpolateSplines(hDz_sim_pos_orig_new, point.X(), point.Y(), point.Z(), 3, 1, 2);
        }
        else{
          //true_endz += RecoFwd_Displacement_Z_Neg->GetBinContent(RecoFwd_Displacement_Z_Neg->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
          //offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Neg, point.X(), point.Y(), point.Z(), 3, 1, 1);
          offset_z = InterpolateSplines(hDz_sim_neg_orig_new, point.X(), point.Y(), point.Z(), 3, 1, 1);
        }
        std::cout<<true_beam_endX<<" "<<true_beam_endY<<" "<<true_beam_endZ<<" "<<offset_z<<std::endl;
        true_endz += offset_z;
 
     return true_endz;

}
double Main::Maker::Sce_Corrected_endZ(double true_beam_endX, double true_beam_endY, double true_beam_endZ){

        TFile f2("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
        TH3F *RecoFwd_Displacement_Z_Neg = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Neg");
        TH3F *RecoFwd_Displacement_Z_Pos = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Pos");
        double true_endz = true_beam_endZ;
        if (true_beam_endX>0){
          true_endz += RecoFwd_Displacement_Z_Pos->GetBinContent(RecoFwd_Displacement_Z_Pos->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
        }
        else{
          true_endz += RecoFwd_Displacement_Z_Neg->GetBinContent(RecoFwd_Displacement_Z_Neg->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
        }
        //std::cout<<"Sce_Corrected_endZ  "<<true_beam_endZ<<"   "<<true_endz<<std::endl; 
       
        return true_endz;

}

void Main::Maker::PrintConfig()
{
  /*LOG_INFO() << "--- Main::Maker::PrintConfig" << std::endl;

  LOG_INFO() << "--- _breakdownPlots: " << _breakdownPlots << std::endl;
  LOG_INFO() << "--- _makePlots " << _makePlots << std::endl;
  LOG_INFO() << "--- _fill_bootstrap_flux " << _fill_bootstrap_flux << std::endl;
  LOG_INFO() << "--- _fill_bootstrap_genie " << _fill_bootstrap_genie << std::endl;
  LOG_INFO() << "--- _target_flux_syst " << _target_flux_syst << std::endl;
  LOG_INFO() << "--- _check_duplicate_events " << _check_duplicate_events << std::endl;

  LOG_INFO() << "--- _beamSpillStarts " << _beamSpillStarts << std::endl;
  LOG_INFO() << "--- _beamSpillEnds " << _beamSpillEnds << std::endl;
  LOG_INFO() << "--- _flashShift " << _flashShift << std::endl;
  LOG_INFO() << "--- _gainCalib " << _gainCalib << std::endl;

  LOG_INFO() << "--- filen " << filen << std::endl;
  LOG_INFO() << "--- evalPOT " << evalPOT << std::endl;
  LOG_INFO() << "--- maxEntries " << maxEntries << std::endl;
  LOG_INFO() << "--- isdata " << isdata << std::endl;

  LOG_INFO() << "--- _pe_cut " << _pe_cut << std::endl;

  LOG_INFO() << "--- targetPOT " << targetPOT << std::endl;
  */
}

void Main::Maker::PrintMaUpMECOff()
{
/*  for (int i = 0; i < 10; i++) {
    std::cout << "**************************** RUNNING WITH MA+1SIGMA AND MEC OFF ****************************" << std::endl;
  }
*/
}

void Main::Maker::PrintReweighKaons()
{
/*  for (int i = 0; i < 10; i++) {
    std::cout << "**************************** RUNNING WITH KAON FLUX SCALED BY " << _kaon_reweigh_factor << " ****************************" << std::endl;
  }*/
}








//____________________________________________________________________________________________________
void Main::Maker::DrawProgressBar(double progress, double barWidth) {
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

//____________________________________________________________________________________________________
void Main::Maker::DrawPOT2(double pot, double target)
{
  //std::string str = "Simulated POT:" + std::to_string(pot);
  
  std::stringstream sstm;
  sstm << "Simulated POT: " << pot;
  std::string str = sstm.str();
  
  TLatex* pot_latex = new TLatex(.10, .96, str.c_str());
  pot_latex->SetTextColor(kGray+2);
  pot_latex->SetNDC();
  pot_latex->SetTextSize(1/30.);
  pot_latex->SetTextAlign(10); //left adjusted
  pot_latex->Draw();
  
  
  std::stringstream sstm2;
  sstm2 << "Scaled to POT: " << target;
  str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}



//____________________________________________________________________________________________________
double Main::Maker::eff_uncertainty(int _n, int _N) {

  double n = (double) _n;
  double N = (double) _N;

  double unc = 1/std::sqrt(N) * std::sqrt((n/N)*(1-n/N));

  return unc;

}


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_genie_pm1_bs, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap_trkmom_genie_pm1_bs.find(channel_namel);
  if (iter == hmap_trkmom_genie_pm1_bs.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  } //second element of the map is the channel
  std::map<std::string,TH1D*> this_map = iter->second;
  //Fill nominal weighted by evt_wgt first
  this_map["nominal"]->Fill(fill_value, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value, wgts_genie.at(i) * evt_wgt);

  }


} 

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,TH2D*>> hmap, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap.find(channel_namel);
  if (iter == hmap.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  }
  std::map<std::string,TH2D*> this_map = iter->second;

  this_map["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

  }


} 


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,UBTH2Poly*>> hmap, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap.find(channel_namel);
  if (iter == hmap.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  }
  std::map<std::string,UBTH2Poly*> this_map = iter->second;

  this_map["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

  }


} 

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,TH2D*> hmap_trkmom_genie_pm1_bs, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {


  hmap_trkmom_genie_pm1_bs["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    hmap_trkmom_genie_pm1_bs[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1, // reco value x (costheta)
                                double fill_value2, // reco value y (momentum)
                                int m, // true bin m (costheta)
                                int n, // true bin n (momentum)
                                double evt_wgt,
                                std::map<std::string,std::vector<std::vector<TH2D*>>> bs_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  bs_reco_per_true["nominal"][m][n]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    bs_reco_per_true[fname.at(i)][m][n]->Fill(fill_value1, fill_value2, wgts.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1, // reco value x (costheta)
                                double fill_value2, // reco value y (momentum)
                                int m, // true bin m (1 number, unrolled)
                                double evt_wgt,
                                std::map<std::string,std::vector<UBTH2Poly*>> bs_poly_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  bs_poly_reco_per_true["nominal"][m]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    bs_poly_reco_per_true[fname.at(i)][m]->Fill(fill_value1, fill_value2, wgts.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(int m, // true bin m (1 number, unrolled)
                                int j, // reco bin i (1 number, unrolled)
                                double evt_wgt,
                                std::map<std::string,std::vector<std::vector<double>>> & bs_poly_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  if (j < 0) j = 0; // Negative bins are overflows, and are all added to entry 0 of the vector

  bs_poly_reco_per_true["nominal"][m][j] += evt_wgt;

  for (size_t i = 0; i < fname.size(); i++) {

    bs_poly_reco_per_true[fname.at(i)][m][j] += wgts.at(i) * evt_wgt;

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}




//___________________________________________________________________________________________________
void Main::Maker::AddPolyBins(UBTH2Poly * h) {

  // std::map<int, std::pair<int, int>> _exclusion_map;
  // _exclusion_map[0] = std::make_pair(2, 3);

  /*h->SetNBinsX(n_bins_double_mucostheta);

  for (int y = 0; y < n_bins_double_mumom; y++) {
    for (int x = 0; x < n_bins_double_mucostheta; x++) {

      auto it = _exclusion_map.find(x);
      if (it != _exclusion_map.end()) {
        if (y == it->second.first) {
          h->AddBin(bins_double_mucostheta[it->first], bins_double_mumom[it->second.first], bins_double_mucostheta[it->first+1], bins_double_mumom[it->second.second+1]);
          continue;
        } else if (y == it->second.second) {
          continue;
        }
      }
      h->AddBin(bins_double_mucostheta[x], bins_double_mumom[y], bins_double_mucostheta[x+1], bins_double_mumom[y+1]);
    }
  }
*/
}



//___________________________________________________________________________________________________
void Main::Maker::AddPolyBins(BootstrapTH2DPoly h) {

  // std::map<int, std::pair<int, int>> _exclusion_map;
  // _exclusion_map[0] = std::make_pair(2, 3);

  // h->SetNBinsX(n_bins_double_mucostheta);

 /* for (int y = 0; y < n_bins_double_mumom; y++) {
    for (int x = 0; x < n_bins_double_mucostheta; x++) {

      auto it = _exclusion_map.find(x);
      if (it != _exclusion_map.end()) {
        if (y == it->second.first) {
          h.AddBin(bins_double_mucostheta[it->first], bins_double_mumom[it->second.first], bins_double_mucostheta[it->first+1], bins_double_mumom[it->second.second+1]);
          continue;
        } else if (y == it->second.second) {
          continue;
        }
      }
      h.AddBin(bins_double_mucostheta[x], bins_double_mumom[y], bins_double_mucostheta[x+1], bins_double_mumom[y+1]);
    }
  }
*/
}



//================================================================

bool Main::Maker::data_beam_PID(const std::vector<int> *pidCandidates){
  auto pid_search = std::find(pidCandidates->begin(), pidCandidates->end(), 211);
  return (pid_search !=pidCandidates->end());
}




bool Main::Maker::isBeamType(int reco_beam_type){
  return (reco_beam_type == 13);
};

bool Main::Maker:: manual_beamPos_data (int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks) {

  double deltaX = data_startX - data_BI_X;
  double deltaY = data_startY - data_BI_Y;
  double cos = data_BI_dirX*data_dirX + data_BI_dirY*data_dirY +
               data_BI_dirZ*data_dirZ;

  if(data_BI_nMomenta != 1 || data_BI_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

};


bool Main::Maker::manual_beamPos_mc(double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ) {
  double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;

};
bool Main::Maker::endAPA3(double reco_beam_endZ){
  return(reco_beam_endZ < cutAPA3_Z);

}

bool Main::Maker::has_pi0shower(const std::vector<double> &track_score){
   //calculate the angle between shower start point - vertex direction and
   //the direction of the shower
   //if the angle< cut value, this shower come from the primary vertex

   return true;
}
//====================================================================================


const std::vector<double> Main::Maker::truncatedMean_xglu(double tmdqdx_test){

   std::vector<double> trunc_mean;
   return trunc_mean;

};

//=====================================================================================
//truncated mean of SIGMA = cutting %
const std::vector<double> Main::Maker::truncatedMean(double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX){

   size_t size = 0;
   std::vector<double> trunc_mean;
   std::vector<double> help_vec;
   truncate_high = 1 - truncate_high; 
   int i_low = 0;
   int i_high = 0;

   //sort the dEdX vecotrs in matrix
   for(auto &&vec : vecs_dEdX){
      size = vec.size();
      help_vec.clear();

      //check dEdX vector isn't empty!
      if(vec.empty()){
         trunc_mean.push_back(-9999.);
         continue;
      }

      else{
         //Sort Vector
         sort(vec.begin(), vec.end());
       
         //Discard upper and lower part of signal
         //rint rounds to integer
         i_low = rint ( size*truncate_low);
         i_high = rint( size*truncate_high);
         
         
         for(int i = i_low; i <= i_high; i++){
               help_vec.push_back(vec[i]);
         };

         //Mean of help vector

         trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());

      }


   };

   return trunc_mean;

};
//===================================================================
double Main::Maker::GetTruncatedMean(const vector<double> tmparr, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac)
{
  //for proton Bragg peak use 0.4-0.95. Seen by CDF of startE using signal proton samples in drawTracking

  if(nsample1>=tmparr.size()){
    return -999;
  }
  vector<double> array;
  for(unsigned int ii=nsample0; ii<=nsample1; ii++){
    array.push_back(tmparr[ii]);
  }
  std::sort(array.begin(), array.end());
  double sum =0.0;
  const int iter0 = array.size()*lowerFrac;
  const int iter1 = array.size()*upperFrac;
  for(int ii=iter0; ii< iter1; ii++){
    sum += array[ii];
  }
  return sum / ( (iter1-iter0)+1E-10);
}


//====================================================================
const std::vector<double> Main::Maker::truncatedMean_libo(double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX){

   size_t size = 0;
   std::vector<double> trunc_mean;
   std::vector<double> help_vec;
   truncate_high = 1 - truncate_high; 
   int i_low = 0;
   int i_high = 0;

   //sort the dEdX vecotrs in matrix
   for(auto &&vec : vecs_dEdX){
      size = vec.size();
      help_vec.clear();

      //check dEdX vector isn't empty!
      if(vec.empty()){
         trunc_mean.push_back(-9999.);
         continue;
      }

      else{
         //Sort Vector
         sort(vec.begin(), vec.end());
       
         //Discard upper and lower part of signal
         //rint rounds to integer
         /*i_low = rint ( size*truncate_low);
         i_high = rint( size*truncate_high);
         
         
         for(int i = i_low; i <= i_high; i++){
               help_vec.push_back(vec[i]);
         };*/


	 double median_par=0;		
	 double RMS_par=1.0*TMath::RMS(vec.begin(), vec.end());
	 int N_par=vec.size(); //vec number of hits
	 std::sort(vec.begin(), vec.end());
	 if(N_par<=0) continue;
	 if(N_par % 2 == 0) {median_par=(vec[((N_par/2)-1)]+vec[N_par/2])/2; }
         else if(N_par == 1) {median_par=vec[N_par];}
         else {median_par = vec[(N_par+1)/2];}

         


	 std::vector<double> TLMean_par; 
         TLMean_par.clear();
         for(int ii=0; ii<int(vec.size()); ii++){
                 if(vec[ii]<median_par+1*RMS_par && vec[ii]>median_par-1*RMS_par)
                 {TLMean_par.push_back(vec[ii]);}
         }


         //Mean of help vector
         trunc_mean.push_back(TMath::Mean(TLMean_par.begin(), TLMean_par.end()));
         //trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());

      }


   };

   return trunc_mean;

};





//===================================================================
bool Main::Maker::secondary_noPion( const std::vector<double> &track_score, 
                                    const std::vector<int> &trackID,
                                    const std::vector<double> &dEdX) {
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -1) && (track_score[i] > cut_trackScore) &&
    /*dEdX[i]>=cut_dEdX_low  &&*/  (dEdX[i] <= cut_dEdX_high)) {
      return false;
    }
  }

  return true;
};
//===============================================================
bool Main::Maker::secondary_noPion_libo( const std::vector<double> &track_score, 
                                    const std::vector<int> &trackID,
                                    const std::vector<double> &dEdX) {
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -1) && (track_score[i] > cut_trackScore) &&
    dEdX[i]>=cut_dEdX_low  &&  (dEdX[i] <= cut_dEdX_high)) {
      return false;
    }
  }

  return true;
};


bool Main::Maker::has_shower_nHits_distance(const std::vector<double> &track_score,
                                    const std::vector<int> &nHits,
                                    const std::vector<double> &distance) {

  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (nHits[i] < cut_nHits_shower_high) && (track_score[i] != -999.) &&
         (distance[i] < cut_daughter_shower_distance_high) &&
         (distance[i] > cut_daughter_shower_distance_low)) {
       return true;
     }
  }

  return false;


}
//--------------------------------------------------------------------------------
bool Main::Maker::has_shower_energy_distance(const std::vector<double> &track_score,
                                    const std::vector<double> &energy,
                                    const std::vector<double> &distance) {

  if(track_score.empty() || energy.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (energy[i] > cut_energy_shower_low) &&
         (energy[i] < cut_energy_shower_high) && (track_score[i] != -999.) &&
         (distance[i] < cut_daughter_shower_distance_high) &&
         (distance[i] > cut_daughter_shower_distance_low)) {
       return true;
     }
  }

  return false;


}



//-------------------------------------------------------------------------------------

bool Main::Maker::has_shower_nHits(const std::vector<double> &track_score,
                                  const std::vector<int> &nHits) {
  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};
//-----------------------------------------------------------------------
bool Main::Maker::has_shower_Eng(const std::vector<double> &track_score,
                                  const std::vector<double> &shweng) {
  if(track_score.empty() || shweng.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (shweng[i] > cut_energy_shower_low) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};
//---------------------------------------------------------------
bool Main::Maker::has_shower_Ang(const std::vector<double> &track_score,
                                  const std::vector<double> &shwang) {
  if(track_score.empty() || shwang.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (shwang[i] < cut_angle_shower_low && shwang[i] > cut_angle_shower_high) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};



//------------------------------------------------------------------
bool Main::Maker::has_shower_Dist(const std::vector<double> &track_score,
                                  const std::vector<double> &distance) {
  if(track_score.empty() || distance.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (distance[i] > cut_daughter_shower_distance_low) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};






//--------------------------------------------------------------------
const std::vector<double> Main::Maker::compute_distanceVertex (double beam_endX,
                                 double beam_endY,
                                 double beam_endZ, 
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ) {
  std::vector<double> distance;
  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return distance;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + 
                   diff_Z_end*diff_Z_end);

    dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start +
                   diff_Z_start*diff_Z_start);
    //std::cout<<"dummy_1 "<<dummy_1<<std::endl;
    if(dummy_1 < dummy_2)
      distance.push_back(dummy_1);
    else 
      distance.push_back(dummy_2);
  }

  return distance;
}
const std::vector<double> Main::Maker::compute_angleVertex   ( double beam_endX,
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
                                 const std::vector<double> &d_endPhi){


  std::vector<double> distance;
  std::vector<double> angle;

  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return angle;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + 
                   diff_Z_end*diff_Z_end);

    dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start +
                   diff_Z_start*diff_Z_start);
    //std::cout<<"dummy_1 "<<dummy_1<<std::endl;

    TVector3 startv;
    TVector3 trkv;
    if(dummy_1 < dummy_2){
      distance.push_back(dummy_1);
      startv.SetXYZ(diff_X_end, diff_Y_end, diff_Z_end); 
      trkv.SetXYZ(sin(d_endTheta[i])*cos(d_endPhi[i]),sin(d_endTheta[i])*sin(d_endPhi[i]) ,cos(d_endTheta[i]));
      angle.push_back(startv.Angle(trkv));
    }
    else{ 
      distance.push_back(dummy_2);
      startv.SetXYZ(diff_X_start, diff_Y_start, diff_Z_start);
      trkv.SetXYZ(sin(d_startTheta[i])*cos(d_startPhi[i]),sin(d_startTheta[i])*sin(d_startPhi[i]) ,cos(d_startTheta[i]));
      angle.push_back(startv.Angle(trkv));
   

    }
  } 

  return angle;



}
//------------------------------------------------------------------------
const std::vector<double> Main::Maker::compute_angleVertex_shower   ( double beam_endX,
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
                                 const std::vector<double> &d_dirZ){


  std::vector<double> distance;
  std::vector<double> angle;

  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return angle;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + 
                   diff_Z_end*diff_Z_end);

    dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start +
                   diff_Z_start*diff_Z_start);
    //std::cout<<"dummy_1 "<<dummy_1<<std::endl;

    TVector3 startv;
    TVector3 trkv;
    if(dummy_1 < dummy_2){
      distance.push_back(dummy_1);
      startv.SetXYZ(diff_X_end, diff_Y_end, diff_Z_end); 
      trkv.SetXYZ(d_dirX[i], d_dirY[i], d_dirZ[i]);
      angle.push_back(startv.Angle(trkv));
    }
    else{ 
      distance.push_back(dummy_2);
      startv.SetXYZ(diff_X_start, diff_Y_start, diff_Z_start);
      trkv.SetXYZ(d_dirX[i], d_dirY[i], d_dirZ[i]);
      angle.push_back(startv.Angle(trkv));
   

    }
  } 

  return angle;



}



//--------------------------------------------------------------------------
void Main::Maker::SetMomThreshCut(double vmom)
{
  momthreshcut = vmom;
}
void Main::Maker::SetTrkScoreCut(double vtrk)
{
  trkscorecut = vtrk;
}


bool Main::Maker::inFV(double x, double y, double z) {
  if (x>-330. && x<330. && y>50. && y<550. && z>50. && z<645.) {return true; }
  else {return false;}
}
void Main::Maker::SetFVCut(bool v){
    FVcuton = v;
}

double Main::Maker::Ecalcmiss(double Esum, double PTmiss, int np) {
   Esum *= 1000; //convert to MeV
   PTmiss *= 1000; //convert to MeV
   double Eexcit = 30.4; //in MeV
   double Mass = 0; // in MeV
   if(np == 0) Mass = 37.2050e3; //Ar40
   else if(np == 1) Mass = 36.2758e3; //Ar39
   else if(np == 2) Mass = 35.3669e3; //Cl38
   else if(np == 3) Mass = 34.4201e3; //S37
   else if(np == 4) Mass = 33.4957e3; //P36
   else if(np == 5) Mass = 32.5706e3; //Si35
   else if(np == 6) Mass = 31.6539e3; //Al34
   else if(np == 7) Mass = 30.7279e3; //Mg33
   else if(np == 8) Mass = 29.8111e3; //Na32
   else if(np == 9) Mass = 28.8918e3; //Ne31
   else if(np >= 10) Mass = 27.9789e3; //F30

   float Ekinrecoil = sqrt(PTmiss*PTmiss + Mass*Mass) - Mass;
   return Esum + Eexcit + Ekinrecoil; // return result in MeV
}



double Main::Maker::thetax(double theta,double phi){
  TVector3 v;
  v.SetMagThetaPhi(1,theta,phi);
  TVector3 x_axis(1,0,0);
  double theta_x = v.Angle(x_axis);
  return theta_x;
}

//========================================================
































void Main::Maker::LoadHist(){

  TFile *infile = TFile::Open("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
  std::cout<<"start to load histograms for SCE correction"<<std::endl;
  //Load in files
  TH3F* hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
  TH3F* hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
  TH3F* hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
  TH3F* hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
  TH3F* hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
  TH3F* hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
        hDz_sim_pos_orig_new = (TH3F*)hDz_sim_pos_orig->Clone(); 
  TH3F* hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
  TH3F* hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
  TH3F* hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
  TH3F* hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
  TH3F* hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
  TH3F* hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
        hDz_sim_neg_orig_new = (TH3F*)hDz_sim_neg_orig->Clone();
  TH3F* hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
  TH3F* hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
  TH3F* hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
  TH3F* hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
  TH3F* hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
  TH3F* hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");
  
  TH3F* hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
  TH3F* hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
  TH3F* hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
  TH3F* hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
  TH3F* hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
  TH3F* hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");
  
  hDx_sim_pos->SetDirectory(0);
  hDy_sim_pos->SetDirectory(0);
  hDz_sim_pos->SetDirectory(0);
  hEx_sim_pos->SetDirectory(0);
  hEy_sim_pos->SetDirectory(0);
  hEz_sim_pos->SetDirectory(0);
  
  hDx_sim_neg->SetDirectory(0);
  hDy_sim_neg->SetDirectory(0);
  hDz_sim_neg->SetDirectory(0);
  hEx_sim_neg->SetDirectory(0);
  hEy_sim_neg->SetDirectory(0);
  hEz_sim_neg->SetDirectory(0);
  
  for(int y = 1; y <= 31; y++){
    for(int z = 1; z <= 37; z++){
      spline_dx_fwd_neg[y-1][z-1] = MakeSpline(hDx_sim_neg,1,y,z,1,1);
      spline_dx_fwd_pos[y-1][z-1] = MakeSpline(hDx_sim_pos,1,y,z,1,2);
      spline_dEx_neg[y-1][z-1] = MakeSpline(hEx_sim_neg,1,y,z,3,1);
      spline_dEx_pos[y-1][z-1] = MakeSpline(hEx_sim_pos,1,y,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int z = 1; z <= 37; z++){
      spline_dy_fwd_neg[x-1][z-1] = MakeSpline(hDy_sim_neg,2,x,z,1,1);
      spline_dy_fwd_pos[x-1][z-1] = MakeSpline(hDy_sim_pos,2,x,z,1,2);
      spline_dEy_neg[x-1][z-1] = MakeSpline(hEy_sim_neg,2,x,z,3,1);
      spline_dEy_pos[x-1][z-1] = MakeSpline(hEy_sim_pos,2,x,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int y = 1; y <= 31; y++){
      spline_dz_fwd_neg[x-1][y-1] = MakeSpline(hDz_sim_neg,3,x,y,1,1);
      spline_dz_fwd_pos[x-1][y-1] = MakeSpline(hDz_sim_pos,3,x,y,1,2);
      spline_dEz_neg[x-1][y-1] = MakeSpline(hEz_sim_neg,3,x,y,3,1);
      spline_dEz_pos[x-1][y-1] = MakeSpline(hEz_sim_pos,3,x,y,3,2);
    }
  }
}

TSpline3* Main::Maker::MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const
{
  TSpline3 *spline = 0;
  
  if(dim1 == 1)
  {
    double a[19];
    double b[19];
    for(int x = 1; x <= 19; x++)
    {
      a[x-1] = spline_hist->GetXaxis()->GetBinCenter(x);
      b[x-1] = spline_hist->GetBinContent(x,dim2_bin,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 2)
  {
    double a[31];
    double b[31];
    for(int y = 1; y <= 31; y++)
    {
      a[y-1] = spline_hist->GetYaxis()->GetBinCenter(y);
      b[y-1] = spline_hist->GetBinContent(dim2_bin,y,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 3)
  {
    double a[37];
    double b[37];
    for(int z = 1; z <= 37; z++)
    {
      a[z-1] = spline_hist->GetZaxis()->GetBinCenter(z);
      b[z-1] = spline_hist->GetBinContent(dim2_bin,dim3_bin,z);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }

  return spline;
}

double Main::Maker::InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const
{
  int bin_x = interp_hist->GetXaxis()->FindBin(xVal);
  int bin_y = interp_hist->GetYaxis()->FindBin(yVal);
  int bin_z = interp_hist->GetZaxis()->FindBin(zVal);

  int bincenter_x = interp_hist->GetXaxis()->GetBinCenter(bin_x);
  int bincenter_y = interp_hist->GetYaxis()->GetBinCenter(bin_y);
  int bincenter_z = interp_hist->GetZaxis()->GetBinCenter(bin_z);

  int max_x = interp_hist->GetNbinsX();
  int max_y = interp_hist->GetNbinsY();
  int max_z = interp_hist->GetNbinsZ();
  
  int low_x;
  int high_x;
  if(bin_x <= 1)
  {
    low_x = 1;
    high_x = 2;
  }
  else if(bin_x >= max_x)
  {
    low_x = max_x-1;
    high_x = max_x;
  }
  else if(xVal > bincenter_x)
  {
    low_x = bin_x;
    high_x = bin_x+1;
  }
  else
  {
    low_x = bin_x-1;
    high_x = bin_x;
  }

  int low_y;
  int high_y;
  if(bin_y <= 1)
  {
    low_y = 1;
    high_y = 2;
  }
  else if(bin_y >= max_y)
  {
    low_y = max_y-1;
    high_y = max_y;
  }
  else if(yVal > bincenter_y)
  {
    low_y = bin_y;
    high_y = bin_y+1;
  }
  else
  {
    low_y = bin_y-1;
    high_y = bin_y;
  }

  int low_z;
  int high_z;
  if(bin_z <= 1)
  {
    low_z = 1;
    high_z = 2;
  }
  else if(bin_z >= max_z)
  {
    low_z = max_z-1;
    high_z = max_z;
  }
  else if(zVal > bincenter_z)
  {
    low_z = bin_z;
    high_z = bin_z+1;
  }
  else
  {
    low_z = bin_z-1;
    high_z = bin_z;
  }

  double interp_val = 0.0;
  
  if(dim == 1)
  {
    double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_neg[high_y-1][high_z-1]->Eval(xVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_pos[high_y-1][high_z-1]->Eval(xVal);
      }
    }

    interp_val = (f_11*(a_2-yVal)*(b_2-zVal) + f_21*(yVal-a_1)*(b_2-zVal) + f_12*(a_2-yVal)*(zVal-b_1) + f_22*(yVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 2)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_neg[high_x-1][high_z-1]->Eval(yVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_pos[high_x-1][high_z-1]->Eval(yVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-zVal) + f_21*(xVal-a_1)*(b_2-zVal) + f_12*(a_2-xVal)*(zVal-b_1) + f_22*(xVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 3)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_neg[high_x-1][high_y-1]->Eval(zVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_pos[high_x-1][high_y-1]->Eval(zVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-yVal) + f_21*(xVal-a_1)*(b_2-yVal) + f_12*(a_2-xVal)*(yVal-b_1) + f_22*(xVal-a_1)*(yVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }

  return interp_val;
}

bool Main::Maker::IsInsideBoundaries(TVector3 const& point) const
{
  return !(
           (TMath::Abs(point.X()) <= 0.0) || (TMath::Abs(point.X()) >= 360.0)
           || (point.Y()             <= 5.2) || (point.Y()             >= 604.0)
           || (point.Z()             <= -0.5) || (point.Z()             >= 695.3)
           );
} 
  
bool Main::Maker::IsTooFarFromBoundaries(TVector3 const& point) const
{
  return (
          (TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X())  >= 360.0)
          || (point.Y()             < -14.8) || (point.Y()              >  624.0)
          || (point.Z()             < -20.5) || (point.Z()              >  715.3)
          );
}

TVector3 Main::Maker::PretendAtBoundary(TVector3 const& point) const
{
  double x = point.X(), y = point.Y(), z = point.Z();
  
  
  if      (TMath::Abs(point.X()) ==    0.0    ) x =                           -0.00001;
  else if (TMath::Abs(point.X()) <	 0.00001) x =   TMath::Sign(point.X(),1)*0.00001; 
  else if (TMath::Abs(point.X()) >=    360.0  ) x = TMath::Sign(point.X(),1)*359.99999;
  
  if      (point.Y() <=   5.2) y =   5.20001;
  else if (point.Y() >= 604.0) y = 603.99999;
  
  if      (point.Z() <=   -0.5) z =   -0.49999;
  else if (point.Z() >= 695.3) z = 695.29999;
  
  return {x, y, z};
}







void Main::Maker::MakeFile() 
{

	clock_t begin = clock();



  system("mkdir -p output/");
  

   

  

  

  //*************************
  //* Starting ROOT application
  //*************************

  
  //TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  //gROOT->ProcessLine(".x rootlogon.C");



  LOG_NORMAL() << "Opening output file with name " << fileoutn << std::endl;
  TFile *file_out = new TFile(fileoutn.c_str(),"RECREATE");
  if ( file_out->IsOpen() ) {
    LOG_NORMAL() << "File opened successfully." << std::endl;
  } else {
    LOG_CRITICAL() << "File not opened (maybe not found?). File: " << fileoutn << std::endl;
    exit(0);
  }
  
  string pattern = filen;
  //string pattern_add = filen_add;
  
  
  
  //*************************
  //* Getting POTs
  //*************************
  
  //TFile* infSCE = new TFile(pattern_add.c_str()); 
  
  //TH3F *RecoFwd_Displacement_Z_Neg = (TH3F*)infSCE->Get("RecoFwd_Displacement_Z_Neg");
  //TH3F *RecoFwd_Displacement_Z_Pos = (TH3F*)infSCE->Get("RecoFwd_Displacement_Z_Pos");

  //LOG_NORMAL()<<"The size of the SCE correction is "<<RecoFwd_Displacement_Z_Pos->GetNbinsX()<<std::endl;  
  
  TChain *chain_ubxsec;
  chain_ubxsec = new TChain("pionana/beamana");
  chain_ubxsec->Add(pattern.c_str());
  //chain_ubxsec->Add(pattern_add.c_str()); 

 
  LOG_NORMAL() << "Using file: " << pattern << endl;
  
  int Nfiles = chain_ubxsec->GetNtrees();
  LOG_NORMAL() << "Number of files: " << Nfiles << endl;
  
  int evts = chain_ubxsec -> GetEntries();
  LOG_NORMAL() << "Number of events used is: " << evts << endl;
   
  UBXSecEvent * t = new UBXSecEvent(chain_ubxsec);
  //ActivateBranches(t);
  LOG_NORMAL() << "Declear the OBJECT of UBXSecEvent"<<std::endl;
  _event_histo_1d = new UBXSecEventHisto1D();
  _event_histo_1d->InitializeBootstraps();

  _event_histo = new UBXSecEventHisto();
  _event_histo->InitializeBootstraps();

  //========================================================================
  //int evts = chain_ubxsec->GetEntries();
  LOG_NORMAL()<<"total number of events is "<<evts<<std::endl;
  //loop over all the events
  int barWidth = 70;

  int Ntotal_beam = 0;
  int Ntotal_tmdqdx = 0;
  int Ntotal_chi2 = 0;
  int Ntotal_shwid = 0;

  int Noriabs = 0;
  int Norichx = 0;
  int Norirea = 0;
  int Noriother = 0;


  int Noriabs_withthresh=0;
  int Norichx_withthresh=0;
  int Norirea_withthresh=0;
  int Noriother_withthresh = 0;

   


  int Nshwcutabs_withthresh = 0;
  int Nshwcutchx_withthresh = 0;
  int Nshwcutrea_withthresh = 0;
  int Nshwcutother_withthresh = 0;


  int Ntmcutabs_withthresh = 0;
  int Ntmcutchx_withthresh = 0;
  int Ntmcutrea_withthresh = 0;
  int Ntmcutother_withthresh = 0;



  int Nchi2abs_withthresh = 0;
  int Nchi2chx_withthresh = 0;
  int Nchi2rea_withthresh = 0;
  int Nchi2other_withthresh = 0;

  int TestGenSig=0;

  int n_bins_mumom = 50;
  int n_bins_mucostheta = 50;
  int n_bins_muphi = 50;
  double bins_mumom[51];
  double bins_mucostheta[51];
  double bins_muphi[51];
  for(int i=0; i<51; i++){
     bins_mumom[i]=0.0 + 1.2/50.0*i;
     bins_mucostheta[i]=-1.0 + 2.0/50.0*i;
     bins_muphi[i]=-TMath::Pi() + 2*TMath::Pi()*i/50.0;
  }

  TH2D *h_pcostheta1st2nd_ptmissing = new TH2D("h_pcostheta1st2nd_ptmissing", "h_pcostheta1st2nd_ptmissing", 50, 0.0, 1.2, 50, -1.0, 1.0);
  TH1D *h_sel_pcostheta1st2nd=new TH1D("h_sel_pcostheta1st2nd", "h_sel_pcostheta1st2nd", 50, -1.0, 1.0);
  TH1D *h_sig_pcostheta1st2nd=new TH1D("h_sig_pcostheta1st2nd", "h_sig_pcostheta1st2nd", 50, -1.0, 1.0);
  TH1D *h_chxbac_pcostheta1st2nd=new TH1D("h_chxbac_pcostheta1st2nd", "h_chxbac_pcostheta1st2nd", 50, -1.0, 1.0);
  TH1D *h_reabac_pcostheta1st2nd=new TH1D("h_reabac_pcostheta1st2nd", "h_reabac_pcostheta1st2nd", 50, -1.0, 1.0);

  TH2D *h_true_reco_mom = new TH2D("h_true_reco_mom", "h_true_reco_mom", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
  TH2D *h_true_reco_costheta = new TH2D("h_true_reco_costheta", "h_true_reco_costheta", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
  TH2D *h_true_reco_phi = new TH2D("h_true_reco_phi", "h_true_reco_phi",  n_bins_muphi, bins_muphi, n_bins_muphi, bins_muphi);



  // True v.s. reco histograms for constructing smearing matrix
  std::map<std::string,TH2D*> bs_geant_pm1_true_reco_mom;
  bs_geant_pm1_true_reco_mom["nominal"] = new TH2D("bs_geant_pm1_true_reco_mom_nominal", ";Proton Momentum (Truth) [GeV]; Proton Momentum (Reco) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
   
  std::map<std::string,TH2D*> bs_geant_pm1_true_reco_costheta;
  bs_geant_pm1_true_reco_costheta["nominal"] = new TH2D("bs_geant_pm1_true_reco_costheta_nominal", ";Proton CosTheta(Truth); Proton CosTheta (Reco) [GeV]", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
 

 
  std::map<std::string, std::map<std::string,TH1D*>> hmap_trkmom_geant_pm1_bs;
  hmap_trkmom_geant_pm1_bs["total"]["nominal"]=new TH1D("h_trkmom_total_geant_pm1_nominal", ";Track Momentum;", 30, 0.0, 1.2);
  hmap_trkmom_geant_pm1_bs["signal"]["nominal"]=new TH1D("h_trkmom_signal_geant_pm1_nominal", ";Track Momentum;", 30, 0.0, 1.2);
  hmap_trkmom_geant_pm1_bs["chxbac"]["nominal"]=new TH1D("h_trkmom_chxbac_geant_pm1_nominal", ";Track Momentum;", 30, 0.0, 1.2);
  hmap_trkmom_geant_pm1_bs["reabac"]["nominal"]=new TH1D("h_trkmom_reabac_geant_pm1_nominal", ";Track Momentum;", 30, 0.0, 1.2);

  std::map<std::string, std::map<std::string,TH1D*>> hmap_trkcostheta_geant_pm1_bs;
  hmap_trkcostheta_geant_pm1_bs["total"]["nominal"] =new TH1D("h_trkcostheta_total_geant_pm1_nominal", ";Track CosTheta;", 30, -1.0, 1.0);
  hmap_trkcostheta_geant_pm1_bs["signal"]["nominal"]=new TH1D("h_trkcostheta_signal_geant_pm1_nominal", ";Track CosTheta;", 30, -1.0, 1.0);
  hmap_trkcostheta_geant_pm1_bs["chxbac"]["nominal"]=new TH1D("h_trkcostheta_chxbac_geant_pm1_nominal", ";Track CosTheta;", 30, -1.0, 1.0);
  hmap_trkcostheta_geant_pm1_bs["reabac"]["nominal"]=new TH1D("h_trkcostheta_reabac_geant_pm1_nominal", ";Track CosTheta;", 30, -1.0, 1.0);
  

  //Efficiency - GEANT pm1sigma
  BootstrapTH1D bs_geant_pm1_eff_mumom_num("bs_geant_pm1_eff_mumom_num", "bs_geant_pm1_eff_mumom_num_title", n_bins_mumom, bins_mumom);
  BootstrapTH1D bs_geant_pm1_eff_mumom_den("bs_geant_pm1_eff_mumom_den", "bs_geant_pm1_eff_mumom_den_title", n_bins_mumom, bins_mumom);

  //Efficiency - GEANT pm1sigma
  BootstrapTH1D bs_geant_pm1_eff_mucostheta_num("bs_geant_pm1_eff_mucostheta_num", "bs_geant_pm1_eff_mucostheta_num_title", n_bins_mucostheta, bins_mucostheta);
  BootstrapTH1D bs_geant_pm1_eff_mucostheta_den("bs_geant_pm1_eff_mucostheta_den", "bs_geant_pm1_eff_mucostheta_den_title", n_bins_mucostheta, bins_mucostheta);



  TH1D* h_Evttot = new TH1D("h_Evttot", "h_Evttot", 1, 0.5, 1.5);
  TH1D *dabsintE = new TH1D("dabsintE","#DeltaE(Reco - True)/True ",100, -3.0,2.0);
  TH1D *dchxintE = new TH1D("dchxintE",   "#DeltaE(Reco - True)/True",100, -3.0,2.0);
  
  TH1D *dsizewire_abs[nslices];
  TH1D *dsizewire_chx[nslices]; 

  TH1D *dslcID[nslices];  

  TH1D *incE[nslices];
  TH1D *pitch[nslices];

  TH1D *incE_gen[nslices];

  TH1D *incE_sel[nslices];

  TH1D *h_energetic_pmom_gen[nslices];
  TH1D *h_energetic_pcostheta_gen[nslices];
  TH1D *h_energetic_pphi_gen[nslices];
  


  TGraph *gr_inc0_slc;

  TGraph *gr_int_slc;
 
  TGraph *gr_intabs_slc;
  TGraph *gr_intchx_slc;

  TGraph *gr_inc_slc;

  TGraphAsymmErrors *gr_incE_slc;


  TGraphAsymmErrors *gr_pitch_slc;


  TGraph *gr_intsel_slc;

  TGraph *gr_intgen_slc;

  TGraph *gr_efficiency_slc;

  TGraph *gr_tempxsec_slc;
  TGraph *gr_tempxsec_chx_slc;





  TGraph *gr_tempxsec_test_slc;
  TGraph *gr_tempxsec_test2_slc;
  TGraph *gr_tempxsec_test3_slc;
  TGraph *gr_tempxsec_test4_slc;

  TGraph *gr_tempxsec_chx_test_slc;
  TGraph *gr_tempxsec_abs_test_slc;

  TGraph *gr_tempxsec_vs_incE_slc;

  TGraph *gr_tempxsec_test_vs_incE_slc;

  TGraph *gr_tempxsec_chx_vs_incE_slc;

  TGraph *gr_selected_tot_slc;


  TGraph *gr_selected_bkg_slc;


   
  std::vector<std::string> fname_geant_pm1;
  std::vector<int> fpion_evt_index;
  fname_geant_pm1.clear();
  fpion_evt_index.clear();
  for(int i= _initial_entry; i<evts; i++){
	chain_ubxsec->GetEntry(i);
        if(isdata==0){
        if(abs(t->true_beam_PDG) != 211 && abs(t->true_beam_PDG) !=13) continue; 
        if(!isBeamType(t->reco_beam_type)) continue;
        if(!manual_beamPos_mc(t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ, 
                              t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                              t->true_beam_startDirX, t->true_beam_startDirY, t->true_beam_startDirZ,
                              t->true_beam_startX, t->true_beam_startY, t->true_beam_startZ)) continue;
        if(!endAPA3(t->reco_beam_endZ) )continue;
        fpion_evt_index.push_back(i);
        }
  }
  
  int totalshower_trkscore=0;
  int totalshower_trknhits=0;
  int totalshower_trkdist=0;
  int totalshower_trkeng=0;
 
  int totalgamma_trkscore=0;
  int totalgamma_trknhits=0;
  int totalgamma_trkdist=0;
  int totalgamma_trkeng=0;
 
  //std::cout<<"start looping over all the events the first pion event is "<<fpion_evt_index[0]<<std::endl;
  std::cout<<"initial entry is "<<_initial_entry<<" total number of events is "<<evts<<std::endl;
  //get the first pion event
  double event_weight = 1.0;

  Evttot = 0;
  for(int ix=0; ix<nslices+1; ix++){
      for(int iy=0; iy<nslices+1; iy++){
        intabs_array_den[ix][iy]=0.0;
        intchx_array_den[ix][iy]=0.0;
        intabs_array_num[ix][iy]=0.0;
        intchx_array_num[ix][iy]=0.0;
      }
  }



  for (int i = 0; i<nslices+1; ++i){
     interaction[i] = 0;
     incident[i] = 0;
     true_incident[i]=0;
     true_interaction[i]=0;
     true_abs[i]=0;
     true_chx[i]=0;
     interaction_gen[i] = 0;
     interaction_sel[i] = 0;

     selected_tot[i]=0;
     selected_bkg[i]=0;


     incident_gen[i] = 0;
     incident_sel[i] = 0;
     if(i<nslices){
     dslcID[i] = new TH1D(Form("dslcID_%d", i),Form("dslcID, %d<wire #<%d", i*nwires_in_slice, (i+1)*nwires_in_slice),30,-14.5,15.5);
     dslcID[i]->SetDirectory(file_out);
     incE[i] = new TH1D(Form("incE_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     //incE[i]->SetDirectory(file_out);
     pitch[i] = new TH1D(Form("pitch_%d",i),Form("Slice thickness, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
     //pitch[i]->SetDirectory(file_out);
     incE_gen[i] = new TH1D(Form("incE_gen_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     //incE_gen[i]->SetDirectory(file_out);
     incE_sel[i] = new TH1D(Form("incE_sel_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     //incE_sel[i]->SetDirectory(file_out);
     h_energetic_pmom_gen[i] = new TH1D(Form("h_energetic_pmom_gen_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinspmom, 0, 1.200);
     h_energetic_pmom_gen[i]->SetDirectory(file_out);
     h_energetic_pcostheta_gen[i] = new TH1D(Form("h_energetic_pcostheta_gen_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinspcostheta, -1.0, 1.0);
     h_energetic_pcostheta_gen[i]->SetDirectory(file_out);
     h_energetic_pphi_gen[i] = new TH1D(Form("h_energetic_pphi_gen_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinspphi, -TMath::Pi(), TMath::Pi());
     h_energetic_pphi_gen[i]->SetDirectory(file_out);
     dsizewire_abs[i] = new TH1D(Form("dsizewire_abs_%d", i), Form("dsizewire_abs_%d", i), 26, -5.5, 20.5); 
     //dsizewire_abs[i]->SetDirectory(file_out);
     dsizewire_chx[i] = new TH1D(Form("dsizewire_chx_%d", i), Form("dsizewire_chx_%d", i), 26, -5.5, 20.5); 
     //dsizewire_chx[i]->SetDirectory(file_out);
     }
  } 




  LoadHist();


  for(int i= _initial_entry; i<evts; i++){
	if(i !=0) DrawProgressBar((double)i/(double)evts, barWidth);
        
        
	chain_ubxsec->GetEntry(i);
        //std::cout<<"i = "<<i<<std::endl;
        if(isdata==0){
        if(i==fpion_evt_index[0] && isdata==0 && _fill_bootstrap_geant){
          //std::cout<<"total g4 reweight factor size is "<<t->g4rw_primary_var->size()<<std::endl;
          for(long unsigned int k=0; k<t->g4rw_primary_var->size(); k++){
             string name = t->g4rw_primary_var->at(k);
             fname_geant_pm1.push_back(name + "_p1");
             fname_geant_pm1.push_back(name + "_m1");
             
          }    
          //Number of events
          for(auto iter : hmap_trkmom_geant_pm1_bs){
             std::string this_name = iter.first;
             std::map<std::string, TH1D*> bs_map = iter.second;
             for(size_t i=0; i<fname_geant_pm1.size(); i++){
                 std::string histo_name = "h_trkmom_"+this_name+"_"+fname_geant_pm1.at(i);
                 std::cout<<"histo_name = "<<histo_name<<std::endl;
                 //double this_bins_mom
                 hmap_trkmom_geant_pm1_bs[this_name][fname_geant_pm1.at(i)]=new TH1D(histo_name.c_str(), ";Track Momentum[GeV];", 30, 0.0, 1.2);
                 histo_name = "h_trkcostheta_"+this_name+"_"+fname_geant_pm1.at(i);

                 hmap_trkcostheta_geant_pm1_bs[this_name][fname_geant_pm1.at(i)]=new TH1D(histo_name.c_str(), ";Track CosTheta;", 30, -1.0, 1.0);
                  
             }
          }
          //Efficiency 
          for(size_t i=0; i<fname_geant_pm1.size(); i++){
               std::string histo_name;
               histo_name ="bs_geant_pm1_true_reco_mom_"+fname_geant_pm1.at(i);
               bs_geant_pm1_true_reco_mom[fname_geant_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Proton Momentum (Truth) [GeV]; Proton Momentum (Reco)[GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);

               histo_name ="bs_geant_pm1_true_reco_costheta_"+fname_geant_pm1.at(i);
               bs_geant_pm1_true_reco_costheta[fname_geant_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Proton Momentum (Truth); Proton CosTheta(Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
       
          }
          bs_geant_pm1_eff_mumom_num.SetWeightNames(fname_geant_pm1);
          bs_geant_pm1_eff_mumom_den.SetWeightNames(fname_geant_pm1);
          
          bs_geant_pm1_eff_mucostheta_num.SetWeightNames(fname_geant_pm1);
          bs_geant_pm1_eff_mucostheta_den.SetWeightNames(fname_geant_pm1);
          
        }  

        }//end of is isdata==0
        //std::cout<<"end of setting all the Weight Names of the bootstraps; is it data? "<<isdata<<std::endl;	
	//========================================================================
	//select muon and pion beam events
	if(isdata==0){
	
          //if(abs(t->true_beam_PDG) != 211 && abs(t->true_beam_PDG) !=13) continue;
          if(abs(t->true_beam_PDG) != 211) continue;
          //Evttot += event_weight; 
          //h_Evttot->Fill(event_weight);
          if(!isBeamType(t->reco_beam_type)) continue;
          if(!manual_beamPos_mc(t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ, 
                              t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                              t->true_beam_startDirX, t->true_beam_startDirY, t->true_beam_startDirZ,
                              t->true_beam_startX, t->true_beam_startY, t->true_beam_startZ)) continue;
          //if(!endAPA3(t->reco_beam_endZ) )continue;
        } 
     
        if(isdata==1){
          if(!data_beam_PID(t->beam_inst_PDG_candidates)) continue;
          //Evttot += event_weight; 
          //h_Evttot->Fill(event_weight);
          if(!manual_beamPos_data(t->event, t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ,
                                t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                                t->beam_inst_X, t->beam_inst_Y, t->beam_inst_dirX, t->beam_inst_dirY,
                                t->beam_inst_dirZ, t->beam_inst_nMomenta, t->beam_inst_nTracks)) continue;
          //if(!endAPA3(t->reco_beam_endZ) )continue;
        }
        h_Evttot->Fill(event_weight);
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        //set the map of wire number and incident particles
        //hit-by-hit wire z position
        for(unsigned int ih=0; ih<t->reco_beam_calo_wire_z->size(); ih++){
            _event_histo->htotinc_reco_beamz->Fill(t->reco_beam_calo_wire_z->at(ih));
            _event_histo->htotinc_reco_beamwire->Fill(t->reco_beam_calo_wire->at(ih));
        }
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        int sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        if(sliceID<0) sliceID=0;
        if(sliceID>=nslices) sliceID = nslices;

        double true_endz = Sce_Corrected_endZ_2nd(t->true_beam_endX, t->true_beam_endY, t->true_beam_endZ);
        std::cout<<"Sce_Corrected "<<t->true_beam_endZ<<"  to   "<<true_endz<<std::endl;
        int true_sliceID = -1;
        string *temp_Ptr0 = t->true_beam_endProcess;
        if (*temp_Ptr0 == "pi+Inelastic" && abs(t->true_beam_PDG) == 211){  //std::cout<<"signal"<<std::endl;

           true_sliceID = int((true_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); //z=0.56 is z coordinate for wire 0
           _event_histo->htotinc_true_beamwire->Fill((true_endz-0.5603500-0.479/2)/0.479);           
           if (true_sliceID <0) true_sliceID = 0;
           if (true_sliceID >= nslices) sliceID = nslices;
           if (true_sliceID <= nslices){
              ++true_interaction[true_sliceID];
              if (t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 && t->true_daughter_nPiMinus ==0){//true absorption
                  ++true_abs[true_sliceID];
                  dabsintE -> Fill((t->reco_beam_interactingEnergy - t->true_beam_interactingEnergy)/t->true_beam_interactingEnergy);
                  ++intabs_array_den[true_sliceID][sliceID];

              }
              if (t->true_daughter_nPi0 > 0 && t->true_daughter_nPiPlus == 0 && t->true_daughter_nPiMinus ==0){
                  ++true_chx[true_sliceID];
                  dchxintE -> Fill((t->reco_beam_interactingEnergy - t->true_beam_interactingEnergy)/t->true_beam_interactingEnergy);
                  ++intchx_array_den[true_sliceID][sliceID];
              } 
           }
           for (int i = 0; i<=true_sliceID; ++i){
              if (i<nslices) ++true_incident[i];
           }
        }//end of selecting pion inclusive events
        
        //LOG_NORMAL()<<"Start to calculate the average energy and thickness with Tingjun's method"<<std::endl;
        //interaction slice ID based on the track end wire number
        //int sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        //double reco_endz = Sce_Corrected_endZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
        //double reco_endz = t->reco_beam_endZ; 
        //sliceID = int((reco_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); 
        if (true_sliceID!=-1){
          dslcID[true_sliceID]->Fill(sliceID - true_sliceID);
        }
        //increment interaction counter
        if(sliceID<=nslices){
            ++interaction[sliceID];
        }
            //increment incident counter
        for (int i = 0; i<=sliceID; ++i) {
           if(i<=nslices){
              ++incident[i];
           }
        }

        //LOG_NORMAL()<<"Start to Fill Vector to store energy and thickness"<<std::endl;

        std::vector<std::vector<double>> vpitch(nslices);
        std::vector<std::vector<double>> vincE(nslices); 

        

        //LOG_NORMAL()<<" size of wire, incE and trkPitch are : "<<t->reco_beam_calo_wire->size()<<"  "
        //            <<t->reco_beam_incidentEnergies->size()<<"  "
        //            <<t->reco_beam_TrkPitch->size()<<std::endl;

        //if (sliceID>=nslices) continue;
        if(t->reco_beam_incidentEnergies->size()>0){
        for (size_t i = 0; i<t->reco_beam_calo_wire->size(); ++i){
          int this_wire = t->reco_beam_calo_wire->at(i);
          int this_sliceID = this_wire/nwires_in_slice;
          //ignore the last slice for pitch and incident energy calculations
          if (this_sliceID>sliceID || sliceID>=nslices) continue;
          
          double this_incE = t->reco_beam_incidentEnergies->at(i);
          double this_pitch = t->reco_beam_TrkPitch->at(i);
          //std::cout<<"i= "<<i<<" this_pitch is "<<this_pitch<< "  this_incE  is "<<this_incE<<std::endl;
          vpitch[this_sliceID].push_back(this_pitch);
          vincE[this_sliceID].push_back(this_incE);
        }
        }
        //LOG_NORMAL()<<"Calculate the pitch in each slice "<<vpitch.size()<<std::endl;
        for (size_t i = 0; i<vpitch.size(); ++i){
          //std::cout<<"vpitch[i] size = "<<vpitch[i].size()<<std::endl;
          if (!vpitch[i].empty()){
          double sum_pitch = 0;
          for (size_t j = 0; j<vpitch[i].size(); ++j){
            //std::cout<<"i= "<<i<<"  j=  "<<j<<" vpitch at ij is "<<vpitch[i][j]<<std::endl;
            sum_pitch += vpitch[i][j];
          }
         
          pitch[i]->Fill(sum_pitch/vpitch[i].size()*nwires_in_slice);
         
          } //else {continue;}
        }
        //LOG_NORMAL()<<"Calculate the energy in each slice "<<std::endl;
        for (size_t i = 0; i<vincE.size(); ++i){
          if (!vincE[i].empty()){
          double sum_incE = 0;
          for (size_t j = 0; j<vincE[i].size(); ++j){
            sum_incE += vincE[i][j];
          }
          incE[i]->Fill(sum_incE/vincE[i].size());
          }
        }
        //LOG_NORMAL()<<"End of calculating the thickness and energy with Tingjun's method"<<std::endl;
 

        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        for(long unsigned int i=0; i<sizeof(slicebins)/sizeof(slicebins[0]); i++){
          slicebins[i]=20*i;
        }
        //ignore the last hit in the calculation of thickness and energy
        //??????????????????????????????????/???????
        
        for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){
            double eincident = 0.0;
            int nincident = 0;
            double slicethickness = 0.0;
            for(int nh=slicebins[is]; nh<slicebins[is+1]; nh++) {
                 if(std::find(t->reco_beam_calo_wire->begin(), t->reco_beam_calo_wire->end(), nh) != t->reco_beam_calo_wire->end()) {
                    if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                       nincident = nincident+1;
                       int realindex = -999;
                       //nh is the wire number find the readindex
                       for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size()-1; kk++){
                           if(t->reco_beam_calo_wire->at(kk)==nh) {realindex = kk;}
                       }
                       if(realindex>=0){ 
                          eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                          slicethickness = slicethickness + t->reco_beam_TrkPitch->at(realindex);
                       }
                    }
                 }//end of if the wire number  
            } // end of loop over all wire in a certain slices [is]
            if(nincident>0){
               eincident = eincident/(nincident*1.0);
               //incidentEnergy = average value of all the hits within the slice
               
               _event_histo->htotinc_reco_beam_incidentEnergy[is]->Fill(eincident);
               _event_histo->htotinc_reco_beam_sliceThickness[is]->Fill(slicethickness);
            }
            if(t->reco_beam_calo_wire->back()>slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){
                if(isdata==0){
                   _event_histo->htotinc_true_beam_trajz[is]->Fill(t->true_beam_traj_Z->back());
                }
                _event_histo->htotinc_reco_beam_endz[is]->Fill(t->reco_beam_calo_endZ);
            }
        }//loop over all the slices
        if(isdata==0){
           _event_histo->htotinc_reco_beamz_reso->Fill((t->reco_beam_calo_endZ-t->true_beam_traj_Z->back()));
           _event_histo->htotinc_reco_beamwire_reso->Fill((t->reco_beam_calo_wire->back()-t->true_beam_slices->back()));
        }



        //========================================================================================       
        //Fill true proton momentum
        unsigned int ngen_proton=0; //number of generated protons
        unsigned int ngen_par=0;  // number of generated particles
        unsigned int ngen_pipm_withthresh=0;        
        if(isdata==0){
	for(unsigned int ind=0; ind<t->true_beam_daughter_PDG->size(); ind++){

           if(abs(t->true_beam_daughter_PDG->at(ind))==2212) 
           {_event_histo_1d->h_orimom_proton->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==2112) 
           {_event_histo_1d->h_orimom_neutron->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==211) 
           {
               _event_histo_1d->h_orimom_pionpm->Fill(t->true_beam_daughter_startP->at(ind));
               if(t->true_beam_daughter_startP->at(ind)>momthreshcut){
                  ngen_pipm_withthresh ++;
               }
              
           }
           else if(abs(t->true_beam_daughter_PDG->at(ind))==111) 
           {_event_histo_1d->h_orimom_pion0->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==321) 
           {_event_histo_1d->h_orimom_kaon->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==11) 
           {_event_histo_1d->h_orimom_electron->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==13) 
           {_event_histo_1d->h_orimom_muon->Fill(t->true_beam_daughter_startP->at(ind));
           }
           else if(abs(t->true_beam_daughter_PDG->at(ind))==22) 
           {_event_histo_1d->h_orimom_photon->Fill(t->true_beam_daughter_startP->at(ind));}
           else {_event_histo_1d->h_orimom_other->Fill(t->true_beam_daughter_startP->at(ind));

           }


           ngen_par++;
           if(abs(t->true_beam_daughter_PDG->at(ind))!=2212) continue;
 		     ngen_proton++;
        }
        }//end of if isdata==0 

        bool isSignal = false;
        bool isChxBKG = false;
        bool isReaBKG = false;
        bool isOtherBKG = false;

        bool isSignal_withthresh = false;
        bool isChxBKG_withthresh = false;
        bool isReaBKG_withthresh = false;
        bool isOtherBKG_withthresh = false;

        unsigned int nproton02=0; //number of reconstructed  protons
        unsigned int npion02=0; //number of reconstructed pions
        if(isdata==0){
        for(unsigned int vcand=0; vcand<t->reco_daughter_allTrack_len->size(); vcand++){
           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(vcand)) == 2212) {nproton02++;}
           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(vcand)) == 211) {

              npion02++;

           }

        }
        string *temp_Ptr=t->true_beam_endProcess;
       
        if(abs(t->true_beam_PDG)==211 && *temp_Ptr == "pi+Inelastic"){

        if(ngen_pipm_withthresh==0 && t->true_daughter_nPi0 == 0){
              _event_histo_1d->h_true_beam_endE_den->Fill(TMath::Sqrt(t->true_beam_endP*t->true_beam_endP + PionMass*PionMass));
              isSignal_withthresh = true;
              Noriabs_withthresh++;
        } else if(t->true_daughter_nPi0 > 0 && ngen_pipm_withthresh ==0){
              isChxBKG_withthresh = true;
              Norichx_withthresh++;
        } else if(t->true_daughter_nPi0 ==0 && ngen_pipm_withthresh>0){
              isReaBKG_withthresh = true;
              Norirea_withthresh++;
        } 
        }//enf of if this is a pion beam event 
        else {
              isOtherBKG_withthresh = true;
              Noriother_withthresh++;
        }
 
        if(abs(t->true_beam_PDG)==211 && *temp_Ptr == "pi+Inelastic"){

        if(t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 
           && t->true_daughter_nPiMinus ==0 ) {
		Noriabs++;
		isSignal = true;
                _event_histo_1d->h_true_sig_trkmult->Fill(nproton02);
        }
        else if(t->true_daughter_nPi0 > 0 && t->true_daughter_nPiPlus == 0
           && t->true_daughter_nPiMinus ==0) {Norichx++;
                isChxBKG = true;
                _event_histo_1d->h_true_chxbac_trkmult->Fill(npion02+nproton02);
        }
        else if(t->true_daughter_nPi0 == 0 && (t->true_daughter_nPiPlus > 0 
           || t->true_daughter_nPiMinus >0)  ) {Norirea++;
                isReaBKG = true;
                _event_histo_1d->h_true_reabac_trkmult->Fill(npion02+nproton02);
        }
        }
        else {  Noriother++;
                isOtherBKG = true;
        }
        }//end of if isdata
        //check all the generated absorption/charge exchange events !!!!!!!!!!!
        //Fill the histograms for the proton momentum and angle
        if(isSignal_withthresh){
              //loop over all the MC particles in signal event
              int temp_genpindex=-999.0; double temp_genpmom=-999.0; double temp_genpcostheta=-999.0; double temp_genpphi=-999.0;
              for(unsigned int tmk=0; tmk<t->true_beam_daughter_startP->size(); tmk++){
                      if(abs(t->true_beam_daughter_PDG->at(tmk)) !=2212) continue;
                      _event_histo->h_true_pphi_test->Fill(TMath::ATan2(t->true_beam_daughter_startPy->at(tmk), t->true_beam_daughter_startPx->at(tmk)));
                      if(t->true_beam_daughter_startP->at(tmk)>temp_genpmom){
                      temp_genpindex = tmk;
                      temp_genpmom=t->true_beam_daughter_startP->at(tmk);
                      temp_genpcostheta=t->true_beam_daughter_startPz->at(tmk)/t->true_beam_daughter_startP->at(tmk);
                      temp_genpphi= TMath::ATan2(t->true_beam_daughter_startPy->at(tmk), t->true_beam_daughter_startPx->at(tmk));
                      }           
              }
                       _event_histo->h_eff_pmom_den->Fill(temp_genpmom);
                       _event_histo->h_eff_pcostheta_den->Fill(temp_genpcostheta);
                       _event_histo->h_eff_pphi_den->Fill(temp_genpphi);
        //-----------------------------------------------------------------------------
          for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

            if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){
                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }
                }
                if(nincident>0){
                       eincident = eincident/(nincident*1.0);
                       _event_histo->hgen_reco_beam_incidentEnergy[is]->Fill(eincident);
                }
            } 
          }//end of loop over all the slices


        //increment interaction counter
        ++interaction_gen[true_sliceID];
        //increment incident counter

        for (int i = 0; i<=sliceID; ++i) ++incident_gen[i];
        if(true_sliceID>=0 && true_sliceID<nslices){
          h_energetic_pmom_gen[true_sliceID]->Fill(temp_genpmom);
          h_energetic_pcostheta_gen[true_sliceID]->Fill(temp_genpcostheta);
          h_energetic_pphi_gen[true_sliceID]->Fill(temp_genpphi);
        }





        //-------------------------------------------------------------------------------
        }//end of if this is a generated signal event
        //================================================================
        //std::cout<<"start to fill the vector with reweight factors"<<t->g4rw_primary_var->size()<<std::endl;
        std::vector<double> wgts_geant_pm1;  //size of wgts_geant_pm1 would be 2*reweight names
        if(isdata==0 && _fill_bootstrap_geant){
        for(long unsigned int i=0; i<t->g4rw_primary_var->size(); i++){
            wgts_geant_pm1.push_back(t->g4rw_primary_plus_sigma_weight->at(i));
            wgts_geant_pm1.push_back(t->g4rw_primary_minus_sigma_weight->at(i));
        }
        }//end of if isdata for filling reweight factors for Geant4

        //std::cout<<"end of filling the vector with reweight factors"<<t->g4rw_primary_var->size()<<std::endl;

        //==================================================================
        if(isdata==0){
        //get the denomator of the energetic protons for momentum, angles
        //TestGenSig is the number of events with at least one proton
        if(t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 && t->true_daughter_nPiMinus ==0 )      {
          int temp_genpindex=-999; double temp_genpmom=-999.0;
          for(unsigned int ntrk=0; ntrk<t->true_beam_daughter_startP->size(); ntrk++){
            if(abs(t->true_beam_daughter_PDG->at(ntrk)) !=2212) continue;
            if(t->true_beam_daughter_startP->at(ntrk)>temp_genpmom){
               temp_genpindex = ntrk;
               temp_genpmom=t->true_beam_daughter_startP->at(ntrk);
            }           
          }
          if(temp_genpindex>=0){
            _event_histo_1d->h_PiAbs_gen_sig_energeticproton_mom->Fill(t->true_beam_daughter_startP->at(temp_genpindex));
            //_event_histo_1d->h_PiAbs_gen_sig_energeticproton_costheta->Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex));
            _event_histo_1d->h_PiAbs_gen_sig_energeticproton_phi->Fill(TMath::ATan2(t->true_beam_daughter_startPy->at(temp_genpindex),t->true_beam_daughter_startPx->at(temp_genpindex)));
           if(_fill_bootstrap_geant){
             bs_geant_pm1_eff_mumom_den.Fill(t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);
             bs_geant_pm1_eff_mucostheta_den.Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);
           }
            TestGenSig++;
          } 
             
        }//end of selected signal before event selection
        
        }//end of if isdata==0
        //std::cout<<"selected the most energetic protons "<<std::endl;
        //==================================================================


        TVector3 truebeamendv3;

        truebeamendv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);


        TVector3 truedaughterstartv3;
        TVector3 truedaughterendv3;
        TVector3 truegranddaughterstartv3;
        TVector3 truegranddaughterendv3;

        //loop over all the grand daughter particles
        //save the parent id of the pions into a vector pargdid
        vector<int> pargdid;
        pargdid.clear();
        if(isdata==0){  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        for(unsigned ingd=0; ingd<t->true_beam_grand_daughter_ID->size(); ingd++){
            if(abs(t->true_beam_grand_daughter_PDG->at(ingd))==111 || abs(t->true_beam_grand_daughter_PDG->at(ingd))==211) {
                 pargdid.push_back(t->true_beam_grand_daughter_parID->at(ingd));
            }
        } 
        //====================================================================
        //----------------------------------------------------------------------------  
        //check how many nucleons 
        //LOG_NORMAL()<<"loop over all the true beam daughter particles"<<std::endl;
        for(unsigned indp=0; indp<t->true_beam_daughter_ID->size(); indp++){


             truedaughterstartv3.SetXYZ(t->true_beam_daughter_startX->at(indp), t->true_beam_daughter_startY->at(indp),t->true_beam_daughter_startZ->at(indp));
             truedaughterendv3.SetXYZ(t->true_beam_daughter_endX->at(indp), t->true_beam_daughter_endY->at(indp),t->true_beam_daughter_endZ->at(indp));
 
             if((truebeamendv3-truedaughterstartv3).Mag()<(truebeamendv3-truedaughterendv3).Mag()) {
               _event_histo_1d->h_daughter_beam_dist->Fill((truebeamendv3-truedaughterstartv3).Mag());
             }
             else {_event_histo_1d->h_daughter_beam_dist->Fill((truebeamendv3-truedaughterendv3).Mag());}



             //tempit: loop over all the parent id of the pions of granddaughter
             std::vector<int>::iterator tempit;

             tempit = std::find(pargdid.begin(), pargdid.end(), t->true_beam_daughter_ID->at(indp));
             if(tempit !=pargdid.end()){
                          if(abs(t->true_beam_daughter_PDG->at(indp))==2212 || abs(t->true_beam_daughter_PDG->at(indp))==2112) {
                           _event_histo_1d->h_gdfromproton->Fill(t->true_beam_daughter_startP->at(indp));
                          }
                          else /*if(abs(t->true_beam_daughter_PDG->at(indp))==)*/ {
                                 _event_histo_1d->h_gdfromother->Fill(t->true_beam_daughter_startP->at(indp));
                                 //std::cout<<"GRANDDAUGHTER PDG CODE is "<<t->true_beam_daughter_PDG->at(indp)<<std::endl;
                          }
      
             } else {continue;}
                          
        } //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }//end of if isdata=0
        //LOG_NORMAL()<<"Start to Fill histograms for chi2 and calculate truncated mean dqdx"<<std::endl;
        //======================================================================
	//loop over all the daughter particles and get the chi2
        for(unsigned int ipfp=0; ipfp<t->reco_daughter_allTrack_ID->size(); ipfp++){ 
           if(t->reco_daughter_PFP_trackScore_collection->at(ipfp)<trkscorecut && t->reco_daughter_PFP_nHits->at(ipfp)>=40) continue;
             _event_histo_1d->h_chi2_phypo_mc->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           if(isdata==0){ 
           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==2212) {
             _event_histo_1d->h_chi2_phypo_proton->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
             if (t->reco_daughter_PFP_true_byHits_parID->at(ipfp) == t->true_beam_ID) {
                  _event_histo_1d->h_nhits_proton->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
             }
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==211) {
              _event_histo_1d->h_nhits_pionpm->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
              _event_histo_1d->h_chi2_phypo_pionpm->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==13) {
              _event_histo_1d->h_nhits_muon->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
              _event_histo_1d->h_chi2_phypo_muon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==11) {
              _event_histo_1d->h_nhits_electron->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
              _event_histo_1d->h_chi2_phypo_electron->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==321) {
              _event_histo_1d->h_nhits_kaon->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
              _event_histo_1d->h_chi2_phypo_kaon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else{
              _event_histo_1d->h_nhits_other->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));

             _event_histo_1d->h_chi2_phypo_other->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           } // end of if isdata ==0
        }//end of loop over all the daughter PDF particles  
        //--------------------------------------------------------------------------------

        //LOG_NORMAL()<<"Loop over all the true beam daughters and fill histograms for start momentum"<<std::endl;
	std::vector<double> trunmeandqdx_test;
	std::vector<double> poop_par;
        //======================================================================
        //loop over all the true beam daughters and fill the histogram for the start momentum
        //of pi0, pipm and protons before performing any cut, which are the denominator 
        //of the track reconstruction efficiency
        if(isdata==0){ 
        for(unsigned int indt=0; indt<t->true_beam_daughter_PDG->size(); indt++){

           if(abs(t->true_beam_daughter_PDG->at(indt)) == 111) { 
              _event_histo_1d->h_mom_gentruepion0->Fill(t->true_beam_daughter_startP->at(indt));
           }
           if(abs(t->true_beam_daughter_PDG->at(indt)) == 211) { 
              _event_histo_1d->h_mom_gentruepionpm->Fill(t->true_beam_daughter_startP->at(indt));
              auto piCSearch = std::find (t->reco_daughter_PFP_true_byHits_ID->begin(), t->reco_daughter_PFP_true_byHits_ID->end(), t->true_beam_daughter_ID->at(indt) );
              if(piCSearch !=t->reco_daughter_PFP_true_byHits_ID->end() ){
                  _event_histo_1d->h_mom_recotruepionpm->Fill(t->true_beam_daughter_startP->at(indt));
              }
           }
           if(abs(t->true_beam_daughter_PDG->at(indt)) !=2112 && abs(t->true_beam_daughter_PDG->at(indt)) !=22 && abs(t->true_beam_daughter_PDG->at(indt)) !=111){
           if(abs(t->true_beam_daughter_PDG->at(indt)) !=2212 && abs(t->true_beam_daughter_PDG->at(indt)) !=211){
              _event_histo_1d->h_mom_gentrueother->Fill(t->true_beam_daughter_startP->at(indt));
             auto otherSearch = std::find(t->reco_daughter_PFP_true_byHits_ID->begin(), t->reco_daughter_PFP_true_byHits_ID->end(), t->true_beam_daughter_ID->at(indt) );
             if(otherSearch !=t->reco_daughter_PFP_true_byHits_ID->end()) {
                _event_histo_1d->h_mom_recotrueother->Fill(t->true_beam_daughter_startP->at(indt));
             } 
           }
           }
           if(abs(t->true_beam_daughter_PDG->at(indt))!= 2212) continue;
	   _event_histo_1d->h_mom_gentruep->Fill(t->true_beam_daughter_startP->at(indt));
           _event_histo_1d->h_thetax_gentruep->Fill(TMath::ACos(t->true_beam_daughter_startPx->at(indt)/t->true_beam_daughter_startP->at(indt)));

           //check if this proton is reaconstructed, then get the initial reconstruction efficiency
           int foundreco=0;
           auto itSearch = std::find( t->reco_daughter_PFP_true_byHits_ID->begin(), t->reco_daughter_PFP_true_byHits_ID->end(), t->true_beam_daughter_ID->at(indt) );
           if( itSearch != t->reco_daughter_PFP_true_byHits_ID->end() ){
                 foundreco++;
                 _event_histo_1d->h_mom_recotruep->Fill(t->true_beam_daughter_startP->at(indt));
           } else { 
            //study the non-reco protons            
            //if(t->true_beam_daughter_startP->at(indt)> 0.3 && t->true_beam_daughter_nHits->at(indt)  > 0 )
            {
                  _event_histo_1d->h_lownhitsp_thetax->Fill(t->true_beam_daughter_startPx->at(indt)/t->true_beam_daughter_startP->at(indt));
                  _event_histo_1d->h_lownhitsp_thetay->Fill(t->true_beam_daughter_startPy->at(indt)/t->true_beam_daughter_startP->at(indt));
                  _event_histo_1d->h_lownhitsp_thetaz->Fill(t->true_beam_daughter_startPz->at(indt)/t->true_beam_daughter_startP->at(indt));

                  _event_histo_1d->h_lownhitsp_startX->Fill(t->true_beam_daughter_startX->at(indt)); 
                  _event_histo_1d->h_lownhitsp_startY->Fill(t->true_beam_daughter_startY->at(indt)); 
                  _event_histo_1d->h_lownhitsp_startZ->Fill(t->true_beam_daughter_startZ->at(indt)); 

                  _event_histo_1d->h_lownhitsp_endX->Fill(t->true_beam_daughter_endX->at(indt)); 
                  _event_histo_1d->h_lownhitsp_endY->Fill(t->true_beam_daughter_endY->at(indt)); 
                  _event_histo_1d->h_lownhitsp_endZ->Fill(t->true_beam_daughter_endZ->at(indt)); 


                  //std::cout<<"libotest "<<t->run<<"   "<<t->subrun<<"   "<<t->event<<"  "
                  //std::cout<<t->true_beam_daughter_endProcess->at(indt)<<"  "
                  //std::cout<<t->true_beam_daughter_startP->at(indt)<<"  "
                  //std::cout<<t->true_beam_daughter_nHits->at(indt)<<std::endl;    
                  _event_histo_1d->h_nonrecop_momvsnhits->Fill(t->true_beam_daughter_startP->at(indt), t->true_beam_daughter_nHits->at(indt));
                  _event_histo_1d->h_mom_nonrecop->Fill(t->true_beam_daughter_startP->at(indt));
            }
           }
        }// end of loop over all the true particles
        }//end of isdata==0
        //======================================================================
        //LOG_NORMAL()<<"Loop over all the reconstructed daughter particles "<<std::endl;
        /*
        * loop over all the reconstructed daughter particles
        * get the photons and fill the parent ID of the photon into 
        * vec_parid vector, which will be checked by the PDG 
        */
        //=========================================================================
        if(isdata==0){
        vector<int> vec_parid; //find out the parID of the photons but no loop
        vec_parid.clear();     //save the photon parent ID into vec_parid

        for(unsigned int mm=0; mm<t->reco_daughter_PFP_true_byHits_PDG->size(); mm++){
             if(t->reco_daughter_PFP_true_byHits_PDG->at(mm) !=22) continue;
              
             std::vector<int>::iterator it;

             it = std::find(vec_parid.begin(), vec_parid.end(), t->reco_daughter_PFP_true_byHits_parID->at(mm));
             if (it != vec_parid.end()) 
             {
               
             } else {
                   vec_parid.push_back(t->reco_daughter_PFP_true_byHits_parID->at(mm));
             }
        }     
        
        vector<int> vec_pionparid;
        vec_pionparid.clear();
        //----------------------------------------------------------------------- 
        double Pgamma = 0.;
        double Egamma = 0.;
        //loop over all the parID from reco photons if there are photons exist in the daughter particles
        if(vec_parid.size()>0){
        for(unsigned int ind_parid=0; ind_parid<vec_parid.size(); ind_parid++){
          // loop over all the true beam daughter ID and check if it is pion 0
          for(unsigned int ind_cd=0; ind_cd<t->true_beam_daughter_ID->size(); ind_cd++){
               if(t->true_beam_daughter_ID->at(ind_cd) == vec_parid.at(ind_parid) && t->true_beam_daughter_PDG->at(ind_cd)==111) {
                        vec_pionparid.push_back(vec_parid.at(ind_parid));
                        _event_histo_1d->h_mom_recotruepion0->Fill(t->true_beam_daughter_startP->at(ind_cd));
	       } else if(t->true_beam_daughter_ID->at(ind_cd) == vec_parid.at(ind_parid)){
                        _event_histo_1d->h_mom_recotruenonpi0->Fill(t->true_beam_daughter_startP->at(ind_cd));
               }
          }
        }
        }
        //---------------------------------------------------------------------------
        // vec pion0 parID produce photons  
        if(vec_pionparid.size()>0){
        for(unsigned int ind_parid=0; ind_parid<vec_pionparid.size(); ind_parid++){
          Pgamma =0.;
          double Pxgamma = 0.0;
          double Pygamma = 0.0; 
          double Pzgamma = 0.0;
          int Ngamma=0;
          for(unsigned int hh=0; hh<t->reco_daughter_PFP_true_byHits_PDG->size(); hh++){
            if(t->reco_daughter_PFP_true_byHits_PDG->at(hh) !=22) continue;
            if(t->reco_daughter_PFP_true_byHits_parID->at(hh) != vec_pionparid.at(ind_parid)) continue;
            Pgamma += t->reco_daughter_PFP_true_byHits_startP->at(hh);
            Egamma += t->reco_daughter_PFP_true_byHits_startE->at(hh);
            Pxgamma += t->reco_daughter_PFP_true_byHits_startPx->at(hh);
            Pygamma += t->reco_daughter_PFP_true_byHits_startPy->at(hh);
            Pzgamma += t->reco_daughter_PFP_true_byHits_startPz->at(hh);
            Ngamma++;
          } // end of all the reco pf objects
          _event_histo_1d->h_ngamma_frompi0->Fill(Ngamma);
          if(Ngamma >=2 && Pgamma>0){
          _event_histo_1d->h_pgamma_frompi0->Fill(TMath::Sqrt(Egamma*Egamma-(Pxgamma*Pxgamma+Pygamma*Pygamma+Pzgamma*Pzgamma)));
                
          }
        } //loop over all the vector of parid
        }
        }//end of if isdata==0
        //======================================================================== 
	
        int nshwcand = 0;

        int ntmdqdxcand_withthresh = 0;

        int nonprotoncand_withthresh = 0;
        TVector3 recobeamendv3, recodauendv3, recodaustartv3, recogranddaustartv3, recogranddauendv3;
        recobeamendv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);







        std::vector<double> daughter_distance3D, daughter_distance3D_shower;
        std::vector<double> daughter_angle3D, daughter_angle3D_shower;

        std::vector<std::vector<double>> *trkdedx_test_ptr=t->reco_daughter_allTrack_calibrated_dEdX_SCE;

        const std::vector<double> reco_daughter_allTrack_truncLibo_dEdX_test = truncatedMean_libo(libo_low, libo_high, *trkdedx_test_ptr);
 


        std::vector<double> *rd_startX_ptr = t->reco_daughter_allTrack_startX;
        std::vector<double> *rd_startY_ptr = t->reco_daughter_allTrack_startY;
        std::vector<double> *rd_startZ_ptr = t->reco_daughter_allTrack_startZ;
   
        std::vector<double> *rd_endX_ptr = t->reco_daughter_allTrack_endX;
        std::vector<double> *rd_endY_ptr = t->reco_daughter_allTrack_endY;
        std::vector<double> *rd_endZ_ptr = t->reco_daughter_allTrack_endZ;
            
        std::vector<double> newrd_startX = *rd_startX_ptr;
        std::vector<double> newrd_startY = *rd_startY_ptr;
        std::vector<double> newrd_startZ = *rd_startZ_ptr;

        std::vector<double> newrd_endX = *rd_endX_ptr;
        std::vector<double> newrd_endY = *rd_endY_ptr;
        std::vector<double> newrd_endZ = *rd_endZ_ptr;

       
            
        daughter_distance3D =   compute_distanceVertex(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ,
                                    newrd_startX, newrd_startY, newrd_startZ,      
                                    newrd_endX, newrd_endY, newrd_endZ);      
        std::vector<double> *rd_theta_ptr = t->reco_daughter_allTrack_Theta;       
        std::vector<double> *rd_phi_ptr = t->reco_daughter_allTrack_Phi;       
  
        std::vector<double> newrd_theta = *rd_theta_ptr;
        std::vector<double> newrd_phi = *rd_phi_ptr;

        daughter_angle3D = compute_angleVertex(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ,
                                    newrd_startX, newrd_startY, newrd_startZ,      
                                    newrd_endX, newrd_endY, newrd_endZ, 
                                    newrd_theta, newrd_phi, newrd_theta, newrd_phi); 




        //calculate the shower reco beam distance
        std::vector<double> *rdshwr_startX_ptr = t->reco_daughter_allShower_startX;
        std::vector<double> *rdshwr_startY_ptr = t->reco_daughter_allShower_startY;
        std::vector<double> *rdshwr_startZ_ptr = t->reco_daughter_allShower_startZ;
   
        /*std::vector<double> *rdshwr_endX_ptr = t->reco_daughter_allShower_endX;
        std::vector<double> *rdshwr_endY_ptr = t->reco_daughter_allShower_endY;
        std::vector<double> *rdshwr_endZ_ptr = t->reco_daughter_allShower_endZ;
        */
        std::vector<double> newrdshwr_startX = *rdshwr_startX_ptr;
        std::vector<double> newrdshwr_startY = *rdshwr_startY_ptr;
        std::vector<double> newrdshwr_startZ = *rdshwr_startZ_ptr;
            
        /*std::vector<double> &newrdshwr_endX = *rdshwr_endX_ptr;
        std::vector<double> &newrdshwr_endY = *rdshwr_endY_ptr;
        std::vector<double> &newrdshwr_endZ = *rdshwr_endZ_ptr;
        */
        std::vector<double> *rdshwr_dirX_ptr = t->reco_daughter_allShower_dirX;
        std::vector<double> *rdshwr_dirY_ptr = t->reco_daughter_allShower_dirY;
        std::vector<double> *rdshwr_dirZ_ptr = t->reco_daughter_allShower_dirZ;

        std::vector<double> newrdshwr_dirX = *rdshwr_dirX_ptr;
        std::vector<double> newrdshwr_dirY = *rdshwr_dirY_ptr;
        std::vector<double> newrdshwr_dirZ = *rdshwr_dirZ_ptr;;

        daughter_distance3D_shower =   compute_distanceVertex(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ,
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ,      
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ);      
         
 
        daughter_angle3D_shower    =   compute_angleVertex_shower(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ,
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ,      
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ,      
                                    newrdshwr_dirX, newrdshwr_dirY, newrdshwr_dirZ); 

       
        
        //LOG_NORMAL()<<"Calculated the shower angle and Shower distance from the vertex"<<std::endl;

        Int_t totalrecogamma = 0;
        //loop over all the daughter particles and get the most energetic shower(score<0.3)
        double temp_energy_shw=-999.0;
        temp_index_shw=-999;
        for(unsigned int gg=0; gg<t->reco_daughter_PFP_trackScore_collection->size(); gg++){
              if(t->reco_daughter_PFP_trackScore_collection->at(gg)>=trkscorecut) continue;
              if(t->reco_daughter_allShower_energy->at(gg)>temp_energy_shw){
                  temp_energy_shw=t->reco_daughter_allShower_energy->at(gg);
                  temp_index_shw=gg;
              }
        }

        
       



        for(unsigned int jj=0; jj<t->reco_daughter_PFP_trackScore_collection->size(); jj++){
            _event_histo_1d->h_trackscore_all->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut){
               totalrecogamma++;
               _event_histo_1d->h_trackscore_shwlike_all->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));

               if(daughter_distance3D_shower[jj]>cut_daughter_shower_distance_low){
                 _event_histo_1d->h_shweng_all->Fill(t->reco_daughter_allShower_energy->at(jj));
                 _event_histo_1d->h_shwang_all->Fill(daughter_angle3D_shower[jj]);
               }
               _event_histo_1d->h_nhits_shwlike_all->Fill(t->reco_daughter_PFP_nHits->at(jj));
               _event_histo_1d->h_shwdis_all->Fill(daughter_distance3D_shower[jj]); 



             

               if(isdata == 0){
                   //-------------------------------------------------------------------------------
                   if(t->reco_daughter_PFP_trackScore_collection->at(jj)!= -999.){
                   totalshower_trkscore++;
                   if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
                       totalgamma_trkscore++;
                   }
                   if(t->reco_daughter_PFP_nHits->at(jj)>cut_energy_shower_low){
                       totalshower_trknhits++;
                       if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22 ||
                          abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                            totalgamma_trknhits++;
                       }
                   }
                   if(daughter_distance3D_shower[jj]>cut_daughter_shower_distance_low){
                       totalshower_trkdist++;
                       if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22 ||
                       abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                           totalgamma_trkdist++;
                       }
                       if(t->reco_daughter_allShower_energy->at(jj)>cut_energy_shower_low){
                           totalshower_trkeng++;
                           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22 ||
                              abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                              totalgamma_trkeng++;
                           }
                       }
                   } 
                   }
                 



                   //------------------------------------------------------------------------------
                   if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
                        _event_histo_1d->h_trackscore_shwlike_photon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        if(daughter_distance3D_shower[jj]>cut_daughter_shower_distance_low){
                           _event_histo_1d->h_shweng_photon->Fill(t->reco_daughter_allShower_energy->at(jj));
                           _event_histo_1d->h_shwang_photon->Fill(daughter_angle3D_shower[jj]);
                        }  
                        _event_histo_1d->h_nhits_shwlike_photon->Fill(t->reco_daughter_PFP_nHits->at(jj));
                        _event_histo_1d->h_shwdis_photon->Fill(daughter_distance3D_shower[jj]);

                   }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){ 
                        _event_histo_1d->h_trackscore_shwlike_electron->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        if(daughter_distance3D_shower[jj]>cut_daughter_shower_distance_low){
                           _event_histo_1d->h_shweng_electron->Fill(t->reco_daughter_allShower_energy->at(jj));
                           _event_histo_1d->h_shwang_electron->Fill(daughter_angle3D_shower[jj]);
                        }
                        _event_histo_1d->h_nhits_shwlike_electron->Fill(t->reco_daughter_PFP_nHits->at(jj));
                        _event_histo_1d->h_shwdis_electron->Fill(daughter_distance3D_shower[jj]);
                   }else {
                        _event_histo_1d->h_trackscore_shwlike_nonphoton->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        if(daughter_distance3D_shower[jj]>cut_daughter_shower_distance_low){
                        _event_histo_1d->h_shweng_nonphoton->Fill(t->reco_daughter_allShower_energy->at(jj));
                        _event_histo_1d->h_shwang_nonphoton->Fill(daughter_angle3D_shower[jj]);
                        }
                        _event_histo_1d->h_nhits_shwlike_nonphoton->Fill(t->reco_daughter_PFP_nHits->at(jj));
                        _event_histo_1d->h_shwdis_nonphoton->Fill(daughter_distance3D_shower[jj]);
                   }
               } 
            }

            //LOG_NORMAL()<<"Filled out track score distributions "<<std::endl;
            if(isdata == 0){
	    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	            _event_histo_1d->h_trackscore_proton->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
                    if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                    _event_histo_1d->h_mom_recotruep_test->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                    _event_histo_1d->h_thetax_recotruep_test->Fill(TMath::ACos(t->reco_daughter_PFP_true_byHits_startPx->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj)));
                    }
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
        	    _event_histo_1d->h_trackscore_pionpm->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));  
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
	            _event_histo_1d->h_trackscore_electron->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));      
                    if(abs(t->reco_daughter_PFP_true_byHits_parPDG->at(jj))==13){
                          _event_histo_1d->h_michele_trackscore->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                    } else if(abs(t->reco_daughter_PFP_true_byHits_parPDG->at(jj))==22){
                          _event_histo_1d->h_gammae_trackscore->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                    } else if(abs(t->reco_daughter_PFP_true_byHits_parPDG->at(jj))==111){
                          _event_histo_1d->h_pione_trackscore->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                    } else{
                          _event_histo_1d->h_othere_trackscore->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                          //std::cout<<"the PDG code of the electron is "<<t->reco_daughter_PFP_true_byHits_parPDG->at(jj)<<std::endl;
                    } 
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
		    _event_histo_1d->h_trackscore_photon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
        	    _event_histo_1d->h_trackscore_muon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==111){
	            _event_histo_1d->h_trackscore_pion0->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else {
        	    _event_histo_1d->h_trackscore_other->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            }
            } //end of if isdata==0
            //---------------------------------------------------------------------

            //LOG_NORMAL()<<"Filled out the hisgrams of track score according to PDG ID"<<std::endl;

            //if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && t->reco_daughter_PFP_nHits->at(jj)>=40) {nshwcand ++; }
            
            //if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && t->reco_daughter_PFP_nHits->at(jj)>=40) continue;

            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut &&
               t->reco_daughter_PFP_trackScore_collection->at(jj) !=-999. &&  
               ((t->reco_daughter_allShower_energy->at(jj)>= cut_energy_shower_low /*&& 
                 t->reco_daughter_allShower_energy->at(jj)<= cut_energy_shower_high*/) ||
                 (daughter_distance3D_shower[jj]>=cut_daughter_shower_distance_low /*&& 
                  daughter_distance3D_shower[jj]<=cut_daughter_shower_distance_high*/)))
               
                {nshwcand ++; 

            }
            
            /*if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && 
               t->reco_daughter_allShower_energy->at(jj)>=cut_energy_shower_low &&
               t->reco_daughter_allShower_energy->at(jj)<=cut_energy_shower_high) continue;
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && 
               daughter_distance3D_shower[jj]>= cut_daughter_shower_distance_low &&
               daughter_distance3D_shower[jj]<= cut_daughter_shower_distance_high) continue;
            */ 
            //LOG_NORMAL()<<"Performing emergy cut shower energy and shower distance cut to the TPC objects before PID"<<std::endl;
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut &&
               t->reco_daughter_PFP_trackScore_collection->at(jj) !=-999. &&  
               ((t->reco_daughter_allShower_energy->at(jj)>= cut_energy_shower_low /*&& 
                 t->reco_daughter_allShower_energy->at(jj)<= cut_energy_shower_high*/) ||
                 (daughter_distance3D_shower[jj]>=cut_daughter_shower_distance_low /*&& 
                  daughter_distance3D_shower[jj]<=cut_daughter_shower_distance_high*/))) continue;
 


            if(isdata==0){//fill histograms for the momentum distributions of the proton/pion true momentum distribution
                          //after track score cut for proton & pion reco effiency calcuation
            if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                      _event_histo_1d->h_mom_trkscoretruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                }
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211) {
                      _event_histo_1d->h_mom_trkscoretruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                }
                   
            }
            }// end of if isdata==0 ----------------------------------------------------------------
            int dedx_vector_size = t->reco_daughter_allTrack_calibrated_dEdX_SCE->at(jj).size();
            double tme=GetTruncatedMean(t->reco_daughter_allTrack_calibrated_dEdX_SCE->at(jj), 2, dedx_vector_size - 5, 0.1, 0.6); 
 
            _event_histo_1d->h_tmdqdx->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));
            _event_histo_1d->h_tmdqdx_new->Fill(tme);
            if(isdata ==0){
	    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	           _event_histo_1d->h_tmdqdx_proton->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
                   _event_histo_1d->h_tmdqdx_proton_new->Fill(tme);
            _event_histo_1d->h_tmdqdxvsrange_proton->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));

	    }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
	           _event_histo_1d->h_tmdqdx_pionpm->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
                   _event_histo_1d->h_tmdqdx_pionpm_new->Fill(tme);
            _event_histo_1d->h_tmdqdxvsrange_pionpm->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));
	    }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
	           _event_histo_1d->h_tmdqdx_photon->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
                   _event_histo_1d->h_tmdqdx_photon_new->Fill(tme);
            _event_histo_1d->h_tmdqdxvsrange_photon->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));
	    }else{
	           _event_histo_1d->h_tmdqdx_other->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
                   _event_histo_1d->h_tmdqdx_other_new->Fill(tme);
            _event_histo_1d->h_tmdqdxvsrange_other->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));
            }		
            } //end of if isdata ==0 

            
            


            //------------------------------------------------------------------
            bool isPCand = false;
            bool isPiCand = false;
            double temp_tmdqdx=reco_daughter_allTrack_truncLibo_dEdX_test.at(jj);

            //==================================================================================
            if(temp_tmdqdx<cut_dEdX_low || (temp_tmdqdx>cut_dEdX_high && temp_tmdqdx<cut_for_proton )){ //fill out chi2 of the transition region
             _event_histo_1d->htmdqdx_chi2_phypo_mc->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            if(isdata==0){
            if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
             _event_histo_1d->htmdqdx_chi2_phypo_proton->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            }
            else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211) {
              _event_histo_1d->htmdqdx_chi2_phypo_pionpm->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            }
            else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13) {
              _event_histo_1d->htmdqdx_chi2_phypo_muon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            }
            else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11) {
              _event_histo_1d->htmdqdx_chi2_phypo_electron->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            }
            else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321) {
              _event_histo_1d->htmdqdx_chi2_phypo_kaon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));}
            else{
             _event_histo_1d->htmdqdx_chi2_phypo_other->Fill(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj));
            }
            } // end of if isdata ==0
            }//end of selecting proton region and transition region

            
            //LOG_NORMAL()<<"Filled out histograms for trunmean dEdx and Chi2"<<std::endl;
            //================================================================================= 
            if(temp_tmdqdx>=cut_dEdX_low && temp_tmdqdx<=cut_dEdX_high 
                  /*&& t->reco_daughter_PFP_trackScore_collection->at(jj)>trkscorecut*/) 
            {ntmdqdxcand_withthresh++;}
            if(temp_tmdqdx>cut_for_proton) {  //proton region 1
               isPCand = true;
               if(isdata==0){
               if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                   _event_histo_1d->h_tmdqdxnhits_proton->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(jj));
                   _event_histo_1d->h_mom_selectedtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                   _event_histo_1d->h_mom_tmdqdxtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
               }
               }//end of if isdata==0
            } 

            if(temp_tmdqdx>cut_dEdX_low && temp_tmdqdx<cut_dEdX_high) {//pion region
               isPiCand = true;
               if(isdata==0){
               if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                     _event_histo_1d->h_mom_selectedtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     _event_histo_1d->h_mom_tmdqdxtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
               }
               }//end of if isdata==0
            }
            if((temp_tmdqdx<cut_for_proton && temp_tmdqdx>cut_dEdX_high) || temp_tmdqdx<cut_dEdX_low) {  //transition region
                if(isdata==0){
                if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                     _event_histo_1d->h_tmdqdxnhits_proton->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(jj));
                     _event_histo_1d->h_mom_tmdqdxtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                     _event_histo_1d->h_mom_tmdqdxtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
                }//selected the primary daughters
                }//end of if isdata==0
                if(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) < 50
                   && t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 0 ) 
                {
                  isPCand = true;
                  if(isdata==0){
                  if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                     if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {              
                        _event_histo_1d->h_mom_selectedtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     }//end of if the tracks with chi2>40
                  }//end of selecting the true daughter particles 
                  }//end of if isdata==0
                } //end of apply chi2 cut
                else if(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 50
                   && t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 0 ) 
                { 
                 isPiCand = true;
                 if(isdata==0){
                  if(t->reco_daughter_PFP_true_byHits_parID->at(jj) ==t->true_beam_ID) {
                     if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                         _event_histo_1d->h_mom_selectedtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     }   
                  }
                 }//end of if isdata==0 
                }
                else {
                    isPiCand = true;
                if(isdata==0) {
                    if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                    if(t->reco_daughter_PFP_true_byHits_startP->at(jj) > 1.0) {
                        _event_histo_1d->h_chi2cutp_thetax->Fill(t->reco_daughter_PFP_true_byHits_startPx->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        _event_histo_1d->h_chi2cutp_thetay->Fill(t->reco_daughter_PFP_true_byHits_startPy->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        _event_histo_1d->h_chi2cutp_thetaz->Fill(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        
                    }}} 
                } //end of if isdata==0
                }
            } //end of performing tmdqdx cut
            //LOG_NORMAL()<<"End of Treat Different Regions of the Truncated Mean dQdx"<<std::endl;
            if(isPCand == false /*&& t->reco_daughter_allTrack_momByRange_alt_proton->at(jj) > momthreshcut*/){
              if(t->reco_daughter_allTrack_Chi2_ndof->at(jj)>10){
                nonprotoncand_withthresh ++;
              }
            }
            //LOG_NORMAL()<<"Count the number of protons with the number of hits greater than 10"<<std::endl;
            if(isPCand == true){
                //-------------------------------------------------------------------------------
                if(isdata==0){
                //std::cout<<"size of dEdx "<<t->reco_daughter_allTrack_calibrated_dEdX_SCE->size()<<std::endl;
                //std::cout<<"size of resRange "<<t->reco_daughter_allTrack_resRange_SCE->size()<<std::endl;
                if(t->reco_daughter_allTrack_calibrated_dEdX_SCE->at(jj).size()>0){ 
                 for(unsigned int ndr=0; ndr<t->reco_daughter_allTrack_calibrated_dEdX_SCE->at(jj).size(); ndr++){
                     _event_histo_1d->h_selproton_dedxRR->Fill(t->reco_daughter_allTrack_resRange_SCE->at(jj).at(ndr),t->reco_daughter_allTrack_calibrated_dEdX_SCE->at(jj).at(ndr));
             
                 }
                }
                //---------------------------------------------------------------------------------
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                    if (t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID) {
                      _event_histo_1d->h_reconhits_proton->Fill(t->reco_daughter_allTrack_Chi2_ndof->at(jj));
                    }
                    _event_histo_1d->h_recomom_selected_proton->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                    _event_histo_1d->h_recomom_selected_pionpm->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                    _event_histo_1d->h_recomom_selected_electron->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
                    _event_histo_1d->h_recomom_selected_muon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321){
                    _event_histo_1d->h_recomom_selected_kaon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else{
                    _event_histo_1d->h_recomom_selected_other->Fill(t->reco_daughter_allTrack_len->at(jj));
                }
                //----------------------------------------------------------------------------
                //std::cout<<"check if this proton is a daughter particles"<<std::endl;

                std::vector<int>::iterator temprecoit;

                temprecoit=std::find(t->true_beam_daughter_ID->begin(), t->true_beam_daughter_ID->end(), t->reco_daughter_PFP_true_byHits_parID->at(jj));
                if(temprecoit != t->true_beam_daughter_ID->end()){
                  recodaustartv3.SetXYZ(t->reco_daughter_allTrack_startX->at(jj),t->reco_daughter_allTrack_startY->at(jj), t->reco_daughter_allTrack_startZ->at(jj));
                  recodauendv3.SetXYZ(t->reco_daughter_allTrack_endX->at(jj),t->reco_daughter_allTrack_endY->at(jj), t->reco_daughter_allTrack_endZ->at(jj));


                  if((recobeamendv3-recodaustartv3).Mag()<(recobeamendv3-recodauendv3).Mag()){
                  _event_histo_1d->h_reco_grand_daughter_beam_dist->Fill((recobeamendv3-recodaustartv3).Mag());
                  }
                  else {_event_histo_1d->h_reco_grand_daughter_beam_dist->Fill((recobeamendv3-recodauendv3).Mag());
                  }

                } else {
                  recogranddaustartv3.SetXYZ(t->reco_daughter_allTrack_startX->at(jj),t->reco_daughter_allTrack_startY->at(jj), t->reco_daughter_allTrack_startZ->at(jj));
                  recogranddauendv3.SetXYZ(t->reco_daughter_allTrack_endX->at(jj),t->reco_daughter_allTrack_endY->at(jj), t->reco_daughter_allTrack_endZ->at(jj));


                  if((recobeamendv3-recogranddaustartv3).Mag()<(recobeamendv3-recogranddauendv3).Mag()){
                   _event_histo_1d->h_reco_daughter_beam_dist->Fill((recobeamendv3-recogranddaustartv3).Mag());

                   _event_histo_1d->h_reco_daughter_deltax->Fill(t->reco_daughter_allTrack_startX->at(jj) - t->reco_beam_endX);
                   _event_histo_1d->h_reco_daughter_deltay->Fill(t->reco_daughter_allTrack_startY->at(jj) - t->reco_beam_endY);
                   _event_histo_1d->h_reco_daughter_deltaz->Fill(t->reco_daughter_allTrack_startZ->at(jj) - t->reco_beam_endZ);




                  }
                  else { _event_histo_1d->h_reco_daughter_beam_dist->Fill((recobeamendv3-recogranddauendv3).Mag());

                  }
                } //     end of this proton id found in the daughter particle ID
                //std::cout<<"---fill histograms for momentum track range, angular resolution of the selected proton"<<std::endl;
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                  _event_histo_1d->h_selproton_momreso->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj),t->reco_daughter_allTrack_momByRange_alt_proton->at(jj));
                  //_event_histo_1d->h_trklen_reso->Fill((t->reco_daughter_allTrack_alt_len->at(jj)-t->reco_daughter_PFP_true_byHits_len->at(jj))/t->reco_daughter_PFP_true_byHits_len->at(jj));
                  _event_histo_1d->h_trklen_reso->Fill((t->reco_daughter_allTrack_len->at(jj)-t->reco_daughter_PFP_true_byHits_len->at(jj))/t->reco_daughter_PFP_true_byHits_len->at(jj));
                  _event_histo_1d->h_selproton_costhetareso->Fill(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj), TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj)));
                  _event_histo_1d->h_totalp_mom->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  if(t->reco_daughter_PFP_true_byHits_endProcess->at(jj) !="protonInelastic") {
                  _event_histo_1d->h_selproton_momreso_new->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj)-1);
                  double recoPx=t->reco_daughter_allTrack_momByRange_alt_proton->at(jj)*TMath::Cos(t->reco_daughter_allTrack_Phi->at(jj))*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj));
                  double recoPy=t->reco_daughter_allTrack_momByRange_alt_proton->at(jj)*TMath::Sin(t->reco_daughter_allTrack_Phi->at(jj))*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj));
                  double recoPz=t->reco_daughter_allTrack_momByRange_alt_proton->at(jj)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj));
                  TVector3 temprecoP, tempbeamP;
                  temprecoP.SetXYZ(recoPx, recoPy, recoPz);
                  tempbeamP.SetXYZ(t->reco_beam_trackEndDirX, t->reco_beam_trackEndDirY, t->reco_beam_trackEndDirZ);
                  TVector3 direction=tempbeamP.Unit();
                  temprecoP.RotateUz(direction);

                  double truePx=t->reco_daughter_PFP_true_byHits_startPx->at(jj);
                  double truePy=t->reco_daughter_PFP_true_byHits_startPy->at(jj);
                  double truePz=t->reco_daughter_PFP_true_byHits_startPz->at(jj);


                  _event_histo_1d->h_selproton_momPxreso_new->Fill(recoPx/truePx-1);
                  _event_histo_1d->h_selproton_momPyreso_new->Fill(recoPy/truePy-1);
                  _event_histo_1d->h_selproton_momPzreso_new->Fill(recoPz/truePz-1);
                 
                  _event_histo_1d->h_selproton_costhetareso_new->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj))/(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj))-1);




                  }
                  if(t->reco_daughter_PFP_true_byHits_endProcess->at(jj)=="protonInelastic"){
                    _event_histo_1d->h_reintp_mom->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }                
                  
                }
                //std::cout<<"end of if this is a proton "<<std::endl;
                }//end of if isdata==0 
                //-----------------------------------------------------------------------------
            } //end of if pcand is true
            //LOG_NORMAL()<<"checked if this is a pion candidate or proton candidate"<<std::endl;
            if(isdata==0 && !isPCand && isPiCand){


                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                    _event_histo_1d->h_recopimom_selected_proton->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                    _event_histo_1d->h_recopimom_selected_pionpm->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                    _event_histo_1d->h_recopimom_selected_electron->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
                    _event_histo_1d->h_recopimom_selected_muon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321){
                    _event_histo_1d->h_recopimom_selected_kaon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else{
                    _event_histo_1d->h_recopimom_selected_other->Fill(t->reco_daughter_allTrack_len->at(jj));
                }

            }


            //--------------------------------------------------------------------
	}// end of loop over all the PFP_true_
        //LOG_NORMAL()<<"Start doing event selection for all the pion events"<<std::endl; 
        Ntotal_beam++;  //total number of events past beam selection

        if(isSignal_withthresh){         _event_histo_1d->hsig_nrecogamma->Fill(totalrecogamma); }
        if(isChxBKG_withthresh){         _event_histo_1d->hchxbac_nrecogamma->Fill(totalrecogamma);}
        if(isReaBKG_withthresh){         _event_histo_1d->hreabac_nrecogamma->Fill(totalrecogamma);}
        if(isOtherBKG_withthresh){          _event_histo_1d->hother_nrecogamma->Fill(totalrecogamma);}
 

        std::vector<int> *daughter_nhits_ptr=t->reco_daughter_PFP_nHits;
        std::vector<int> nhitsref=*daughter_nhits_ptr;

        std::vector<double> &daughter_shwdis_ptr=daughter_distance3D_shower;
        //std::vector<double> shwdisref=&daughter_shwdis_ptr;

        std::vector<double> &daughter_shwang_ptr=daughter_angle3D_shower;
        //std::vector<double> shwangref=&daughter_shwang_ptr;

        std::vector<double> *daughter_trkscore_ptr=t->reco_daughter_PFP_trackScore_collection;
        std::vector<double> trkscoreref=*daughter_trkscore_ptr;

        std::vector<double> *daughter_shweng_ptr=t->reco_daughter_allShower_energy;
        std::vector<double> shwengref= *daughter_shweng_ptr;


        std::vector<int> *daughter_ID_ptr=t->reco_daughter_allTrack_ID;
        std::vector<int> trkIDref=*daughter_ID_ptr;

        std::vector<std::vector<double>> *trkdedx_ptr=t->reco_daughter_allTrack_calibrated_dEdX_SCE;
        std::vector<std::vector<double>> trkdedxref=*trkdedx_ptr;


        const std::vector<double> reco_daughter_allTrack_truncLibo_dEdX = truncatedMean_libo(libo_low, libo_high, *trkdedx_ptr);
 

        //if(!secondary_noPion_libo(*daughter_trkscore_ptr,
        //     *daughter_ID_ptr,
        //     reco_daughter_allTrack_truncLibo_dEdX)) continue;



        if(ntmdqdxcand_withthresh> 0) continue;

        Ntotal_tmdqdx++;
        if(isSignal_withthresh) {Ntmcutabs_withthresh++;
               _event_histo_1d->hsig_tmdqdxcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isChxBKG_withthresh) {Ntmcutchx_withthresh++;
               _event_histo_1d->hchxbac_tmdqdxcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }  
        if(isReaBKG_withthresh) {Ntmcutrea_withthresh++;
               _event_histo_1d->hreabac_tmdqdxcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isOtherBKG_withthresh) {Ntmcutother_withthresh++;
        }

        _event_histo_1d->h_tmdqdxcut_mult->Fill(t->reco_daughter_allTrack_ID->size());

        //If there are any more cuts, add them here
        if(nonprotoncand_withthresh> 0) continue; 
        Ntotal_chi2++;

        if(isSignal_withthresh) {Nchi2abs_withthresh++;
               _event_histo_1d->hsig_chi2cut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isChxBKG_withthresh) {Nchi2chx_withthresh++;
               _event_histo_1d->hchxbac_chi2cut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }  
        if(isReaBKG_withthresh) {Nchi2rea_withthresh++;
               _event_histo_1d->hreabac_chi2cut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isOtherBKG_withthresh) {Nchi2other_withthresh++;}
        _event_histo_1d->h_chi2cut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        


        
        /*
        */
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if(temp_index_shw>=0){
             _event_histo_1d->h_energetic_shower_eng_all->Fill(t->reco_daughter_allShower_energy->at(temp_index_shw));
             _event_histo_1d->h_energetic_shower_ang_all->Fill(daughter_angle3D_shower[temp_index_shw]);
             _event_histo_1d->h_energetic_shower_dist_all->Fill(daughter_distance3D_shower[temp_index_shw]);

             if(isSignal_withthresh){ 
                 _event_histo_1d->h_energetic_shower_eng_abs->Fill(t->reco_daughter_allShower_energy->at(temp_index_shw));
                 _event_histo_1d->h_energetic_shower_ang_abs->Fill(daughter_angle3D_shower[temp_index_shw]);
                 _event_histo_1d->h_energetic_shower_dist_abs->Fill(daughter_distance3D_shower[temp_index_shw]);
             } 
             if(isChxBKG_withthresh){
                 _event_histo_1d->h_energetic_shower_eng_chx->Fill(t->reco_daughter_allShower_energy->at(temp_index_shw));
                 _event_histo_1d->h_energetic_shower_ang_chx->Fill(daughter_angle3D_shower[temp_index_shw]);
                 _event_histo_1d->h_energetic_shower_dist_chx->Fill(daughter_distance3D_shower[temp_index_shw]);

             }
             if(isReaBKG_withthresh){
                 _event_histo_1d->h_energetic_shower_eng_rea->Fill(t->reco_daughter_allShower_energy->at(temp_index_shw));
                 _event_histo_1d->h_energetic_shower_ang_rea->Fill(daughter_angle3D_shower[temp_index_shw]);
                 _event_histo_1d->h_energetic_shower_dist_rea->Fill(daughter_distance3D_shower[temp_index_shw]);

             }
        }

        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        if(sel_abs) if(has_shower_Dist(trkscoreref, daughter_shwdis_ptr)) continue;
        if(sel_abs) if(has_shower_Eng(trkscoreref, shwengref)) continue;
        //if(sel_abs) if(has_shower_Ang( trkscoreref, daughter_shwang_ptr)) continue;
        if(sel_abs) if(temp_index_shw>=0 && t->reco_daughter_allShower_energy->at(temp_index_shw)>100) continue;
        if(sel_abs) if(temp_index_shw>=0 && daughter_distance3D_shower[temp_index_shw]>7) continue;


        //if(sel_chx) if(trkscoreref.size()==0) continue;
        if(sel_chx) if(!has_shower_Dist(trkscoreref, daughter_shwdis_ptr)) continue;
        if(sel_chx) if(!has_shower_Eng(trkscoreref, shwengref)) continue;
        
        //if(sel_chx) if(temp_index_shw>=0 && t->reco_daughter_allShower_energy->at(temp_index_shw)<=100) continue;
        //if(sel_chx) if(temp_index_shw>=0 && daughter_distance3D_shower[temp_index_shw]<=6) continue;
        //if(sel_chx) if(temp_index_shw<0) continue;
        
        /*if(sel_abs) if(has_shower_nHits_distance(trkscoreref, nhitsref, daughter_distance3D_shower)) continue;
        if(sel_chx) if(trkscoreref.size()==0) continue;
        if(sel_chx) if(!has_shower_nHits_distance(trkscoreref, nhitsref, daughter_distance3D_shower)) continue;
        */
        /*if(sel_abs) if(has_shower_energy_distance(trkscoreref, shwengref, daughter_shwdis_ptr)) continue;
        if(sel_chx) if(trkscoreref.size()==0) continue;
        if(sel_chx) if(!has_shower_energy_distance(trkscoreref, shwengref, daughter_shwdis_ptr)) continue;
        */  
 

 
        Ntotal_shwid++;

        _event_histo_1d->h_shwcut_mult->Fill(t->reco_daughter_allTrack_ID->size());



        if(isSignal_withthresh) {Nshwcutabs_withthresh++;
        _event_histo_1d->hsig_shwcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isChxBKG_withthresh) {Nshwcutchx_withthresh++;
        _event_histo_1d->hchxbac_shwcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }  
        if(isReaBKG_withthresh) {Nshwcutrea_withthresh++;
        _event_histo_1d->hreabac_shwcut_mult->Fill(t->reco_daughter_allTrack_ID->size());
        }
        if(isOtherBKG_withthresh) {Nshwcutother_withthresh++;
        }
 
        

        double true_Ptmissing = 0.0; double Emissing_Qsubtracted=0.0;
        double reco_Ptmissing = 0.0; double reco_Pz = 0.0; double Emissing_Qsubtracted_reco=0.0;
        double true_Pmissing = 0.0; double reco_Pmissing=0.0;
        double Eexcit = 0.0304;
        double total_truepPx=0.; double total_truepPy=0.; double total_truepPz=0.; 
        double total_pPx = 0.; double total_pPy = 0.; double total_pPz=0.;
        double total_truepKE=0.0;
        double total_recopKE=0.0; 



        //loop over all the tracks , calculate total P,KE
        //and select the most energetic proton candidates
        double temp_pmom=-999.0; double temp_pcostheta=-999.0; double temp_pphi=-999.0; double temp_pcosthetax=-999.;
        double temp_plen = -999.0; int temp_pindex = -999;

        Int_t totalgd=0;

        for(unsigned int ntrk=0; ntrk<t->reco_daughter_allTrack_len->size(); ntrk++){
             _event_histo_1d->h_PiAbs_sel_pmom->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk));
             //rotate vectors to get the angle with respect to the beam direction


             _event_histo_1d->h_PiAbs_sel_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk)));
             _event_histo_1d->h_PiAbs_sel_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));
             
             _event_histo_1d->h_PiAbs_sel_pcosthetax->Fill(TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk))));
             _event_histo_1d->h_PiAbs_sel_pphi->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));
             if(isdata==0){ 
               std::vector<int>::iterator selrecoit;
               selrecoit=std::find(t->true_beam_daughter_ID->begin(), t->true_beam_daughter_ID->end(), t->reco_daughter_PFP_true_byHits_parID->at(ntrk));
 
   
               if(t->reco_daughter_PFP_true_byHits_origin->at(ntrk)==2){
                   _event_histo_1d->h_seltrk_ptheta_cosmic->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));
                   _event_histo_1d->h_seltrk_pphi_cosmic->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));
               } //select the beam daughter particles
               else if (selrecoit != t->true_beam_daughter_ID->end()){
                     totalgd++;
                     _event_histo_1d->h_seltrk_ptheta_beam_granddaughter->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));
                     _event_histo_1d->h_seltrk_pphi_beam_granddaughter->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));
                     _event_histo_1d->h_seltrk_distvtx_granddaughter->Fill(daughter_distance3D.at(ntrk));
                     _event_histo_1d->h_seltrk_angle_granddaughter->Fill(daughter_angle3D.at(ntrk));
                     //std::cout<<"Event Number is "<<t->event<<"  run number is "<<t->run<<" Particle PDG = "<<t->reco_daughter_PFP_true_byHits_PDG->at(ntrk)<<" Parent ID= "<<t->reco_daughter_PFP_true_byHits_parID->at(ntrk)<<" Parent PDG = "<<t->reco_daughter_PFP_true_byHits_parPDG->at(ntrk)<<std::endl;
               } else { //select the beam grand daughter particles
                     _event_histo_1d->h_seltrk_ptheta_beam_daughter->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));
                     _event_histo_1d->h_seltrk_pphi_beam_daughter->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));
                     _event_histo_1d->h_seltrk_distvtx_daughter->Fill(daughter_distance3D.at(ntrk));
                     _event_histo_1d->h_seltrk_angle_daughter->Fill(daughter_angle3D.at(ntrk));
               }
               total_truepPx += t->reco_daughter_PFP_true_byHits_startPx->at(ntrk);
               total_truepPy += t->reco_daughter_PFP_true_byHits_startPy->at(ntrk);
               total_truepKE += t->reco_daughter_PFP_true_byHits_startE->at(ntrk) - ProtonMass;
             }
             total_recopKE += TMath::Sqrt(t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk)*t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk) + ProtonMass*ProtonMass) - ProtonMass;
             

             total_pPx += t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPy += t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPz += t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));


             
             if(t->reco_daughter_allTrack_len->at(ntrk) > temp_plen){
                                temp_pindex = ntrk;
                                temp_plen=t->reco_daughter_allTrack_len->at(ntrk);
                                temp_pmom = t->reco_daughter_allTrack_momByRange_alt_proton->at(ntrk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(ntrk);
                                temp_pcosthetax=TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk)));
             }
            

        } //end of loop over all the tracks and get the information of most energetic 
        


        //==========================================================================================
        _event_histo_1d->h_sel_gdvsd->Fill(t->reco_daughter_allTrack_len->size(),totalgd);

        if(isdata==0 &&(isSignal_withthresh || isChxBKG_withthresh || isReaBKG_withthresh) && temp_pindex>-999.){
        if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "total", fname_geant_pm1, wgts_geant_pm1);
        if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "total", fname_geant_pm1, wgts_geant_pm1);
        } 

        true_Ptmissing = TMath::Sqrt(total_truepPx*total_truepPx + total_truepPy*total_truepPy);

        reco_Ptmissing = TMath::Sqrt(total_pPx*total_pPx + total_pPy*total_pPy);
        reco_Pz = total_pPz; 

        if(isdata==0){
           Emissing_Qsubtracted = TMath::Sqrt(t->true_beam_startP*t->true_beam_startP + PionMass*PionMass) - total_truepKE + Eexcit;
           true_Pmissing = TMath::Sqrt((t->true_beam_endPx - total_truepPx)*(t->true_beam_endPx - total_truepPx) +
                                       (t->true_beam_endPy - total_truepPy)*(t->true_beam_endPy - total_truepPy) +
                                       (t->true_beam_endPz - total_truepPz)*(t->true_beam_endPz - total_truepPz));
        }

        Emissing_Qsubtracted_reco = t->reco_beam_interactingEnergy/1000. - total_recopKE + Eexcit; 
        //LOG_NORMAL()<<"Reco Energy of Beam is "<<t->reco_beam_interactingEnergy<<std::endl;
        double reco_beam_px = 0.0; double reco_beam_py = 0.0; double reco_beam_pz = 0.0;

        

        reco_beam_px = TMath::Sqrt(t->reco_beam_interactingEnergy/1000.*t->reco_beam_interactingEnergy/1000.-PionMass*PionMass)*t->reco_beam_trackEndDirX; 
        reco_beam_py = TMath::Sqrt(t->reco_beam_interactingEnergy/1000.*t->reco_beam_interactingEnergy/1000.-PionMass*PionMass)*t->reco_beam_trackEndDirY;
        reco_beam_pz = TMath::Sqrt(t->reco_beam_interactingEnergy/1000.*t->reco_beam_interactingEnergy/1000.-PionMass*PionMass)*t->reco_beam_trackEndDirZ;

        reco_Pmissing = TMath::Sqrt((reco_beam_px - total_pPx)*(reco_beam_px - total_pPx) + 
                                    (reco_beam_py - total_pPy)*(reco_beam_py - total_pPy) + 
                                    (reco_beam_pz - total_pPz)*(reco_beam_pz - total_pPz));
        if(t->reco_daughter_allTrack_ID->size()>0){
           _event_histo_1d->h_sel_Pmissing->Fill(reco_Pmissing);
           _event_histo_1d->h_sel_Emissing->Fill(Emissing_Qsubtracted_reco);
           _event_histo_1d->h_sel_Ptmissing->Fill(reco_Ptmissing);
           _event_histo_1d->h_sel_Plongit->Fill(reco_Pz);
        }

        //std::cout<<"end of loop over all the tracks and calculated total px py pz "<<std::endl;;
        double p1st_px=0.;
        double p2nd_px=0.;
        double p1st_py=0.;
        double p2nd_py=0.;
        double p1st_pz=0.;
        double p2nd_pz=0.;
        double E1st=0.;
        double E2nd=0.;
        double pcostheta1st2nd=0.;
        if(t->reco_daughter_allTrack_len->size()==2){
                  if(t->reco_daughter_allTrack_momByRange_alt_proton->at(0) > 
                               t->reco_daughter_allTrack_momByRange_alt_proton->at(1)){
                               p1st_px=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(0))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(0));
                               p2nd_px=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(1))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(1));
                               p1st_py=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(0))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(0));
                               p2nd_py=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(1))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(1));
                               p1st_pz=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(0));
                               p2nd_pz=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(1));
                               E1st = TMath::Sqrt(ProtonMass*ProtonMass+t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*t->reco_daughter_allTrack_momByRange_alt_proton->at(0));
                               E2nd = TMath::Sqrt(ProtonMass*ProtonMass+t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*t->reco_daughter_allTrack_momByRange_alt_proton->at(1));
  
 
                  }else{
                               p2nd_px=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(0))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(0));
                               p1st_px=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(1))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(1));
                               p2nd_py=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(0))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(0));
                               p1st_py=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(1))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(1));
                               p2nd_pz=t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(0));
                               p1st_pz=t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(1));
                               E2nd = TMath::Sqrt(ProtonMass*ProtonMass+t->reco_daughter_allTrack_momByRange_alt_proton->at(0)*t->reco_daughter_allTrack_momByRange_alt_proton->at(0));
                               E1st = TMath::Sqrt(ProtonMass*ProtonMass+t->reco_daughter_allTrack_momByRange_alt_proton->at(1)*t->reco_daughter_allTrack_momByRange_alt_proton->at(1));
 


                  }
                  TLorentzVector p1st;
                  p1st.SetPxPyPzE(p1st_px, p1st_py, p1st_pz, E1st);
                  TLorentzVector p2nd;
                  p2nd.SetPxPyPzE(p2nd_px, p2nd_py, p2nd_pz, E2nd);

                  TVector3 p1stv3;
                  TVector3 p2ndv3;
                  p1stv3.SetXYZ(p1st_px, p1st_py, p1st_pz);
                  p2ndv3.SetXYZ(p2nd_px, p2nd_py, p2nd_pz);
   
                  //boost p1st and p2nd to the CM frame
                  TVector3 bv;
                  bv.SetXYZ(p1st_px+p2nd_px, p1st_py+p2nd_py, p1st_pz+p2nd_pz);
                  p1st.Boost(bv);
                  pcostheta1st2nd = TMath::Cos(p1stv3.Angle(p2ndv3));
                  h_sel_pcostheta1st2nd->Fill(pcostheta1st2nd);
                  if(isdata==0 && isSignal_withthresh) {h_sig_pcostheta1st2nd->Fill(pcostheta1st2nd); }
                  if(isdata==0 && isChxBKG_withthresh) {h_chxbac_pcostheta1st2nd->Fill(pcostheta1st2nd); }
                  if(isdata==0 && isReaBKG_withthresh) {h_reabac_pcostheta1st2nd->Fill(pcostheta1st2nd); }
 
                  h_pcostheta1st2nd_ptmissing->Fill(reco_Ptmissing, pcostheta1st2nd);

        } //end of calculating kinematic variables for the 2p events
        for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

            if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){
                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }
                }
                if(nincident>0){
                       eincident = eincident/(nincident*1.0);
                       _event_histo->hsel_reco_beam_incidentEnergy[is]->Fill(eincident);
                }
            } 
        }

        /*int*/ sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        /*double*/// reco_endz = Sce_Corrected_endZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
        //reco_endz = t->reco_beam_endZ;
        //sliceID = int((reco_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); 

        //if (sliceID>=nslices) continue;  //ignore the events with the slice ID >  nslices (24);

        //increment interaction counter
        ++selected_tot[sliceID];
 


        //-----------------------------------------------------------------------------------


        if(isdata==0 && isSignal_withthresh) {
                 _event_histo_1d->h_true_beam_endE_num->Fill(TMath::Sqrt(t->true_beam_endP*t->true_beam_endP + PionMass*PionMass));
 
                 //fill histograms for multiplicities of different 
                 _event_histo_1d->h_PiAbs_sig_pmult->Fill(t->reco_daughter_allTrack_ID->size());
                 _event_histo_1d->h_PiAbs_sig_nmult->Fill(t->true_daughter_nNeutron);
                 _event_histo_1d->h_PiAbs_sig_truevsreco_pmult->Fill(t->true_daughter_nProton, t->reco_daughter_allTrack_ID->size());
                 _event_histo_1d->h_sig_pvsnmult->Fill(t->true_daughter_nProton, t->true_daughter_nNeutron);


                 if( t->reco_daughter_allTrack_ID->size() == 3){
                       //std::cout<<"checking events with 3 reco protons "<<t->true_daughter_nProton<<"     "<<t->true_daughter_nNeutron<<std::endl;
                 }

                 //true to get the true particle multiplicity
                 if(t->reco_daughter_allTrack_ID->size()>0){

                       _event_histo_1d->h_sig_Pmissing->Fill(reco_Pmissing);
                       _event_histo_1d->h_sig_Emissing->Fill(Emissing_Qsubtracted_reco);
                       _event_histo_1d->h_sig_Ptmissing->Fill(reco_Ptmissing);
                       _event_histo_1d->h_sig_Plongit->Fill(reco_Pz);
                       _event_histo_1d->h_sig_Ptmissing_vs_pcand->Fill(reco_Ptmissing, t->reco_daughter_allTrack_ID->size());
                       
                 }



                 //loop over all the selected particles and find out the most energetic proton
                 temp_plen = -999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0; temp_pindex = -999;

                 for(unsigned int tmk=0; tmk<t->reco_daughter_allTrack_ID->size(); tmk++){
                      _event_histo_1d->h_PiAbs_sig_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                      _event_histo_1d->h_PiAbs_sig_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                      _event_histo_1d->h_PiAbs_sig_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                      _event_histo_1d->h_PiAbs_sig_pmom->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk)); 


                      if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_pindex=tmk;
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                       }

                 }//loop over all the proton candidates




                 if(temp_pindex>-999){ //fill histograms with at least one proton
                   if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "signal", fname_geant_pm1, wgts_geant_pm1);
                   if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "signal", fname_geant_pm1, wgts_geant_pm1);
                       

                 }


                 int temp_genpindex=-999.0; double temp_genpmom=-999.0; double temp_genpcostheta=-999.0; double temp_genpphi=-999.0;
                 for(unsigned int tmk=0; tmk<t->true_beam_daughter_startP->size(); tmk++){
                      if(abs(t->true_beam_daughter_PDG->at(tmk)) !=2212) continue;
                      if(t->true_beam_daughter_startP->at(tmk)>temp_genpmom){
                      temp_genpindex = tmk;
                      temp_genpmom=t->true_beam_daughter_startP->at(tmk);
                      temp_genpcostheta=t->true_beam_daughter_startPz->at(tmk)/t->true_beam_daughter_startP->at(tmk);
                      temp_genpphi= TMath::ATan2(t->true_beam_daughter_startPy->at(tmk), t->true_beam_daughter_startPx->at(tmk));
                      }           
                 }
 
                 //std::cout<<"got the most energetic information of pindex costheta and mom"<<std::endl;
                 if(temp_genpindex> -999){ 
                    if(!isdata && _fill_bootstrap_geant) bs_geant_pm1_eff_mumom_num.Fill(temp_genpmom, event_weight, wgts_geant_pm1);        
                    if(!isdata && _fill_bootstrap_geant) bs_geant_pm1_eff_mucostheta_num.Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);        


                    if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_genpmom, temp_pmom, event_weight, bs_geant_pm1_true_reco_mom, fname_geant_pm1, wgts_geant_pm1);
                    if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_genpcostheta, temp_pcostheta, event_weight, bs_geant_pm1_true_reco_costheta, fname_geant_pm1, wgts_geant_pm1);
                    if(!isdata ) {
                       _event_histo->h_eff_pmom_num->Fill(temp_genpmom);
                       _event_histo->h_eff_pcostheta_num->Fill(temp_genpcostheta);
                       _event_histo->h_eff_pphi_num->Fill(temp_genpphi);
                    } 
                    h_true_reco_mom->Fill(temp_genpmom, temp_pmom, event_weight);
                    h_true_reco_costheta->Fill(temp_genpcostheta, temp_pcostheta, event_weight);
                    h_true_reco_phi->Fill(temp_genpphi, temp_pphi, event_weight);

                 }

                 double temp_pmom2=-999.0; 
                 int temp_selpindex=-999;
                 
                 double totalPxtest=0.0;
                 double totalPytest=0.0;
                 double totalPxtest_onlyproton=0.0;
                 double totalPytest_onlyproton=0.0;
              
                 double totalKE_proton=0.;
                 double totalKE_neutron=0.;
                 for(unsigned int tmks=0; tmks<t->true_beam_daughter_startP->size(); tmks++){
                      totalPxtest += t->true_beam_daughter_startPx->at(tmks);
                      totalPytest += t->true_beam_daughter_startPy->at(tmks);
                      if(abs(t->true_beam_daughter_PDG->at(tmks) ) ==2212) {

                          totalPxtest_onlyproton += t->true_beam_daughter_startPx->at(tmks);
                          totalPytest_onlyproton += t->true_beam_daughter_startPy->at(tmks);
 
                          totalKE_proton += TMath::Sqrt(t->true_beam_daughter_startP->at(tmks)*t->true_beam_daughter_startP->at(tmks)+ProtonMass*ProtonMass) - ProtonMass;
                          _event_histo_1d->h_PiAbs_sig_proton_mom->Fill(t->true_beam_daughter_startP->at(tmks));
                          _event_histo_1d->h_PiAbs_sig_proton_costheta->Fill(t->true_beam_daughter_startPz->at(tmks)/t->true_beam_daughter_startP->at(tmks));
                          
                      }

                      if(abs(t->true_beam_daughter_PDG->at(tmks) ) ==2112) {
                          totalKE_neutron +=TMath::Sqrt(t->true_beam_daughter_startP->at(tmks)*t->true_beam_daughter_startP->at(tmks)+NeutronMass*NeutronMass) - NeutronMass;
                          _event_histo_1d->h_PiAbs_sig_neutron_mom->Fill(t->true_beam_daughter_startP->at(tmks));
                          _event_histo_1d->h_PiAbs_sig_neutron_costheta->Fill(t->true_beam_daughter_startPz->at(tmks)/t->true_beam_daughter_startP->at(tmks));
                      }

                      if(abs(t->true_beam_daughter_PDG->at(tmks)) !=2212) continue;
                      if(t->true_beam_daughter_startP->at(tmks)>temp_pmom2){
                          temp_selpindex = tmks;
                          temp_pmom2=t->true_beam_daughter_startP->at(tmks);
                      }
                 }

                 _event_histo_1d->h_PiAbs_sig_true_totalKE_protonvsneutron->Fill(totalKE_proton, totalKE_neutron);

                 _event_histo_1d->h_PiAbs_sig_true_Ptmissing->Fill(true_Ptmissing);
                 _event_histo_1d->h_PiAbs_sig_true_Ptmissing_onlyproton->Fill(TMath::Sqrt(totalPxtest_onlyproton*totalPxtest_onlyproton+totalPytest_onlyproton*totalPytest_onlyproton));


                 if(temp_selpindex<0) continue;

                   _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_mom->Fill(t->true_beam_daughter_startP->at(temp_selpindex),temp_pmom);
                   _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_costheta->Fill(t->true_beam_daughter_startPz->at(temp_selpindex)/t->true_beam_daughter_startP->at(temp_selpindex),temp_pcostheta);
                 if(t->true_beam_daughter_startPy->at(temp_selpindex)>0){
                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_phi->Fill(TMath::ACos(t->true_beam_daughter_startPx->at(temp_selpindex)/TMath::Sqrt(t->true_beam_daughter_startPx->at(temp_selpindex)*t->true_beam_daughter_startPx->at(temp_selpindex)+t->true_beam_daughter_startPy->at(temp_selpindex)*t->true_beam_daughter_startPy->at(temp_selpindex))),temp_pphi);
                 } else {
                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_phi->Fill(-TMath::ACos(t->true_beam_daughter_startPx->at(temp_selpindex)/TMath::Sqrt(t->true_beam_daughter_startPx->at(temp_selpindex)*t->true_beam_daughter_startPx->at(temp_selpindex)+t->true_beam_daughter_startPy->at(temp_selpindex)*t->true_beam_daughter_startPy->at(temp_selpindex))),temp_pphi);

                 }

            for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

              if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){

                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                   if(t->reco_beam_calo_wire->size()>0 && t->reco_beam_incidentEnergies->size()>0){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }//end of if the wire number found in the vector of tracks
                  }   
                } //end of loop over all the wires before the end point of a beam track
                if(nincident>0){
                     eincident = eincident/(nincident*1.0);
                     _event_histo->hsig_reco_beam_incidentEnergy[is]->Fill(eincident);
                }
              }
            } 


        //increment interaction counter
        ++interaction_sel[true_sliceID];
   
        //increment incident counter
        for (int i = 0; i<=sliceID; ++i) ++incident_sel[i];




        }//end of if this is a signal event

 




        else if(isdata==0 && isChxBKG_withthresh) {
             _event_histo_1d->h_PiAbs_chxbac_pmult->Fill(t->reco_daughter_allTrack_ID->size());
             int Nreco_pion0=0;
             temp_plen = -999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;
             for(unsigned int tmk=0; tmk<t->reco_daughter_allTrack_len->size(); tmk++){
                 _event_histo_1d->h_PiAbs_chxbac_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_chxbac_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_chxbac_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                 _event_histo_1d->h_PiAbs_chxbac_pmom->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk)); 



                 if(t->reco_daughter_allTrack_ID->size()>0){
                      _event_histo_1d->h_chxbac_Pmissing->Fill(reco_Pmissing);
                      _event_histo_1d->h_chxbac_Emissing->Fill(Emissing_Qsubtracted_reco);
                      _event_histo_1d->h_chxbac_Ptmissing->Fill(reco_Ptmissing);  
                      _event_histo_1d->h_chxbac_Plongit->Fill(reco_Pz);
                 }
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk)) == 22 
                 || abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk)) == 11){
                     Nreco_pion0++;
                 }
                 //get the most energetic proton momentum and angle
                 if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles

             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "chxbac", fname_geant_pm1, wgts_geant_pm1);
             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "chxbac", fname_geant_pm1, wgts_geant_pm1);




             _event_histo_1d->h_reco_photon_chx->Fill(Nreco_pion0);
 
        for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

            if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){

                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                   if(t->reco_beam_calo_wire->size()>0 && t->reco_beam_incidentEnergies->size()>0){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }//end of if the wire number found in the vector of tracks
 
                   }
                }//end of loop over all the wire numbers before the end of beam track
                if(nincident>0){
                   eincident = eincident/(nincident*1.0);
                   _event_histo->hchxbac_reco_beam_incidentEnergy[is]->Fill(eincident);
                } 
            } 
        }
        int sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        //double reco_endz = Sce_Corrected_endZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
        //double reco_endz = t->reco_beam_endZ;
        //sliceID = int((reco_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); 

        //if (sliceID>=nslices) continue;  //ignore the events with the slice ID >  nslices (24);

        //increment interaction counter
        ++selected_bkg[sliceID];
  
        }//end of if this is a chxbac event
        else if(isdata==0 && isReaBKG_withthresh) {
             _event_histo_1d->h_PiAbs_reabac_pmult->Fill(t->reco_daughter_allTrack_ID->size());
             if(t->reco_daughter_allTrack_ID->size()>0){
                      _event_histo_1d->h_reabac_Ptmissing->Fill(reco_Ptmissing);
                      _event_histo_1d->h_reabac_Plongit->Fill(reco_Pz); 
                      _event_histo_1d->h_reabac_Pmissing->Fill(reco_Pmissing); 
                      _event_histo_1d->h_reabac_Emissing->Fill(Emissing_Qsubtracted_reco);
             }
             int Nreco_pionpm=0;
             int Nreco_pion0_rea=0;
             temp_plen=-999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;
             for(unsigned int tmk=0; tmk<t->reco_daughter_PFP_true_byHits_PDG->size(); tmk++){
                 _event_histo_1d->h_PiAbs_reabac_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_reabac_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_reabac_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                 _event_histo_1d->h_PiAbs_reabac_pmom->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk)); 
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk)) == 22 
                 || abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk)) == 11 ){
                     Nreco_pion0_rea++;
                     //std::cout<<"Parent PDG code of this photon is "<<t->reco_daughter_PFP_true_byHits_parPDG->at(tmk)<<std::endl;
                 }

                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk) == 211)){
                     Nreco_pionpm++;
                 }
                 //get the most energetic proton momentum and angle
                 /*
                 */
 
                 if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles

             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "reabac", fname_geant_pm1, wgts_geant_pm1);
             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "reabac", fname_geant_pm1, wgts_geant_pm1);



             _event_histo_1d->h_reco_pionpm_rea->Fill(Nreco_pionpm + Nreco_pion0_rea);
             _event_histo_1d->h_reco_photon_rea->Fill(Nreco_pion0_rea);


        for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

            if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){

                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                   if(t->reco_beam_calo_wire->size()>0 && t->reco_beam_incidentEnergies->size()>0){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }//end of if the wire number found in the vector of tracks
 
                   }
                }//end of loop over all the wire numbers before the end of beam track
                if(nincident>0){
                   eincident = eincident/(nincident*1.0);
                   _event_histo->hreabac_reco_beam_incidentEnergy[is]->Fill(eincident);
                } 
            } 
        }
        int sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        //double reco_endz = Sce_Corrected_endZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
        //double reco_endz = t->reco_beam_endZ;
        //sliceID = int((reco_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); 

        //if (sliceID>=nslices) continue;  //ignore the events with the slice ID >  nslices (24);

        //increment interaction counter
        ++selected_bkg[sliceID];
  
        }//end of if this is a reabac event
        else if(isdata==0)  { 
             std::cout<<"this is nonpion beam event"<<std::endl;
             for(unsigned int tmk=0; tmk<t->reco_daughter_PFP_true_byHits_PDG->size(); tmk++){
                 _event_histo_1d->h_PiAbs_other_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_other_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_other_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                 _event_histo_1d->h_PiAbs_other_pmom->Fill(t->reco_daughter_allTrack_momByRange_alt_proton->at(tmk)); 
             } 
             if(t->reco_daughter_allTrack_ID->size()>0){
                      _event_histo_1d->h_other_Pmissing->Fill(reco_Pmissing);
                      _event_histo_1d->h_other_Emissing->Fill(Emissing_Qsubtracted_reco);
                      _event_histo_1d->h_other_Ptmissing->Fill(reco_Ptmissing);  
                      _event_histo_1d->h_other_Plongit->Fill(reco_Pz);
             }
 
        for(long unsigned int is=0; is<sizeof(slicebins)/sizeof(slicebins[0])-1; is++){

            if(t->reco_beam_calo_wire->back()>=slicebins[is] && t->reco_beam_calo_wire->back()<slicebins[is+1]){
                double eincident = 0.0;
                int nincident = t->reco_beam_calo_wire->back()-slicebins[is]+1;
                for(int ic=slicebins[is]; ic<=t->reco_beam_calo_wire->back()-1; ic++){
                   if(t->reco_beam_calo_wire->size()>0 && t->reco_beam_incidentEnergies->size()>0){
                    //ic is the wire number over here
                    if(std::find(t->reco_beam_calo_wire->begin(),t->reco_beam_calo_wire->end(), ic) !=t->reco_beam_calo_wire->end()){
                      if(t->reco_beam_incidentEnergies->size()>0 && t->reco_beam_TrkPitch->size()>0){                    
                         int realindex = -999;
                         //nh is the wire number find the readindex
                         for(long unsigned int kk=0; kk<t->reco_beam_calo_wire->size(); kk++){
                            if(t->reco_beam_calo_wire->at(kk)==ic) {realindex = kk;}
                         }
                         if(realindex >=0){
                            eincident = eincident + t->reco_beam_incidentEnergies->at(realindex);
                         }
                      } 
                       
                    }//end of if the wire number found in the vector of tracks
 
                   }
                }//end of loop over all the wire numbers before the end of beam track
                if(nincident>0){
                   eincident = eincident/(nincident*1.0);
                   _event_histo->hotherbac_reco_beam_incidentEnergy[is]->Fill(eincident);
                } 
 
            } 
        }
        int sliceID = t->reco_beam_calo_wire->back()/nwires_in_slice;
        //double reco_endz = Sce_Corrected_endZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
        //double reco_endz = t->reco_beam_endZ;
        //sliceID = int((reco_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); 

        //if (sliceID>=nslices) continue;  //ignore the events with the slice ID >  nslices (24);

        //increment interaction counter
        ++selected_bkg[sliceID];
        }//end of if this is an other background event
        //-----------------------------------------------------------------------------------
        //------------------------------------------------------------
        //LOG_NORMAL()<<"end of checking signal background"<<std::endl;

	//-----------------------------------------------------------  

       

  } //end of loop over all the events


	



  //=========================================================================
  LOG_NORMAL() << "Saving 1D Event Histo." << std::endl;

  file_out->WriteObject(_event_histo_1d, "UBXSecEventHisto1D");
  file_out->WriteObject(_event_histo, "UBXSecEventHisto");

  file_out->WriteObject(&hmap_trkmom_geant_pm1_bs, "hmap_trkmom_geant_pm1_bs");

  file_out->WriteObject(&bs_geant_pm1_eff_mumom_num, "bs_geant_pm1_eff_mumom_num");
  file_out->WriteObject(&bs_geant_pm1_eff_mumom_den, "bs_geant_pm1_eff_mumom_den");

  file_out->WriteObject(&hmap_trkcostheta_geant_pm1_bs, "hmap_trkcostheta_geant_pm1_bs");

  file_out->WriteObject(&bs_geant_pm1_eff_mucostheta_num, "bs_geant_pm1_eff_mucostheta_num");
  file_out->WriteObject(&bs_geant_pm1_eff_mucostheta_den, "bs_geant_pm1_eff_mucostheta_den");


  file_out->WriteObject(&bs_geant_pm1_true_reco_mom, "bs_geant_pm1_true_reco_mom");
  file_out->WriteObject(&bs_geant_pm1_true_reco_costheta, "bs_geant_pm1_true_reco_costheta");


 
  LOG_NORMAL() << "1D Event Histo saved." << std::endl;
  
 
  double slcid[nslices];
  double avg_incE[nslices];
  double avg_pitch[nslices];

  double efficiency_tj[nslices];
  double efficiency_tj_err[nslices];
 
  double tempxsec[nslices];
  double tempxsec_chx[nslices];

  double tempxsec_pmom[nslices][nbinspmom];
  double tempxsec_pcostheta[nslices][nbinspcostheta];
  double tempxsec_pphi[nslices][nbinspphi];


  double tempxsec_test[nslices];
  double tempxsec_test2[nslices];
  double tempxsec_test3[nslices];
  double tempxsec_test4[nslices];

  double tempxsec_chx_test[nslices];
  double tempxsec_abs_test[nslices];
  double bkgsuberr[nslices];

  double rms_incE[nslices];
  double rms_pitch[nslices];
  double unverr[nslices];
  TEfficiency* pEff = 0;






  for (int i = 0; i<nslices; ++i){
     slcid[i] = i;
     
     unverr[i] = 0;
     
     avg_incE[i] = incE[i]->GetMean();

     avg_pitch[i] = pitch[i]->GetMean();

     rms_incE[i] = incE[i]->GetRMS()/TMath::Sqrt(incE[i]->GetEntries());;

     rms_pitch[i] = pitch[i]->GetRMS()/TMath::Sqrt(pitch[i]->GetEntries());

     efficiency_tj[i]=interaction_sel[i]/interaction_gen[i];

     efficiency_tj_err[i]=efficiency_tj[i]*TMath::Sqrt(1/interaction_sel[i]+1/interaction_gen[i]);

     if( incident[i]>0 && avg_pitch[i]>0 && efficiency_tj[i]>0 ){


     tempxsec_test4[i]  =MAr/(Density*NA*avg_pitch[i])*(selected_tot[i]-selected_bkg[i])/incident[i]/efficiency_tj[i];
   

     tempxsec_test3[i]  =MAr/(Density*NA*avg_pitch[i])*interaction_sel[i]/incident[i];


     tempxsec_test2[i]  =MAr/(Density*NA*avg_pitch[i])*interaction_gen[i]/incident[i];


     tempxsec_test[i]     =MAr/(Density*NA*avg_pitch[i])*true_abs[i]/incident[i];

     tempxsec[i]          =MAr/(Density*NA*avg_pitch[i])*true_abs[i]/true_incident[i];

     for(int nb=0; nb<nbinspmom; nb++){
          tempxsec_pmom[i][nb]=(1/1e-27)*MAr/(Density*NA*avg_pitch[i])*h_energetic_pmom_gen[i]->GetBinContent(nb+1)/true_incident[i];
     }
     for(int nb=0; nb<nbinspcostheta; nb++){
          tempxsec_pcostheta[i][nb]=(1/1e-27)*MAr/(Density*NA*avg_pitch[i])*h_energetic_pcostheta_gen[i]->GetBinContent(nb+1)/true_incident[i];
     }
     for(int nb=0; nb<nbinspphi; nb++){
          tempxsec_pphi[i][nb]=(1/1e-27)*MAr/(Density*NA*avg_pitch[i])*h_energetic_pphi_gen[i]->GetEntries()/true_incident[i];
     }
          
     tempxsec_chx[i]      =MAr/(Density*NA*avg_pitch[i])*true_chx[i]/true_incident[i];         

     double avgthickness = 12;
     tempxsec_chx_test[i] = true_chx[i]/true_incident[i];
     tempxsec_abs_test[i] = true_abs[i]/true_incident[i];
     bkgsuberr[i] = TMath::Sqrt(selected_tot[i]+selected_bkg[i]);
   
          



     } else {tempxsec[i]=0.0; tempxsec_chx[i]=0.0; tempxsec_chx_test[i]=0.0; tempxsec_abs_test[i] = 0.0; tempxsec_test[i]=0.0; tempxsec_test2[i]=0.0; tempxsec_test3[i]=0.0; tempxsec_test4[i]=0.0;}
           
     tempxsec[i] = tempxsec[i]/1e-27;  
     tempxsec_chx[i] = tempxsec_chx[i]/1e-27;

     tempxsec_test[i] = tempxsec_test[i]/1e-27;
     tempxsec_test2[i] = tempxsec_test2[i]/1e-27;
     tempxsec_test3[i] = tempxsec_test3[i]/1e-27;
     tempxsec_test4[i] = tempxsec_test4[i]/1e-27;
     

  }//end of loop over all the nslices





  file_out->Write();

  LOG_NORMAL() << "All saved." << std::endl;
  
  file_out->Close();

  LOG_NORMAL() << "Output file closed." << std::endl;


  string outfilename0 = "/dune/app/users/jiang/dunetpc_analysis/PionAbsChxAna/UBAna/Main/mac/output_sliceIDmat.root";
  TFile output_sliceIDmat(outfilename0.c_str(), "RECREATE");
  /*h_true_reco_mom->Write();
  h_true_reco_costheta->Write();
  h_true_reco_phi->Write();

  h_Evttot->Write();
  dabsintE->Write();
  dchxintE->Write();
  for(int ig=0; ig<nslices; ig++){
      dsizewire_abs[ig]->Write();
      dsizewire_chx[ig]->Write();
      dslcID[ig]->Write();
      incE[ig]->Write();
      pitch[ig]->Write();
  }
  */ 
  std::cout<<"nslices+1 = "<<nslices+1<<std::endl;
  for(int k=0; k<nslices+1; k++){
      for(int l=0; l<nslices+1; l++){ 
         sliceIDmat_den->SetBinContent(k+1, l+1, intabs_array_den[k][l]); 
      }
  }

  //sliceIDmat_den->Write();


  output_sliceIDmat.Close();
  LOG_NORMAL()<<"output sliceIDmat file closed "<<std::endl; 
  string outfilename="/dune/app/users/jiang/dunetpc_analysis/PionAbsChxAna/UBAna/Main/mac/output_graphs.root";

  TFile output_newsf(outfilename.c_str(), "RECREATE");

  gr_inc0_slc = new TGraph(nslices, slcid, incident);

  gr_int_slc = new TGraph(nslices, slcid, true_interaction);
  
  gr_intabs_slc = new TGraph(nslices, slcid, true_abs);
  gr_intchx_slc = new TGraph(nslices, slcid, true_chx);

  gr_inc_slc = new TGraph(nslices, slcid, true_incident);

  gr_incE_slc = new TGraphAsymmErrors(nslices, slcid, avg_incE, unverr, unverr, rms_incE, rms_incE);


  gr_pitch_slc = new TGraphAsymmErrors(nslices, slcid, avg_pitch, unverr, unverr, rms_pitch, rms_pitch);


  gr_intsel_slc = new TGraph(nslices, slcid, interaction_sel);

  gr_intgen_slc = new TGraph(nslices, slcid, interaction_gen);

  gr_efficiency_slc = new TGraph(nslices, slcid, efficiency_tj);

  gr_tempxsec_slc = new TGraph(nslices, slcid, tempxsec);
  gr_tempxsec_chx_slc = new TGraph(nslices, slcid, tempxsec_chx);

  gr_tempxsec_chx_test_slc = new TGraph(nslices, slcid, tempxsec_chx_test);
  gr_tempxsec_abs_test_slc = new TGraph(nslices, slcid, tempxsec_abs_test);


  TGraph *gr_tempxsec_pmom_slc[nslices];

  TGraph *gr_tempxsec_pcostheta_slc[nslices];

  for(int idx=0; idx<nslices; idx++){
       double pxmom[nbinspmom];
       double xsecmom[nbinspmom];
       double pxcostheta[nbinspcostheta];
       double xseccostheta[nbinspcostheta];
       double binwidthpmom=1.2/nbinspmom;
       double binwidthpcostheta=1.2/nbinspcostheta;

       for(int jk=0; jk<nbinspmom; jk++){
          pxmom[jk] = (jk+0.5)*binwidthpmom;
          xsecmom[jk]=tempxsec_pmom[idx][jk];
       } 
       for(int jk=0; jk<nbinspcostheta; jk++){
          pxcostheta[jk] = (jk+0.5)*binwidthpcostheta;
          xseccostheta[jk]=tempxsec_pcostheta[idx][jk];
       } 

       /*TGraph *gr_libo_xsecvspmom = new TGraph(nbinspmom, pxmom,  xsecmom);
       gr_libo_xsecvspmom->Write(Form("gr_tempxsec_pmom_slc_%d",idx));

       TGraph *gr_libo_xsecvspcostheta = new TGraph(nbinspcostheta, pxcostheta,  xseccostheta);
       gr_libo_xsecvspcostheta->Write(Form("gr_tempxsec_pcostheta_slc_%d",idx));
       */

       gr_tempxsec_pmom_slc[idx]=new TGraph(nbinspmom, pxmom,  xsecmom);
       gr_tempxsec_pmom_slc[idx]->Write(Form("gr_tempxsec_pmom_%d",idx));
       gr_tempxsec_pcostheta_slc[idx]=new TGraph(nbinspcostheta, pxcostheta, xseccostheta);
       gr_tempxsec_pcostheta_slc[idx]->Write(Form("gr_tempxsec_pcostheta_%d",idx));

  }
  gr_tempxsec_test_slc = new TGraph(nslices, slcid, tempxsec_test);
  gr_tempxsec_test2_slc = new TGraph(nslices, slcid, tempxsec_test2);
  gr_tempxsec_test3_slc = new TGraph(nslices, slcid, tempxsec_test3);
  gr_tempxsec_test4_slc = new TGraph(nslices, slcid, tempxsec_test4);

  gr_tempxsec_vs_incE_slc = new TGraph(nslices, avg_incE, tempxsec);

  gr_tempxsec_chx_vs_incE_slc = new TGraph(nslices, avg_incE, tempxsec_chx);

  gr_tempxsec_test_vs_incE_slc = new TGraph(nslices, avg_incE, tempxsec_test);

  gr_selected_tot_slc = new TGraph(nslices, slcid, selected_tot);


  gr_selected_bkg_slc = new TGraph(nslices, slcid, selected_bkg);


  gr_inc0_slc->Write("gr_inc0_slc");
   
  gr_int_slc->Write("gr_int_slc");
  gr_intabs_slc->Write("gr_intabs_slc");
  gr_intchx_slc->Write("gr_intchx_slc");
  gr_intsel_slc->Write("gr_intsel_slc");
  gr_intgen_slc->Write("gr_intgen_slc");


  gr_inc_slc->Write("gr_inc_slc");
  gr_incE_slc->Write("gr_incE_slc");
  gr_pitch_slc->Write("gr_pitch_slc");

  gr_efficiency_slc->Write("gr_efficiency_slc");

  gr_tempxsec_slc->Write("gr_tempxsec_slc");
  gr_tempxsec_chx_slc->Write("gr_tempxsec_chx_slc");

  gr_tempxsec_chx_test_slc->Write("gr_tempxsec_chx_test_slc");
  gr_tempxsec_abs_test_slc->Write("gr_tempxsec_abs_test_slc");

  gr_tempxsec_vs_incE_slc->SetName("gr_tempxsec_vs_incE_slc");
  gr_tempxsec_vs_incE_slc->Write();

  gr_tempxsec_test_vs_incE_slc->SetName("gr_test_tempxsec_vs_incE_slc");
  gr_tempxsec_test_vs_incE_slc->Write();

  gr_tempxsec_chx_vs_incE_slc->SetName("gr_chx_tempxsec_vs_incE_slc");
  gr_tempxsec_chx_vs_incE_slc->Write();
   
  gr_tempxsec_test_slc->Write("gr_tempxsec_test_slc");
  gr_tempxsec_test2_slc->Write("gr_tempxsec_test2_slc");
  gr_tempxsec_test3_slc->Write("gr_tempxsec_test3_slc");
  gr_tempxsec_test4_slc->Write("gr_tempxsec_test4_slc");

 
  gr_selected_tot_slc->Write("gr_selected_tot_slc");
  gr_selected_bkg_slc->Write("gr_selected_bkg_slc");

  


  output_newsf.Write();





  //=============================================================================================
  TH1D  *selected_signal_percut = new TH1D("selected_signal_percut", "selected_percut", 4, 0, 4);

  TH1D  *selected_percut = new TH1D("selected_percut", "selected_percut", 4, 0, 4);

  TH1D  *generated_signal_percut = new TH1D("generated_signal_percut", "generated_signal_percut", 4, 0, 4);

  TH1D  *acceptance_percut = new TH1D("acceptance_percut", "acceptance_percut", 4, 0, 4);

  //selected events after each cut
  /*_event_histo_1d->*/selected_percut->SetBinContent(1, Noriabs_withthresh+Norichx_withthresh+Norirea_withthresh);
  /*_event_histo_1d->*/selected_percut->SetBinContent(2, Ntmcutabs_withthresh+Ntmcutchx_withthresh+Ntmcutrea_withthresh);
  /*_event_histo_1d->*/selected_percut->SetBinContent(3, Nchi2abs_withthresh+Nchi2chx_withthresh+Nchi2rea_withthresh);
  /*_event_histo_1d->*/selected_percut->SetBinContent(4, Nshwcutabs_withthresh+Nshwcutchx_withthresh+Nshwcutrea_withthresh);


  //selected signal events after each cut
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(1, Noriabs_withthresh);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(2, Ntmcutabs_withthresh);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(3, Nchi2abs_withthresh);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(4, Nshwcutabs_withthresh);


  //generated signal events
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(1, Noriabs_withthresh);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(2, Noriabs_withthresh);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(3, Noriabs_withthresh);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(4, Noriabs_withthresh);


  double acceptance_tmcut = (double)Ntmcutabs_withthresh/(Ntmcutabs_withthresh+Ntmcutchx_withthresh+Ntmcutrea_withthresh);
  acceptance_tmcut = acceptance_tmcut*Ntmcutabs_withthresh/Noriabs_withthresh;

  double acceptance_chi2cut = (double)Nchi2abs_withthresh/(Nchi2abs_withthresh+Nchi2chx_withthresh+Nchi2rea_withthresh);
  acceptance_chi2cut = acceptance_chi2cut*Nchi2abs_withthresh/Noriabs_withthresh;

  double acceptance_shwcut = (double)Nshwcutabs_withthresh/(Nshwcutabs_withthresh+Nshwcutchx_withthresh+Nshwcutrea_withthresh);
  acceptance_shwcut = acceptance_shwcut*Nshwcutabs_withthresh/Noriabs_withthresh;

  


  acceptance_percut->SetBinContent(1, (double)Noriabs_withthresh/(Noriabs_withthresh+Norichx_withthresh+Norirea_withthresh)*(Noriabs_withthresh/Noriabs_withthresh) );
  acceptance_percut->SetBinContent(2, acceptance_tmcut);
  acceptance_percut->SetBinContent(3, acceptance_chi2cut);
  acceptance_percut->SetBinContent(4, acceptance_shwcut );
  
  


  std::vector<std::string> cut_names = {"initial", "TrunMeandEdX", "Chi2", "TrackScore"};

  for(int i=0; i<4; i++){
      std::cout<<cut_names.at(i)<<" & "<</*_event_histo_1d->*/selected_signal_percut->GetBinContent(i+1)<<std::endl;
  }


  TCanvas * canvas_eff_pur_graph_percut = new TCanvas();

  canvas_eff_pur_graph_percut->SetLeftMargin(0.05157593);
  canvas_eff_pur_graph_percut->SetRightMargin(0.1475645);
  canvas_eff_pur_graph_percut->SetTopMargin(0.04210526);
  canvas_eff_pur_graph_percut->SetBottomMargin(0.1578947);

  TH1F *h = new TH1F("h","",4, 0, 4);
  h->SetMaximum(1);
  h->GetXaxis()->SetBinLabel(1,"Initial");
  h->GetXaxis()->SetBinLabel(2,"TrunMeandQdx");
  h->GetXaxis()->SetBinLabel(3,"Chi2");
  h->GetXaxis()->SetBinLabel(4,"TrackScore");

  h->GetXaxis()->SetLabelOffset(0.009);
  h->GetXaxis()->SetLabelSize(0.06);

  h->Draw();  


  TEfficiency* pEff_percut = new TEfficiency(*selected_signal_percut,*generated_signal_percut);
  pEff_percut->SetTitle("EfficiencyPerCut;Cut index;Efficiency");
  pEff_percut->SetLineColor(kGreen+3);
  pEff_percut->SetMarkerColor(kGreen+3);
  pEff_percut->SetMarkerStyle(20);
  pEff_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pEff_percut_graph = pEff_percut->CreateGraph();
  for (int i = 0; i < 4; i++) {
    pEff_percut_graph->SetPointEXhigh(i, 0.);
    pEff_percut_graph->SetPointEXlow(i, 0.);
  }
  auto axis = pEff_percut_graph->GetYaxis();
  axis->SetLimits(0.,1.); 

  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(1,"Initial");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(2,"TrunMeandQdx");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(3,"Chi2");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(4,"TrackScore");
 
  pEff_percut_graph->Draw("PL"); 



	  
  TEfficiency* pPur_percut = new TEfficiency(/*_event_histo_1d->*/*selected_signal_percut, /*_event_histo_1d->*/*selected_percut);
  pPur_percut->SetTitle("Purity Per Cut; Cut index; Purity");
  pPur_percut->SetLineColor(kRed+3);
  pPur_percut->SetMarkerColor(kRed+3);
  pPur_percut->SetMarkerStyle(20);
  pPur_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pPur_percut_graph = pPur_percut->CreateGraph();
  
  for (int i = 0; i < 4; i++) {
    pPur_percut_graph->SetPointEXhigh(i, 0.);
    pPur_percut_graph->SetPointEXlow(i, 0.);
  }

  pPur_percut_graph->Draw("PL");





  TLegend* l = new TLegend(0.4842407,0.8168421,0.777937,0.9221053,NULL,"brNDC");
  l->AddEntry(pEff_percut_graph,"Efficiency");
  l->AddEntry(pPur_percut_graph,"Purity");
  //leg->AddEntry(gr3,"Neutrino MCFlash","l");
  //  
  l->Draw();
 
  canvas_eff_pur_graph_percut->SaveAs("eff_pur_test.png");

  TCanvas * canvas_acceptance_graph_percut = new TCanvas();
  TH1F *hh = new TH1F("hh","",4, 0, 4);
  hh->SetMaximum(1);
  hh->GetXaxis()->SetBinLabel(1,"Initial");
  hh->GetXaxis()->SetBinLabel(2,"TrackScore");
  hh->GetXaxis()->SetBinLabel(3,"TrunMeandQdx");
  hh->GetXaxis()->SetBinLabel(4,"Chi2");

  hh->GetXaxis()->SetLabelOffset(0.009);
  hh->GetXaxis()->SetLabelSize(0.06);

  hh->Draw();  
 



  Int_t nn=4;
  Double_t x[nn], y[nn];
  for (int i = 0; i < 4; i++) {
    x[i]=i+0.5;
    y[i]=acceptance_percut->GetBinContent(i+1);
    //std::cout<<"xnn = "<<x[i]<<" ynn = "<<y[i]<<std::endl;
  }

 

  TGraph *gr_percut= new TGraph(nn, x, y);
  gr_percut->Draw("PL");
  gr_percut->SetTitle("Acceptance Per Cut; Cut index; Acceptance");
  gr_percut->SetLineColor(kBlue);
  gr_percut->SetMarkerColor(kBlue);
  gr_percut->SetMarkerStyle(20);
  gr_percut->SetMarkerSize(0.8);
  gr_percut->Draw("PL");
  canvas_acceptance_graph_percut->SaveAs("acceptance_test.png");
	  

  std::cout<<"Noriabs = "<<Noriabs<<std::endl;
  std::cout<<"Norichx = "<<Norichx<<std::endl;
  std::cout<<"Norirea = "<<Norirea<<std::endl;
  std::cout<<"Noriother = "<<Noriother<<std::endl;
  std::cout<<"TestGenSig= "<<TestGenSig<<std::endl;

  std::cout<<"Total pion events is "<<fpion_evt_index.size()<<std::endl;
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
  std::cout<<"Noriabs_withthresh = "<<Noriabs_withthresh<<std::endl;
  std::cout<<"Norichx_withthresh = "<<Norichx_withthresh<<std::endl;
  std::cout<<"Norirea_withthresh = "<<Norirea_withthresh<<std::endl;
  std::cout<<"Noriother_withthresh = "<<Noriother_withthresh<<std::endl;
 
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  std::cout<<"Ntmcutabs_withthresh = "<<Ntmcutabs_withthresh<<std::endl;
  std::cout<<"Ntmcutchx_withthresh = "<<Ntmcutchx_withthresh<<std::endl;
  std::cout<<"Ntmcutrea_withthresh = "<<Ntmcutrea_withthresh<<std::endl;
  std::cout<<"Ntmcutother_withthresh = "<<Ntmcutother_withthresh<<std::endl;

  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  std::cout<<"Nchi2abs_withthresh = "<<Nchi2abs_withthresh<<std::endl;
  std::cout<<"Nchi2chx_withthresh = "<<Nchi2chx_withthresh<<std::endl;
  std::cout<<"Nchi2rea_withthresh = "<<Nchi2rea_withthresh<<std::endl;
  std::cout<<"Nchi2other_withthresh = "<<Nchi2other_withthresh<<std::endl;

  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  std::cout<<"Nshwcutabs_withthresh = "<<Nshwcutabs_withthresh<<std::endl;
  std::cout<<"Nshwcutchx_withthresh = "<<Nshwcutchx_withthresh<<std::endl;
  std::cout<<"Nshwcutrea_withthresh = "<<Nshwcutrea_withthresh<<std::endl;
  std::cout<<"Nshwcutother_withthresh = "<<Nshwcutother_withthresh<<std::endl;

  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;


  std::cout<<"Ntotal beam event is "<<Ntotal_beam<<std::endl; 
  std::cout<<"Ntotal tmdqdx event is "<<Ntotal_tmdqdx<<std::endl; 
  std::cout<<"Ntotal ch2 event is "<<Ntotal_chi2<<std::endl; 
  std::cout<<"Ntotal shwid event is "<<Ntotal_shwid<<std::endl; 
 
  //=========================================================================
  std::cout<<"total number of showers (trkscore) "<<totalshower_trkscore<<std::endl;
  std::cout<<"total number of showers (trknhits) "<<totalshower_trknhits<<std::endl;
  std::cout<<"total number of showers (trkdist) "<<totalshower_trkdist<<std::endl;
  std::cout<<"total number of showers (trkeng) "<<totalshower_trkeng<<std::endl;

  std::cout<<"total number of gammas (trkscore) "<<totalgamma_trkscore<<std::endl;
  std::cout<<"total number of gammas (trknhits) "<<totalgamma_trknhits<<std::endl;
  std::cout<<"total number of gammas (trkdist) "<<totalgamma_trkdist<<std::endl;
  std::cout<<"total number of gammas (trkeng) "<<totalgamma_trkeng<<std::endl;


  //==========================================================================

  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << std::endl << std::endl;
  LOG_NORMAL() << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  //rootapp->Run();
  //rootapp->Terminate(0);
  
  return;
}

#endif
