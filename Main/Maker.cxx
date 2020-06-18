#ifndef __MAIN_MAKER_CXX__
#define __MAIN_MAKER_CXX__

#include "Maker.h"

using namespace Base;

void Main::Maker::SetInputFile(std::string in)
{
  filen = in;
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

  h->SetNBinsX(n_bins_double_mucostheta);

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
}



//___________________________________________________________________________________________________
void Main::Maker::AddPolyBins(BootstrapTH2DPoly h) {

  // std::map<int, std::pair<int, int>> _exclusion_map;
  // _exclusion_map[0] = std::make_pair(2, 3);

  // h->SetNBinsX(n_bins_double_mucostheta);

  for (int y = 0; y < n_bins_double_mumom; y++) {
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
const std::vector<double> Main::Maker::truncatedMean_xglu(double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX){

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
  string pattern_add = filen_add;
  
  
  
  //*************************
  //* Getting POTs
  //*************************
  
  
  
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

  _event_histo_1d = new UBXSecEventHisto1D();
  _event_histo_1d->InitializeBootstraps();
  //========================================================================
  //int evts = chain_ubxsec->GetEntries();
  std::cout<<"total number of events is "<<evts<<std::endl;
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

  int n_bins_mumom = 10;
  int n_bins_mucostheta = 10;
  double bins_mumom[11];
  double bins_mucostheta[11];
  for(int i=0; i<11; i++){
     bins_mumom[i]=1.2/10.0*i;
     bins_mucostheta[i]=-1.0+2.0/10.0*i;
  }
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





  //std::cout<<"libo test before looping over all the events"<<std::endl;
  
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
  

  //get the first pion event
  double event_weight = 1.0;
  for(int i= _initial_entry; i<evts; i++){
	if(i !=0) DrawProgressBar((double)i/(double)evts, barWidth);
        

	chain_ubxsec->GetEntry(i);
        if(isdata==0){
        if(i==fpion_evt_index[0] && isdata==0 && _fill_bootstrap_geant){
          
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
	
	//========================================================================
	//select muon and pion beam events
	if(isdata==0){
	
          if(abs(t->true_beam_PDG) != 211 && abs(t->true_beam_PDG) !=13) continue; 
          if(!isBeamType(t->reco_beam_type)) continue;
          if(!manual_beamPos_mc(t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ, 
                              t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                              t->true_beam_startDirX, t->true_beam_startDirY, t->true_beam_startDirZ,
                              t->true_beam_startX, t->true_beam_startY, t->true_beam_startZ)) continue;
          if(!endAPA3(t->reco_beam_endZ) )continue;
        } 
        if(isdata==1){
          if(!data_beam_PID(t->data_BI_PDG_candidates)) continue;
          if(!manual_beamPos_data(t->event, t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ,
                                t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                                t->data_BI_X, t->data_BI_Y, t->data_BI_dirX, t->data_BI_dirY,
                                t->data_BI_dirZ, t->data_BI_nMomenta, t->data_BI_nTracks)) continue;
          if(!endAPA3(t->reco_beam_endZ) )continue;
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
        } else if( ngen_pipm_withthresh>0){
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
        else if(t->true_daughter_nPi0 > 0 ) {Norichx++;
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
        //================================================================
        std::vector<double> wgts_geant_pm1;  //size of wgts_geant_pm1 would be 2*reweight names
        if(isdata==0){
        for(long unsigned int i=0; i<t->g4rw_primary_var->size(); i++){
            wgts_geant_pm1.push_back(t->g4rw_primary_plus_sigma_weight->at(i));
            wgts_geant_pm1.push_back(t->g4rw_primary_minus_sigma_weight->at(i));
        }
        }//end of if isdata for filling reweight factors for Geant4
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
           bs_geant_pm1_eff_mumom_den.Fill(t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);
           bs_geant_pm1_eff_mucostheta_den.Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);
            TestGenSig++;
          } 
             
        }//end of selected signal before event selection
        
        }//end of if isdata==0
        
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
        //loop over all the true beam daughter particles
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
        /*
        * loop over all the reconstructed daughter particles
        * get the photons and fill the parent ID of the photon into 
        * vec_parid vector, which will be checked by the PDG 
        */
        //=========================================================================
        if(isdata==0){
        vector<int> vec_parid; //find out the parID of the photons but no loop
        vec_parid.clear();

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
	
        int nonprotoncand = 0;
        int nshwcand = 0;

	int ntmdqdxcand = 0;
        int ntmdqdxcand_withthresh = 0;

        int nonprotoncand_withthresh = 0;
        TVector3 recobeamendv3, recodauendv3, recodaustartv3, recogranddaustartv3, recogranddauendv3;
        recobeamendv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);







        std::vector<double> daughter_distance3D, daughter_distance3D_shower;
        std::vector<double> daughter_angle3D;

        std::vector<std::vector<double>> *trkdedx_test_ptr=t->reco_daughter_allTrack_calibrated_dEdX_SCE;

        const std::vector<double> reco_daughter_allTrack_truncLibo_dEdX_test = truncatedMean_xglu(libo_low, libo_high, *trkdedx_test_ptr);
 


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
        daughter_distance3D_shower =   compute_distanceVertex(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ,
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ,      
                                    newrdshwr_startX, newrdshwr_startY, newrdshwr_startZ);      
         
 


	//for(unsigned int jj=0; jj<t->reco_daughter_PFP_true_byHits_PDG->size(); jj++){
        for(unsigned int jj=0; jj<t->reco_daughter_PFP_trackScore_collection->size(); jj++){
            _event_histo_1d->h_trackscore_all->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut){
               _event_histo_1d->h_trackscore_shwlike_all->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
               _event_histo_1d->h_nhits_shwlike_all->Fill(t->reco_daughter_PFP_nHits->at(jj));
               if(isdata == 0){
                   if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
                        _event_histo_1d->h_trackscore_shwlike_photon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        _event_histo_1d->h_nhits_shwlike_photon->Fill(t->reco_daughter_PFP_nHits->at(jj));

                   }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){ 
                        _event_histo_1d->h_trackscore_shwlike_electron->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        _event_histo_1d->h_nhits_shwlike_electron->Fill(t->reco_daughter_PFP_nHits->at(jj));
                         
                   }else {
                        _event_histo_1d->h_trackscore_shwlike_nonphoton->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));
                        _event_histo_1d->h_nhits_shwlike_nonphoton->Fill(t->reco_daughter_PFP_nHits->at(jj));
              
                   }
               } 
            }
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

            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && t->reco_daughter_PFP_nHits->at(jj)>=40) {nshwcand ++; }
            
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<trkscorecut && t->reco_daughter_PFP_nHits->at(jj)>=40) continue;
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

             
 
            _event_histo_1d->h_tmdqdx->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));
            if(isdata ==0){
	    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	           _event_histo_1d->h_tmdqdx_proton->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
            _event_histo_1d->h_tmdqdxvsrange_proton->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));

	    }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
	           _event_histo_1d->h_tmdqdx_pionpm->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj));        
            _event_histo_1d->h_tmdqdxvsrange_pionpm->Fill(reco_daughter_allTrack_truncLibo_dEdX_test.at(jj), t->reco_daughter_allTrack_len->at(jj));
	    }
            } //end of if isdata ==0 


            


            //------------------------------------------------------------------
            bool isPCand = false;
            bool isPiCand = false;
            double temp_tmdqdx=reco_daughter_allTrack_truncLibo_dEdX_test.at(jj);

            //==================================================================================
            if(temp_tmdqdx<cut_dEdX_low || (temp_tmdqdx>cut_dEdX_high && temp_tmdqdx<3.8 )){ //fill out chi2 of the transition region
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
            //================================================================================= 
            if(temp_tmdqdx>=cut_dEdX_low && temp_tmdqdx<=cut_dEdX_high 
                  && t->reco_daughter_PFP_trackScore_collection->at(jj)>trkscorecut) 
            {ntmdqdxcand_withthresh++;}
            if(temp_tmdqdx>3.8) {
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
                  ntmdqdxcand++;
               if(isdata==0){
               if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                     _event_histo_1d->h_mom_selectedtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     _event_histo_1d->h_mom_tmdqdxtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
               }
               }//end of if isdata==0
            }
            if((temp_tmdqdx<3.8 && temp_tmdqdx>cut_dEdX_high) || temp_tmdqdx<cut_dEdX_low) {  //transition region
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
            if(isPCand == false /*&& t->reco_daughter_allTrack_momByRange_proton->at(jj) > momthreshcut*/){
              if(t->reco_daughter_allTrack_Chi2_ndof->at(jj)>10){
                nonprotoncand_withthresh ++;
              }
            }
            if(isPCand == true){
                //-------------------------------------------------------------------------------
                if(isdata==0){
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
                //check if this proton is a daughter particles

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
                //-------fill histograms for momentum track range, angular resolution of the selected proton
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                  _event_histo_1d->h_selproton_momreso->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj),t->reco_daughter_allTrack_momByRange_proton->at(jj));
                 _event_histo_1d->h_trklen_reso->Fill((t->reco_daughter_allTrack_len->at(jj)-t->reco_daughter_PFP_true_byHits_len->at(jj))/t->reco_daughter_PFP_true_byHits_len->at(jj));
                  _event_histo_1d->h_selproton_costhetareso->Fill(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj), TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj)));

                  _event_histo_1d->h_totalp_mom->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));

                  if(t->reco_daughter_PFP_true_byHits_endProcess->at(jj) !="protonInelastic") {
                  _event_histo_1d->h_selproton_momreso_new->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj)-1);
                  _event_histo_1d->h_selproton_momPxreso_new->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Cos(t->reco_daughter_allTrack_Phi->at(jj))*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj))/t->reco_daughter_PFP_true_byHits_startPx->at(jj)-1);
                  _event_histo_1d->h_selproton_momPyreso_new->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Sin(t->reco_daughter_allTrack_Phi->at(jj))*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj))/t->reco_daughter_PFP_true_byHits_startPy->at(jj)-1);
                  _event_histo_1d->h_selproton_momPzreso_new->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj))/t->reco_daughter_PFP_true_byHits_startPz->at(jj)-1);
                 
                  _event_histo_1d->h_selproton_costhetareso_new->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj))/(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj))-1);
                  }
                  if(t->reco_daughter_PFP_true_byHits_endProcess->at(jj)=="protonInelastic"){
                    _event_histo_1d->h_reintp_mom->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }                

                }
                }//end of if isdata==0 
                //-----------------------------------------------------------------------------
            } //end of if pcand is true
            else if(isPCand == false ) {
              if(t->reco_daughter_allTrack_Chi2_ndof->at(jj)>10){
                nonprotoncand++;
              }
            }
            //check if this is a pion candidate and fill out track range of the selected pion candidate
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
         
        Ntotal_beam++;  //total number of events past beam selection

        std::vector<int> *daughter_nhits_ptr=t->reco_daughter_PFP_nHits;
        std::vector<int> nhitsref=*daughter_nhits_ptr;
        std::vector<double> *daughter_trkscore_ptr=t->reco_daughter_PFP_trackScore_collection;
        std::vector<double> trkscoreref=*daughter_trkscore_ptr;

        std::vector<int> *daughter_ID_ptr=t->reco_daughter_allTrack_ID;
        std::vector<int> trkIDref=*daughter_ID_ptr;

        std::vector<std::vector<double>> *trkdedx_ptr=t->reco_daughter_allTrack_calibrated_dEdX_SCE;
        std::vector<std::vector<double>> trkdedxref=*trkdedx_ptr;


        const std::vector<double> reco_daughter_allTrack_truncLibo_dEdX = truncatedMean_xglu(libo_low, libo_high, *trkdedx_ptr);
 

        if(!secondary_noPion_libo(*daughter_trkscore_ptr,
             *daughter_ID_ptr,
             reco_daughter_allTrack_truncLibo_dEdX)) continue;



        //if(ntmdqdxcand_withthresh> 0) continue;

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
        



        if(sel_abs) if(has_shower_nHits(trkscoreref, nhitsref)) continue;


        if(sel_chx) if(trkscoreref.size()==0) continue;
        if(sel_chx) if(!has_shower_nHits(trkscoreref, nhitsref)) continue;
        //if(sel_chx) if(nshwcand>0) continue;

 
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
 
             _event_histo_1d->h_PiAbs_sel_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(ntrk));
             _event_histo_1d->h_PiAbs_sel_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk)));
             _event_histo_1d->h_PiAbs_sel_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));
             
             _event_histo_1d->h_PiAbs_sel_pcosthetax->Fill(TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk))));
             _event_histo_1d->h_PiAbs_sel_pphi->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));

             if(isdata==0){ 
               //check if the parID can be found in the true beam daughters
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
                     std::cout<<"Event Number is "<<t->event<<"  run number is "<<t->run<<" Particle PDG = "<<t->reco_daughter_PFP_true_byHits_PDG->at(ntrk)<<" Parent ID= "<<t->reco_daughter_PFP_true_byHits_parID->at(ntrk)<<" Parent PDG = "<<t->reco_daughter_PFP_true_byHits_parPDG->at(ntrk)<<std::endl;
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
             
             total_recopKE += TMath::Sqrt(t->reco_daughter_allTrack_momByRange_proton->at(ntrk)*t->reco_daughter_allTrack_momByRange_proton->at(ntrk) + ProtonMass*ProtonMass) - ProtonMass;
             total_pPx += t->reco_daughter_allTrack_momByRange_proton->at(ntrk)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPy += t->reco_daughter_allTrack_momByRange_proton->at(ntrk)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPz += t->reco_daughter_allTrack_momByRange_proton->at(ntrk)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));

             if(t->reco_daughter_allTrack_len->at(ntrk) > temp_plen){
                                temp_pindex = ntrk;
                                temp_plen=t->reco_daughter_allTrack_len->at(ntrk);
                                temp_pmom = t->reco_daughter_allTrack_momByRange_proton->at(ntrk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(ntrk);
                                temp_pcosthetax=TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk)));
             }

        }

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
        //std::cout<<"Reco Energy of Beam is "<<t->reco_beam_interactingEnergy<<std::endl;
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
                      _event_histo_1d->h_PiAbs_sig_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 


                      if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_pindex=tmk;
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=t->reco_daughter_allTrack_momByRange_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                       }

                 }//loop over all the proton candidates




                 if(temp_pindex>-999){ //fill histograms with at least one proton
                   if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "signal", fname_geant_pm1, wgts_geant_pm1);
                   if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "signal", fname_geant_pm1, wgts_geant_pm1);
                 }


                 int temp_genpindex=-999; double temp_genpmom=-999.0; double temp_genpcostheta=-999.0;
                 for(unsigned int tmk=0; tmk<t->true_beam_daughter_startP->size(); tmk++){
                      if(abs(t->true_beam_daughter_PDG->at(tmk)) !=2212) continue;
                      if(t->true_beam_daughter_startP->at(tmk)>temp_genpmom){
                      temp_genpindex = tmk;
                      temp_genpmom=t->true_beam_daughter_startP->at(tmk);
                      temp_genpcostheta=t->true_beam_daughter_startPz->at(tmk)/t->true_beam_daughter_startP->at(tmk);
                      }           
                 }
 
                 //std::cout<<"got the most energetic information of pindex costheta and mom"<<std::endl;
                 if(temp_genpindex> -999){ 
                    if(!isdata && _fill_bootstrap_geant) bs_geant_pm1_eff_mumom_num.Fill(temp_genpmom, event_weight, wgts_geant_pm1);        
                    if(!isdata && _fill_bootstrap_geant) bs_geant_pm1_eff_mucostheta_num.Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex), event_weight, wgts_geant_pm1);        


                    if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_genpmom, temp_pmom, event_weight, bs_geant_pm1_true_reco_mom, fname_geant_pm1, wgts_geant_pm1);
                    if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_genpcostheta, temp_pcostheta, event_weight, bs_geant_pm1_true_reco_costheta, fname_geant_pm1, wgts_geant_pm1);
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
        }

 




        else if(isdata==0 && isChxBKG_withthresh) {
             _event_histo_1d->h_PiAbs_chxbac_pmult->Fill(t->reco_daughter_allTrack_ID->size());
             int Nreco_pion0=0;
             temp_plen = -999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;
             for(unsigned int tmk=0; tmk<t->reco_daughter_allTrack_len->size(); tmk++){
                 _event_histo_1d->h_PiAbs_chxbac_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_chxbac_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_chxbac_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                 _event_histo_1d->h_PiAbs_chxbac_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 



                 if(t->reco_daughter_allTrack_ID->size()>0){
                      _event_histo_1d->h_chxbac_Pmissing->Fill(reco_Pmissing);
                      _event_histo_1d->h_chxbac_Emissing->Fill(Emissing_Qsubtracted_reco);
                      _event_histo_1d->h_chxbac_Ptmissing->Fill(reco_Ptmissing);  
                      _event_histo_1d->h_chxbac_Plongit->Fill(reco_Pz);
                 }
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk) == 22)){
                     Nreco_pion0++;
                 }
                 //get the most energetic proton momentum and angle
                 if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=t->reco_daughter_allTrack_momByRange_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles

             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "chxbac", fname_geant_pm1, wgts_geant_pm1);
             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "chxbac", fname_geant_pm1, wgts_geant_pm1);




             _event_histo_1d->h_reco_photon_chx->Fill(Nreco_pion0);
 

        }
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
                 _event_histo_1d->h_PiAbs_reabac_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk) == 22)){
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
                                temp_pmom=t->reco_daughter_allTrack_momByRange_proton->at(tmk);
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles

             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pmom, event_weight, hmap_trkmom_geant_pm1_bs, "reabac", fname_geant_pm1, wgts_geant_pm1);
             if(!isdata && _fill_bootstrap_geant) FillBootstrap(temp_pcostheta, event_weight, hmap_trkcostheta_geant_pm1_bs, "reabac", fname_geant_pm1, wgts_geant_pm1);



             _event_histo_1d->h_reco_pionpm_rea->Fill(Nreco_pionpm);
             _event_histo_1d->h_reco_photon_rea->Fill(Nreco_pion0_rea);
        }
        else if(isdata==0)  { 
             std::cout<<"this is nonpion beam event"<<std::endl;
             for(unsigned int tmk=0; tmk<t->reco_daughter_PFP_true_byHits_PDG->size(); tmk++){
                 _event_histo_1d->h_PiAbs_other_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_other_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_other_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));
                 _event_histo_1d->h_PiAbs_other_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 
             } 
             if(t->reco_daughter_allTrack_ID->size()>0){
                      _event_histo_1d->h_other_Pmissing->Fill(reco_Pmissing);
                      _event_histo_1d->h_other_Emissing->Fill(Emissing_Qsubtracted_reco);
                      _event_histo_1d->h_other_Ptmissing->Fill(reco_Ptmissing);  
                      _event_histo_1d->h_other_Plongit->Fill(reco_Pz);
             }
 
        }
        //-----------------------------------------------------------------------------------
        //------------------------------------------------------------
        //std::cout<<"end of checking signal background"<<std::endl;

	//-----------------------------------------------------------  

       

  } //end of loop over all the events


	



  //=========================================================================
  LOG_NORMAL() << "Saving 1D Event Histo." << std::endl;
  
  file_out->WriteObject(_event_histo_1d, "UBXSecEventHisto1D");

  file_out->WriteObject(&hmap_trkmom_geant_pm1_bs, "hmap_trkmom_geant_pm1_bs");

  file_out->WriteObject(&bs_geant_pm1_eff_mumom_num, "bs_geant_pm1_eff_mumom_num");
  file_out->WriteObject(&bs_geant_pm1_eff_mumom_den, "bs_geant_pm1_eff_mumom_den");

  file_out->WriteObject(&hmap_trkcostheta_geant_pm1_bs, "hmap_trkcostheta_geant_pm1_bs");

  file_out->WriteObject(&bs_geant_pm1_eff_mucostheta_num, "bs_geant_pm1_eff_mucostheta_num");
  file_out->WriteObject(&bs_geant_pm1_eff_mucostheta_den, "bs_geant_pm1_eff_mucostheta_den");


  file_out->WriteObject(&bs_geant_pm1_true_reco_mom, "bs_geant_pm1_true_reco_mom");
  file_out->WriteObject(&bs_geant_pm1_true_reco_costheta, "bs_geant_pm1_true_reco_costheta");

  LOG_NORMAL() << "1D Event Histo saved." << std::endl;

 


  




  file_out->Write();

  LOG_NORMAL() << "All saved." << std::endl;
  
  file_out->Close();

  LOG_NORMAL() << "Output file closed." << std::endl;
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
