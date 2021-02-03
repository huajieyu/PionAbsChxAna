#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <iostream>
#include <vector>

//For MC from Owen Goodwins studies
double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

double cutAPA3_Z = 226.;

bool manual_beamPos_mc(double beam_startX, double beam_startY,
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

bool endAPA3(double reco_beam_endZ){
  return(reco_beam_endZ < cutAPA3_Z);
}


void ana::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ana.C
//      root> ana t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   const int nwires_in_slice = 20;
   const int nslices = 480/nwires_in_slice;
   Int_t nbinse=12; 
   Int_t nbinsthickness = 100;

   TFile f("hist.root","recreate");
   double interaction[nslices];
   double incident[nslices];
   TH1D *incE[nslices];
   TH1D *pitch[nslices];
   for (int i = 0; i<nslices; ++i){
     interaction[i] = 0;
     incident[i] = 0;
     incE[i] = new TH1D(Form("incE_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     pitch[i] = new TH1D(Form("pitch_%d",i),Form("Slice thickness, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
   }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(abs(true_beam_PDG) != 211 && abs(true_beam_PDG) !=13) continue;
      if(!manual_beamPos_mc(reco_beam_startX, reco_beam_startY, reco_beam_startZ, 
                            reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ,
                            true_beam_startDirX, true_beam_startDirY, true_beam_startDirZ,
                            true_beam_startX, true_beam_startY, true_beam_startZ)) continue;
      if(!endAPA3(reco_beam_endZ) )continue;

      //interaction slice ID based on the track end wire number
      int sliceID = reco_beam_calo_wire->back()/nwires_in_slice;
      if (sliceID>=nslices) continue;

      //increment interaction counter
      ++interaction[sliceID];
      //increment incident counter
      for (int i = 0; i<=sliceID; ++i) ++incident[i];

      std::vector<std::vector<double>> vpitch(nslices);
      std::vector<std::vector<double>> vincE(nslices); 
      for (size_t i = 0; i<reco_beam_calo_wire->size(); ++i){
        int this_wire = (*reco_beam_calo_wire)[i];
        int this_sliceID = this_wire/nwires_in_slice;
        //ignore the last slice for pitch and incident energy calculations
        if (this_sliceID>sliceID) continue;

        double this_incE = (*reco_beam_incidentEnergies)[i];
        double this_pitch = (*reco_beam_TrkPitch)[i];
        //std::cout<<this_pitch<<std::endl;
        vpitch[this_sliceID].push_back(this_pitch);
        vincE[this_sliceID].push_back(this_incE);
      }
      for (size_t i = 0; i<vpitch.size(); ++i){
        if (!vpitch[i].empty()){
          double sum_pitch = 0;
          for (size_t j = 0; j<vpitch[i].size(); ++j){
            //std::cout<<vpitch[i][j]<<std::endl;
            sum_pitch += vpitch[i][j];
          }
          //std::cout<<sum_pitch<<" "<<vpitch[i].size()<<std::endl;
          pitch[i]->Fill(sum_pitch/vpitch[i].size()*nwires_in_slice);
        }
      }
      for (size_t i = 0; i<vincE.size(); ++i){
        if (!vincE[i].empty()){
          double sum_incE = 0;
          for (size_t j = 0; j<vincE[i].size(); ++j){
            sum_incE += vincE[i][j];
          }
          incE[i]->Fill(sum_incE/vincE[i].size());
        }
      }
   }

   double slcid[nslices];
   double avg_incE[nslices];
   double avg_pitch[nslices];
   for (int i = 0; i<nslices; ++i){
     slcid[i] = i;
     avg_incE[i] = incE[i]->GetMean();
     avg_pitch[i] = pitch[i]->GetMean();
   }

   TGraph *gr_int_slc = new TGraph(nslices, slcid, interaction);
   TGraph *gr_inc_slc = new TGraph(nslices, slcid, incident);
   TGraph *gr_incE_slc = new TGraph(nslices, slcid, avg_incE);
   TGraph *gr_pitch_slc = new TGraph(nslices, slcid, avg_pitch);
   
   f.Write();
   gr_int_slc->Write("gr_int_slc");
   gr_inc_slc->Write("gr_inc_slc");
   gr_incE_slc->Write("gr_incE_slc");
   gr_pitch_slc->Write("gr_pitch_slc");
   f.Close();

}
