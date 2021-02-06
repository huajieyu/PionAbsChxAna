#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3F.h>
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

   TFile f2("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
   TH3F *RecoFwd_Displacement_Z_Neg = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Neg");
   TH3F *RecoFwd_Displacement_Z_Pos = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Pos");

   const int nwires_in_slice = 20;
   const int nslices = 480/nwires_in_slice;
   Int_t nbinse=12; 
   Int_t nbinsthickness = 100;

   double NA=6.02214076e23;
   double MAr=35.95; //gmol
   double Density = 1.39; // g/cm^3

   TFile f("hist.root","recreate");
   double interaction[nslices];
   double true_interaction[nslices];
   double true_abs[nslices];
   double signal[nslices];
   double incident[nslices];
   double true_incident[nslices];
   TH1D *incE[nslices];
   TH1D *pitch[nslices];
   for (int i = 0; i<nslices; ++i){
     interaction[i] = 0;
     true_interaction[i] = 0;
     true_abs[i] = 0;
     signal[i] = 0;
     incident[i] = 0;
     true_incident[i] = 0;
     incE[i] = new TH1D(Form("incE_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     pitch[i] = new TH1D(Form("pitch_%d",i),Form("Slice thickness, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
   }

   TH1D *dslcID = new TH1D("dslcID","reco slice ID - true slice ID",20,-10,10);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(abs(true_beam_PDG) != 211) continue;
      //if(abs(true_beam_PDG) != 211 && abs(true_beam_PDG) !=13) continue;
      if(!manual_beamPos_mc(reco_beam_startX, reco_beam_startY, reco_beam_startZ, 
                            reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ,
                            true_beam_startDirX, true_beam_startDirY, true_beam_startDirZ,
                            true_beam_startX, true_beam_startY, true_beam_startZ)) continue;
      //if(!endAPA3(reco_beam_endZ) )continue;
      //std::cout<<*true_beam_endProcess<<" "<<true_beam_endX<<" "<<true_beam_endY<<" "<<true_beam_endZ<<std::endl;

      int true_sliceID = -1;

      if ((*true_beam_endProcess) == "pi+Inelastic" && abs(true_beam_PDG) == 211){  //std::cout<<"signal"<<std::endl;
        double true_endz = true_beam_endZ;
        if (true_beam_endX>0){
          true_endz += RecoFwd_Displacement_Z_Pos->GetBinContent(RecoFwd_Displacement_Z_Pos->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
        }
        else{
          true_endz += RecoFwd_Displacement_Z_Neg->GetBinContent(RecoFwd_Displacement_Z_Neg->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
        }
        true_sliceID = int((true_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); //z=0.56 is z coordinate for wire 0
        if (true_sliceID <0) true_sliceID = 0;
        if (true_sliceID < nslices){
          ++true_interaction[true_sliceID];
          if (true_daughter_nPi0 == 0 && true_daughter_nPiPlus == 0){//true absorption
            ++true_abs[true_sliceID];
          }
        }
        for (int i = 0; i<=true_sliceID; ++i){
          if (i<nslices) ++true_incident[i];
        }
      }

      //interaction slice ID based on the track end wire number
      int sliceID = reco_beam_calo_wire->back()/nwires_in_slice;
      //Add code to determine if it is absorption
      bool isAbs = false;

      if (sliceID>=nslices) continue;

      if (true_sliceID!=-1){
        dslcID->Fill(sliceID - true_sliceID);
        if (sliceID == true_sliceID){
          ++signal[sliceID];
        }
      }

      //increment interaction counter
      ++interaction[sliceID];
      //increment incident counter
      for (int i = 0; i<=sliceID; ++i){
        ++incident[i];
      }

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
   double eff_int[nslices];
   double eff_inc[nslices];
   double pur[nslices];
   double xs[nslices];
   double xs_abs[nslices];
   for (int i = 0; i<nslices; ++i){
     slcid[i] = i;
     avg_incE[i] = incE[i]->GetMean();
     avg_pitch[i] = pitch[i]->GetMean();
     if (interaction[i]){
       eff_int[i] = signal[i]/true_interaction[i];
       eff_inc[i] = true_incident[i]/incident[i];
       pur[i] = signal[i]/interaction[i];
     }
     if (avg_pitch[i]&&true_incident[i]){
       xs[i] = MAr/(Density*NA*avg_pitch[i])*true_interaction[i]/true_incident[i]*1e27;
       xs_abs[i] = MAr/(Density*NA*avg_pitch[i])*true_abs[i]/true_incident[i]*1e27;
     }
     else xs[i] = 0;
     std::cout<<MAr<<" "<<Density<<" "<<NA<<" "<<avg_pitch[i]<<" "<<signal[i]<<" "<<incident[i]<<" "<<xs[i]<<std::endl;
   }

   TGraph *gr_int_slc = new TGraph(nslices, slcid, interaction);
   TGraph *gr_trueint_slc = new TGraph(nslices, slcid, true_interaction);
   TGraph *gr_inc_slc = new TGraph(nslices, slcid, incident);
   TGraph *gr_incE_slc = new TGraph(nslices, slcid, avg_incE);
   TGraph *gr_pitch_slc = new TGraph(nslices, slcid, avg_pitch);
   TGraph *gr_effint_slc = new TGraph(nslices, slcid, eff_int);
   TGraph *gr_effinc_slc = new TGraph(nslices, slcid, eff_inc);
   TGraph *gr_pur_slc = new TGraph(nslices, slcid, pur);
   TGraph *gr_xs_incE = new TGraph(nslices, avg_incE, xs);
   TGraph *gr_xsabs_incE = new TGraph(nslices, avg_incE, xs_abs);
   f.Write();
   gr_int_slc->Write("gr_int_slc");
   gr_trueint_slc->Write("gr_trueint_slc");
   gr_inc_slc->Write("gr_inc_slc");
   gr_incE_slc->Write("gr_incE_slc");
   gr_pitch_slc->Write("gr_pitch_slc");
   gr_effint_slc->Write("gr_effint_slc");
   gr_effinc_slc->Write("gr_effinc_slc");
   gr_pur_slc->Write("gr_pur_slc");
   gr_xs_incE->Write("gr_xs_incE");
   gr_xsabs_incE->Write("gr_xsabs_incE");
   f.Close();

}
