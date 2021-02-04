#ifndef __MAIN_PIABSANALYSIS_CXX__
#define __MAIN_PIABSANALYSIS_CXX__

#include "PiAbsanalysis.h"
namespace Main{
void PiAbsanalysis::SetBNBCosmicFile(std::string f) {
   
    mc_bnbcosmic_file_name = f;

}

void PiAbsanalysis::SetBNBONFile(std::string f) {
   
    bnbon_file_name = f;

}
void PiAbsanalysis::stackHists(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac){
  histarray_data[0]->SetLineColor(kBlack);
  histarray_data[0]->SetLineWidth(2);
  histarray_data[0]->SetLineStyle(1);

  histarray_sig[0]-> SetFillColor(kRed); 
  histarray_sig[0]-> SetLineWidth(1); 
  histarray_sig[0]->Scale(normfac);

  histarray_bac[0]-> SetFillColor(kBlue+1);
  histarray_bac[0]-> SetLineWidth(1);
  histarray_bac[0]->Scale(normfac);

  histarray_bac[1]-> SetFillColor(kGreen+2);
  histarray_bac[1]-> SetLineWidth(1);
  histarray_bac[1]->Scale(normfac);

  histarray_bac[2]-> SetFillColor(kOrange+1);
  histarray_bac[2]-> SetLineWidth(1);
  histarray_bac[2]->Scale(normfac);

  //stack->Add(histarray_sig[0]); // signal
  stack->Add(histarray_bac[0]); // 
  stack->Add(histarray_bac[1]); // 
  stack->Add(histarray_bac[2]); // 
  stack->Add(histarray_sig[0]); // signal
}

void PiAbsanalysis::subtractBKG(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac){
  histarray_data[0]->SetLineColor(kBlack);
  histarray_data[0]->SetLineWidth(2);
  histarray_data[0]->SetLineStyle(1);

  histarray_sig[0]-> SetFillColor(2); 
  histarray_sig[0]-> SetLineWidth(1); 
  histarray_sig[0]->Scale(normfac);

  histarray_bac[0]-> SetFillColor(kBlue+1);
  histarray_bac[0]-> SetLineWidth(1);
  histarray_bac[0]->Scale(normfac);

  histarray_bac[1]-> SetFillColor(kGreen+2);
  histarray_bac[1]-> SetLineWidth(1);
  histarray_bac[1]->Scale(normfac);

  histarray_bac[2]-> SetFillColor(kOrange+1);
  histarray_bac[2]-> SetLineWidth(1);
  histarray_bac[2]->Scale(normfac);

  histarray_data[0]->Add(histarray_bac[0], -1);
  histarray_data[0]->Add(histarray_bac[1], -1);
  histarray_data[0]->Add(histarray_bac[2], -1);
  //(interacted data-interacted MCbkg) / incident data
  histarray_data[0]->Divide(histarray_data[1]);

  stack->Add(histarray_data[0]); // 
  
}






void PiAbsanalysis::DoPiAbsAnalyze(){

  clock_t begin = clock();

  std::string analyser_outdir = std::getenv("MYSW_OUTDIR");
  analyser_outdir += "output_data_mc/";
  std::string mkdir_command = "mkdir -p ";
  mkdir_command += analyser_outdir;
  system(mkdir_command.c_str());

  LOG_NORMAL() << "Created output folder with name " << analyser_outdir << "." << std::endl;
  
  //gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->SetBatch(1);

  // *************************************
  // Opening files
  // *************************************
  TFile* mc_bnbcosmic_file = TFile::Open(mc_bnbcosmic_file_name.c_str(), "READ");
  TFile* bnbon_file = TFile::Open(bnbon_file_name.c_str(), "READ");

  LOG_NORMAL()<<"Test of the script for xsec calculation with slicing method"<<std::endl;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  TH1D *h_Evttot_mc = (TH1D*)mc_bnbcosmic_file->Get("h_Evttot");
  TH1D *h_Evttot_data = (TH1D*)bnbon_file->Get("h_Evttot");

  Double_t Nevtmc= h_Evttot_mc->GetBinContent(1);
  Double_t Nevtdata = h_Evttot_data->GetBinContent(1);

  Double_t norm_fac = Nevtdata/(Nevtmc*1.0);

  TH1D *htotinc_reco_beamz;
  htotinc_reco_beamz = (TH1D*)mc_bnbcosmic_file->Get("htotinc_reco_beamz");
  TH1D *htotinc_reco_beamwire;
  htotinc_reco_beamwire = (TH1D*)mc_bnbcosmic_file->Get("htotinc_reco_beamwire");
  

  TH1D *htotinc_reco_beamz_reso;
  htotinc_reco_beamz_reso = (TH1D*)mc_bnbcosmic_file->Get("htotinc_reco_beamz_reso");
 
  TH1D *htotinc_reco_beamwire_reso;
  htotinc_reco_beamwire_reso = (TH1D*)mc_bnbcosmic_file->Get("htotinc_reco_beamwire_reso");

  
  TH1D* htotinc_true_beam_trajz[25];
  for(int i=0; i<25; i++){
       htotinc_true_beam_trajz[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("htotinc_true_beam_trajz_%d",i));
  }
  TH1D* htotinc_reco_beam_endz[25];
  for(int i=0; i<25; i++){
       htotinc_reco_beam_endz[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("htotinc_reco_beam_endz_%d",i));
  }
  
  std::cout<<"libo test -1"<<std::endl;
  
  TH1D* hgen_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hgen_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hgen_reco_beam_incidentEnergy_%d",i));
  }
  TH1D* hsel_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hsel_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hsel_reco_beam_incidentEnergy_%d",i));
  }
  std::cout<<"libo test 0"<<std::endl;
  TH1D* hsig_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hsig_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hsig_reco_beam_incidentEnergy_%d",i));
  }
  TH1D* hchxbac_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hchxbac_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hchxbac_reco_beam_incidentEnergy_%d",i));
  }
  TH1D* hreabac_reco_beam_incidentEnergy[25];

  for(int i=0; i<25; i++){
       hreabac_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hreabac_reco_beam_incidentEnergy_%d",i));
  }
  TH1D* hotherbac_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hotherbac_reco_beam_incidentEnergy[i] = (TH1D*)mc_bnbcosmic_file->Get(Form("hotherbac_reco_beam_incidentEnergy_%d",i));
  }




  TH1D* htotinc_data_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       htotinc_data_reco_beam_incidentEnergy[i] = (TH1D*)bnbon_file->Get(Form("htotinc_reco_beam_incidentEnergy_%d",i));
  }

  TH1D* htotinc_data_reco_beam_sliceThickness[25];
  for(int i=0; i<25; i++){
       htotinc_data_reco_beam_sliceThickness[i] = (TH1D*)bnbon_file->Get(Form("htotinc_reco_beam_sliceThickness_%d",i));
  }


  TH1D* hsel_data_reco_beam_incidentEnergy[25];
  for(int i=0; i<25; i++){
       hsel_data_reco_beam_incidentEnergy[i] = (TH1D*)bnbon_file->Get(Form("hsel_reco_beam_incidentEnergy_%d",i));
  }


   
  TCanvas *cs = new TCanvas("cs", "cs", 900, 700);
 
  cs->cd(); 
  htotinc_reco_beamz->SetLineWidth(2);
  htotinc_reco_beamz->SetLineColor(kBlue);
  htotinc_reco_beamz->SetTitle("");
  htotinc_reco_beamz->GetXaxis()->SetTitle("Beam Reco End Z (cm)");
  htotinc_reco_beamz->GetYaxis()->SetTitle("Nevts");
  htotinc_reco_beamz->SetMarkerStyle(20);
  htotinc_reco_beamz->Draw("P");
  cs->SaveAs("h_beam_reco_endz.png");
  cs->Update();

  TCanvas *cs1 = new TCanvas("cs1", "cs1", 900, 700);
  htotinc_reco_beamwire->SetLineWidth(2);
  htotinc_reco_beamwire->SetLineColor(kBlue);
  htotinc_reco_beamwire->SetTitle("");
  htotinc_reco_beamwire->GetXaxis()->SetTitle("Beam Reco End Wire#");
  htotinc_reco_beamwire->GetYaxis()->SetTitle("Nevts");
  htotinc_reco_beamwire->SetMarkerStyle(20);
  int rwbinmax = htotinc_reco_beamwire->GetMaximumBin();
  double ymax = htotinc_reco_beamwire->GetBinContent(rwbinmax);
  htotinc_reco_beamwire->SetMaximum(ymax);
  htotinc_reco_beamwire->Draw("P");
  TLine *tl[26];
  for(int i=0; i<26; i++){
     tl[i] = new TLine(i*20,0 , i*20, ymax);
     tl[i] ->SetLineColor(kRed);
     tl[i] ->SetLineWidth(2);
     tl[i] ->SetLineStyle(2);
     tl[i]->Draw("same");
  }

  cs1->SaveAs("h_beam_reco_wire.png");

  
  TCanvas *cs2 = new TCanvas("cs2", "cs2", 900, 700);
  htotinc_reco_beamz_reso->SetLineWidth(2);
  htotinc_reco_beamz_reso->SetLineColor(kBlue);
  htotinc_reco_beamz_reso->SetTitle("");
  htotinc_reco_beamz_reso->GetXaxis()->SetTitle("Beam End Z reso (cm)");
  htotinc_reco_beamz_reso->GetYaxis()->SetTitle("Nevts");
  htotinc_reco_beamz_reso->SetMarkerStyle(20);
  htotinc_reco_beamz_reso->Draw("hist");
  cs2->SaveAs("h_beamz_reco_reso.png");

  TCanvas *cs3 = new TCanvas("cs3", "cs3", 900, 700);
  htotinc_reco_beamwire_reso->SetLineWidth(2);
  htotinc_reco_beamwire_reso->SetLineColor(kBlue);
  htotinc_reco_beamwire_reso->SetTitle("");
  htotinc_reco_beamwire_reso->GetXaxis()->SetTitle("Beam End wire reso");
  htotinc_reco_beamwire_reso->GetYaxis()->SetTitle("Nevts");
  htotinc_reco_beamwire_reso->SetMarkerStyle(20);
  htotinc_reco_beamwire_reso->Draw("hist");
  cs3->SaveAs("h_beamwire_reco_reso.png");


  Double_t slicebins[26];
  for(long unsigned int i=0; i<sizeof(slicebins)/sizeof(slicebins[0]); i++){
          slicebins[i]=20*i;
  }
  //TLine *tlmin;
  //TLine *tlmax;
  TCanvas *css[44];
  //double ymaxvalue[4]={500, 130, 80, 80 };
  TLegend *legendz= new TLegend(0.5, 0.7, 0.85, 0.85);
  legendz->SetBorderSize(0);
  legendz->SetFillStyle(0);

  for(int i=0; i<25; i++){
     css[i]=new TCanvas(Form("css_%d",i), Form("css_%d",i), 900, 700);

     css[i]->cd();
     htotinc_true_beam_trajz[i]->SetTitle("");
     htotinc_true_beam_trajz[i]->GetXaxis()->SetTitle("True Beam End Z (cm)");
     htotinc_true_beam_trajz[i]->GetYaxis()->SetTitle("Nevts");

     int binmax = htotinc_true_beam_trajz[i]->GetMaximumBin();
     double ymaxvalue = 3.0*htotinc_true_beam_trajz[i]->GetBinContent(binmax);
     htotinc_true_beam_trajz[i]->SetMaximum(ymaxvalue);
     htotinc_true_beam_trajz[i]->SetLineColor(kBlue);
     htotinc_true_beam_trajz[i]->SetLineWidth(2);
     htotinc_true_beam_trajz[i]->Draw("hist");     
     htotinc_reco_beam_endz[i]->SetLineColor(kRed);
     htotinc_reco_beam_endz[i]->SetLineWidth(2);
     htotinc_reco_beam_endz[i]->Draw("hist,same");

     //Draw vertical lines at each boundary
     //tlmin = new TLine(slicebins[i], 0, slicebins[i], ymaxvalue[i]);
     //tlmin->SetLineWidth(2);
     //tlmin->SetLineColor(kRed);
     //tlmin->Draw("same");

     //tlmax = new TLine(slicebins[i+1], 0, slicebins[i+1], ymaxvalue[i]);
     //tlmax->SetLineWidth(2);
     //tlmax->SetLineColor(kRed);
     //tlmax->Draw("same");
     if(i==1){
        legendz->AddEntry(htotinc_true_beam_trajz[i], "MC");
        legendz->AddEntry(htotinc_reco_beam_endz[i],"protoDUNE data");
     } 
     legendz->Draw("same");
     std::string filename = "h_true_beam_trajz"; 
     filename = filename + to_string(slicebins[i])+"_MeV_to_"+to_string(slicebins[i+1])+"_MeV";
     filename = filename + ".png";
     css[i]->SaveAs(filename.c_str());
  }

  TCanvas *css_incidente[44];
  //double ymaxvalue[4]={500, 130, 80, 80 };
 
  double Aver_incidentEnergy[25]; 
  for(int i=0; i<25; i++){
     css_incidente[i]=new TCanvas(Form("css_incidente_%d",i), Form("css_incidente_%d",i), 900, 700);

     css_incidente[i]->cd();
     htotinc_data_reco_beam_incidentEnergy[i]->SetTitle("");
     htotinc_data_reco_beam_incidentEnergy[i]->GetXaxis()->SetTitle("Incident Energy [MeV]");
     htotinc_data_reco_beam_incidentEnergy[i]->GetYaxis()->SetTitle("Nevts");

     int binmax = htotinc_data_reco_beam_incidentEnergy[i]->GetMaximumBin();
     double ymaxvalue = 2.0*htotinc_data_reco_beam_incidentEnergy[i]->GetBinContent(binmax);
     htotinc_data_reco_beam_incidentEnergy[i]->SetMaximum(ymaxvalue);
     htotinc_data_reco_beam_incidentEnergy[i]->SetLineColor(kBlue);
     htotinc_data_reco_beam_incidentEnergy[i]->SetLineWidth(2);
     htotinc_data_reco_beam_incidentEnergy[i]->Draw("hist");     

     double temp_tot=0.;
     double temp_area=0.;
     int temp_Nbins = htotinc_data_reco_beam_incidentEnergy[i]->GetNbinsX();
     for(int t=0; t<temp_Nbins; t++){
         if(htotinc_data_reco_beam_incidentEnergy[i]->GetBinContent(t)>0){
         temp_tot +=htotinc_data_reco_beam_incidentEnergy[i]->GetBinContent(t);
         temp_area +=htotinc_data_reco_beam_incidentEnergy[i]->GetBinCenter(t)*htotinc_data_reco_beam_incidentEnergy[i]->GetBinContent(t);
         }
     } 
     if(temp_tot>0){ 
        Aver_incidentEnergy[i]=temp_area/temp_tot;
     }  else {Aver_incidentEnergy[i]=0.0;}
     std::cout<<"i = "<<i<<"  "<<Aver_incidentEnergy[i]<<std::endl;

     std::string filename = "htotinc_data_reco_beam_incidentE"; 
     filename = filename + to_string(slicebins[i])+"_"+to_string(slicebins[i+1])+"_";
     filename = filename + ".png";
     css_incidente[i]->SaveAs(filename.c_str());
  }
 


  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  TCanvas *css_sele[44];
  //double ymaxvalue[4]={500, 130, 80, 80 };

  for(int i=0; i<25; i++){
     css_sele[i]=new TCanvas(Form("css_sele_%d",i), Form("css_sele_%d",i), 900, 700);

     css_sele[i]->cd();
     hsel_data_reco_beam_incidentEnergy[i]->SetTitle("");
     hsel_data_reco_beam_incidentEnergy[i]->GetXaxis()->SetTitle("Incident Energy [MeV]");
     hsel_data_reco_beam_incidentEnergy[i]->GetYaxis()->SetTitle("Nevts");

     int binmax = hsel_data_reco_beam_incidentEnergy[i]->GetMaximumBin();
     double ymaxvalue = 2.0*hsel_data_reco_beam_incidentEnergy[i]->GetBinContent(binmax);
     hsel_data_reco_beam_incidentEnergy[i]->SetMaximum(ymaxvalue);
     hsel_data_reco_beam_incidentEnergy[i]->SetLineColor(kBlue);
     hsel_data_reco_beam_incidentEnergy[i]->SetLineWidth(2);
     hsel_data_reco_beam_incidentEnergy[i]->Draw("hist");     

     std::string filename = "hsel_data_reco_beam_incidentE"; 
     filename = filename + to_string(slicebins[i])+"_"+to_string(slicebins[i+1])+"_";
     filename = filename + ".png";
     css_sele[i]->SaveAs(filename.c_str());
  }
  //--------------------------------------------------------------------------------
  TCanvas *css_thickness[44];
  //double ymaxvalue[4]={500, 130, 80, 80 };

  for(int i=0; i<25; i++){
     css_thickness[i]=new TCanvas(Form("css_thickness_%d",i), Form("css_thickness_%d",i), 900, 700);

     css_thickness[i]->cd();
     htotinc_data_reco_beam_sliceThickness[i]->SetTitle("");
     htotinc_data_reco_beam_sliceThickness[i]->GetXaxis()->SetTitle("Thickness (cm)");
     htotinc_data_reco_beam_sliceThickness[i]->GetYaxis()->SetTitle("Nevts");

     int binmax = htotinc_data_reco_beam_sliceThickness[i]->GetMaximumBin();
     double ymaxvalue = 2.0*htotinc_data_reco_beam_sliceThickness[i]->GetBinContent(binmax);
     htotinc_data_reco_beam_sliceThickness[i]->SetMaximum(ymaxvalue);
     htotinc_data_reco_beam_sliceThickness[i]->SetLineColor(kBlue);
     htotinc_data_reco_beam_sliceThickness[i]->SetLineWidth(2);
     htotinc_data_reco_beam_sliceThickness[i]->Draw("hist");     

     std::string filename = "htotinc_data_reco_beam_thickness"; 
     filename = filename + to_string(slicebins[i])+"_"+to_string(slicebins[i+1])+"_";
     filename = filename + ".png";
     css_thickness[i]->SaveAs(filename.c_str());
  }
 
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  double N_incident[25];
  double N_incident_err[25];
  double N_interact[25];
  double N_interact_err[25];
  double N_background[25];
  double N_background_err[25];
  double Aver_thickness[25];
  double sliceID[25];
  double sliceID_err[25];
  double sliceeff[25];
  double sliceeff_err[25];
  double tempxsec[25];
  for(int j=0; j<25; j++){
        N_incident[j]=0;
        N_interact[j]=0;
        N_background[j]=0.0;
        N_incident_err[j]=0;
        N_interact_err[j]=0;
        N_background_err[j]=0.0;
        Aver_thickness[j]=0.0;
        sliceID[j]=j;
        sliceID_err[j]=0;
        sliceeff[j]=0; 
        sliceeff_err[j]=0;
        tempxsec[j]=0;
  }

  

 
  TCanvas *css_xsec[44];
  //double ymaxvalue[4]={500, 130, 80, 80 };
  int NtotalData=0;
  int NtotalData_beam=0;
  for(int i=0; i<25; i++){
     css_xsec[i]=new TCanvas(Form("css_xsec_%d",i), Form("css_xsec_%d",i), 900, 700);

     css_xsec[i]->cd();

     TH1D* histdata[2];
     histdata[0]=(TH1D*)hsel_data_reco_beam_incidentEnergy[i]->Clone();
     histdata[1]=(TH1D*)htotinc_data_reco_beam_incidentEnergy[i]->Clone();
     NtotalData += histdata[0]->GetEntries();
     NtotalData_beam += histdata[1]->GetEntries();
     TH1D* histsig[1];
     histsig[0]=(TH1D*)hsig_reco_beam_incidentEnergy[i]->Clone();

     TH1D* histbac[3];
     histbac[0]=(TH1D*)hchxbac_reco_beam_incidentEnergy[i]->Clone();
     histbac[1]=(TH1D*)hreabac_reco_beam_incidentEnergy[i]->Clone();
     histbac[2]=(TH1D*)hotherbac_reco_beam_incidentEnergy[i]->Clone();
     //histdata[1]->Draw("hist");
     //histdata[0]->Draw("hist");

     N_incident[i]=histdata[1]->GetEntries();
     N_incident_err[i]=TMath::Sqrt(N_incident[i]);
     N_interact[i]=histdata[0]->GetEntries();
     N_interact_err[i]=TMath::Sqrt(N_interact[i]);

     

     N_background[i]=norm_fac*(hchxbac_reco_beam_incidentEnergy[i]->GetEntries()
                              +hreabac_reco_beam_incidentEnergy[i]->GetEntries()
                              +hotherbac_reco_beam_incidentEnergy[i]->GetEntries());

     N_background_err[i]=TMath::Sqrt(N_background[i]);

     sliceeff[i] = hsig_reco_beam_incidentEnergy[i]->GetEntries()/hgen_reco_beam_incidentEnergy[i]->GetEntries();
     sliceeff_err[i] = sliceeff[i]*TMath::Sqrt(1/hsig_reco_beam_incidentEnergy[i]->GetEntries()
                                              +1/hgen_reco_beam_incidentEnergy[i]->GetEntries());

     //std::cout<<"libo hhhh"<<std::endl; 
     THStack *hs = new THStack("hs", "");
     stackHists(hs, histsig, histbac, histdata, norm_fac);
     //subtractBKG(hs, histsig, histbac, histdata, norm_fac);
     int binmax = histdata[0]->GetMaximumBin();
     double ymaxvalue = 3.0*histdata[0]->GetBinContent(binmax);
     histdata[0]->SetMarkerStyle(20);
     histdata[0]->SetLineWidth(2);
     histdata[0]->SetMaximum(ymaxvalue);
     histdata[0]->Draw("same");

     //hs->SetMaximum(ymaxvalue);
     hs->Draw("HIST,SAME");     
     histdata[0]->Draw("same");

     //average value of each slice thickness
     

     double temp_tot=0.;
     double temp_area=0.;
     int temp_Nbins = htotinc_data_reco_beam_sliceThickness[i]->GetNbinsX();
     for(int t=0; t<temp_Nbins; t++){
         if(htotinc_data_reco_beam_sliceThickness[i]->GetBinContent(t)>0){
         temp_tot +=htotinc_data_reco_beam_sliceThickness[i]->GetBinContent(t);
         temp_area +=htotinc_data_reco_beam_sliceThickness[i]->GetBinCenter(t)*htotinc_data_reco_beam_sliceThickness[i]->GetBinContent(t);
         }
     }
     if(temp_tot>0){ 
        Aver_thickness[i]=temp_area/temp_tot;
     } else {Aver_thickness[i]=0.0;}
     std::cout<<" average of thick ness of "<<i<<" is "<<Aver_thickness[i]<<std::endl; 
     if(N_interact[i]>N_background[i]) {
        
        tempxsec[i]=MAr/(Density*NA*Aver_thickness[i])*(N_interact[i]-N_background[i])/N_incident[i]/sliceeff[i];
         
        tempxsec[i]=tempxsec[i]/1e-27;

     } else {tempxsec[i]=0;}


     TLegend *legj=new TLegend(0.6, 0.6, 0.85, 0.85);
     legj->SetBorderSize(0);
     legj->SetFillStyle(0);
     legj->AddEntry(histsig[0], "Absorption");
     legj->AddEntry(histbac[0], "Charge Exchange");
     legj->AddEntry(histbac[1], "Reaction");
     legj->AddEntry(histbac[2], "Other BKG");
     legj->AddEntry(histdata[0], "Data");
     legj->Draw("same");
     std::string filename = "h_xsec_incidentE"; 
     filename = filename + to_string(slicebins[i])+"_"+to_string(slicebins[i+1]);
     filename = filename + ".png";
     css_xsec[i]->SaveAs(filename.c_str());
   }
   
   std::cout<<"total number of beam data is "<<NtotalData_beam<<std::endl; 
   std::cout<<"total number of selected data is "<<NtotalData<<std::endl; 
   double Ntotalbkg=0.0;
   for(int i=0; i<25; i++){
       Ntotalbkg +=N_background[i];
   }
   std::cout<<"total number of background are "<<Ntotalbkg<<std::endl;
   TCanvas *cgr = new TCanvas("cgr","cgr", 900, 700);
   TGraph *gr=new TGraph(25, sliceID, tempxsec);
   gr->GetXaxis()->SetTitle("slice ID");
   gr->GetYaxis()->SetTitle("#sigma(mb)");
   gr->SetTitle("");
   gr->SetMaximum(300); 
   gr->Draw("AP*");
   cgr->SaveAs("h_temp_xsec.png");

   TCanvas *cgr1 = new TCanvas("cgr1","cgr1", 900, 700);
   TGraphAsymmErrors *gr1=new TGraphAsymmErrors(25, sliceID, N_incident, sliceID_err, sliceID_err, N_incident_err, N_incident_err);
   gr1->GetXaxis()->SetTitle("slice ID");
   gr1->GetYaxis()->SetTitle("Nincident");
   gr1->SetMaximum(4e4);
   gr1->SetTitle(""); 
   gr1->Draw("AP*");
   cgr1->SaveAs("h_temp_Nincident.png");


   TCanvas *cgr2 = new TCanvas("cgr2","cgr2", 900, 700);
   TGraphAsymmErrors *gr2=new TGraphAsymmErrors(25, sliceID, N_interact, sliceID_err, sliceID_err, N_interact_err, N_interact_err);
   gr2->GetXaxis()->SetTitle("slice ID");
   gr2->GetYaxis()->SetTitle("Ninteract)");
   //gr2->SetMaximum(100.);
   gr2->SetMinimum(0.0);
   gr2->SetTitle(""); 
   gr2->Draw("AP*");
   cgr2->SaveAs("h_temp_Ninteract.png");


   TCanvas *cgr3 = new TCanvas("cgr3","cgr3", 900, 700);
   TGraphAsymmErrors *gr3=new TGraphAsymmErrors(25, sliceID, N_background, sliceID_err, sliceID_err, N_background_err, N_background_err);
   gr3->GetXaxis()->SetTitle("slice ID");
   gr3->GetYaxis()->SetTitle("Nbackground");
   //gr3->SetMaximum(70);
   gr3->SetTitle(""); 
   gr3->Draw("AP*");
   cgr3->SaveAs("h_temp_Nbackground.png");
   std::cout<<"total size of slice ID is "<<sizeof(sliceID)/sizeof(sliceID[0])<<std::endl;
   std::cout<<"total size of slice ID is "<<sizeof(Aver_incidentEnergy)/sizeof(Aver_incidentEnergy[0])<<std::endl;
   for(int i=0; i<25; i++){
      std::cout<<"slice ID and incident Energy is "<<sliceID[i]<<"    "<<Aver_incidentEnergy[i]<<std::endl;
   }
   TCanvas *cgr4 = new TCanvas("cgr4","cgr4", 900, 700);
   TGraph *gr4=new TGraph(25, sliceID, Aver_incidentEnergy);
   gr4->GetXaxis()->SetTitle("slice ID");
   gr4->GetYaxis()->SetTitle("incident Energy");
   //gr4->SetMaximum(1500);
   gr4->SetTitle(""); 
   gr4->Draw("AP*");
   cgr4->SaveAs("h_temp_incidentEnergy.png");
 
   TCanvas *cgr5 = new TCanvas("cgr5","cgr5", 900, 700);
   TGraph *gr5=new TGraph(25, Aver_incidentEnergy, tempxsec);
   gr5->GetYaxis()->SetTitle("#sigma (mb)");
   gr5->GetXaxis()->SetTitle("incident Energy");
   gr5->SetMaximum(700);
   gr5->SetTitle(""); 
   gr5->Draw("AP*");
   cgr5->SaveAs("h_temp_xsec_vs_incidentEnergy.png");
      
   TCanvas *cgr6 = new TCanvas("cgr6","cgr6", 900, 700);
   TGraph *gr6=new TGraph(25, sliceID, Aver_thickness);
   gr6->GetXaxis()->SetTitle("slice ID");
   gr6->GetYaxis()->SetTitle("Thickness");
   gr6->SetMaximum(20);
   gr6->SetTitle(""); 
   gr6->Draw("AP*");
   cgr6->SaveAs("h_temp_thichness.png");

   TCanvas *cgr7 = new TCanvas("cgr7","cgr7", 900, 700);
   TGraphAsymmErrors *gr7 = new TGraphAsymmErrors(25, sliceID, sliceeff, sliceID_err, sliceID_err,  sliceeff_err, sliceeff_err);
   gr7->GetXaxis()->SetTitle("slice ID");
   gr7->GetYaxis()->SetTitle("Efficiency");
   gr7->SetMaximum(1.0);
   gr7->SetTitle(""); 
   gr7->Draw("AP*");
   cgr7->SaveAs("h_temp_sliceeff.png");



  // Computing time<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
            
  //rootapp->Run();
  //rootapp->Terminate(0);
                 
  return;




}
} //end of namespace Main
#endif
