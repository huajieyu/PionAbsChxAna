void stackHists(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac){
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
void subtractBKG(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac){
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


void plot_reso(){
  //loadStyle();
 
   

  int tune=1;
  int cosmicCut=1;

  TFile *input0;
  TFile *input1;
  //Avogadro constant NA
  double NA=6.02214076e23;
  double MAr=35.95; //gmol
  double Density = 1.39; // g/cm^3


  if(tune==1){
  if (cosmicCut){
    //input0 = new TFile("output_test_data.root");
    input0 = new TFile("output_graphs.root");
    input1 = new TFile("output_graphs.root");
  }
  } 

  TFile *input2;
  input2 = new TFile("/dune/app/users/jiang/geant4reweight-dev/npi0_xsec.root");
  TGraph *abs_KE = (TGraph*)input2->Get("abs_KE");

  //----------------------------------------------------------------------
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(1111);
    //gStyle->SetPadTickX(1);
    //gStyle->SetPadTickY(1);
    gStyle->SetPadColor(kWhite);
    gStyle->SetStatY(0.90);
    gStyle->SetStatX(0.90);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.2);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelFont(62,"X");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleFont(62,"X");
    gStyle->SetTitleOffset(0.85,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelFont(62,"Y");
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleFont(62,"Y");
    gStyle->SetTitleOffset(1.0,"Y");
    gStyle->SetTitleX(0.22);
    gStyle->SetTitleY(0.98);
    gStyle->SetTitleSize(0.04,"t");
    //gStyle->SetTitleTextColor(kRed);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    //gStyle->SetTitleFontSize(0);
    //gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    //gStyle->SetGridStyle();
  //----------------------------------------------------------------------
  gROOT->SetBatch(1);

   TCanvas *cgr = new TCanvas("cgr","cgr", 900, 700);
   TGraph *gr=(TGraph*)input0->Get("gr_tempxsec_slc");
   gr->GetXaxis()->SetTitle("slice ID");
   gr->GetYaxis()->SetTitle("#sigma(mb)");
   gr->SetTitle("");
   gr->SetMaximum(300); 
   gr->Draw("AP*");
   cgr->SaveAs("h_temp_xsec.png");
   
   TCanvas *cgr1 = new TCanvas("cgr1","cgr1", 900, 700);
   TGraph *gr1=(TGraph*)input0->Get("gr_inc_slc");
   gr1->GetXaxis()->SetTitle("slice ID");
   gr1->GetYaxis()->SetTitle("Nincident");
   gr1->SetMaximum(4e4);
   gr1->SetTitle(""); 
   gr1->Draw("AP*");
   cgr1->SaveAs("h_temp_Nincident.png");


   TCanvas *cgr2 = new TCanvas("cgr2","cgr2", 900, 700);
   TGraph *gr2=(TGraph*)input0->Get("gr_int_slc");
   gr2->GetXaxis()->SetTitle("slice ID");
   gr2->GetYaxis()->SetTitle("Ninteract)");
   //gr2->SetMaximum(100.);
   gr2->SetMinimum(0.0);
   gr2->SetTitle(""); 
   gr2->Draw("AP*");
   cgr2->SaveAs("h_temp_Ninteract.png");


   TCanvas *cgr3 = new TCanvas("cgr3","cgr3", 900, 700);
   TGraph *gr3=(TGraph*)input0->Get("gr_selected_bkg_slc");
   gr3->GetXaxis()->SetTitle("slice ID");
   gr3->GetYaxis()->SetTitle("Nbackground");
   //gr3->SetMaximum(70);
   gr3->SetTitle(""); 
   gr3->Draw("AP*");
   cgr3->SaveAs("h_temp_Nbackground.png");

   TCanvas *cgr4 = new TCanvas("cgr4","cgr4", 900, 700);
   TGraph *gr4=(TGraph*)input0->Get("gr_incE_slc");
   gr4->GetXaxis()->SetTitle("slice ID");
   gr4->GetYaxis()->SetTitle("incident Energy");
   //gr4->SetMaximum(1500);
   gr4->SetTitle(""); 
   gr4->Draw("AP*");
   cgr4->SaveAs("h_temp_incidentEnergy.png");
   
   TCanvas *cgr5 = new TCanvas("cgr5","cgr5", 900, 700);
   TGraph *gr5=(TGraph*)input0->Get("gr_tempxsec_vs_incE_slc");
   gr5->GetYaxis()->SetTitle("#sigma (mb)");
   gr5->GetXaxis()->SetTitle("incident Energy");
   gr5->SetMaximum(700);
   gr5->SetTitle(""); 
   gr5->Draw("AP*");
   cgr5->SaveAs("h_temp_xsec_vs_incidentEnergy.png");
      
   TCanvas *cgr6 = new TCanvas("cgr6","cgr6", 900, 700);
   TGraph *gr6=(TGraph*)input0->Get("gr_pitch_slc");
   gr6->GetXaxis()->SetTitle("slice ID");
   gr6->GetYaxis()->SetTitle("Thickness");
   gr6->SetMaximum(20);
   gr6->SetTitle(""); 
   gr6->Draw("AP*");
   cgr6->SaveAs("h_temp_thichness.png");

   TCanvas *cgr7 = new TCanvas("cgr7","cgr7", 900, 700);
   TGraph *gr7=(TGraph*)input0->Get("gr_efficiency_slc");
   gr7->GetXaxis()->SetTitle("slice ID");
   gr7->GetYaxis()->SetTitle("Efficiency");
   gr7->SetMaximum(1.0);
   gr7->SetTitle(""); 
   gr7->Draw("AP*");
   cgr7->SaveAs("h_temp_efficiency.png");

   Double_t xsec_temp[16];
   Double_t incE_temp[16];
   Double_t a[24];
   Double_t b[24];
   Double_t c[24];
   Double_t d[24];

   for(int i=0; i<gr4->GetN(); i++){
     gr4->GetPoint(i,a[i], b[i]);
     gr->GetPoint(i, c[i], d[i]);
   }
   for(int j=0; j<16; j++){
         xsec_temp[j]=d[j+4];
         incE_temp[j]=b[j+4];
   }


   TGraph *temp_xs = new TGraph(16, incE_temp, xsec_temp);

   TCanvas *cgr8 = new TCanvas("cgr8","cgr8", 900, 700);
   TMultiGraph *mg = new TMultiGraph();
   temp_xs->SetMarkerColor(kRed);
   temp_xs->SetMarkerStyle(20);
   abs_KE->SetLineWidth(2);
   mg->Add(abs_KE);
   mg->Add(temp_xs);
   TLegend *legxs = new TLegend(0.5, 0.7, 0.85, 0.85);
   legxs->SetBorderSize(0);
   legxs->SetFillStyle(0);
   legxs->AddEntry(abs_KE, "G4 Prediction");
   legxs->AddEntry(temp_xs, "Measured (MC)");
   mg->Draw("ap");
   legxs->Draw("same");
   gPad->Update();
   mg->GetXaxis()->SetTitle("Beam Kinetic Energy[MeV]");
   mg->GetYaxis()->SetTitle("#sigma(mb)");
   gPad->Modified();
   cgr8->SaveAs("h_tempxsec_vs_incE_G4.png");
 } //========================================================
