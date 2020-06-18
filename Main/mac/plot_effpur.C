void plot_effpur(){


gROOT->ForceStyle();
gStyle->SetOptStat(0);
gStyle->SetTitle("");
//gStyle->SetCanvasColor(-1);
//gStyle->SetPadColor(-1);
//gStyle->SetFrameFillColor(-1);
//gStyle->SetHistFillColor(-1);
//gStyle->SetTitleFillColor(-1);
//gStyle->SetFillColor(-1);
//gStyle->SetFillStyle(4000);
gStyle->SetStatStyle(0);
gStyle->SetTitleStyle(0);
gStyle->SetCanvasBorderSize(0);
gStyle->SetFrameBorderSize(0);
gStyle->SetLegendBorderSize(0);
gStyle->SetStatBorderSize(0);
gStyle->SetTitleBorderSize(0);
gStyle->SetTitleX(0.7);
gStyle->SetTitleY(0.7);



TFile *inputfile = new TFile("output_test.root");
TFile *inputfile_data = new TFile("output_test_data.root");


TH1D *h_mom_recotruep;
TH1D *h_mom_gentruep;
TH1D *h_mom_selectedtruep;
TH1D *h_mom_recotruep_test;

TH1D *h_mom_trkscoretruep;
TH1D *h_mom_tmdqdxtruep;

h_mom_recotruep = (TH1D*)inputfile->Get("h_mom_recotruep");
h_mom_recotruep_test = (TH1D*)inputfile->Get("h_mom_recotruep_test");
h_mom_gentruep = (TH1D*)inputfile->Get("h_mom_gentruep");
h_mom_selectedtruep = (TH1D*)inputfile->Get("h_mom_selectedtruep");

h_mom_trkscoretruep = (TH1D*)inputfile->Get("h_mom_trkscoretruep");
h_mom_tmdqdxtruep = (TH1D*)inputfile->Get("h_mom_tmdqdxtruep");

TH1D *h_mom_selectedtruepipm;
TH1D *h_mom_trkscoretruepipm;
TH1D *h_mom_tmdqdxtruepipm;
h_mom_selectedtruepipm = (TH1D*)inputfile->Get("h_mom_selectedtruepipm");
h_mom_trkscoretruepipm = (TH1D*)inputfile->Get("h_mom_trkscoretruepipm");
h_mom_tmdqdxtruepipm = (TH1D*)inputfile->Get("h_mom_tmdqdxtruepipm");



TCanvas *cv=new TCanvas("cv","cv", 900,600);
  TEfficiency *pEff_nocut0 = new TEfficiency(*h_mom_recotruep_test, *h_mom_gentruep);
  pEff_nocut0->SetLineColor(kMagenta+1); 
  pEff_nocut0->SetMarkerColor(kMagenta+1);
  pEff_nocut0->SetLineWidth(2);
  pEff_nocut0->SetMarkerStyle(20);
  pEff_nocut0->SetMarkerSize(0.5);
  //pEff_nocut0->GetPaintedGraph()->GetXaxis()->SetTitleSize(1.0);
  pEff_nocut0->SetTitle(";True Momentum [GeV];No. of Tracks");
  pEff_nocut0->Draw();
   
cv->SaveAs("h_beforecut_eff_proton.png");
TH1D *h_mom_recotruepionpm;
TH1D *h_mom_gentruepionpm;

h_mom_recotruepionpm = (TH1D*)inputfile->Get("h_mom_recotruepionpm");
h_mom_gentruepionpm = (TH1D*)inputfile->Get("h_mom_gentruepionpm");

  TEfficiency *pEff_nocut1 = new TEfficiency(*h_mom_recotruepionpm, *h_mom_gentruepionpm);
  pEff_nocut1->SetLineColor(kMagenta+1); 
  pEff_nocut1->SetMarkerColor(kMagenta+1);
  pEff_nocut1->SetLineWidth(2);
  pEff_nocut1->SetMarkerStyle(20);
  pEff_nocut1->SetMarkerSize(0.5);
  pEff_nocut1->SetTitle(";True Momentum [GeV];No. of Tracks");
  pEff_nocut1->Draw();
   
cv->SaveAs("h_beforecut_eff_pionpm.png");
//----------------------------------------------------------------







//--------------------------------------------------------------
TH1D *h_mom_recotruepion0;
TH1D *h_mom_gentruepion0;

h_mom_recotruepion0 = (TH1D*)inputfile->Get("h_mom_recotruepion0");
h_mom_gentruepion0 = (TH1D*)inputfile->Get("h_mom_gentruepion0");

  TEfficiency *pEff_nocut2 = new TEfficiency(*h_mom_recotruepion0, *h_mom_gentruepion0);
  pEff_nocut2->SetLineColor(kMagenta+1); 
  pEff_nocut2->SetMarkerColor(kMagenta+1);
  pEff_nocut2->SetLineWidth(2);
  pEff_nocut2->SetMarkerStyle(20);
  pEff_nocut2->SetMarkerSize(0.5);
  pEff_nocut2->SetTitle(";True Momentum [GeV];No. of Tracks");
  pEff_nocut2->Draw();
   
cv->SaveAs("h_beforecut_eff_pion0.png");
//---------------------------------------------------------------------
TCanvas *cpp = new TCanvas("cpp", "cpp", 900, 600);
h_mom_gentruepionpm->SetLineColor(kGreen+1);
h_mom_gentruepionpm->SetFillColor(kGreen+1);
h_mom_gentruepionpm->SetFillStyle(3003);
h_mom_gentruepionpm->GetXaxis()->SetTitle("True Momentum [GeV]");
h_mom_gentruepionpm->GetYaxis()->SetTitle("No. of Tracks");
h_mom_gentruepionpm->SetTitle("");
h_mom_recotruepionpm->SetLineColor(kGreen+1);
h_mom_recotruepionpm->SetFillColor(kGreen+1);
h_mom_recotruepionpm->SetFillStyle(3001);
h_mom_recotruepionpm->SetTitle("");
h_mom_gentruepionpm->Draw("hist");
h_mom_recotruepionpm->Draw("hist,same");

TLegend *lp=new TLegend(0.4, 0.7, 0.85, 0.85);
lp->SetBorderSize(0);
lp->AddEntry(h_mom_gentruepionpm, "Generated #pi^{#pm}");
lp->AddEntry(h_mom_recotruepionpm, "Reconstructed #pi^{#pm}");
lp->Draw("same");

cpp->SaveAs("h_genvstruepipm_truemomentum.png");


TCanvas *cp0 = new TCanvas("cp0", "cp0", 900, 600);
h_mom_gentruepion0->SetLineColor(kGreen+1);
h_mom_gentruepion0->SetFillColor(kGreen+1);
h_mom_gentruepion0->SetFillStyle(3003);
h_mom_gentruepion0->GetXaxis()->SetTitle("True Momentum [GeV]");
h_mom_gentruepion0->GetYaxis()->SetTitle("No. of Tracks");
h_mom_gentruepion0->SetTitle("");
h_mom_recotruepion0->SetLineColor(kGreen+1);
h_mom_recotruepion0->SetFillColor(kGreen+1);
h_mom_recotruepion0->SetFillStyle(3001);
h_mom_recotruepion0->SetTitle("");
h_mom_gentruepion0->Draw("hist");
h_mom_recotruepion0->Draw("hist,same");

TLegend *lp0=new TLegend(0.4, 0.7, 0.85, 0.85);
lp0->SetBorderSize(0);
lp0->AddEntry(h_mom_gentruepion0, "Generated #pi^{0}");
lp0->AddEntry(h_mom_recotruepion0, "Reconstructed #pi^{0}");
lp0->Draw("same");

cp0->SaveAs("h_genvstruepi0_truemomentum.png");



  TH1D *h_recomom_selected_proton;
  TH1D *h_recomom_selected_neutron;
  TH1D *h_recomom_selected_pionpm;
  TH1D *h_recomom_selected_pion0;
  TH1D *h_recomom_selected_kaon;
  TH1D *h_recomom_selected_electron;
  TH1D *h_recomom_selected_muon;
  TH1D *h_recomom_selected_photon;
  TH1D *h_recomom_selected_other;

  h_recomom_selected_proton=(TH1D*)inputfile->Get("h_recomom_selected_proton");
  //h_recomom_selected_neutron=(TH1D*)inputfile->Get("h_recomom_selected_neutron");
  h_recomom_selected_pionpm=(TH1D*)inputfile->Get("h_recomom_selected_pionpm");
  //h_recomom_selected_pion0=(TH1D*)inputfile->Get("h_recomom_selected_pion0");
  h_recomom_selected_kaon=(TH1D*)inputfile->Get("h_recomom_selected_kaon");
  h_recomom_selected_electron=(TH1D*)inputfile->Get("h_recomom_selected_electron");
  h_recomom_selected_muon=(TH1D*)inputfile->Get("h_recomom_selected_muon");
  //h_recomom_selected_photon=(TH1D*)inputfile->Get("h_recomom_selected_photon");
  h_recomom_selected_other=(TH1D*)inputfile->Get("h_recomom_selected_other");

  TCanvas *crecomom_selected = new TCanvas("crecomom_selected", "crecomom_selected", 900, 600);
  THStack *hs_recomom_selected = new THStack("hs_pmom", "");
  h_recomom_selected_proton->SetLineColor(kMagenta);
  h_recomom_selected_proton->SetFillColor(kMagenta);
  //h_recomom_selected_neutron->SetLineColor(kBlack);
  //h_recomom_selected_neutron->SetFillColor(kBlack);
  h_recomom_selected_pionpm->SetLineColor(kBlue);
  h_recomom_selected_pionpm->SetFillColor(kBlue);
  //h_recomom_selected_pion0->SetLineColor(kGreen);
  //h_recomom_selected_pion0->SetFillColor(kGreen);

  h_recomom_selected_kaon->SetLineColor(kOrange-3);
  h_recomom_selected_kaon->SetFillColor(kOrange-3);
  h_recomom_selected_electron->SetLineColor(kCyan);
  h_recomom_selected_electron->SetFillColor(kCyan);
  h_recomom_selected_muon->SetLineColor(kRed);
  h_recomom_selected_muon->SetFillColor(kRed);
  //h_recomom_selected_photon->SetLineColor(kGreen+1);
  //h_recomom_selected_photon->SetFillColor(kGreen+1);
  h_recomom_selected_other->SetLineColor(kYellow);
  h_recomom_selected_other->SetFillColor(kYellow);



  hs_recomom_selected->Add(h_recomom_selected_proton);
  //hs_recomom_selected->Add(h_recomom_selected_neutron);
  hs_recomom_selected->Add(h_recomom_selected_pionpm);
  //hs_recomom_selected->Add(h_recomom_selected_pion0);
  hs_recomom_selected->Add(h_recomom_selected_kaon);
  hs_recomom_selected->Add(h_recomom_selected_electron);
  hs_recomom_selected->Add(h_recomom_selected_muon);
  //hs_recomom_selected->Add(h_recomom_selected_photon);
  hs_recomom_selected->Add(h_recomom_selected_other);

  hs_recomom_selected->Draw("hist");

  hs_recomom_selected->GetXaxis()->SetTitle("Track Range [cm]");
  hs_recomom_selected->GetYaxis()->SetTitle("No. of Tracks");
  gPad->Modified();
  //leg2->Draw("same");

  //crecomom_selected->SaveAs("h_recomom_selected_alltrack.png");  
  //-----------------------------------------------------------------------

  double total_entries= h_recomom_selected_proton->GetEntries()+h_recomom_selected_pionpm->GetEntries()+
                        h_recomom_selected_kaon->GetEntries()+h_recomom_selected_electron->GetEntries()+
                        h_recomom_selected_muon->GetEntries()+h_recomom_selected_other->GetEntries();

  std::cout<<"over all purity is "<<h_recomom_selected_proton->GetEntries()/total_entries<<std::endl;;   

  //===========================================================================
}
