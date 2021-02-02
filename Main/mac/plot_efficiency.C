void plot_efficiency(){


gROOT->ForceStyle();
gROOT->SetBatch(1);
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
//gStyle->SetLabelOffset(1.2);
//gStyle->SetLabelFont(72);
gStyle->SetTitleSize(0.05, "X");
gStyle->SetTitleSize(0.05, "Y");
gStyle->SetTitleOffset(0.85, "X");
gStyle->SetTitleOffset(0.85, "Y");

/*TLatex* prelim = new TLatex(0.9,0.93, "protoDUNE-SP: Preliminary");
prelim->SetTextFont(62);
prelim->SetTextColor(kGray+2);
prelim->SetNDC();
prelim->SetTextSize(1/30.);
prelim->SetTextAlign(32);
*/

gROOT->SetBatch(1);

TFile *inputfile = new TFile("output_test.root");
TFile *inputfile_data = new TFile("output_test_data.root");

TH1D *h_PiAbs_gen_sig_energeticproton_mom;
h_PiAbs_gen_sig_energeticproton_mom = (TH1D*)inputfile->Get("h_PiAbs_gen_sig_energeticproton_mom");
h_PiAbs_gen_sig_energeticproton_mom->Rebin(2);


TH1D *h_PiAbs_gen_sig_energeticproton_costheta;
h_PiAbs_gen_sig_energeticproton_costheta = (TH1D*)inputfile->Get("h_PiAbs_gen_sig_energeticproton_costheta");
h_PiAbs_gen_sig_energeticproton_costheta->Rebin(2);


TH1D *h_thetax_gentruep;
h_thetax_gentruep = (TH1D*)inputfile->Get("h_thetax_gentruep");

TH1D *h_thetax_recotruep_test;
h_thetax_recotruep_test = (TH1D*)inputfile->Get("h_thetax_recotruep_test");


TLatex* prelim = new TLatex(0.9,0.93, "protoDUNE-SP: Preliminary");
prelim->SetTextFont(62);
prelim->SetTextColor(kGray+2);
prelim->SetNDC();
prelim->SetTextSize(1/30.);
prelim->SetTextAlign(32);

/*
TCanvas *c_thetax=new TCanvas("c_thetax", "c_thetax", 900, 700);

TEfficiency *pEff_thetax = new TEfficiency(*h_thetax_recotruep_test, *h_thetax_gentruep);
pEff_thetax->SetLineColor(kGreen+1); 
pEff_thetax->SetMarkerColor(kGreen+1);
pEff_thetax->SetLineWidth(2);
pEff_thetax->SetMarkerStyle(20);
pEff_thetax->SetMarkerSize(0.5);
pEff_thetax->SetTitle(";True #theta_{x};No. of Tracks");
pEff_thetax->Draw();
c_thetax->SaveAs("h_efficiency_thetax_energetic.png");
*/





/*TCanvas *c_denvsnum_thetax = new TCanvas("c_denvsnum_thetax", "c_denvsnum_thetax", 900, 700);
h_PiAbs_gen_sig_energeticproton_thetax->SetFillStyle(3002);
h_PiAbs_gen_sig_energeticproton_thetax->SetFillColor(kGreen+1);

h_PiAbs_gen_sig_energeticproton_thetax->SetTitle(";True #theta_{x} of Leading Proton[Rad]; No. of Events");
h_PiAbs_sig_energeticproton_thetax->SetFillStyle(3001);
h_PiAbs_sig_energeticproton_thetax->SetFillColor(kGreen+1);

h_PiAbs_gen_sig_energeticproton_thetax->Draw("hist");
h_PiAbs_sig_energeticproton_thetax->Draw("hist,same");
*/

  TH1D *h_PiAbs_sig_energeticproton_reco_mom;
  TH1D *h_PiAbs_sig_energeticproton_reco_costheta;
  TH1D *h_PiAbs_sig_energeticproton_reco_phi;

  TH1D *h_PiAbs_chxbac_energeticproton_mom;
  TH1D *h_PiAbs_chxbac_energeticproton_costheta;
  TH1D *h_PiAbs_chxbac_energeticproton_phi;

  TH1D *h_PiAbs_reabac_energeticproton_mom;
  TH1D *h_PiAbs_reabac_energeticproton_costheta;
  TH1D *h_PiAbs_reabac_energeticproton_phi;

  TH1D *h_PiAbs_seldata_energeticproton_mom;
  TH1D *h_PiAbs_seldata_energeticproton_costheta;
  TH1D *h_PiAbs_seldata_energeticproton_phi;



  //===================================================================
  //===================================================================
  //===============================================================
  //===================================================================

TH1D* h_track_trackscore_data;
h_track_trackscore_data=(TH1D*)inputfile_data->Get("h_trackscore_all");
 
TH1D* h_track_trackscore[7];

h_track_trackscore[0] = (TH1D*)inputfile->Get("h_trackscore_pionpm");
h_track_trackscore[1] = (TH1D*)inputfile->Get("h_trackscore_pion0");
h_track_trackscore[2] = (TH1D*)inputfile->Get("h_trackscore_proton");
h_track_trackscore[3] = (TH1D*)inputfile->Get("h_trackscore_electron");
h_track_trackscore[4] = (TH1D*)inputfile->Get("h_trackscore_muon");
h_track_trackscore[5] = (TH1D*)inputfile->Get("h_trackscore_other");
h_track_trackscore[6] = (TH1D*)inputfile->Get("h_trackscore_photon");

Double_t normfac=0.0;
for(int i=0; i<7; i++){
   normfac += h_track_trackscore[i]->Integral();
}
std::cout<<"normfac test 0 = "<<normfac<<std::endl;


normfac= normfac/h_track_trackscore_data->Integral();

std::cout<<"entries of data is "<<h_track_trackscore_data->Integral()<<std::endl;
std::cout<<"normfac test 0 = "<<normfac<<std::endl;
h_track_trackscore_data->Scale(normfac);

std::cout<<"total scaled data entries is "<<h_track_trackscore_data->Integral()<<std::endl;

TCanvas *cc= new TCanvas("cc", "cc", 900, 700);
THStack *hs_new = new THStack("hs_new", "");

h_track_trackscore[0]->SetLineColor(kBlue);
h_track_trackscore[0]->SetFillColor(kBlue);
h_track_trackscore[1]->SetLineColor(kGreen);
h_track_trackscore[1]->SetFillColor(kGreen);
h_track_trackscore[2]->SetLineColor(kMagenta);
h_track_trackscore[2]->SetFillColor(kMagenta);
h_track_trackscore[3]->SetLineColor(kCyan);
h_track_trackscore[3]->SetFillColor(kCyan);
h_track_trackscore[4]->SetLineColor(kRed);
h_track_trackscore[4]->SetFillColor(kRed);
h_track_trackscore[5]->SetLineColor(kYellow);
h_track_trackscore[5]->SetFillColor(kYellow);
h_track_trackscore[6]->SetLineColor(kGreen+3);
h_track_trackscore[6]->SetFillColor(kGreen+3);


TLegend *leg=new TLegend(0.35, 0.65, 0.65, 0.85);
leg->SetBorderSize(0);
leg->SetNColumns(2);
leg->AddEntry(h_track_trackscore[0], "#pi^{#pm}");
//leg->AddEntry(h_track_trackscore[1], "#pi0");
leg->AddEntry(h_track_trackscore[2], "proton");
leg->AddEntry(h_track_trackscore[3], "electron");
leg->AddEntry(h_track_trackscore[4], "muon");
leg->AddEntry(h_track_trackscore[5], "other");
leg->AddEntry(h_track_trackscore[6], "photon");

hs_new->Add(h_track_trackscore[2]);
hs_new->Add(h_track_trackscore[0]);
hs_new->Add(h_track_trackscore[3]);
hs_new->Add(h_track_trackscore[4]);
hs_new->Add(h_track_trackscore[5]);
hs_new->Add(h_track_trackscore[6]);
hs_new->Add(h_track_trackscore[1]);

hs_new->Draw("hist");
h_track_trackscore_data->SetLineWidth(3);
h_track_trackscore_data->Draw("same");
leg->AddEntry(h_track_trackscore_data, "Data");
leg->Draw("same");
prelim->Draw("same");
hs_new->GetXaxis()->SetTitle("Track Score of Track");
hs_new->GetYaxis()->SetTitle("No. of Tracks");
hs_new->SetTitle("");
gPad->Modified();

cc->SaveAs("h_trackscore_data+mc.png");
//==================================================================

TCanvas *cc_trkscr =new TCanvas("cc_trkscr", "cc_trkscr", 900, 700);
TH1D *h0=(TH1D*)inputfile->Get("h_trackscore_pionpm");
TH1D *h1=(TH1D*)inputfile->Get("h_trackscore_pion0");
TH1D *h2=(TH1D*)inputfile->Get("h_trackscore_proton");
TH1D *h3=(TH1D*)inputfile->Get("h_trackscore_electron");
TH1D *h4=(TH1D*)inputfile->Get("h_trackscore_muon");
TH1D *h5=(TH1D*)inputfile->Get("h_trackscore_other");
TH1D *h6=(TH1D*)inputfile->Get("h_trackscore_photon");
h0->Rebin(5); h1->Rebin(5);
h2->Rebin(5); h3->Rebin(5);
h4->Rebin(5); h5->Rebin(5);
h6->Rebin(5);

h0->Add(h1);
h0->Add(h2);
//h0->Add(h3);
h0->Add(h4);
h0->Add(h5);
h0->Add(h6);
h6->Add(h3);
h6->Divide(h0);
h6->SetLineColor(kRed);
h6->SetLineWidth(3);
h6->GetXaxis()->SetTitle("Track Score");
h6->GetYaxis()->SetTitle("Fraction of #gamma,e");
h6->SetTitle("");
h6->SetFillColor(kWhite);
h6->Draw("hist");
cc_trkscr->SaveAs("h_trackScore_ratio.png");

//=====================================================================
TH1D *h_tmdqdx;
h_tmdqdx=(TH1D*)inputfile->Get("h_tmdqdx");

TH1D *h_tmdqdx_data;
h_tmdqdx_data=(TH1D*)inputfile_data->Get("h_tmdqdx");

TCanvas *c_tmdqdx = new TCanvas("c_tmdqdx", "c_tmdqdx", 900, 700);
h_tmdqdx_data->Scale(h_tmdqdx->Integral()/h_tmdqdx_data->Integral());
h_tmdqdx_data->SetLineWidth(3);
h_tmdqdx->SetLineColor(kRed);
h_tmdqdx->SetLineWidth(3);
h_tmdqdx_data->SetTitle(";Truncated Mean dQdx[ADC/cm]; No. of Tracks");
h_tmdqdx_data->Draw();
h_tmdqdx->Draw("same");

TLegend *leg_dm=new TLegend(0.7, 0.7, 0.85, 0.85);
leg_dm->AddEntry(h_tmdqdx, "MC");
leg_dm->AddEntry(h_tmdqdx_data, "Data");
leg_dm->Draw("same");
prelim->Draw("same");
c_tmdqdx->SaveAs("h_tmdqdx_data_vs_MC.png");


TH1D *h_proton_chi2;
TH1D *h_pion_chi2;
TH1D *h_muon_chi2;
TH1D *h_electron_chi2;
TH1D *h_kaon_chi2;
TH1D *h_other_chi2;



h_proton_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_proton");
h_proton_chi2 ->GetXaxis()->SetTitle("#chi^{2}/dof");
h_proton_chi2 ->GetYaxis()->SetTitle("No. of Tracks");
h_proton_chi2->SetFillColor(kMagenta);
h_proton_chi2->SetLineWidth(4);
h_proton_chi2->SetTitle("");

h_pion_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_pionpm");
h_pion_chi2 ->SetFillColor(kBlue);
h_pion_chi2 ->SetLineWidth(4);

h_muon_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_muon");
h_muon_chi2 -> SetFillColor(kRed);
h_electron_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_electron");
h_electron_chi2->SetFillColor(kCyan);
h_kaon_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_kaon");
h_kaon_chi2->SetFillColor(kGreen+2);
h_other_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_other");
h_other_chi2->SetFillColor(kYellow);


h_proton_chi2->Rebin(2);
h_pion_chi2->Rebin(2);
h_muon_chi2->Rebin(2);
h_electron_chi2->Rebin(2);
h_kaon_chi2->Rebin(2);
h_other_chi2->Rebin(2);

TH1D *h_chi2_phypo_data;
h_chi2_phypo_data=(TH1D*)inputfile_data->Get("h_chi2_phypo_mc");

h_chi2_phypo_data->Rebin(2);

TCanvas *cchi2 = new TCanvas("cchi2", "cchi2", 900, 600);
THStack *hs_chi2 = new THStack("hs_chi2", "");


hs_chi2->Add(h_proton_chi2);
hs_chi2->Add(h_pion_chi2);
hs_chi2->Add(h_muon_chi2);
hs_chi2->Add(h_electron_chi2);
hs_chi2->Add(h_kaon_chi2);
hs_chi2->Add(h_other_chi2);
hs_chi2->Draw("hist");

hs_chi2->GetXaxis()->SetTitle("#chi^{2}/ndof");
hs_chi2->GetYaxis()->SetTitle("No. of Tracks");
hs_chi2->SetTitle("");
gPad->Modified();


h_chi2_phypo_data->Scale((h_proton_chi2->Integral()+h_pion_chi2->Integral()+
                h_muon_chi2->Integral()+h_electron_chi2->Integral()+h_kaon_chi2->Integral()+h_other_chi2->Integral()

)/h_chi2_phypo_data->Integral());
h_chi2_phypo_data->SetLineWidth(3);
h_chi2_phypo_data->Draw("same");

TLegend *leg2=new TLegend(0.55, 0.45, 0.85, 0.85);
leg2->SetBorderSize(0);
leg2->SetNColumns(1);
leg2->AddEntry(h_proton_chi2, "Proton");
leg2->AddEntry(h_pion_chi2, "#pi^{+}/#pi^{-}");
leg2->AddEntry(h_muon_chi2, "Muon");
leg2->AddEntry(h_electron_chi2, "Electron");
leg2->AddEntry(h_kaon_chi2, "Kaon");
leg2->AddEntry(h_other_chi2, "Other");
leg2->AddEntry(h_chi2_phypo_data, "Data");
leg2->Draw("same");
prelim->Draw("same");
cchi2->SaveAs("chi2_data_vs_mc.png");


} 
