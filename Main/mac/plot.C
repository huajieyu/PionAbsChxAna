void plot(){

gROOT->SetBatch(1);
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
TFile *inputfile_nothresh = new TFile("output_test_nothresh.root");
TH1D  *h_shwlike_nhits[3];
h_shwlike_nhits[0] = (TH1D*)inputfile->Get("h_nhits_shwlike_photon");
h_shwlike_nhits[1] = (TH1D*)inputfile->Get("h_nhits_shwlike_electron");
h_shwlike_nhits[2] = (TH1D*)inputfile->Get("h_nhits_shwlike_nonphoton");

for(int i=0; i<3; i++){
  h_shwlike_nhits[i]->Rebin(2);
}
TCanvas *c_shwlike=new TCanvas("c_shwlike", "c_shwlike", 900,700);

THStack *hs_shwlike = new THStack("hs_shwlike", "");

h_shwlike_nhits[0]->SetLineColor(kGreen+1);
h_shwlike_nhits[0]->SetFillColor(kGreen+1);
h_shwlike_nhits[1]->SetLineColor(kCyan);
h_shwlike_nhits[1]->SetFillColor(kCyan);
h_shwlike_nhits[2]->SetLineColor(kMagenta);
h_shwlike_nhits[2]->SetFillColor(kMagenta);

hs_shwlike->Add(h_shwlike_nhits[0]);
hs_shwlike->Add(h_shwlike_nhits[1]);
hs_shwlike->Add(h_shwlike_nhits[2]);
TLegend* legshw=new TLegend(0.7, 0.7, 0.85, 0.85);
legshw->SetBorderSize(0);
legshw->AddEntry(h_shwlike_nhits[0], "photon");
legshw->AddEntry(h_shwlike_nhits[1], "electron");
legshw->AddEntry(h_shwlike_nhits[2], "Other");

hs_shwlike->Draw("hist");
legshw->Draw("same");
c_shwlike->SaveAs("h_nhits_showerlike.png");


TH1D* h_track_trackscore[7];

h_track_trackscore[0] = (TH1D*)inputfile->Get("h_trackscore_pionpm");
h_track_trackscore[1] = (TH1D*)inputfile->Get("h_trackscore_pion0");
h_track_trackscore[2] = (TH1D*)inputfile->Get("h_trackscore_proton");
h_track_trackscore[3] = (TH1D*)inputfile->Get("h_trackscore_electron");
h_track_trackscore[4] = (TH1D*)inputfile->Get("h_trackscore_muon");
h_track_trackscore[5] = (TH1D*)inputfile->Get("h_trackscore_other");
h_track_trackscore[6] = (TH1D*)inputfile->Get("h_trackscore_photon");


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

TH1D *h_mom_recotrueother;
TH1D *h_mom_gentrueother;

h_mom_recotrueother = (TH1D*)inputfile->Get("h_mom_recotrueother");
h_mom_gentrueother = (TH1D*)inputfile->Get("h_mom_gentrueother");

  TEfficiency *pEff_nocut11 = new TEfficiency(*h_mom_recotrueother, *h_mom_gentrueother);
  pEff_nocut11->SetLineColor(kMagenta+1); 
  pEff_nocut11->SetMarkerColor(kMagenta+1);
  pEff_nocut11->SetLineWidth(2);
  pEff_nocut11->SetMarkerStyle(20);
  pEff_nocut11->SetMarkerSize(0.5);
  pEff_nocut11->SetTitle(";True Momentum [GeV];No. of Tracks");
  pEff_nocut11->Draw();
   
cv->SaveAs("h_beforecut_eff_other.png");

//----------------------------------------------------------------
TH1D *h_true_beam_endE_num;
TH1D *h_true_beam_endE_den;

TH1D *h_true_beam_endE_num_nothresh;
TH1D *h_true_beam_endE_den_nothresh;
h_true_beam_endE_num = (TH1D*)inputfile->Get("h_true_beam_endE_num");
h_true_beam_endE_den = (TH1D*)inputfile->Get("h_true_beam_endE_den");

h_true_beam_endE_num_nothresh = (TH1D*)inputfile_nothresh->Get("h_true_beam_endE_num");
h_true_beam_endE_den_nothresh = (TH1D*)inputfile_nothresh->Get("h_true_beam_endE_den");

TCanvas *c_endEff=new TCanvas("c_endEff", "c_endEff", 900, 700);
TEfficiency *pEff_nothresh = new TEfficiency(*h_true_beam_endE_num_nothresh, *h_true_beam_endE_den_nothresh);
pEff_nothresh->SetLineColor(kRed);
pEff_nothresh->SetMarkerColor(kRed);
pEff_nothresh->SetLineWidth(2);
pEff_nothresh->SetMarkerStyle(20);
pEff_nothresh->SetMarkerSize(0.5);
pEff_nothresh->SetTitle(";True Beam End Energy[GeV];Efficiency");

TEfficiency *pEff_withthresh = new TEfficiency(*h_true_beam_endE_num, *h_true_beam_endE_den);
pEff_withthresh->SetLineColor(kBlue);
pEff_withthresh->SetMarkerColor(kBlue);
pEff_withthresh->SetLineWidth(2);
pEff_withthresh->SetMarkerStyle(20);
pEff_withthresh->SetMarkerSize(0.5);
pEff_withthresh->SetTitle(";True Beam End Energy[GeV];Efficiency");
pEff_nothresh->Draw("");

pEff_withthresh->Draw("same");

TLegend *legEff=new TLegend(0.45, 0.7, 0.85, 0.85);
legEff->AddEntry(pEff_nothresh, "Efficiency (nothresh)");
legEff->AddEntry(pEff_withthresh, "Efficiency (withthresh)");
legEff->Draw("same");

gPad->Update();
auto graph_libo = pEff_nothresh->GetPaintedGraph();
graph_libo->SetMinimum(0.0);
graph_libo->SetMaximum(1.0);
gPad->Update();

c_endEff->SaveAs("h_efficiency_true_beam_end_point_energy.png");

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
  pEff_nocut2->SetTitle(";True Momentum [GeV];Efficiency");
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




//---------------------------------------------------------------------

TCanvas *cc= new TCanvas("cc", "cc", 900, 600);
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
leg->AddEntry(h_track_trackscore[1], "#pi0");
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


hs_new->Draw("hist");
leg->Draw("same");
hs_new->GetXaxis()->SetTitle("Track Score of Track");
hs_new->GetYaxis()->SetTitle("No. of Tracks");
hs_new->SetTitle("");
gPad->Modified();

cc->SaveAs("h_trackscore.png");

//=============================================================
TH1D *h_proton_tmdqdxchi2;
TH1D *h_pion_tmdqdxchi2;
TH1D *h_muon_tmdqdxchi2;
TH1D *h_electron_tmdqdxchi2;
TH1D *h_kaon_tmdqdxchi2;
TH1D *h_other_tmdqdxchi2;
TH1D *h_mc_tmdqdxchi2;


h_proton_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_proton");
h_proton_tmdqdxchi2 ->GetXaxis()->SetTitle("#chi^{2}/dof");
h_proton_tmdqdxchi2 ->GetYaxis()->SetTitle("No. of Tracks");
h_proton_tmdqdxchi2->SetFillColor(kMagenta);
h_proton_tmdqdxchi2->SetLineWidth(4);
h_proton_tmdqdxchi2->SetTitle("");

h_pion_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_pionpm");
h_pion_tmdqdxchi2 ->SetFillColor(kBlue);
h_pion_tmdqdxchi2 ->SetLineWidth(4);

h_muon_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_muon");
h_muon_tmdqdxchi2 -> SetFillColor(kRed);
h_electron_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_electron");
h_electron_tmdqdxchi2->SetFillColor(kCyan);
h_kaon_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_kaon");

h_other_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_other");
h_other_tmdqdxchi2->SetFillColor(kYellow);

h_mc_tmdqdxchi2 = (TH1D*)inputfile->Get("htmdqdx_chi2_phypo_mc");
h_mc_tmdqdxchi2->SetLineWidth(2);


TCanvas *ctmdqdxchi2 = new TCanvas("ctmdqdxchi2", "ctmdqdxchi2", 900, 600);
THStack *hs_tmdqdxchi2 = new THStack("hs_tmdqdxchi2", "");


hs_tmdqdxchi2->Add(h_proton_tmdqdxchi2);
hs_tmdqdxchi2->Add(h_pion_tmdqdxchi2);
hs_tmdqdxchi2->Add(h_muon_tmdqdxchi2);
hs_tmdqdxchi2->Add(h_electron_tmdqdxchi2);
//hs_tmdqdxchi2->Add(h_kaon_tmdqdxchi2);
hs_tmdqdxchi2->Add(h_other_tmdqdxchi2);
hs_tmdqdxchi2->Draw("hist");

h_mc_tmdqdxchi2->Draw("same");

hs_tmdqdxchi2->GetXaxis()->SetTitle("#chi^2/dof");
hs_tmdqdxchi2->GetYaxis()->SetTitle("No. of Tracks");
//cchi2->SaveAs("h_chi2_pionvsproton.png");
gPad->Modified();
ctmdqdxchi2->SaveAs("htmdqdx_chi2_pionvsproton.png");



//================================================================
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

h_other_chi2 = (TH1D*)inputfile->Get("h_chi2_phypo_other");
h_other_chi2->SetFillColor(kYellow);

TCanvas *cchi2 = new TCanvas("cchi2", "cchi2", 900, 600);
THStack *hs_chi2 = new THStack("hs_chi2", "");


hs_chi2->Add(h_proton_chi2);
hs_chi2->Add(h_pion_chi2);
hs_chi2->Add(h_muon_chi2);
hs_chi2->Add(h_electron_chi2);
//hs_chi2->Add(h_kaon_chi2);
hs_chi2->Add(h_other_chi2);
hs_chi2->Draw("hist");

TLegend *leg2=new TLegend(0.55, 0.55, 0.85, 0.85);
leg2->SetBorderSize(0);
leg2->SetNColumns(1);
leg2->AddEntry(h_proton_chi2, "Proton");
leg2->AddEntry(h_pion_chi2, "#pi^{+}/#pi^{-}");
leg2->AddEntry(h_muon_chi2, "Muon");
leg2->AddEntry(h_electron_chi2, "Electron");
//leg2->AddEntry(h_kaon_chi2, "Kaon");
leg2->AddEntry(h_other_chi2, "Other");
leg2->Draw("same");

hs_chi2->GetXaxis()->SetTitle("#chi^2/dof");
hs_chi2->GetYaxis()->SetTitle("No. of Tracks");
//cchi2->SaveAs("h_chi2_pionvsproton.png");
gPad->Modified();
cchi2->SaveAs("h_chi2_pionvsproton.png");
//======================================================

TH1D *h_proton_nhits;
TH1D *h_pion_nhits;
TH1D *h_muon_nhits;
TH1D *h_electron_nhits;
TH1D *h_kaon_nhits;
TH1D *h_other_nhits;

TH1D *h_proton_reconhits;

h_proton_reconhits = (TH1D*)inputfile->Get("h_reconhits_proton");


h_proton_nhits = (TH1D*)inputfile->Get("h_nhits_proton");
h_proton_nhits ->GetXaxis()->SetTitle("#chi^{2}/dof");
h_proton_nhits ->GetYaxis()->SetTitle("No. of Tracks");
h_proton_nhits->SetFillColor(kMagenta);
h_proton_nhits->SetLineWidth(4);
h_proton_nhits->SetTitle("");

h_pion_nhits = (TH1D*)inputfile->Get("h_nhits_pionpm");
h_pion_nhits ->SetFillColor(kBlue);
h_pion_nhits ->SetLineWidth(4);

h_muon_nhits = (TH1D*)inputfile->Get("h_nhits_muon");
h_muon_nhits -> SetFillColor(kRed);
h_electron_nhits = (TH1D*)inputfile->Get("h_nhits_electron");
h_electron_nhits->SetFillColor(kCyan);
h_kaon_nhits = (TH1D*)inputfile->Get("h_nhits_kaon");

h_other_nhits = (TH1D*)inputfile->Get("h_nhits_other");
h_other_nhits->SetFillColor(kYellow);

TCanvas *cnhits = new TCanvas("cnhits", "cnhits", 900, 600);
THStack *hs_nhits = new THStack("hs_nhits", "");


hs_nhits->Add(h_proton_nhits);
hs_nhits->Add(h_pion_nhits);
hs_nhits->Add(h_muon_nhits);
hs_nhits->Add(h_electron_nhits);
//hs_chi2->Add(h_kaon_chi2);
hs_nhits->Add(h_other_nhits);
hs_nhits->Draw("hist");

leg2->Draw("same");

hs_nhits->GetXaxis()->SetTitle("Number of Hits");
hs_nhits->GetYaxis()->SetTitle("No. of Tracks");
//cnhits->SaveAs("h_nhits_pionvsproton.png");
gPad->Modified();
cnhits->SaveAs("h_nhits_pionvsproton.png");
//==========================================================
h_proton_reconhits->Rebin(2);
h_proton_nhits->Rebin(2);
TH1D *h_proton_tmdqdxnhits;

h_proton_tmdqdxnhits = (TH1D*)inputfile->Get("h_tmdqdxnhits_proton");
h_proton_tmdqdxnhits->Rebin(2);


TCanvas *ceffnhits = new TCanvas("ceffnhits", "ceffnhits",900,600);
TEfficiency* pEffnhits = 0;
pEffnhits = new TEfficiency(*h_proton_reconhits, *h_proton_nhits);
pEffnhits->SetLineColor(kMagenta+1);
pEffnhits->SetMarkerColor(kMagenta+1);
pEffnhits->SetMarkerStyle(2);
pEffnhits->SetMarkerSize(0.5);
pEffnhits->SetLineWidth(2);
pEffnhits->SetTitle(";Number of Hits; Efficiency");
pEffnhits->Draw();
gPad->Update();
auto graph= pEffnhits->GetPaintedGraph();
graph->SetMaximum(1.2);
graph->SetMinimum(0.0);
gPad->Update();


TEfficiency* pEffnhits2 = 0;
pEffnhits2 = new TEfficiency(*h_proton_tmdqdxnhits, *h_proton_nhits);
pEffnhits2->SetLineColor(kGreen+1);
pEffnhits2->SetMarkerColor(kGreen+1);
pEffnhits2->SetMarkerStyle(2);
pEffnhits2->SetMarkerSize(0.5);
pEffnhits2->SetLineWidth(2);
pEffnhits2->SetTitle(";Number of Hits; Efficiency");
pEffnhits2->Draw("same");


TLegend* legforeff=new TLegend(0.15, 0.7, 0.455, 0.85);
legforeff->SetBorderSize(0);
legforeff->AddEntry(pEffnhits2, "TrunMean dQdx");
legforeff->AddEntry(pEffnhits, "Chi2/ndof");
legforeff->Draw("same");

ceffnhits->SaveAs("h_eff_nhits.png");

//========================================================
TH1D *h_pion_tmdqdx;
TH1D *h_proton_tmdqdx;

h_pion_tmdqdx = (TH1D*)inputfile->Get("h_tmdqdx_pionpm");
h_proton_tmdqdx = (TH1D*)inputfile->Get("h_tmdqdx_proton");


TCanvas *cdqdx= new TCanvas("cdqdx", "cdqdx", 900, 1200);
//upper plot will be in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.0, 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetGridx();         // Vertical grid
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pa

h_proton_tmdqdx->GetXaxis()->SetTitle("Truncated Mean dQdx [ADC/cm]");
h_proton_tmdqdx->GetYaxis()->SetTitle("No. of Tracks");
h_proton_tmdqdx->SetLineWidth(4);
h_proton_tmdqdx->SetLineColor(kBlue);
h_proton_tmdqdx->SetTitle("");

h_pion_tmdqdx->GetXaxis()->SetTitle("Truncated Mean dQdx [ADC/cm]");
h_pion_tmdqdx->GetYaxis()->SetTitle("No. of Tracks");
h_pion_tmdqdx->SetLineWidth(4);
h_pion_tmdqdx->SetLineColor(kRed);
h_pion_tmdqdx->SetTitle("");


h_pion_tmdqdx->Draw("hist");
h_proton_tmdqdx->Draw("hist,same");

TLegend *legpp= new TLegend(0.7, 0.6, 0.85, 0.85);

legpp->AddEntry(h_proton_tmdqdx, "Proton");
legpp->AddEntry(h_pion_tmdqdx, "#pi^{+}/#pi^{-}");

legpp->Draw("same");
cdqdx->cd(); // Go back to the main canvas before defining pad 2

TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();       // pad2 becomes the current pa

TH1D *h1=(TH1D*)h_pion_tmdqdx->Clone();
TH1D *h2=(TH1D*)h_proton_tmdqdx->Clone();
h1->Rebin(2);
h2->Rebin(2);
h1->Add(h2);
h2->Divide(h1);
h2->GetYaxis()->SetTitle("Proton/(#pi^{#pm}+Proton)");
h2->GetYaxis()->SetTitleSize(24);
h2->GetYaxis()->SetTitleFont(43);
h2->GetYaxis()->SetTitleOffset(1.3);
h2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h2->GetYaxis()->SetLabelSize(11);
h2->GetXaxis()->SetTitleSize(0.1);
h2->GetXaxis()->SetLabelFont(43);
h2->GetXaxis()->SetLabelSize(11);
h2->Draw("hist");

cdqdx->SaveAs("h_trunmeandqdx.png");

/*

*/
TH1D *h_PiAbs_other_pmom;
TH1D *h_PiAbs_sig_pmom;
TH1D *h_PiAbs_chxbac_pmom;
TH1D *h_PiAbs_reabac_pmom;
TH1D *h_PiAbs_sel_pmom_data;

h_PiAbs_other_pmom=(TH1D*)inputfile->Get("h_PiAbs_other_pmom");

h_PiAbs_sig_pmom=(TH1D*)inputfile->Get("h_PiAbs_sig_pmom");
h_PiAbs_chxbac_pmom=(TH1D*)inputfile->Get("h_PiAbs_chxbac_pmom");
h_PiAbs_reabac_pmom=(TH1D*)inputfile->Get("h_PiAbs_reabac_pmom");
h_PiAbs_sel_pmom_data=(TH1D*)inputfile_data->Get("h_PiAbs_sel_pmom");


h_PiAbs_sig_pmom->SetFillColor(kMagenta);
h_PiAbs_sig_pmom->SetLineColor(kMagenta);
h_PiAbs_chxbac_pmom->SetFillColor(kGreen+1);
h_PiAbs_chxbac_pmom->SetLineColor(kGreen+1);
h_PiAbs_reabac_pmom->SetFillColor(kBlue);
h_PiAbs_reabac_pmom->SetLineColor(kBlue);
h_PiAbs_other_pmom->SetFillColor(kRed);
h_PiAbs_other_pmom->SetLineColor(kRed);
Double_t scale_fac=(h_PiAbs_sig_pmom->Integral()+h_PiAbs_chxbac_pmom->Integral()+h_PiAbs_reabac_pmom->Integral()+h_PiAbs_other_pmom->Integral())/h_PiAbs_sel_pmom_data->Integral();
h_PiAbs_sel_pmom_data->Scale(scale_fac);

h_PiAbs_sig_pmom->Rebin(2);
h_PiAbs_chxbac_pmom->Rebin(2);
h_PiAbs_reabac_pmom->Rebin(2);
h_PiAbs_other_pmom->Rebin(2);
h_PiAbs_sel_pmom_data->Rebin(2);

TCanvas *cpmom= new TCanvas("cpmom", "cpmom", 900, 600);

THStack *hs_pmom = new THStack("hs_pmom", "");
hs_pmom->Add(h_PiAbs_sig_pmom);
hs_pmom->Add(h_PiAbs_chxbac_pmom);
hs_pmom->Add(h_PiAbs_reabac_pmom);
hs_pmom->Add(h_PiAbs_other_pmom);
hs_pmom->Draw("hist");
h_PiAbs_sel_pmom_data->SetLineWidth(4);
h_PiAbs_sel_pmom_data->Draw("same");
hs_pmom->GetXaxis()->SetTitle("Track Momentum [GeV]");
hs_pmom->GetYaxis()->SetTitle("No. of Events");
hs_pmom->SetMaximum(100);
gPad->Modified();

TLegend *pmultleg = new TLegend(0.55, 0.6, 0.85, 0.85);
pmultleg->AddEntry(h_PiAbs_sig_pmom, "Signal");
pmultleg->AddEntry(h_PiAbs_chxbac_pmom, "BKG-Charge Exchange");
pmultleg->AddEntry(h_PiAbs_reabac_pmom, "BKG-Pion Reaction");
pmultleg->AddEntry(h_PiAbs_other_pmom, "Other"); 
pmultleg->AddEntry(h_PiAbs_sel_pmom_data, "Data");
pmultleg->Draw("same");
cpmom->SaveAs("h_proton_momentum.png");

TH1D *h_PiAbs_other_pcostheta;
TH1D *h_PiAbs_sig_pcostheta;
TH1D *h_PiAbs_chxbac_pcostheta;
TH1D *h_PiAbs_reabac_pcostheta;
TH1D *h_PiAbs_sel_pcostheta_data;

h_PiAbs_other_pcostheta=(TH1D*)inputfile->Get("h_PiAbs_other_pcostheta");
h_PiAbs_sig_pcostheta=(TH1D*)inputfile->Get("h_PiAbs_sig_pcostheta");
h_PiAbs_chxbac_pcostheta=(TH1D*)inputfile->Get("h_PiAbs_chxbac_pcostheta");
h_PiAbs_reabac_pcostheta=(TH1D*)inputfile->Get("h_PiAbs_reabac_pcostheta");
h_PiAbs_sel_pcostheta_data=(TH1D*)inputfile_data->Get("h_PiAbs_sel_pcostheta");

h_PiAbs_sel_pcostheta_data->Scale(scale_fac);

h_PiAbs_other_pcostheta->Rebin(2);
h_PiAbs_sig_pcostheta->Rebin(2);
h_PiAbs_chxbac_pcostheta->Rebin(2);
h_PiAbs_reabac_pcostheta->Rebin(2);
h_PiAbs_sel_pcostheta_data->Rebin(2);


h_PiAbs_sig_pcostheta->SetFillColor(kMagenta);
h_PiAbs_sig_pcostheta->SetLineColor(kMagenta);
h_PiAbs_chxbac_pcostheta->SetFillColor(kGreen+1);
h_PiAbs_chxbac_pcostheta->SetLineColor(kGreen+1);
h_PiAbs_reabac_pcostheta->SetFillColor(kBlue);
h_PiAbs_reabac_pcostheta->SetLineColor(kBlue);
h_PiAbs_other_pcostheta->SetFillColor(kRed);
h_PiAbs_other_pcostheta->SetLineColor(kRed);

TCanvas *cpcostheta= new TCanvas("cpcostheta", "cpcostheta", 900, 600);

THStack *hs_pcostheta = new THStack("hs_pcostheta", "");
hs_pcostheta->Add(h_PiAbs_sig_pcostheta);
hs_pcostheta->Add(h_PiAbs_chxbac_pcostheta);
hs_pcostheta->Add(h_PiAbs_reabac_pcostheta);
hs_pcostheta->Add(h_PiAbs_other_pcostheta);
hs_pcostheta->Draw("hist");
h_PiAbs_sel_pcostheta_data->SetLineWidth(4);
h_PiAbs_sel_pcostheta_data->Draw("same");

hs_pcostheta->GetXaxis()->SetTitle("cos#theta");
hs_pcostheta->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cpcostheta->SaveAs("h_proton_costheta.png");
//=====================================================================
TH1D *h_PiAbs_other_ptheta;
TH1D *h_PiAbs_sig_ptheta;
TH1D *h_PiAbs_chxbac_ptheta;
TH1D *h_PiAbs_reabac_ptheta;
TH1D *h_PiAbs_sel_ptheta_data;


h_PiAbs_other_ptheta=(TH1D*)inputfile->Get("h_PiAbs_other_ptheta");
h_PiAbs_sig_ptheta=(TH1D*)inputfile->Get("h_PiAbs_sig_ptheta");
h_PiAbs_chxbac_ptheta=(TH1D*)inputfile->Get("h_PiAbs_chxbac_ptheta");
h_PiAbs_reabac_ptheta=(TH1D*)inputfile->Get("h_PiAbs_reabac_ptheta");
h_PiAbs_sel_ptheta_data=(TH1D*)inputfile_data->Get("h_PiAbs_sel_ptheta");

h_PiAbs_sel_ptheta_data->Scale(scale_fac);

h_PiAbs_other_ptheta->Rebin(2);
h_PiAbs_sig_ptheta->Rebin(2);
h_PiAbs_chxbac_ptheta->Rebin(2);
h_PiAbs_reabac_ptheta->Rebin(2);
h_PiAbs_sel_ptheta_data->Rebin(2);

h_PiAbs_sig_ptheta->SetFillColor(kMagenta);
h_PiAbs_sig_ptheta->SetLineColor(kMagenta);
h_PiAbs_chxbac_ptheta->SetFillColor(kGreen+1);
h_PiAbs_chxbac_ptheta->SetLineColor(kGreen+1);
h_PiAbs_reabac_ptheta->SetFillColor(kBlue);
h_PiAbs_reabac_ptheta->SetLineColor(kBlue);
h_PiAbs_other_ptheta->SetFillColor(kRed);
h_PiAbs_other_ptheta->SetLineColor(kRed);


TCanvas *cptheta= new TCanvas("cptheta", "cptheta", 900, 600);

THStack *hs_ptheta = new THStack("hs_ptheta", "");
hs_ptheta->Add(h_PiAbs_sig_ptheta);
hs_ptheta->Add(h_PiAbs_chxbac_ptheta);
hs_ptheta->Add(h_PiAbs_reabac_ptheta);
hs_ptheta->Add(h_PiAbs_other_ptheta);

hs_ptheta->Draw("hist");
h_PiAbs_sel_ptheta_data->SetLineWidth(4);
h_PiAbs_sel_ptheta_data->Draw("same");
hs_ptheta->GetXaxis()->SetTitle("#theta");
hs_ptheta->GetYaxis()->SetTitle("No. of Events");
hs_ptheta->SetMaximum(100);
gPad->Modified();
pmultleg->Draw("same");
cptheta->SaveAs("h_proton_theta.png");




//====================================================================
TH1D *h_PiAbs_other_pphi;
TH1D *h_PiAbs_sig_pphi;
TH1D *h_PiAbs_chxbac_pphi;
TH1D *h_PiAbs_reabac_pphi;
TH1D *h_PiAbs_sel_pphi_data;

h_PiAbs_other_pphi=(TH1D*)inputfile->Get("h_PiAbs_other_pphi");
h_PiAbs_sig_pphi=(TH1D*)inputfile->Get("h_PiAbs_sig_pphi");
h_PiAbs_chxbac_pphi=(TH1D*)inputfile->Get("h_PiAbs_chxbac_pphi");
h_PiAbs_reabac_pphi=(TH1D*)inputfile->Get("h_PiAbs_reabac_pphi");
h_PiAbs_sel_pphi_data=(TH1D*)inputfile_data->Get("h_PiAbs_sel_pphi");

h_PiAbs_sel_pphi_data->Scale(scale_fac);


h_PiAbs_other_pphi->Rebin(2);
h_PiAbs_sig_pphi->Rebin(2);
h_PiAbs_chxbac_pphi->Rebin(2);
h_PiAbs_reabac_pphi->Rebin(2);
h_PiAbs_sel_pphi_data->Rebin(2);

h_PiAbs_sig_pphi->SetFillColor(kMagenta);
h_PiAbs_sig_pphi->SetLineColor(kMagenta);
h_PiAbs_chxbac_pphi->SetFillColor(kGreen+1);
h_PiAbs_chxbac_pphi->SetLineColor(kGreen+1);
h_PiAbs_reabac_pphi->SetFillColor(kBlue);
h_PiAbs_reabac_pphi->SetLineColor(kBlue);
h_PiAbs_other_pphi->SetFillColor(kRed);
h_PiAbs_other_pphi->SetLineColor(kRed);


TCanvas *cpphi= new TCanvas("cpphi", "cpphi", 900, 600);

THStack *hs_pphi = new THStack("hs_pphi", "");
hs_pphi->Add(h_PiAbs_sig_pphi);
hs_pphi->Add(h_PiAbs_chxbac_pphi);
hs_pphi->Add(h_PiAbs_reabac_pphi);
hs_pphi->Add(h_PiAbs_other_pphi);
hs_pphi->Draw("hist");
h_PiAbs_sel_pphi_data->SetLineWidth(4);
h_PiAbs_sel_pphi_data->Draw("same");
hs_pphi->SetMaximum(100);
hs_pphi->GetXaxis()->SetTitle("Azimuthal Angle #phi [Rad]");
hs_pphi->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cpphi->SaveAs("h_proton_phi.png");
//=========================================================
TH1D *h_other_Emissing;
TH1D *h_sig_Emissing;
TH1D *h_chxbac_Emissing;
TH1D *h_reabac_Emissing;
TH1D *h_sel_Emissing_data;

h_other_Emissing=(TH1D*)inputfile->Get("h_other_Emissing");
h_sig_Emissing=(TH1D*)inputfile->Get("h_sig_Emissing");
h_chxbac_Emissing=(TH1D*)inputfile->Get("h_chxbac_Emissing");
h_reabac_Emissing=(TH1D*)inputfile->Get("h_reabac_Emissing");
h_sel_Emissing_data=(TH1D*)inputfile_data->Get("h_sel_Emissing");


scale_fac=(h_sig_Emissing->Integral()+h_chxbac_Emissing->Integral()+h_reabac_Emissing->Integral()+h_other_Emissing->Integral())/h_sel_Emissing_data->Integral();

h_sel_Emissing_data->Scale(scale_fac);

h_sig_Emissing->SetFillColor(kMagenta);
h_sig_Emissing->SetLineColor(kMagenta);
h_chxbac_Emissing->SetFillColor(kGreen+1);
h_chxbac_Emissing->SetLineColor(kGreen+1);
h_reabac_Emissing->SetFillColor(kBlue);
h_reabac_Emissing->SetLineColor(kBlue);
h_other_Emissing->SetFillColor(kRed);
h_other_Emissing->SetLineColor(kRed);


TCanvas *cEmissing= new TCanvas("cEmissing", "cEmissing", 900, 600);

THStack *hs_Emissing = new THStack("hs_Emissing", "");
hs_Emissing->Add(h_sig_Emissing);
hs_Emissing->Add(h_chxbac_Emissing);
hs_Emissing->Add(h_reabac_Emissing);
hs_Emissing->Add(h_other_Emissing);
hs_Emissing->Draw("hist");
h_sel_Emissing_data->SetLineWidth(4);
h_sel_Emissing_data->Draw("same");
hs_Emissing->SetMaximum(50);
hs_Emissing->GetXaxis()->SetTitle("Missing Energy[GeV]");
hs_Emissing->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cEmissing->SaveAs("h_Emissing.png");
//====================================================
TH1D *h_other_Pmissing;
TH1D *h_sig_Pmissing;
TH1D *h_chxbac_Pmissing;
TH1D *h_reabac_Pmissing;
TH1D *h_sel_Pmissing_data;

h_other_Pmissing=(TH1D*)inputfile->Get("h_other_Pmissing");
h_sig_Pmissing=(TH1D*)inputfile->Get("h_sig_Pmissing");
h_chxbac_Pmissing=(TH1D*)inputfile->Get("h_chxbac_Pmissing");
h_reabac_Pmissing=(TH1D*)inputfile->Get("h_reabac_Pmissing");
h_sel_Pmissing_data=(TH1D*)inputfile_data->Get("h_sel_Pmissing");

h_sel_Pmissing_data->Scale(scale_fac);

h_sig_Pmissing->SetFillColor(kMagenta);
h_sig_Pmissing->SetLineColor(kMagenta);
h_chxbac_Pmissing->SetFillColor(kGreen+1);
h_chxbac_Pmissing->SetLineColor(kGreen+1);
h_reabac_Pmissing->SetFillColor(kBlue);
h_reabac_Pmissing->SetLineColor(kBlue);
h_other_Pmissing->SetFillColor(kRed);
h_other_Pmissing->SetLineColor(kRed);


TCanvas *cPmissing= new TCanvas("cPmissing", "cPmissing", 900, 600);

THStack *hs_Pmissing = new THStack("hs_Pmissing", "");
hs_Pmissing->Add(h_sig_Pmissing);
hs_Pmissing->Add(h_chxbac_Pmissing);
hs_Pmissing->Add(h_reabac_Pmissing);
hs_Pmissing->Add(h_other_Pmissing);
hs_Pmissing->Draw("hist");
h_sel_Pmissing_data->SetLineWidth(4);
h_sel_Pmissing_data->Draw("same");
hs_Pmissing->SetMaximum(50);
hs_Pmissing->GetXaxis()->SetTitle("Missing Momentum[GeV]");
hs_Pmissing->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cPmissing->SaveAs("h_Pmissing.png");
//=========================================================
TH1D *h_other_Ptmissing;
TH1D *h_sig_Ptmissing;
TH1D *h_chxbac_Ptmissing;
TH1D *h_reabac_Ptmissing;
TH1D *h_sel_Ptmissing_data;

h_other_Ptmissing=(TH1D*)inputfile->Get("h_other_Ptmissing");
h_sig_Ptmissing=(TH1D*)inputfile->Get("h_sig_Ptmissing");
h_chxbac_Ptmissing=(TH1D*)inputfile->Get("h_chxbac_Ptmissing");
h_reabac_Ptmissing=(TH1D*)inputfile->Get("h_reabac_Ptmissing");
h_sel_Ptmissing_data=(TH1D*)inputfile_data->Get("h_sel_Ptmissing");

h_sel_Ptmissing_data->Scale(scale_fac);

h_sig_Ptmissing->SetFillColor(kMagenta);
h_sig_Ptmissing->SetLineColor(kMagenta);
h_chxbac_Ptmissing->SetFillColor(kGreen+1);
h_chxbac_Ptmissing->SetLineColor(kGreen+1);
h_reabac_Ptmissing->SetFillColor(kBlue);
h_reabac_Ptmissing->SetLineColor(kBlue);
h_other_Ptmissing->SetFillColor(kRed);
h_other_Ptmissing->SetLineColor(kRed);


TCanvas *cPtmissing= new TCanvas("cPtmissing", "cPtmissing", 900, 600);

THStack *hs_Ptmissing = new THStack("hs_Ptmissing", "");
hs_Ptmissing->Add(h_sig_Ptmissing);
hs_Ptmissing->Add(h_chxbac_Ptmissing);
hs_Ptmissing->Add(h_reabac_Ptmissing);
hs_Ptmissing->Add(h_other_Ptmissing);
hs_Ptmissing->Draw("hist");
h_sel_Ptmissing_data->SetLineWidth(4);
h_sel_Ptmissing_data->Draw("same");
hs_Ptmissing->SetMaximum(50);
hs_Ptmissing->GetXaxis()->SetTitle("Transverse Missing P[GeV]");
hs_Ptmissing->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cPtmissing->SaveAs("h_Ptmissing.png");
//====================================================
TH1D *h_other_Plongit;
TH1D *h_sig_Plongit;
TH1D *h_chxbac_Plongit;
TH1D *h_reabac_Plongit;
TH1D *h_sel_Plongit_data;

h_other_Plongit=(TH1D*)inputfile->Get("h_other_Plongit");
h_sig_Plongit=(TH1D*)inputfile->Get("h_sig_Plongit");
h_chxbac_Plongit=(TH1D*)inputfile->Get("h_chxbac_Plongit");
h_reabac_Plongit=(TH1D*)inputfile->Get("h_reabac_Plongit");
h_sel_Plongit_data=(TH1D*)inputfile_data->Get("h_sel_Plongit");

h_sel_Plongit_data->Scale(scale_fac);

h_sig_Plongit->SetFillColor(kMagenta);
h_sig_Plongit->SetLineColor(kMagenta);
h_chxbac_Plongit->SetFillColor(kGreen+1);
h_chxbac_Plongit->SetLineColor(kGreen+1);
h_reabac_Plongit->SetFillColor(kBlue);
h_reabac_Plongit->SetLineColor(kBlue);
h_other_Plongit->SetFillColor(kRed);
h_other_Plongit->SetLineColor(kRed);


TCanvas *cPlongit= new TCanvas("cPlongit", "cPlongit", 900, 600);

THStack *hs_Plongit = new THStack("hs_Plongit", "");
hs_Plongit->Add(h_sig_Plongit);
hs_Plongit->Add(h_chxbac_Plongit);
hs_Plongit->Add(h_reabac_Plongit);
hs_Plongit->Add(h_other_Plongit);
hs_Plongit->Draw("hist");
h_sel_Plongit_data->SetLineWidth(4);
h_sel_Plongit_data->Draw("same");
hs_Plongit->SetMaximum(30);
hs_Plongit->GetXaxis()->SetTitle("Longitudinal Momentum[GeV]");
hs_Plongit->GetYaxis()->SetTitle("No. of Events");
gPad->Modified();
pmultleg->Draw("same");
cPlongit->SaveAs("h_Plongit.png");



//========================================================
TCanvas *ckk = new TCanvas("ckk", "ckk", 900, 600);
h_mom_gentruep->SetLineColor(kGreen+1);
h_mom_gentruep->SetFillColor(kGreen+1);
h_mom_gentruep->SetFillStyle(3003);
h_mom_gentruep->GetXaxis()->SetTitle("True Momentum [GeV]");
h_mom_gentruep->GetYaxis()->SetTitle("No. of Tracks");
h_mom_gentruep->SetTitle("");
h_mom_recotruep_test->SetLineColor(kGreen+1);
h_mom_recotruep_test->SetFillColor(kGreen+1);
h_mom_recotruep_test->SetFillStyle(3001);
h_mom_recotruep_test->SetTitle("");
h_mom_gentruep->Draw("hist");
h_mom_recotruep_test->Draw("hist,same");

TLegend *lg=new TLegend(0.4, 0.7, 0.85, 0.85);
lg->SetBorderSize(0);
lg->AddEntry(h_mom_gentruep, "Generated Protons");
lg->AddEntry(h_mom_recotruep_test, "Reconstructed Protons");
lg->Draw("same");

ckk->SaveAs("h_genvstruep_truemomentum.png");
//=============================================================

 


TCanvas *c=new TCanvas("c", "c", 900, 600);


  //h_mom_recotruep->Divide(h_mom_gentruep);
  //h_mom_recotruep->Draw();
  TEfficiency *pEff_nocut = new TEfficiency(*h_mom_recotruep_test, *h_mom_gentruep);
  pEff_nocut->SetLineColor(kMagenta+1); 
  pEff_nocut->SetMarkerColor(kMagenta+1);
  pEff_nocut->SetLineWidth(2);
  pEff_nocut->SetMarkerStyle(20);
  pEff_nocut->SetMarkerSize(0.5);
  //pEff_nocut->SetMaximum(1.4);
  pEff_nocut->SetTitle(";True Momentum [GeV];No. of Tracks");
  pEff_nocut->Draw();
  
  TEfficiency *pEff_trkscore = new TEfficiency(*h_mom_trkscoretruep, *h_mom_gentruep);
  pEff_trkscore->SetLineColor(kRed); 
  pEff_trkscore->SetMarkerColor(kRed);
  pEff_trkscore->SetLineWidth(2);
  pEff_trkscore->SetMarkerStyle(20);
  pEff_trkscore->SetMarkerSize(0.5);
  pEff_trkscore->Draw("same");
  
  TEfficiency *pEff_tmdqdx = new TEfficiency(*h_mom_tmdqdxtruep, *h_mom_gentruep);
  pEff_tmdqdx->SetLineColor(kBlue); 
  pEff_tmdqdx->SetMarkerColor(kBlue);
  pEff_tmdqdx->SetLineWidth(2);
  pEff_tmdqdx->SetMarkerStyle(20);
  pEff_tmdqdx->SetMarkerSize(0.5);
  pEff_tmdqdx->Draw("same");
  
  TEfficiency *pEff_sel = new TEfficiency(*h_mom_selectedtruep, *h_mom_gentruep);
  pEff_sel->SetLineColor(kGreen+1); 
  pEff_sel->SetMarkerColor(kGreen+1);
  pEff_sel->SetLineWidth(2);
  pEff_sel->SetMarkerStyle(20);
  pEff_sel->SetMarkerSize(0.5);
  pEff_sel->Draw("same");

  std::cout<<"total gen tracks "<<h_mom_gentruep->GetEntries()<<std::endl;
  std::cout<<"total sel tracks "<<h_mom_selectedtruep->GetEntries()<<std::endl; 

 
  TLegend *effleg = new TLegend(0.1, 0.6, 0.4, 0.85);
  effleg->SetBorderSize(0);
  effleg->SetFillStyle(0);
  effleg->AddEntry(pEff_nocut, "No Cut");
  effleg->AddEntry(pEff_trkscore, "TrkScore Cut");
  effleg->AddEntry(pEff_tmdqdx, "TrunMean dEdX Cut");
  effleg->AddEntry(pEff_sel, "Chi2 Cut");
 
  effleg->Draw("same");
  c->SaveAs("h_all_efficiency.png");   
  //-------------------------------------------------------------
  

 TCanvas *cpi=new TCanvas("cpi", "cpi", 900, 600);


  //h_mom_recotruep->Divide(h_mom_gentruep);
  //h_mom_recotruep->Draw();
  TEfficiency *pipmEff_nocut = new TEfficiency(*h_mom_recotruepionpm, *h_mom_gentruepionpm);
  pipmEff_nocut->SetLineColor(kMagenta+1); 
  pipmEff_nocut->SetMarkerColor(kMagenta+1);
  pipmEff_nocut->SetLineWidth(2);
  pipmEff_nocut->SetMarkerStyle(20);
  pipmEff_nocut->SetMarkerSize(0.5);
  //pipmEff_nocut->GetYaxis()->SetRange(0, 1.4);
  pipmEff_nocut->SetTitle(";True Momentum [GeV];No. of Tracks");
  pipmEff_nocut->Draw();
  
  TEfficiency *pipmEff_trkscore = new TEfficiency(*h_mom_trkscoretruepipm, *h_mom_gentruepionpm);
  pipmEff_trkscore->SetLineColor(kRed); 
  pipmEff_trkscore->SetMarkerColor(kRed);
  pipmEff_trkscore->SetLineWidth(2);
  pipmEff_trkscore->SetMarkerStyle(20);
  pipmEff_trkscore->SetMarkerSize(0.5);
  pipmEff_trkscore->Draw("same");
  
  TEfficiency *pipmEff_tmdqdx = new TEfficiency(*h_mom_tmdqdxtruepipm, *h_mom_gentruepionpm);
  pipmEff_tmdqdx->SetLineColor(kBlue); 
  pipmEff_tmdqdx->SetMarkerColor(kBlue);
  pipmEff_tmdqdx->SetLineWidth(2);
  pipmEff_tmdqdx->SetMarkerStyle(20);
  pipmEff_tmdqdx->SetMarkerSize(0.5);
  pipmEff_tmdqdx->Draw("same");
  
  TEfficiency *pipmEff_sel = new TEfficiency(*h_mom_selectedtruepipm, *h_mom_gentruepionpm);
  pipmEff_sel->SetLineColor(kGreen+1); 
  pipmEff_sel->SetMarkerColor(kGreen+1);
  pipmEff_sel->SetLineWidth(2);
  pipmEff_sel->SetMarkerStyle(20);
  pipmEff_sel->SetMarkerSize(0.5);
  pipmEff_sel->Draw("same");
 
  //pipmEff_nocut->GetYaxis()->SetUserRange(0, 1.4);
 
  effleg->Draw("same");
  gPad->Update();
  cpi->SaveAs("h_all_efficiency_chargedpion.png");   







  //------------------------------------------------------------
  TH1D *h_orimom_proton;
  TH1D *h_orimom_neutron;
  TH1D *h_orimom_pionpm;
  TH1D *h_orimom_pion0;
  TH1D *h_orimom_kaon;
  TH1D *h_orimom_electron;
  TH1D *h_orimom_muon;
  TH1D *h_orimom_photon;
  TH1D *h_orimom_other;

  h_orimom_proton=(TH1D*)inputfile->Get("h_orimom_proton");
  h_orimom_neutron=(TH1D*)inputfile->Get("h_orimom_neutron");
  h_orimom_pionpm=(TH1D*)inputfile->Get("h_orimom_pionpm");
  h_orimom_pion0=(TH1D*)inputfile->Get("h_orimom_pion0");
  h_orimom_kaon=(TH1D*)inputfile->Get("h_orimom_kaon");
  h_orimom_electron=(TH1D*)inputfile->Get("h_orimom_electron");
  h_orimom_muon=(TH1D*)inputfile->Get("h_orimom_muon");
  h_orimom_photon=(TH1D*)inputfile->Get("h_orimom_photon");
  h_orimom_other=(TH1D*)inputfile->Get("h_orimom_other");

  TCanvas *corimom = new TCanvas("corimom", "corimom", 900, 600);
  THStack *hs_orimom = new THStack("hs_pmom", "");
  h_orimom_proton->SetLineColor(kMagenta);
  h_orimom_proton->SetFillColor(kMagenta);
  h_orimom_neutron->SetLineColor(kBlack);
  h_orimom_neutron->SetFillColor(kBlack);
  h_orimom_pionpm->SetLineColor(kBlue);
  h_orimom_pionpm->SetFillColor(kBlue);
  h_orimom_pion0->SetLineColor(kGreen);
  h_orimom_pion0->SetFillColor(kGreen);

  h_orimom_kaon->SetLineColor(kOrange-3);
  h_orimom_kaon->SetFillColor(kOrange-3);
  h_orimom_electron->SetLineColor(kCyan);
  h_orimom_electron->SetFillColor(kCyan);
  h_orimom_muon->SetLineColor(kRed);
  h_orimom_muon->SetFillColor(kRed);
  h_orimom_photon->SetLineColor(kGreen+1);
  h_orimom_photon->SetFillColor(kGreen+1);
  h_orimom_other->SetLineColor(kYellow);
  h_orimom_other->SetFillColor(kYellow);



  hs_orimom->Add(h_orimom_proton);
  //hs_orimom->Add(h_orimom_neutron);
  hs_orimom->Add(h_orimom_pionpm);
  //hs_orimom->Add(h_orimom_pion0);
  hs_orimom->Add(h_orimom_kaon);
  hs_orimom->Add(h_orimom_electron);
  hs_orimom->Add(h_orimom_muon);
  //hs_orimom->Add(h_orimom_photon);
  hs_orimom->Add(h_orimom_other);

  hs_orimom->Draw("hist");
  gPad->Modified(); gPad->Update();
  hs_orimom->GetXaxis()->SetTitle("True Momentum [GeV]");
  hs_orimom->GetYaxis()->SetTitle("No. of Particles");
  gPad->Modified(); gPad->Update();
  TLegend *legorimom=new TLegend(0.6, 0.5, 0.85, 0.85);
  legorimom->SetBorderSize(0);
  legorimom->AddEntry("h_orimom_proton", "Proton");
  legorimom->AddEntry("h_orimom_pionpm", "#pi^{+}/#pi^{-}");
  legorimom->AddEntry("h_orimom_kaon", "Kaon");
  legorimom->AddEntry("h_orimom_electron", "Electron");
  legorimom->AddEntry("h_orimom_muon", "Muon");
  //legorimom->AddEntry("h_orimom_photon", "#gamma");
  legorimom->AddEntry("h_orimom_other", "Other");

  legorimom->Draw("same");
  corimom->SaveAs("h_Ori_Mom_charged_cand.png");
  //==========================================
  


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
  leg2->Draw("same");

  crecomom_selected->SaveAs("h_recomom_selected_alltrack.png");  
  //-----------------------------------------------------------------------
  TH1D *h_recopimom_selected_proton;
  TH1D *h_recopimom_selected_neutron;
  TH1D *h_recopimom_selected_pionpm;
  TH1D *h_recopimom_selected_pion0;
  TH1D *h_recopimom_selected_kaon;
  TH1D *h_recopimom_selected_electron;
  TH1D *h_recopimom_selected_muon;
  TH1D *h_recopimom_selected_photon;
  TH1D *h_recopimom_selected_other;

  h_recopimom_selected_proton=(TH1D*)inputfile->Get("h_recopimom_selected_proton");
  //h_recopimom_selected_neutron=(TH1D*)inputfile->Get("h_recopimom_selected_neutron");
  h_recopimom_selected_pionpm=(TH1D*)inputfile->Get("h_recopimom_selected_pionpm");
  //h_recopimom_selected_pion0=(TH1D*)inputfile->Get("h_recopimom_selected_pion0");
  h_recopimom_selected_kaon=(TH1D*)inputfile->Get("h_recopimom_selected_kaon");
  h_recopimom_selected_electron=(TH1D*)inputfile->Get("h_recopimom_selected_electron");
  h_recopimom_selected_muon=(TH1D*)inputfile->Get("h_recopimom_selected_muon");
  //h_recopimom_selected_photon=(TH1D*)inputfile->Get("h_recopimom_selected_photon");
  h_recopimom_selected_other=(TH1D*)inputfile->Get("h_recopimom_selected_other");

  TCanvas *crecopimom_selected = new TCanvas("crecopimom_selected", "crecopimom_selected", 900, 600);
  THStack *hs_recopimom_selected = new THStack("hs_pmom", "");
  h_recopimom_selected_proton->SetLineColor(kMagenta);
  h_recopimom_selected_proton->SetFillColor(kMagenta);
  //h_recopimom_selected_neutron->SetLineColor(kBlack);
  //h_recopimom_selected_neutron->SetFillColor(kBlack);
  h_recopimom_selected_pionpm->SetLineColor(kBlue);
  h_recopimom_selected_pionpm->SetFillColor(kBlue);
  //h_recopimom_selected_pion0->SetLineColor(kGreen);
  //h_recopimom_selected_pion0->SetFillColor(kGreen);

  h_recopimom_selected_kaon->SetLineColor(kOrange-3);
  h_recopimom_selected_kaon->SetFillColor(kOrange-3);
  h_recopimom_selected_electron->SetLineColor(kCyan);
  h_recopimom_selected_electron->SetFillColor(kCyan);
  h_recopimom_selected_muon->SetLineColor(kRed);
  h_recopimom_selected_muon->SetFillColor(kRed);
  //h_recopimom_selected_photon->SetLineColor(kGreen+1);
  //h_recopimom_selected_photon->SetFillColor(kGreen+1);
  h_recopimom_selected_other->SetLineColor(kYellow);
  h_recopimom_selected_other->SetFillColor(kYellow);



  hs_recopimom_selected->Add(h_recopimom_selected_proton);
  //hs_recopimom_selected->Add(h_recopimom_selected_neutron);
  hs_recopimom_selected->Add(h_recopimom_selected_pionpm);
  //hs_recopimom_selected->Add(h_recopimom_selected_pion0);
  hs_recopimom_selected->Add(h_recopimom_selected_kaon);
  hs_recopimom_selected->Add(h_recopimom_selected_electron);
  hs_recopimom_selected->Add(h_recopimom_selected_muon);
  //hs_recopimom_selected->Add(h_recopimom_selected_photon);
  hs_recopimom_selected->Add(h_recopimom_selected_other);

  hs_recopimom_selected->Draw("hist");

  hs_recopimom_selected->GetXaxis()->SetTitle("Track Range [cm]");
  hs_recopimom_selected->GetYaxis()->SetTitle("No. of Tracks");
  gPad->Modified();
  leg2->Draw("same");

  crecopimom_selected->SaveAs("h_recopimom_selected_alltrack.png");  

 


  

  //--------------------------------------------------------------------------

  //==========================================================================
  



  //===========================================================================
}
