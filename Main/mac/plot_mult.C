void plot_mult(){


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



TH1D* h_shwcut_mult[5];

h_shwcut_mult[0] = (TH1D*)inputfile->Get("h_shwcut_mult");
h_shwcut_mult[1] = (TH1D*)inputfile->Get("hsig_shwcut_mult");
h_shwcut_mult[2] = (TH1D*)inputfile->Get("hchxbac_shwcut_mult");
h_shwcut_mult[3] = (TH1D*)inputfile->Get("hreabac_shwcut_mult");
h_shwcut_mult[4] = (TH1D*)inputfile_data->Get("h_shwcut_mult");

TH1D* h_tmdqdxcut_mult[5];

h_tmdqdxcut_mult[0] = (TH1D*)inputfile->Get("h_tmdqdxcut_mult");
h_tmdqdxcut_mult[1] = (TH1D*)inputfile->Get("hsig_tmdqdxcut_mult");
h_tmdqdxcut_mult[2] = (TH1D*)inputfile->Get("hchxbac_tmdqdxcut_mult");
h_tmdqdxcut_mult[3] = (TH1D*)inputfile->Get("hreabac_tmdqdxcut_mult");
h_tmdqdxcut_mult[4] = (TH1D*)inputfile_data->Get("h_tmdqdxcut_mult");

TH1D* h_chi2cut_mult[5];

h_chi2cut_mult[0] = (TH1D*)inputfile->Get("h_chi2cut_mult");
h_chi2cut_mult[1] = (TH1D*)inputfile->Get("hsig_chi2cut_mult");
h_chi2cut_mult[2] = (TH1D*)inputfile->Get("hchxbac_chi2cut_mult");
h_chi2cut_mult[3] = (TH1D*)inputfile->Get("hreabac_chi2cut_mult");
h_chi2cut_mult[4] = (TH1D*)inputfile_data->Get("h_chi2cut_mult");



//---------------------------------------------------------------------

TCanvas *cc= new TCanvas("cc", "cc", 900, 600);
THStack *hs_new = new THStack("hs_new", "");

h_shwcut_mult[0]->SetLineColor(kBlue);
//h_shwcut_mult[0]->SetFillColor(kBlue);
h_shwcut_mult[0]->SetLineWidth(2);
h_shwcut_mult[1]->SetLineColor(kRed);
h_shwcut_mult[1]->SetFillColor(kRed);
h_shwcut_mult[2]->SetLineColor(kMagenta);
h_shwcut_mult[2]->SetFillColor(kMagenta);
h_shwcut_mult[3]->SetLineColor(kGreen+1);
h_shwcut_mult[3]->SetFillColor(kGreen+1);


TLegend *leg=new TLegend(0.55, 0.65, 0.85, 0.85);
leg->SetBorderSize(0);
leg->SetNColumns(2);
leg->AddEntry(h_shwcut_mult[0], "Total MC");
leg->AddEntry(h_shwcut_mult[1], "Signal");
leg->AddEntry(h_shwcut_mult[2], "ChxBKG");
leg->AddEntry(h_shwcut_mult[3], "ReaBKG");

hs_new->Add(h_shwcut_mult[1]);
hs_new->Add(h_shwcut_mult[2]);
hs_new->Add(h_shwcut_mult[3]);


hs_new->Draw("hist");
h_shwcut_mult[0]->Draw("same");
leg->Draw("same");
hs_new->GetXaxis()->SetTitle("Track Multiplicity");
hs_new->GetYaxis()->SetTitle("No. of Events");
hs_new->SetTitle("");
gPad->Modified();

cc->SaveAs("h_shwcut_mult.png");
//======================================================
TCanvas *cc_datavsmc = new TCanvas("cc_datavsmc", "cc_datavsmc", 900,600);
h_shwcut_mult[4]->SetLineColor(kRed);
h_shwcut_mult[4]->SetLineWidth(2);
Double_t fac_num=h_shwcut_mult[4]->Integral();
Double_t fac_denom=h_shwcut_mult[0]->Integral();
h_shwcut_mult[0]->Scale(fac_num/fac_denom);

h_shwcut_mult[0]->SetTitle("");
h_shwcut_mult[0]->GetXaxis()->SetTitle("Track Multiplicity");
h_shwcut_mult[0]->GetYaxis()->SetTitle("No. of Events");

h_shwcut_mult[0]->Draw("hist");
h_shwcut_mult[4]->Draw("hist,same");

TLegend *legmd=new TLegend(0.55, 0.65, 0.85, 0.85);
legmd->AddEntry(h_shwcut_mult[0], "MC");
legmd->AddEntry(h_shwcut_mult[4], "Data");
legmd->Draw("same");
cc_datavsmc->SaveAs("h_shwcut_datavsmc.png");


//===================================================
TCanvas *cc2= new TCanvas("cc2", "cc2", 900, 600);
THStack *hs_new2 = new THStack("hs_new2", "");

h_tmdqdxcut_mult[0]->SetLineColor(kBlue);
//h_tmdqdxcut_mult[0]->SetFillColor(kBlue);
h_tmdqdxcut_mult[0]->SetLineWidth(2);
h_tmdqdxcut_mult[1]->SetLineColor(kRed);
h_tmdqdxcut_mult[1]->SetFillColor(kRed);
h_tmdqdxcut_mult[2]->SetLineColor(kMagenta);
h_tmdqdxcut_mult[2]->SetFillColor(kMagenta);
h_tmdqdxcut_mult[3]->SetLineColor(kGreen+1);
h_tmdqdxcut_mult[3]->SetFillColor(kGreen+1);



hs_new2->Add(h_tmdqdxcut_mult[1]);
hs_new2->Add(h_tmdqdxcut_mult[2]);
hs_new2->Add(h_tmdqdxcut_mult[3]);


hs_new2->Draw("hist");
h_tmdqdxcut_mult[0]->Draw("same");
leg->Draw("same");
hs_new2->GetXaxis()->SetTitle("Track Multiplicity");
hs_new2->GetYaxis()->SetTitle("No. of Events");
hs_new2->SetTitle("");
gPad->Modified();

cc2->SaveAs("h_tmdqdxcut_mult.png");
//=======================================================
TCanvas *cc2_datavsmc = new TCanvas("cc2_datavsmc", "cc2_datavsmc", 900,600);
h_tmdqdxcut_mult[4]->SetLineColor(kRed);
h_tmdqdxcut_mult[4]->SetLineWidth(2);
fac_num=h_tmdqdxcut_mult[4]->Integral();
fac_denom=h_tmdqdxcut_mult[0]->Integral();
h_tmdqdxcut_mult[0]->Scale(fac_num/fac_denom);

h_tmdqdxcut_mult[0]->SetTitle("");
h_tmdqdxcut_mult[0]->GetXaxis()->SetTitle("Track Multiplicity");
h_tmdqdxcut_mult[0]->GetYaxis()->SetTitle("No. of Events");

h_tmdqdxcut_mult[0]->Draw("hist");
h_tmdqdxcut_mult[4]->Draw("hist,same");
legmd->Draw("same");
cc2_datavsmc->SaveAs("h_tmdqdxcut_datavsmc.png");






TCanvas *cc3= new TCanvas("cc3", "cc3", 900, 600);
THStack *hs_new3 = new THStack("hs_new3", "");

h_chi2cut_mult[0]->SetLineColor(kBlue);
//h_chi2cut_mult[0]->SetFillColor(kBlue);
h_chi2cut_mult[0]->SetLineWidth(2);
h_chi2cut_mult[1]->SetLineColor(kRed);
h_chi2cut_mult[1]->SetFillColor(kRed);
h_chi2cut_mult[2]->SetLineColor(kMagenta);
h_chi2cut_mult[2]->SetFillColor(kMagenta);
h_chi2cut_mult[3]->SetLineColor(kGreen+1);
h_chi2cut_mult[3]->SetFillColor(kGreen+1);



hs_new3->Add(h_chi2cut_mult[1]);
hs_new3->Add(h_chi2cut_mult[2]);
hs_new3->Add(h_chi2cut_mult[3]);


hs_new3->Draw("hist");
h_chi2cut_mult[0]->Draw("same");
leg->Draw("same");
hs_new3->GetXaxis()->SetTitle("Track Multiplicity");
hs_new3->GetYaxis()->SetTitle("No. of Events");
hs_new3->SetTitle("");
gPad->Modified();

cc3->SaveAs("h_chi2cut_mult.png");
//==================================================
TCanvas *cc3_datavsmc = new TCanvas("cc3_datavsmc", "cc3_datavsmc", 900,600);
h_chi2cut_mult[4]->SetLineColor(kRed);
h_chi2cut_mult[4]->SetLineWidth(2);
fac_num=h_chi2cut_mult[4]->Integral();
fac_denom=h_chi2cut_mult[0]->Integral();
h_chi2cut_mult[0]->Scale(fac_num/fac_denom);

h_chi2cut_mult[0]->SetTitle("");
h_chi2cut_mult[0]->GetXaxis()->SetTitle("Track Multiplicity");
h_chi2cut_mult[0]->GetYaxis()->SetTitle("No. of Events");

h_chi2cut_mult[0]->Draw("hist");
h_chi2cut_mult[4]->Draw("hist,same");

legmd->Draw("same");
cc3_datavsmc->SaveAs("h_chi2cut_datavsmc.png");





}
