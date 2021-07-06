void plot_tmva(){
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

gStyle->SetTitleSize(0.06, "X");
gStyle->SetTitleSize(0.06, "Y");
gStyle->SetTitleOffset(0.85, "X");
gStyle->SetTitleOffset(0.85, "Y");

TLatex* prelim = new TLatex(0.9,0.93, "protoDUNE-SP: Preliminary");
prelim->SetTextFont(62);
prelim->SetTextColor(kGray+2);
prelim->SetNDC();
prelim->SetTextSize(1/30.);
prelim->SetTextAlign(32);









TFile *inputfile_tot=new TFile("TMVApp_data_tot.root");
TFile *inputfile_sig=new TFile("TMVApp_data_sig.root");
TFile *inputfile_bkg=new TFile("TMVApp_data_bkg.root");
TH1D *htot=(TH1D*)inputfile_tot->Get("BDTG_data_minus");
TH1D *hsig=(TH1D*)inputfile_sig->Get("BDTG_data_minus");
TH1D *hbkg=(TH1D*)inputfile_bkg->Get("BDTG_data_minus");


TCanvas *c=new TCanvas("c", "c", 900, 700);
gPad->SetBottomMargin(0.12);
gPad->SetLeftMargin(0.12);
htot->SetLineWidth(3);
hsig->SetLineWidth(3);
hbkg->SetLineWidth(3);
hsig->SetLineColor(kMagenta);
hbkg->SetLineColor(kGreen+2);
hsig->SetFillColor(kMagenta);
hbkg->SetFillColor(kGreen+2);
hsig->SetFillStyle(3004);
hbkg->SetFillStyle(3005);

double intsig=0;
double intbkg=0;
for(int i=0; i<hsig->GetNbinsX(); i++){
	if(hsig->GetBinCenter(i+1)<0) continue;
	intsig +=hsig->GetBinContent(i+1);
	intbkg +=hbkg->GetBinContent(i+1);
}
std::cout<<"Total number of signal is "<<intsig<<std::endl;
std::cout<<"Total number of background is "<<intbkg<<std::endl;
std::cout<<"Purity is "<<(double)intsig/(intsig+intbkg)<<std::endl;
TLegend *leg=new TLegend(0.4, 0.7, 0.7, 0.85);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(hsig, "Signal");
leg->AddEntry(hbkg, "Background");
THStack *hs=new THStack("hs", "");
hs->Add(hsig);
hs->Add(hbkg);
hs->SetTitle(";BDTG Score; Nevts");
hs->Draw("hist");
leg->Draw("same");
prelim->Draw("same");

c->SaveAs("h_BDTScore.png");
}
