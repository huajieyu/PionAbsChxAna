//#include "plot.h"

void GetMuondEdxR(){
  //rr, dedx = array.array('f'), array.array('f')
  std::vector <float> rr, dedx;
  rr.push_back(0.9833);
  rr.push_back(1.786);
  rr.push_back(3.321);
  rr.push_back(6.598);
  rr.push_back(10.58);
  rr.push_back(30.84);
  rr.push_back(42.50);
  rr.push_back(67.50);
  rr.push_back(106.3);
  rr.push_back(172.5);
  rr.push_back(238.5);
  rr.push_back(493.4);
  rr.push_back(616.3);
  rr.push_back(855.2);
  rr.push_back(1202);
  rr.push_back(1758);
  rr.push_back(2297);
  dedx.push_back(5.687);
  dedx.push_back(4.461);
  dedx.push_back(3.502);
  dedx.push_back(2.731);
  dedx.push_back(2.340);
  dedx.push_back(1.771);
  dedx.push_back(1.670);
  dedx.push_back(1.570);
  dedx.push_back(1.519);
  dedx.push_back(1.510);
  dedx.push_back(1.526);
  dedx.push_back(1.610);
  dedx.push_back(1.645);
  dedx.push_back(1.700);
  dedx.push_back(1.761);
  dedx.push_back(1.829);
  dedx.push_back(1.877);
  for (int i=0; i< rr.size(); i++){
      rr[i]/=1.396;
      dedx[i]*=1.396;
   }
  Int_t n2=rr.size();
  float rr2[n2], dedx2[n2];
  for(int i2=0; i2<n2; i2++){
    rr2[i2]=rr[i2];
    dedx2[i2]=dedx[i2];
  }
  TGraph *gr = new TGraph(n2, rr2, dedx2);
  gr->SetLineWidth(2);
  gr->SetLineColor(3);
  gr->Draw("same");
}

void plot_dEdx(){
  //loadStyle();
 
   

  int tune=1;
  int cosmicCut=1;

  TFile *input0;
  TFile *input1;
  if(tune==1){
  if (cosmicCut){
    input0 = new TFile("output_test_data.root");
    input1 = new TFile("output_test.root");
  }
  } 

  TFile *f0=new TFile("TheoreticalPredictions_dedxresrange.root");


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
    //gStyle->SetPadGridX(kTRUE);
    //gStyle->SetPadGridY(kTRUE);
    //gStyle->SetGridStyle();
  //----------------------------------------------------------------------

  gROOT->SetBatch(1);
 TGraph* graph=(TGraph*)f0->Get("proton");
 graph->SetLineColor(2);
 graph->SetLineWidth(2);
 
 TGraph* graphmu=(TGraph*)f0->Get("muon");
 graphmu->SetLineColor(3);
 graphmu->SetLineWidth(2);

 TF1 *f_dEdxVsRR = new TF1("f_dEdxVsRR","17 * x^(-0.42)",0,100);
 TF1 *f_dEdxVsRRmu = new TF1("f_dEdxVsRRmu","8 * x^(-0.37)",0,100);
 f_dEdxVsRRmu->SetLineColor(3);
 TF1 *f_dEdxVsRRMIP = new TF1("f_dEdxVsRRMIP","2.12",0,100);
 f_dEdxVsRRMIP->SetLineColor(6);


  cout<<"libo test 0"<<endl;
  TH2D *selproton_dEdx_vs_resrange_pcont_BNB; 
  selproton_dEdx_vs_resrange_pcont_BNB=(TH2D*)input0->Get("h_selproton_dedxRR");

  TH2D *selproton_dEdx_vs_resrange_pcont_MC; 
  selproton_dEdx_vs_resrange_pcont_MC=(TH2D*)input1->Get("h_selproton_dedxRR");

  TCanvas *cs = new TCanvas("cs");
 
  cs->cd(); 
  selproton_dEdx_vs_resrange_pcont_BNB->SetTitle("");
  selproton_dEdx_vs_resrange_pcont_BNB->GetXaxis()->SetTitle("Track Residual Range[cm]");
  selproton_dEdx_vs_resrange_pcont_BNB->GetYaxis()->SetTitle("dEdx[MeV/cm]");
  selproton_dEdx_vs_resrange_pcont_BNB->Draw("colz");

  TLegend* leg = new TLegend(0.65, 0.65, .9, .9);
  leg->AddEntry(graphmu, "Muon Expectation", "l");
  leg->AddEntry(graph, "Proton Expectation", "l");
  leg->AddEntry(f_dEdxVsRRMIP,"MIP expectation", "l");
  //graph->Draw("same");
  f_dEdxVsRR->Draw("same");
  //f_dEdxVsRRmu->Draw("same");
  GetMuondEdxR();
  f_dEdxVsRRMIP->Draw("same");
  leg->Draw("same");
  cs->Print("selproton_dEdx_vs_resrange_pcont_BNB.png");
  cs->Update();
  
  //88888888888888888888888888888888888888888888888888888888888888888
  selproton_dEdx_vs_resrange_pcont_MC->SetTitle("");
  selproton_dEdx_vs_resrange_pcont_MC->GetXaxis()->SetTitle("Track Residual Range[cm]");
  selproton_dEdx_vs_resrange_pcont_MC->GetYaxis()->SetTitle("dEdx[MeV/cm]");
  selproton_dEdx_vs_resrange_pcont_MC->Draw("colz");

  //graph->Draw("same");
  f_dEdxVsRR->Draw("same");
  //f_dEdxVsRRmu->Draw("same");
  GetMuondEdxR();
  f_dEdxVsRRMIP->Draw("same");
  leg->Draw("same");
  cs->Print("selproton_dEdx_vs_resrange_pcont_MC.png");
  cs->Update();
  
  TH2D* h_selproton_momreso = (TH2D*)input1->Get("h_selproton_momreso"); 
  TCanvas *cd = new TCanvas("cd","cd", 900, 700);
  h_selproton_momreso->GetXaxis()->SetTitle("True Momentum[GeV]");
  h_selproton_momreso->GetYaxis()->SetTitle("Reco Momentum[GeV]");
  h_selproton_momreso->Draw("colz");
  cd->SaveAs("h_true_reco_proton_mom.png");

  TH2D* h_selproton_costhetareso = (TH2D*)input1->Get("h_selproton_costhetareso"); 
  TCanvas *cdd = new TCanvas("cdd","cdd", 900, 700);
  h_selproton_costhetareso->GetXaxis()->SetTitle("True cos#theta");
  h_selproton_costhetareso->GetYaxis()->SetTitle("Reco cos#thetaGeV]");
  h_selproton_costhetareso->Draw("colz");
  cdd->SaveAs("h_true_reco_proton_costheta.png");




 
  //88888888888888888888888888888888888888888888888888888888888888888888
  /*
  */
  //============================================================================ 
} //========================================================
