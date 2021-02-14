{
  gStyle->SetOptStat(0);
  TFile f1("histnew.root");
  //TFile f2("npi0_xsec.root");
  TFile f2("/dune/data2/users/calcuttj/GeantReweight/2_2_21/exclusive_xsec.root");
  TGraphErrors *gr_xs_incE = (TGraphErrors*)f1.Get("gr_xs_incE");
  TGraphErrors *gr_xscex_incE = (TGraphErrors*)f1.Get("gr_xscex_incE");
  TGraphErrors *gr_xsabs_incE = (TGraphErrors*)f1.Get("gr_xsabs_incE");

  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  TGraph *abs_KE = (TGraph*)f2.Get("abs_KE");
  TGraph *cex_KE = (TGraph*)f2.Get("cex_KE");

  TCanvas *c1 = new TCanvas("c1","c1");
  total_inel_KE->SetTitle("Inelastic;KE (MeV); #sigma (mb)");
  total_inel_KE->GetXaxis()->SetRangeUser(0,1000);
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("ac");
  gr_xs_incE->Draw("p");

  TCanvas *c2 = new TCanvas("c2","c2");
  abs_KE->Draw("ac");
  abs_KE->SetTitle("Absorption;KE (MeV); #sigma (mb)");
  abs_KE->GetXaxis()->SetRangeUser(0,1000);
  abs_KE->SetLineColor(2);
  gr_xsabs_incE->Draw("p");

  TCanvas *c3 = new TCanvas("c3","c3");
  cex_KE->Draw("ac");
  cex_KE->SetTitle("Charge exchange;KE (MeV); #sigma (mb)");
  cex_KE->GetXaxis()->SetRangeUser(0,1000);
  cex_KE->SetLineColor(2);
  cex_KE->GetYaxis()->SetRangeUser(0,150);
  gr_xscex_incE->Draw("p");

  const int nwires_in_slice = 20;
  const int nslices = 480/nwires_in_slice;

  TH1D *incE[nslices];
  TH1D *pitch[nslices];
  TH1D *dEdx[nslices];
  
  for (int i = 0; i<nslices; ++i){
    incE[i] = (TH1D*)f1.Get(Form("incE_%d",i));
    pitch[i] = (TH1D*)f1.Get(Form("pitch_%d",i));
    dEdx[i] = (TH1D*)f1.Get(Form("dEdx_%d",i));
  }    

  TCanvas *c4 = new TCanvas("c4","c4");
  pitch[5]->SetTitle(";Slice thickness (cm);N pions");
  pitch[5]->Draw();
  pitch[10]->SetLineColor(2);
  pitch[10]->Draw("same");
  pitch[15]->SetLineColor(3);
  pitch[15]->Draw("same");
  TLegend *leg4 = new TLegend(0.15,0.7,0.6,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(pitch[5],"Slice 5, 100<=wire #<120","l");
  leg4->AddEntry(pitch[10],"Slice 10, 200<=wire #<220","l");
  leg4->AddEntry(pitch[15],"Slice 15, 300<=wire #<320","l");
  leg4->Draw();

  TCanvas *c5 = new TCanvas("c5","c5");
  TGraph *gr_pitch_slc = (TGraph*)f1.Get("gr_pitch_slc");
  gr_pitch_slc->SetTitle("");
  gr_pitch_slc->GetXaxis()->SetTitle("Slice ID");
  gr_pitch_slc->GetYaxis()->SetTitle("Slice thickness (cm)");
  gr_pitch_slc->Draw("ap");

  TCanvas *c6 = new TCanvas("c6","c6");
  incE[5]->SetTitle(";Pion KE (MeV);N pions");
  incE[5]->Draw();
  incE[10]->SetLineColor(2);
  incE[10]->Draw("same");
  incE[15]->SetLineColor(3);
  incE[15]->Draw("same");
  TLegend *leg6 = new TLegend(0.15,0.7,0.6,0.9);
  leg6->SetFillStyle(0);
  leg6->AddEntry(incE[5],"Slice 5, 100<=wire #<120","l");
  leg6->AddEntry(incE[10],"Slice 10, 200<=wire #<220","l");
  leg6->AddEntry(incE[15],"Slice 15, 300<=wire #<320","l");
  leg6->Draw();

  TCanvas *c7 = new TCanvas("c7","c7");
  TGraph *gr_incE_slc = (TGraph*)f1.Get("gr_incE_slc");
  gr_incE_slc->SetTitle("");
  gr_incE_slc->GetXaxis()->SetTitle("Slice ID");
  gr_incE_slc->GetYaxis()->SetTitle("Pion KE (MeV)");
  gr_incE_slc->Draw("ap");

  TCanvas *c8 = new TCanvas("c8","c8");
  TGraph *gr_trueint_slc = (TGraph*)f1.Get("gr_trueint_slc");
  gr_trueint_slc->SetTitle("");
  gr_trueint_slc->GetXaxis()->SetTitle("Slice ID");
  gr_trueint_slc->GetYaxis()->SetTitle("N interactions");
  gr_trueint_slc->Draw("ap");

  TCanvas *c9 = new TCanvas("c9","c9");
  TGraph *gr_trueinc_slc = (TGraph*)f1.Get("gr_trueinc_slc");
  gr_trueinc_slc->SetTitle("");
  gr_trueinc_slc->GetXaxis()->SetTitle("Slice ID");
  gr_trueinc_slc->GetYaxis()->SetTitle("N incidents");
  gr_trueinc_slc->Draw("ap");

  TCanvas *c10 = new TCanvas("c10","c10");
  TGraphErrors *gr_dEdx_slc = (TGraphErrors*)f1.Get("gr_dEdx_slc");
  gr_dEdx_slc->SetTitle("");
  gr_dEdx_slc->GetXaxis()->SetTitle("Slice ID");
  gr_dEdx_slc->GetYaxis()->SetTitle("Calibrated dE/dx (MeV/cm)");
  gr_dEdx_slc->Draw("ap");


  c1->Print("inc.png");
  c2->Print("abs.png");
  c3->Print("cex.png");
  c1->Print("inc.pdf");
  c2->Print("abs.pdf");
  c3->Print("cex.pdf");

  c4->Print("slicet.pdf");
  c5->Print("slicet_id.pdf");
  c6->Print("pionKE.pdf");
  c7->Print("pionKE_id.pdf");
  c8->Print("int_id.pdf");
  c9->Print("inc_id.pdf");
  c10->Print("dEdx_id.pdf");
}
