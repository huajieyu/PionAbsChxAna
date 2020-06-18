#ifndef __BASE_REWEIGHTINGPLOTTER_CXX__
#define __BASE_REWEIGHTINGPLOTTER_CXX__

#include "ReweightingPlotter.h"

namespace Base {

	void ReweightingPlotter::SetEventBootstrapMap(std::map<std::string, BootstrapTH1D> map_bs)
  {
    _map_bs = map_bs;

    _configured = true;
  }

  void ReweightingPlotter::SetEfficiencyBootstraps(BootstrapTH1D eff_num, BootstrapTH1D eff_den)
  {
    _bs_eff_num = eff_num;
    _bs_eff_den = eff_den;
  }

  void ReweightingPlotter::SetXSecBootstrap(BootstrapTH1D bs_xsec) 
  {

    _bs_xsec = bs_xsec;

  }




  void ReweightingPlotter::MakePlots(int variable, bool normalised, bool makeLaTeX) 
  {

    if (!_configured) {
      std::cout << _name << "Not configured." << std::endl;
      throw std::exception();
    }

    TH1D *histo; // the nominal histogram
    std::map<std::string, std::vector<TH1D>> histo_map; // a map from function name to 2 TH1D (+ and - 1 sigma)
    std::map<std::string, std::vector<TH1D>> histo_map_den; // (only for eff plot) a map from function name to 2 TH1D (+ and - 1 sigma)
    TH1D *histo_p1; // for each function, this will reoresent the p1 histo
    TH1D *histo_m1; // for each function, this will reoresent the m1 histo

    double efficiency_nominal;
    double efficiency_p1;
    double efficiency_m1;

    if (variable == 2 || variable == 3) {

      BootstrapTH1D bs = _map_bs["total"];

      histo = (TH1D*) (bs.GetNominal().Clone("histo_hateroot"));

      histo_map = bs.UnpackPMHisto();

    }

    if (variable == 0 || variable == 1) {

      // Construct nominal histogram
      //histo->Reset();
      histo = (TH1D*) _bs_eff_num.GetNominal().Clone("nominal_eff");
      TH1D * temp_eff_den = (TH1D*)_bs_eff_den.GetNominal().Clone("temp_eff_den");
      histo->Divide(temp_eff_den);
      efficiency_nominal = _bs_eff_num.GetNominal().Integral() / _bs_eff_den.GetNominal().Integral();

      // Numerator for efficiency calculation 
      histo_map = _bs_eff_num.UnpackPMHisto();
      histo_map_den = _bs_eff_den.UnpackPMHisto();

    }



    // Make a directory to store the plots
    if (variable == 0) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 1) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 2) system("mkdir ./EvtWgtEventPlots");
    if (variable == 3) system("mkdir ./EvtWgtEventPlots");
  
    // Avoid root to dislay the canvases
    //gROOT->SetBatch(kTRUE);
  
    // Opening a text file to write the integrals
    std::ofstream outfile;
    if (variable == 0) outfile.open("./EvtWgtEfficiencyPlots/IntegralsPmom.txt");
    if (variable == 1) outfile.open("./EvtWgtEfficiencyPlots/IntegralsPcostheta.txt");
    if (variable == 2) outfile.open("./EvtWgtEventPlots/IntegralsPmom.txt");
    if (variable == 3) outfile.open("./EvtWgtEventPlots/IntegralsPcostheta.txt");
  
  
  
  
    // Open LaTeX file to write the table
    std::ofstream latexFile;
    std::ofstream latexFile3;
    if(makeLaTeX) {
      if (variable == 0) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPmom.tex");
      if (variable == 1) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPcostheta.tex");
      if (variable == 2) latexFile.open("./EvtWgtEventPlots/evtwgtEventPmom.tex");
      if (variable == 3) latexFile.open("./EvtWgtEventPlots/evtwgtEventPcostheta.tex");
      latexFile << "\\begin{table}[]" << std::endl;
      latexFile << "\\caption{}" << std::endl;
      latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << std::endl;
      latexFile << "\\label{tab:}" << std::endl;
      latexFile << "\\centering" << std::endl;
      latexFile << "\\begin{tabular}{ccc}" << std::endl;
      latexFile << "\\toprule" << std::endl;
      if (variable == 0 || variable == 1) latexFile << "  &  Efficiency  &  Difference (\\%) \\\\" << std::endl;
      else latexFile << "  &  Integral  &  Difference (\\%) \\\\" << std::endl;
      latexFile << "\\midrule" << std::endl;

      if (variable == 0) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPmom_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 1) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPcostheta_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 2) latexFile3.open("./EvtWgtEventPlots/evtwgtEventPmom_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 3) latexFile3.open("./EvtWgtEventPlots/evtwgtEventPcostheta_figures.tex", std::ofstream::out | std::ofstream::trunc);
      latexFile3 << "\\begin{figure}[t]" << std::endl;
      latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
      latexFile3 << "\\centering" << std::endl;
    }


    bool is_first = true;

    int function_counter = 0;

    for(auto iter : histo_map) {

      function_counter++;

      std::string function_name = iter.first;

      //std::cout << "This is function name: " << function_name << std::endl;


      if (variable == 2 || variable == 3) {

        histo_p1 = (TH1D*) iter.second.at(1).Clone("histo_p1_hateroot");
        histo_m1 = (TH1D*) iter.second.at(0).Clone("histo_m1_hateroot");

      }

      if (variable == 0 || variable == 1) {

        //histo_p1->Reset();
        histo_p1 = (TH1D*) iter.second.at(1).Clone("p1");
        TH1D * temp_eff_den_p1 = (TH1D*)histo_map_den[function_name].at(1).Clone("temp_eff_den_p1");
        histo_p1->Divide(temp_eff_den_p1);

        efficiency_p1 = iter.second.at(1).Integral() / histo_map_den[function_name].at(1).Integral();

        //histo_m1->Reset();
        histo_m1 = (TH1D*) iter.second.at(0).Clone("m1");
        TH1D * temp_eff_den_m1 = (TH1D*)histo_map_den[function_name].at(0).Clone("temp_eff_den_m1");
        histo_m1->Divide(temp_eff_den_m1);

        efficiency_m1 = iter.second.at(0).Integral() / histo_map_den[function_name].at(0).Integral();

      }



      TString SaveName;
      if(variable == 0 || variable == 2) SaveName = "Pmom_"+function_name;
      if(variable == 1 || variable == 3) SaveName = "Pcostheta_"+function_name;
      TString LegName  = GetLegendName(function_name);

      if(normalised) SaveName += "_normalised";


      double histo_Int = histo->Integral();
      double histo_p1_Int = histo_p1->Integral();
      double histo_m1_Int = histo_m1->Integral();


      if (normalised) {
        histo->Scale(1./histo_Int);
        histo_p1->Scale(1./histo_p1_Int);
        histo_m1->Scale(1./histo_m1_Int);
      }

      // Calculate integrals
      outfile << function_name << std::endl;
      outfile << "Integral Nominal:    " << histo->Integral() << std::endl;
      outfile << "Integral nominal_p1: " << histo_p1->Integral() << std::endl;
      outfile << "Integral nominal_m1: " << histo_m1->Integral() << std::endl;

      outfile << "Difference w.r.t. Nominal (%):" << std::endl;
      outfile << "nominal_p1: " << (histo_p1->Integral()-histo->Integral())/(histo->Integral())*100. << std::endl;
      outfile << "nominal_m1: " << (histo_m1->Integral()-histo->Integral())/(histo->Integral())*100. << std::endl;
      outfile << "--------------------------------------" << std::endl << std::endl;


      if(makeLaTeX) {
        if ( (variable == 2) || (variable == 3) ) {
          if (is_first) latexFile << "Nominal" << " & " << histo->Integral() << " & 0" << "\\\\" << std::endl;
          latexFile << "\\midrule" << std::endl;
          latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ & " << histo_p1->Integral() << " & " << (histo_p1->Integral()-histo->Integral())/(histo->Integral())*100. << "\\\\" << std::endl;
          latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ & " << histo_m1->Integral() << " & " << (histo_m1->Integral()-histo->Integral())/(histo->Integral())*100. << "\\\\" << std::endl;
        }

        if (variable == 0 || variable == 1) {  // save the value of the efficiency to the LaTeX file

          double eff_nom = efficiency_nominal * 100.;
          double eff_p1  = efficiency_p1      * 100.;
          double eff_m1  = efficiency_m1      * 100.;

          if (is_first) latexFile << "Nominal" << " & " << std::setprecision(4) << eff_nom << " & 0" << "\\\\" << std::endl;
          latexFile << "\\midrule" << std::endl;
        
          latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ & "
          << std::setprecision(4) << eff_p1 << " & "
          << std::setprecision(4) << (eff_p1-eff_nom)/eff_nom*100. << "\\\\" << std::endl;
        
          latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ & "
          << std::setprecision(4) << eff_m1 << " & "
          << std::setprecision(4) << (eff_m1-eff_nom)/eff_nom*100. << "\\\\" << std::endl;
        }

      }



      // Define the Canvas
      TCanvas *c = new TCanvas("c", "canvas", 800, 800);



      // Upper plot will be in pad1
      TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
      pad1->SetBottomMargin(0); // Upper and lower plot are joined
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();               // pad1 becomes the current pad
      if (variable == 0 || variable == 1) histo_p1->SetMaximum(1.);
      histo_p1->SetStats(0);          // No statistics on upper plot
      histo_m1->SetStats(0);          // No statistics on upper plot
      histo->SetStats(0);          // No statistics on upper plot
      histo_p1->Draw("histo");               // Draw h1
      histo->Draw("histo same");         // Draw h2 on top of h1
      histo_m1->Draw("histo same");
    
      //uBooNESimulation();
    
      if (normalised) {
        // TLatex
        double x;
        x = 0.87;
        if (variable == 3) x = 0.46;
        double y = 0.52;
        double size = 28;
        int color = 1;
        int font = 43;
        int align = 32;
        TLatex *latex = new TLatex( x, y, "Area Normalised" );
        latex->SetNDC();
        latex->SetTextSize(size);
        latex->SetTextColor(color);
        latex->SetTextFont(font);
        latex->SetTextAlign(align);
      
        latex->Draw();
      
      }

      // Do not draw the Y axis label on the upper plot and redraw a small
      // axis instead, in order to avoid the first label (0) to be clipped.
      histo->GetYaxis()->SetLabelSize(0.);
      TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
      axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      axis->SetLabelSize(15);
      axis->Draw();
      
      // Legend for the upper plot
      TLegend* leg;
      if (variable == 0 || variable == 2) leg = new TLegend(0.65,0.6,.85,0.87);
      if (variable == 1 || variable == 3) leg = new TLegend(0.216792,0.5723502,0.4172932,0.843318,NULL,"brNDC");
      leg->SetTextFont(42);
      leg->SetBorderSize(0);
      //leg->SetHeader("");
      //leg->SetTextFont(42);
      //leg->AddEntry(histo, "BBA2005 (Nominal)");
      //leg->AddEntry(histo_p1, "Dipole");
      leg->AddEntry(histo, "Nominal");
      leg->AddEntry(histo_p1, LegName + " + 1#sigma");
      leg->AddEntry(histo_m1, LegName + " - 1#sigma");
      leg->Draw();
    
      // lower plot will be in pad
      c->cd();          // Go back to the main canvas before defining pad2
      TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.5);
      pad2->SetGridx(); // vertical grid
      //pad2->SetGridy(); // orizontal grid
      pad2->Draw();
      pad2->cd();       // pad2 becomes the current pad
    
    
      // Define the first ratio plot
      TH1D *ratio_p1 = (TH1D*)histo_p1->Clone("ratio_p1");
      //ratio_p1->SetMinimum(0.92);  // Define Y ..
      //ratio_p1->SetMaximum(1.08); // .. range
      //ratio_p1->Sumw2();
      ratio_p1->SetStats(0);      // No statistics on lower plot
      ratio_p1->Add(histo, -1.0);
      ratio_p1->Divide(histo);
      //ratio_p1->Add(histo, -1.0);
      ratio_p1->SetLineWidth(2);
      ratio_p1->SetLineColor(kRed+1);
      //ratio_p1->Draw("hist");       // Draw the ratio plot
    
      // Define the second ratio plot
      TH1D *ratio_m1 = (TH1D*)histo_m1->Clone("ratio_m1");
      //ratio_m1->SetMinimum(0.9);  // Define Y ..
      //ratio_m1->SetMaximum(1.1); // .. range
      //ratio_m1->Sumw2();
      ratio_m1->SetStats(0);      // No statistics on lower plot
      ratio_m1->Add(histo, -1.0);
      ratio_m1->Divide(histo);
      //ratio_m1->Add(histo, -1.0);
      ratio_m1->SetLineWidth(2);
      ratio_m1->SetLineColor(kGreen+2);
      //ratio_m1->Draw("hist same");       // Draw the ratio plot
    
    
    
      // Try to set the Y range for the ratio plots
      double max = ratio_p1->GetMaximum();
      double min = ratio_p1->GetMinimum();
      if (ratio_m1->GetMaximum() > max) max = ratio_m1->GetMaximum();
      if (ratio_m1->GetMinimum() < min) min = ratio_m1->GetMinimum();
      //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
      //ratio_p1->SetMinimum(max);  // Define Y ..
      //ratio_p1->SetMaximum(min); // .. range
      //ratio_m1->SetMinimum(max);  // Define Y ..
      //ratio_m1->SetMaximum(min); // .. range

    
      // Draw the ratio plot
      //ratio_p1->Draw("hist");
      //ratio_m1->Draw("hist same");
    
      THStack *hs = new THStack("hs","");
      hs->Add(ratio_p1);
      hs->Add(ratio_m1);
      hs->SetMaximum(hs->GetMaximum("nostack")+0.01);
      hs->SetMinimum(hs->GetMinimum("nostack")-0.01);



      hs->Draw("NOSTACK histo");
    
    
      //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
    
       histo_p1->SetMinimum(0.0001); // Otherwise 0 label overlaps (need to do it after the THStack, otherwise sets the minimum)


      //**********************
      //
      // Settings
      //
      //**********************

      // h1 settings
      histo->SetLineColor(kBlack);
      histo->SetLineWidth(2);

      // Y axis h1 plot settings
      histo_p1->GetYaxis()->SetTitle("Selected Events");
      if (variable == 0 || variable == 1) histo_p1->GetYaxis()->SetTitle("Efficiency");
      histo_p1->GetYaxis()->CenterTitle();
      histo_p1->GetYaxis()->SetTitleSize(25);
      histo_p1->GetYaxis()->SetTitleFont(43);
      histo_p1->GetYaxis()->SetTitleOffset(1.55);

      // h2 settings
      histo_p1->SetLineColor(kRed+1);
      histo_p1->SetLineWidth(2);

      // h3 settings
      histo_m1->SetLineColor(kGreen+2);
      histo_m1->SetLineWidth(2);

      // Ratio plot (ratio_p1) settings
      ratio_p1->SetTitle(""); // Remove the ratio title
    
    
      hs->GetYaxis()->SetTitle("Ratio");
      hs->GetYaxis()->CenterTitle();
      hs->GetYaxis()->SetNdivisions(505);
      hs->GetYaxis()->SetTitleSize(25);
      hs->GetYaxis()->SetTitleFont(43);
      hs->GetYaxis()->SetTitleOffset(1.2);
      hs->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      hs->GetYaxis()->SetLabelSize(15);
    
    
      if(variable == 0 || variable == 2) hs->GetXaxis()->SetTitle("P_{proton} [GeV]");
      if(variable == 1 || variable == 3) hs->GetXaxis()->SetTitle("cos#theta_{proton}");
      hs->GetXaxis()->CenterTitle();
      hs->GetXaxis()->SetTitleSize(25);
      hs->GetXaxis()->SetTitleFont(43);
      hs->GetXaxis()->SetTitleOffset(3.5);
      hs->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      hs->GetXaxis()->SetLabelSize(20);
    
      // Draw linea at 1 in ratio plot
      TLine *line;
      if (variable == 0 || variable == 2) line = new TLine(0,1,2.5,1);
      if (variable == 1 || variable == 3) line = new TLine(-1,1,1,1);
      line->SetLineColor(kBlack);
      line->SetLineStyle(9); // dashed
      line->Draw();
    
      if (variable == 0 || variable == 1) {
        c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".C");
         c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".pdf");
      }
      if (variable == 2 || variable == 3) {
        c->Print("./EvtWgtEventPlots/" + SaveName + ".C");
        c->Print("./EvtWgtEventPlots/" + SaveName + ".pdf");
      }

      if (makeLaTeX) {
        latexFile3 << "\\subfloat[][]" << std::endl;
        latexFile3 << "   {\\includegraphics[width=.35\\textwidth]{images/EvtWgtEfficiencyPlots/" << SaveName << "}" << std::endl;
        latexFile3 << "   \\label{fig:" << "EvtWgtEfficiencyPlots_" << SaveName << "}} \\quad" << std::endl;
      }

      if (makeLaTeX && function_counter % 12 == 0) {
        latexFile3 << "\\caption{Efficiency Plots}" << std::endl;
        latexFile3 << "\\label{fig:EvtWgtEfficiencyPlots}" << std::endl;
        latexFile3 << "\\end{adjustwidth}" << std::endl;
        latexFile3 << "\\end{figure}" << std::endl;

        latexFile3 << "\\begin{figure}[t]" << std::endl;
        latexFile3 << "\\ContinuedFloat" << std::endl;
        latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
        latexFile3 << "\\centering" << std::endl;

      }


      is_first = false;

    } // end loop functions
  
    if(makeLaTeX) {

      latexFile << "\\bottomrule" << std::endl;
      latexFile << "\\end{tabular}" << std::endl;
      latexFile << "\\end{table}" << std::endl;


      latexFile3 << "\\caption{Efficiency Plots}" << std::endl;
      latexFile3 << "\\label{fig:EvtWgtXsecDiffPlots}" << std::endl;
      latexFile3 << "\\end{adjustwidth}" << std::endl;
      latexFile3 << "\\end{figure}" << std::endl;
    }



  }

  void ReweightingPlotter::MakeBackgroundPlots(int variable, bool normalised, bool makeLaTeX) {
    

  
  // Create output directory
  system("mkdir ./EvtWgtBackgroundPlots");
  system("mkdir ./EvtWgtBackgroundPlots_reducedLegend"); 
 
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  
  double scaleFactor = 1.;//5.3e19/1.22e20;
  std::cout << "scaleFactor = " <<scaleFactor<<std::endl;
  std::cout << "Pion Absorption:   " << _map_bs["signal"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "Charge Exchage:  " << _map_bs["chxbac"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "Pion Reaction:    " << _map_bs["reabac"].GetNominal().Integral() * scaleFactor << std::endl;

  std::ofstream latexFile;
  std::ofstream latexFile3;





  double numu_nominal = _map_bs["signal"].GetNominal().Integral();
  double anumu_nominal = _map_bs["chxbac"].GetNominal().Integral();
  double nue_nominal = _map_bs["reabac"].GetNominal().Integral();

  if(makeLaTeX) {
    if (variable == 0) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPmom.tex");
    if (variable == 1) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPcostheta.tex");
    latexFile << "\\begin{table}[]" << std::endl;
    latexFile << "\\begin{adjustwidth}{-2.1cm}{-1cm}" << std::endl;
    latexFile << "\\caption{}" << std::endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << std::endl;
    latexFile << "\\label{tab:}" << std::endl;
    latexFile << "\\centering" << std::endl;
    latexFile << "\\tiny" << std::endl;
    latexFile << "\\begin{tabular}{c   c c  c  }" << std::endl;
    latexFile << "\\toprule" << std::endl;
    latexFile << "  &  \\multicolumn{2}{ c }{PiAbs} &  \\multicolumn{2}{ c }{$PiChx}   & \\multicolumn{2}{ c }{PiRea}   \\\\" << std::endl;
    latexFile << "  &  Events  &  Diff. (\\%) &  Events  &  Diff. (\\%) & Events  &  Diff. (\\%)  \\\\" << std::endl;
    latexFile << "\\midrule" << std::endl;
    latexFile << "Nominal & " << _map_bs["signal"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["chxbac"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["reabac"].GetNominal().Integral() << " & 0 " << "\\\\" << std::endl;

    if (variable == 0) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPmom_figures.tex", std::ofstream::out | std::ofstream::trunc);
    if (variable == 1) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPcostheta_figures.tex", std::ofstream::out | std::ofstream::trunc);
    latexFile3 << "\\begin{figure}[t]" << std::endl;
    latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
    latexFile3 << "\\centering" << std::endl;
  }
  
  
  
  
  
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  
  std::map<std::string, std::vector<TH1D>> map_bs_signal = _map_bs["signal"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_anumu = _map_bs["chxbac"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_nue = _map_bs["reabac"].UnpackPMHisto();

  int function_counter = 0;

  for (auto iter : map_bs_signal) {

    function_counter++;

    std::string function_name = iter.first;

    auto ii = map_bs_anumu.find(function_name);
    std::vector<TH1D> anumu_v = ii->second;
    const TH1D* anumu_nom = &_map_bs["chxbac"].GetNominal();

    ii = map_bs_nue.find(function_name);
    std::vector<TH1D> nue_v = ii->second;
    const TH1D* nue_nom = &_map_bs["reabac"].GetNominal();

    
    TString SaveName;
    if (variable == 0) SaveName = "Pmom_" + function_name;
    if (variable == 1) SaveName = "PCosTheta_" + function_name;
    
    c->SetLogy();                                 // Log scale
    //c->SetGridy();                                // Horizontal grid
    
    
    TH1D* histo_numu;
    TH1D* histo_numu_p1;
    TH1D* histo_numu_m1;
    TH1D* histo_anumu;
    TH1D* histo_anumu_p1;
    TH1D* histo_anumu_m1;
    TH1D* histo_nue;
    TH1D* histo_nue_p1;
    TH1D* histo_nue_m1;


    
    if (variable == 0 || variable == 1) {
      histo_numu        = (TH1D*)_map_bs["signal"].GetNominal().Clone("histo_numu");
      histo_numu_p1     = (TH1D*)iter.second.at(1).Clone("histo_numu_p1");
      histo_numu_m1     = (TH1D*)iter.second.at(0).Clone("histo_numu_m1");
      
      histo_anumu       = (TH1D*) (anumu_nom->Clone("histo_anumu"));
      histo_anumu_p1    = (TH1D*) (anumu_v.at(1).Clone("histo_anumu_p1"));
      histo_anumu_m1    = (TH1D*) (anumu_v.at(0).Clone("histo_anumu_m1"));
      
      histo_nue         = (TH1D*) (nue_nom->Clone("histo_nue"));
      histo_nue_p1      = (TH1D*) (nue_v.at(1).Clone("histo_nue_p1"));
      histo_nue_m1      = (TH1D*) (nue_v.at(0).Clone("histo_nue_m1"));
      
    }
    else {
      std::cout << "Invalid option. Exit." << std::endl;
      exit(0);
    }
    
    
    // Settings
    histo_numu->SetStats(0);          // No statistics on upper plot
    //if (variable == 1) histo_numu->SetMinimum(1);
    //if (variable == 1) histo_numu->SetMaximum(1e5);

    if (variable == 0) histo_numu->GetXaxis()->SetTitle("Reconstructed p_{proton} [GeV]");
    if (variable == 1) histo_numu->GetXaxis()->SetTitle("Reconstructed cos#theta_{proton}");
    histo_numu->GetXaxis()->CenterTitle();
    histo_numu->GetXaxis()->SetTitleSize(25);
    histo_numu->GetXaxis()->SetTitleFont(43);
    histo_numu->GetXaxis()->SetTitleOffset(1.45);
    
    histo_numu->GetYaxis()->SetTitle("Events");
    histo_numu->GetYaxis()->CenterTitle();
    histo_numu->GetYaxis()->SetTitleSize(25);
    histo_numu->GetYaxis()->SetTitleFont(43);
    histo_numu->GetYaxis()->SetTitleOffset(1.55);
    
    double lineWidth = 1;
    
    // Nominal
    histo_numu   ->Draw("histo");
    histo_anumu  ->Draw("histo same");
    histo_nue    ->Draw("histo same");

    histo_numu   ->SetLineColor(kRed+1);
    histo_anumu  ->SetLineColor(kOrange+1);
    histo_nue    ->SetLineColor(kViolet+1);
    
    histo_numu   ->SetLineWidth(lineWidth+1);
    histo_anumu  ->SetLineWidth(lineWidth+1);
    histo_nue    ->SetLineWidth(lineWidth+1);

    // +- 1 sigma
    histo_numu_p1   ->Draw("histo same");
    histo_anumu_p1  ->Draw("histo same");
    histo_nue_p1    ->Draw("histo same");

    histo_numu_m1   ->Draw("histo same");
    histo_anumu_m1  ->Draw("histo same");
    histo_nue_m1    ->Draw("histo same");

    histo_numu_p1  ->SetLineColor(kRed+1);
    histo_anumu_p1 ->SetLineColor(kOrange+1);
    histo_nue_p1   ->SetLineColor(kViolet+1);

    histo_numu_m1  ->SetLineColor(kRed+1);
    histo_anumu_m1 ->SetLineColor(kOrange+1);
    histo_nue_m1   ->SetLineColor(kViolet+1);

    histo_numu_p1  ->SetLineStyle(7);
    histo_anumu_p1 ->SetLineStyle(7);
    histo_nue_p1   ->SetLineStyle(7);

    histo_numu_m1  ->SetLineStyle(3);
    histo_anumu_m1 ->SetLineStyle(3);
    histo_nue_m1   ->SetLineStyle(3);

    histo_numu_p1  ->SetLineWidth(lineWidth);
    histo_anumu_p1 ->SetLineWidth(lineWidth);
    histo_nue_p1   ->SetLineWidth(lineWidth);

    histo_numu_m1  ->SetLineWidth(lineWidth);
    histo_anumu_m1 ->SetLineWidth(lineWidth);
    histo_nue_m1   ->SetLineWidth(lineWidth);

    
    // Legend
    TLegend* leg;
    if (variable == 0) leg = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1) leg = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    
    leg->AddEntry(histo_numu,  "Pion Absorption (nominal)");
    leg->AddEntry(histo_anumu, "Charge Exchange (nominal)");
    leg->AddEntry(histo_nue,   "Pion Reaction (nominal)");

    leg->AddEntry(histo_numu_p1,   "PiAbs (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_anumu_p1,  "PiChx (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_nue_p1,    "PiRea (" + GetLegendName(function_name) + " + 1 #sigma)");

    leg->AddEntry(histo_numu_m1,   "PiAbs (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_anumu_m1,  "PiChx (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_nue_m1,    "PiRea (" + GetLegendName(function_name) + " - 1 #sigma)");

    leg->Draw();
    
    
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".pdf");
    

    // Reduced Legend Plots
    histo_numu->GetXaxis()->SetTitleOffset(1.10);
    histo_numu->GetYaxis()->SetTitleOffset(1.30);
    histo_numu->GetXaxis()->SetTitleSize(30);
    histo_numu->GetYaxis()->SetTitleSize(30);

    histo_numu  ->Draw("histo");
    histo_anumu ->Draw("histo same");
    histo_nue   ->Draw("histo same");

    histo_numu_p1  ->Draw("histo same");
    histo_anumu_p1 ->Draw("histo same");
    histo_nue_p1   ->Draw("histo same");

    histo_numu_m1  ->Draw("histo same");
    histo_anumu_m1 ->Draw("histo same");
    histo_nue_m1   ->Draw("histo same");

    TLegend* leg2;
    if (variable == 0) leg2 = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1) leg2 = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);

    leg2->AddEntry(histo_numu,  "Signal");
    leg2->AddEntry(histo_anumu, "Charge Exchange");
    leg2->AddEntry(histo_nue,   "Pion Reaction");

    leg2->Draw();

    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".pdf");

    if (makeLaTeX) {
        latexFile3 << "\\subfloat[][]" << std::endl;
        latexFile3 << "   {\\includegraphics[width=.35\\textwidth]{images/EvtWgtBackgroundPlots/" << SaveName << "}" << std::endl;
        latexFile3 << "   \\label{fig:" << "EvtWgtEfficiencyPlots_" << SaveName << "}} \\quad" << std::endl;
    }

    if (makeLaTeX && function_counter % 12 == 0) {
      latexFile3 << "\\caption{Signal and Background Plots}" << std::endl;
      latexFile3 << "\\label{fig:EvtWgtBackgroundPlotss}" << std::endl;
      latexFile3 << "\\end{adjustwidth}" << std::endl;
      latexFile3 << "\\end{figure}" << std::endl;

      latexFile3 << "\\begin{figure}[t]" << std::endl;
      latexFile3 << "\\ContinuedFloat" << std::endl;
      latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
      latexFile3 << "\\centering" << std::endl;

    }




    double numu_p1;
    double numu_m1;
    double anumu_p1;
    double anumu_m1;
    double nue_p1;
    double nue_m1;

    if (variable == 0) {
      numu_p1   = histo_numu_p1   ->Integral();
      numu_m1   = histo_numu_m1   ->Integral();
      anumu_p1  = histo_anumu_p1  ->Integral();
      anumu_m1  = histo_anumu_m1  ->Integral();
      nue_p1    = histo_nue_p1    ->Integral();
      nue_m1    = histo_nue_p1    ->Integral();
    }
     
    
    if(makeLaTeX) {
      double numu_diff_p1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      double numu_diff_m1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      
      double anumu_diff_p1 = (anumu_p1 -anumu_nominal)/anumu_nominal * 100.;
      double anumu_diff_m1 = (anumu_m1 -anumu_nominal)/anumu_nominal * 100.;

      double nue_diff_p1 = (nue_p1 -nue_nominal)/nue_nominal * 100.;
      double nue_diff_m1 = (nue_m1 -nue_nominal)/nue_nominal * 100.;
      

      latexFile << "\\midrule" << std::endl;
      latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_p1  << " & " << std::setprecision(2) << numu_diff_p1
      << " & " << std::setprecision(5) << anumu_p1 << " & " << std::setprecision(2) << anumu_diff_p1
      << " & " << std::setprecision(5) << nue_p1   << " & " << std::setprecision(2) << nue_diff_p1
      << "\\\\" << std::endl;
      
      latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_m1  << " & " << std::setprecision(2) << numu_diff_m1
      << " & " << std::setprecision(5) << anumu_m1 << " & " << std::setprecision(2) << anumu_diff_m1
      << " & " << std::setprecision(5) << nue_m1   << " & " << std::setprecision(2) << nue_diff_m1
      << "\\\\" << std::endl;
    }
  }
  
  
  
  
  if(makeLaTeX) {
    latexFile << "\\bottomrule" << std::endl;
    latexFile << "\\end{tabular}" << std::endl;
    latexFile << "\\end{adjustwidth}" << std::endl;
    latexFile << "\\end{table}" << std::endl;

    latexFile3 << "\\caption{Signal and Background Plots}" << std::endl;
    latexFile3 << "\\label{fig:EvtWgtBackgroundPlots}" << std::endl;
    latexFile3 << "\\end{adjustwidth}" << std::endl;
    latexFile3 << "\\end{figure}" << std::endl;
  }
 

  }


  void ReweightingPlotter::MakeXsecDiffPlots(bool makeLaTeX) {




  }

  TString ReweightingPlotter::GetLegendName(std::string fName) {
    
    TString legName = "null";
    if (fName.find("fReac") != std::string::npos) legName = "fReac";
    if (fName.find("fAbs1") != std::string::npos) legName = "fAbs1";
    if (fName.find("fAbs2") != std::string::npos) legName = "fAbs2";
    if (fName.find("fCex") != std::string::npos) legName = "fCex";
  
    return legName;
  }
  
}

#endif
