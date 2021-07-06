#include <cstdlib>

#include <vector>

#include <iostream>

#include <map>

#include <string>

#include <fstream>

#include "TFile.h"

#include "TTree.h"

#include "TString.h"

#include "TSystem.h"

#include "TROOT.h"

#include "TStopwatch.h"

#include "TGraph.h"

#include "TCanvas.h"

#include "TMarker.h"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Tools.h"

#include "TMVA/Reader.h"

#include "TMVA/MethodCuts.h"

#endif

using namespace TMVA;

using namespace std;

void TMVAClassificationApplication_data( TString myMethodList = "" ) 

{   

#ifdef __CINT__

  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT

#endif

  //---------------------------------------------------------------

  // This loads the library

  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested

  std::map<std::string,int> Use;

  // --- Boosted Decision Trees

  //Use["BDTG"] = 1; // uses Gradient Boost
  Use["DNN_CPU"] = 1;
  std::cout << std::endl;

  std::cout << "==> Start TMVAClassificationApplication" << std::endl;

  // Select methods (don't look at this code - not of interest)

  if (myMethodList != "") {

    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );

    for (UInt_t i=0; i<mlist.size(); i++) {

      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {

	std::cout << "Method \"" << regMethod 

		  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;

	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {

	  std::cout << it->first << " ";

	}

	std::cout << std::endl;

	return;

      }

      Use[regMethod] = 1;

    }

  }

  // --- Create the Reader object
  //TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
  TMVA::Reader *reader = new TMVA::Reader();
  // Create a set of variables and declare them to the reader
  
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

  Int_t event;

  UInt_t run;

  Int_t subrun;

  Float_t ShwEng, ShwDist, TrackScore, Chi2;
  
  reader->AddVariable("ShwEng", &ShwEng);    

  reader->AddVariable("ShwDist", &ShwDist);  

  reader->AddVariable("TrackScore", &TrackScore);  

  reader->AddVariable("Chi2", &Chi2);  
  /*
  reader->AddVariable("count_calorimetric_pion_tree", &count_calorimetric_pion_tree);
  */
  ofstream output_data("BDTG_data_minus.txt");

  ofstream output_data_08("BDTG_data_minus_08.txt");
  
  // Book method(s)
  std::cout<<"Try to book method of BDTG"<<std::endl;

  string inputxmlname = "/dune/app/users/jiangl/dunetpc_analysis/newduneana_May/PionAbsChxAna/Main/mac/dataset/weights/TMVAClassification_DNN_CPU.weights.xml";
  //string inputxmlname = "/dune/app/users/jiangl/dunetpc_analysis/newduneana_May/PionAbsChxAna/Main/mac/dataset/weights/TMVAClassification_BDTG.weights.xml";
  string methodname = "DNN_CPU";
  //string methodname = "BDTG";

  reader->BookMVA(methodname, inputxmlname); // WEIGHT FILE IS GREEN FILE .XML
  //reader->BookMVA("BDTG", "/dune/app/users/jiangl/dunetpc_analysis/newduneana_May/PionAbsChxAna/Main/mac/dataset/weights/TMVAClassification_BDTG.weights.xml" ); // WEIGHT FILE IS GREEN FILE .XML
  std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;
  // Book output histograms

  UInt_t nbin = 40;

  //TH1F   *histBdtG(0);
  TH1F    *histBdtG;
  Float_t lowerv=0.0;
  Float_t upperv=1.0;
  if (Use[methodname]) histBdtG = new TH1F("BDTG_data_minus","BDTG_data_minus",nbin, lowerv, upperv);

  // Prepare input tree (this must be replaced by your data source)

  // in this example, there is a toy tree with signal and one with background events

  // we'll later on use only the "signal" events for the test in this example.

  //   

  TFile *input = TFile::Open("../output_test.root");

  std::cout << "--- Select signal sample" << std::endl;

  TTree* theTree = (TTree*)input->Get("TreeBackground");

  std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl; 

  theTree->SetBranchAddress("event", &event);

  theTree->SetBranchAddress("run", &run);

  theTree->SetBranchAddress("subrun", &subrun);

  theTree->SetBranchAddress("ShwEng", &ShwEng);

  theTree->SetBranchAddress("ShwDist", &ShwDist);

  theTree->SetBranchAddress("TrackScore", &TrackScore);

  theTree->SetBranchAddress("Chi2", &Chi2);
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  /*
  theTree->SetBranchAddress("trk_charge_minos_tree", &trk_charge_minos_tree);

  theTree->SetBranchAddress("Minimum_dedx_tree", &Minimum_dedx_tree);

  theTree->SetBranchAddress("nvertices_tree", &nvertices_tree);

  theTree->SetBranchAddress("ratio_cut1_charge1_tree", &ratio_cut1_charge1_tree);

  theTree->SetBranchAddress("ratio_cut2_charge1_tree", &ratio_cut2_charge1_tree);

  theTree->SetBranchAddress("count_calorimetric_pion_tree", &count_calorimetric_pion_tree);

  theTree->SetBranchAddress("vtx_x_tree", &vtx_x_tree);

  theTree->SetBranchAddress("vtx_y_tree", &vtx_y_tree);

  theTree->SetBranchAddress("vtx_z_tree", &vtx_z_tree);

  theTree->SetBranchAddress("reco_trk_around_vertex_tree", &reco_trk_around_vertex_tree);

  theTree->SetBranchAddress("trk_index_minos_tree", &trk_index_minos_tree);

  theTree->SetBranchAddress("trk_mom_minos_tree", &trk_mom_minos_tree);

  theTree->SetBranchAddress("trkstartdcosx_tree", &trkstartdcosx_tree);

  theTree->SetBranchAddress("trkstartdcosy_tree", &trkstartdcosy_tree);

  theTree->SetBranchAddress("trkstartdcosz_tree", &trkstartdcosz_tree);

  theTree->SetBranchAddress("trk_pi_code_tree", &trk_pi_code_tree);

  theTree->SetBranchAddress("trkstartdcosx_pi_tree", &trkstartdcosx_pi_tree);

  theTree->SetBranchAddress("trkstartdcosy_pi_tree", &trkstartdcosy_pi_tree);

  theTree->SetBranchAddress("trkstartdcosz_pi_tree", &trkstartdcosz_pi_tree);

  theTree->SetBranchAddress("pi_startx_tree", &pi_startx_tree);

  theTree->SetBranchAddress("pi_starty_tree", &pi_starty_tree);

  theTree->SetBranchAddress("pi_startz_tree", &pi_startz_tree);

  theTree->SetBranchAddress("pi_endx_tree", &pi_endx_tree);

  theTree->SetBranchAddress("pi_endy_tree", &pi_endy_tree);

  theTree->SetBranchAddress("pi_endz_tree", &pi_endz_tree);
  */

  std::cout<<"Start to loop over all the events and get the values of each variable"<<std::endl;
  for (Long64_t ievt = 0; ievt < theTree->GetEntries(); ievt++){

    theTree->GetEntry(ievt);
    //Float_t value = reader->EvaluateMVA("BDTG method");
    Float_t value = reader->EvaluateMVA(methodname);
    std::cout<<"evaluate the score from BDTG method, value is  "<<value<<std::endl;

    if(Use[methodname]){
      
      histBdtG->Fill(value);

      output_data << value << " "

		  << event << " "

		  << run << " "

		  << subrun << " "
		  /*	
		  << trk_charge_minos_tree << " "

		  << trk_mom_minos_tree << " "

		  << trkstartdcosx_tree << " "

		  << trkstartdcosy_tree << " "

		  << trkstartdcosz_tree << " "

		  << trkstartdcosx_pi_tree << " "

		  << trkstartdcosy_pi_tree << " "

		  << trkstartdcosz_pi_tree << " "
		  */
		  << "\n";
		  	
      if(value >= 0.8){

	output_data_08 << event << " "

		       << run << " "

		       << subrun << " "

		       << value << "  "

		       <<"\n";

      }               

    } 

  }

  std::cout << "--- End of event loop: "; 

  // --- Write histograms

  TFile *target  = new TFile( "TMVApp_data.root","RECREATE" );

  if (Use[methodname])   histBdtG   ->Write(); 

  std::cout << "--- Created root file: \"TMVApp_data.root\" containing the MVA output histograms" << std::endl;

  

  delete reader;

    

  std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;

  TCanvas *c1 = new TCanvas("c1","c1");

  histBdtG->Draw();

} 



