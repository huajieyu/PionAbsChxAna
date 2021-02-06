/**
 * \file PiAbsanalysis.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class PiAbsanalysis
 *
 * @author jiangl
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_PIABSANALYSIS_H__
#define __MAIN_PIABSANALYSIS_H__

#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <math.h>

#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TH2Poly.h>

#include "../DataTypes/UBTH2Poly.h"
#include "../DataTypes/BootstrapTH2DPoly.h"

#include "../DataTypes/UBXSecEventHisto1D.h"
#include "../DataTypes/UBXSecEventHisto.h"

#include "UBXSecEvent.h"
#include "../DataTypes/BootstrapTH1D.h"
#include "../DataTypes/BootstrapTH2D.h"
#include "../Base/PlottingTools.h"
#include "../Base/CrossSectionCalculator1D.h"
#include "../Base/MigrationMatrix2D.h"
#include "../Base/MigrationMatrix4D.h"
#include "../Base/MigrationMatrix4DPoly.h"
#include "../Base/CrossSectionCalculator2D.h"
#include "../Base/ReweightingPlotter.h"
#include "../Base/CovarianceCalculator2D.h"
#include "../Base/CrossSectionBootstrapCalculator1D.h"
#include "../Base/CrossSectionBootstrapCalculator2D.h"
#include "../Base/CrossSectionBootstrapCalculator2DPoly.h"
#include "../Base/CrossSectionCalculator2DPoly.h"
#include "../Base/UncertaintyPlotter.h"

#include "../Base/LoggerFeature.h"


using namespace DataTypes;
using namespace Base;







namespace Main {

  /**
     \class PiAbsanalysis
     User defined class PiAbsanalysis ... these comments are used to generate
     doxygen documentation!
  */
  class PiAbsanalysis : public LoggerFeature{
    
  public:
    
    /// Default constructor
    PiAbsanalysis(std::string name = "PiAbsanalysis") 
    : LoggerFeature(name) {}
     
    /// Default destructor
    ~PiAbsanalysis(){}

    void DoPiAbsAnalyze();
    void SetBNBCosmicFile(std::string f);
    void SetBNBONFile(std::string f);
    void stackHists(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac);
    void subtractBKG(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac);

  private:

    std::string mc_bnbcosmic_file_name     = "ubxsecana_output.root";
    std::string bnbon_file_name            = "ubxsecana_output.root";
    double NA=6.02214076e23;
    double MAr=35.95; //gmol
    double Density = 1.39; // g/cm^3


     
  };
}

#endif
/** @} */ // end of doxygen group 

