import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.PiAbsanalysis()


analyser.SetBNBCosmicFile("/dune/app/users/jiang/dunetpc_analysis/PionAbsChxAna/UBAna/Main/mac/output_test.root") # Tune 1 - full stat - full syst
analyser.SetBNBONFile("/dune/app/users/jiang/dunetpc_analysis/PionAbsChxAna/UBAna/Main/mac/output_test.root")    
#analyser.SetPrefix("cv")
# analyser.SetPrefix("cv_tune3")
# analyser.SetPrefix("cv_cosmicscaled_overlay_test")
# analyser.SetPrefix("cv_nodirt")
# analyser.SetPrefix("kaonup")
# analyser.SetPrefix("kaondown")
#analyser.SetFluxCorrectionWeight(1.028)

#analyser.ImportAlternativeMC("xsec_file_cv_tune3.root")

#analyser.SetBeamOffSubtraction(False)
#analyser.SetBreakdownPlots(False)

#extra_unc = math.sqrt(0.02*0.02) # POT counting
#analyser.SetExtraUncertainty(extra_unc)

#analyser.ImportDetectorSystematics(True)

#analyser.ImportCosmicSystematics(True)

#analyser.ImportDirtSystematics(True)

#analyser.DoGenieSystematics(False)
#analyser.ImportGenieSystematics(True)

#analyser.DoExtraSystematics(False)
#analyser.ImportExtraSystematics(True)

#analyser.DoMCStatSystematics(False)
#analyser.ImportMCStatSystematics(True)


#analyser.DoFluxSystematics(False)
#analyser.ImportFluxSystematics(True)
#analyser.SetExtraFluxUncertainty(0.)
#analyser.SetTargetFluxSystematic("total"); # Other options: "FluxUnisim", "kminus", "kplus", "kzero", "piminus", "piplus"

analyser.DoPiAbsAnalyze();
input("Please press enter to exit.")
#raw_input("Please press enter to exit.")
