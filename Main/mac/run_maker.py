import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()




maker.SetInputFile("/dune/data/users/calcuttj/pduneana_Prod4a_1GeV_5_14_21.root");
maker.SetOutputFile("/dune/app/users/jiang/dunetpc_analysis/newduneana_May/PionAbsChxAna/Main/mac/output_test.root") # Run1 After Neutrino
#maker.SetOutputFile("/dune/app/users/jiang/dunetpc_analysis/UBAna/Main/mac/output_test_data.root") # Run1 After Neutrino

maker.SetEntries(-1)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(5.0)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(False)

maker.SetSignalTypeAbs(True)
maker.SetSignalTypeChx(False)

maker.SetFVCut(True);
maker.FillBootstrapGeant(False);
maker.SetMomThreshCut(0.0);
maker.SetTrkScoreCut(0.4);
maker.SetFVCut(True);



maker.PrintConfig()

maker.MakeFile()




