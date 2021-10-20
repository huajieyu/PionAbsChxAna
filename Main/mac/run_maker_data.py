import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()




#maker.SetInputFile("/dune/data/users/calcuttj/pionana_mc_1GeV_5_29_20.root");
#maker.SetInputFile("/dune/data/users/calcuttj/pionana_5387_5_28_20.root") # Run1 After Neutrino

#maker.SetInputFile("/dune/data/users/calcuttj/pionana_5387_6_15_20.root");
#maker.SetInputFile("/dune/data/users/calcuttj/pionana_5387_6_19_20.root");
#maker.SetInputFile("/dune/data/users/calcuttj/pionana_Prod3_5387_1GeV_7_30_20.root");
#maker.SetInputFile("/dune/data/users/calcuttj/pionana_Prod3_5387_1GeV_9_20_20.root");
#maker.SetInputFile("/dune/data/users/fstocker/ANALYSIS/prod4_files/pionana_Prod4_5387_1GeV_12_16_20.root");
#maker.SetInputFile("/dune/data/users/calcuttj/pduneana_Prod4_1GeV_5387_5_12_21.root");
maker.SetInputFile("/pnfs/dune/scratch/users/yinrui/pduneana_Prod4_ori/6_30_21_5387/output_sce_1GeV/pduneana_data.root");
maker.SetOutputFile("/dune/app/users/mmatt15/PionAbsChxAna/Main/mac/output_test_data.root") # Run1 After Neutrino

maker.SetEntries(-1)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(5.0)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.SetSignalTypeAbs(True)
maker.SetSignalTypeChx(False)

maker.SetFVCut(True);
maker.FillBootstrapGeant(True);
maker.SetMomThreshCut(0.15);
maker.SetTrkScoreCut(0.3);
maker.SetFVCut(True);



maker.PrintConfig()

maker.MakeFile()




