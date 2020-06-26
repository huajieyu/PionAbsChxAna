#ifndef __DATATYPES_UBXSECEVENTHISTO1D_CXX__
#define __DATATYPES_UBXSECEVENTHISTO1D_CXX__

#include "UBXSecEventHisto1D.h"


namespace DataTypes {


    void UBXSecEventHisto1D::InitializeBootstraps()
    {
    	    h_orimom_proton = new TH1D("h_orimom_proton", "h_orimom_proton", 100, 0.0001, 0.8);
	    h_orimom_neutron = new TH1D("h_orimom_neutron", "h_orimom_neutron", 100, 0.0001, 0.8);
	    h_orimom_pionpm = new TH1D("h_orimom_pionpm", "h_orimom_pionpm", 100, 0.0001, 0.8);
	    h_orimom_pion0 = new TH1D("h_orimom_pion0", "h_orimom_pion0", 100, 0.0001, 0.8);
	    h_orimom_kaon = new TH1D("h_orimom_kaon", "h_orimom_kaon", 100, 0.0001, 0.8);
	    h_orimom_electron = new TH1D("h_orimom_electron", "h_orimom_electron", 100, 0.0001, 0.8);
	    h_orimom_muon = new TH1D("h_orimom_muon", "h_orimom_muon", 100, 0.0001, 0.8);
	    h_orimom_photon = new TH1D("h_orimom_photon", "h_orimom_photon", 100, 0.0001, 0.8);
	    h_orimom_other = new TH1D("h_orimom_other", "h_orimom_other", 100, 0.0001, 0.8);



	    h_chi2_phypo_proton= new TH1D("h_chi2_phypo_proton", "h_chi2_phypo_proton", 100, 0, 500);
	    h_chi2_phypo_pionpm= new TH1D("h_chi2_phypo_pionpm", "h_chi2_phypo_pionpm", 100, 0, 500);
	    h_chi2_phypo_muon= new TH1D("h_chi2_phypo_muon", "h_chi2_phypo_muon", 100, 0, 500);
	    h_chi2_phypo_electron= new TH1D("h_chi2_phypo_electron", "h_chi2_phypo_electron", 100, 0, 500);
	    h_chi2_phypo_kaon= new TH1D("h_chi2_phypo_kaon", "h_chi2_phypo_kaon", 100, 0, 500);
	    h_chi2_phypo_other= new TH1D("h_chi2_phypo_other", "h_chi2_phypo_other", 100, 0, 500);	    
	    h_chi2_phypo_mc= new TH1D("h_chi2_phypo_mc", "h_chi2_phypo_mc", 100, 0, 500);	    


	    htmdqdx_chi2_phypo_proton= new TH1D("htmdqdx_chi2_phypo_proton", "htmdqdx_chi2_phypo_proton", 100, 0, 500);
	    htmdqdx_chi2_phypo_pionpm= new TH1D("htmdqdx_chi2_phypo_pionpm", "htmdqdx_chi2_phypo_pionpm", 100, 0, 500);
	    htmdqdx_chi2_phypo_muon= new TH1D("htmdqdx_chi2_phypo_muon", "htmdqdx_chi2_phypo_muon", 100, 0, 500);
	    htmdqdx_chi2_phypo_electron= new TH1D("htmdqdx_chi2_phypo_electron", "htmdqdx_chi2_phypo_electron", 100, 0, 500);
	    htmdqdx_chi2_phypo_kaon= new TH1D("htmdqdx_chi2_phypo_kaon", "htmdqdx_chi2_phypo_kaon", 100, 0, 500);
	    htmdqdx_chi2_phypo_other= new TH1D("htmdqdx_chi2_phypo_other", "htmdqdx_chi2_phypo_other", 100, 0, 500);	    
	    htmdqdx_chi2_phypo_mc= new TH1D("htmdqdx_chi2_phypo_mc", "htmdqdx_chi2_phypo_mc", 100, 0, 500);	    




	    h_nhits_proton= new TH1D("h_nhits_proton", "h_nhits_proton", 60, 0, 60);
	    h_reconhits_proton= new TH1D("h_reconhits_proton", "h_reconhits_proton", 60, 0, 60);
	    h_tmdqdxnhits_proton= new TH1D("h_tmdqdxnhits_proton", "h_tmdqdxnhits_proton", 60, 0, 60);

	    h_nhits_pionpm= new TH1D("h_nhits_pionpm", "h_nhits_pionpm", 60, 0, 60);
	    h_nhits_muon= new TH1D("h_nhits_muon", "h_nhits_muon", 60, 0, 60);
	    h_nhits_electron= new TH1D("h_nhits_electron", "h_nhits_electron", 60, 0, 60);
	    h_nhits_kaon= new TH1D("h_nhits_kaon", "h_nhits_kaon", 60, 0, 60);
	    h_nhits_other= new TH1D("h_nhits_other", "h_nhits_other", 60, 0, 60);	    
	

	    h_recomom_selected_proton= new TH1D("h_recomom_selected_proton", "h_recomom_selected_proton", 100, 0, 100);
	    h_recomom_selected_pionpm= new TH1D("h_recomom_selected_pionpm", "h_recomom_selected_pionpm", 100, 0, 100);
	    h_recomom_selected_muon= new TH1D("h_recomom_selected_muon", "h_recomom_selected_muon", 100, 0, 100);
	    h_recomom_selected_electron= new TH1D("h_recomom_selected_electron", "h_recomom_selected_electron", 100, 0, 100);
	    h_recomom_selected_kaon= new TH1D("h_recomom_selected_kaon", "h_recomom_selected_kaon", 100, 0, 100);
	    h_recomom_selected_other= new TH1D("h_recomom_selected_other", "h_recomom_selected_other", 100, 0, 100);	    

	    h_recopimom_selected_proton= new TH1D("h_recopimom_selected_proton", "h_recopimom_selected_proton", 100, 0, 100);
	    h_recopimom_selected_pionpm= new TH1D("h_recopimom_selected_pionpm", "h_recopimom_selected_pionpm", 100, 0, 100);
	    h_recopimom_selected_muon= new TH1D("h_recopimom_selected_muon", "h_recopimom_selected_muon", 100, 0, 100);
	    h_recopimom_selected_electron= new TH1D("h_recopimom_selected_electron", "h_recopimom_selected_electron", 100, 0, 100);
	    h_recopimom_selected_kaon= new TH1D("h_recopimom_selected_kaon", "h_recopimom_selected_kaon", 100, 0, 100);
	    h_recopimom_selected_other= new TH1D("h_recopimom_selected_other", "h_recopimom_selected_other", 100, 0, 100);	    


	    h_mom_gentruep = new TH1D("h_mom_gentruep", "h_mom_gentruep", 50, 0.0,2.0);  		
	    h_thetax_gentruep = new TH1D("h_thetax_gentruep", "h_thetax_gentruep", 50, -3.14,3.14);  		


	    h_mom_recotruep = new TH1D("h_mom_recotruep", "h_mom_recotruep", 50, 0.0, 2.0);  		
	    h_mom_recotruep_test = new TH1D("h_mom_recotruep_test", "h_mom_recotruep_test", 50, 0.0, 2.0);  		
	    h_thetax_recotruep_test = new TH1D("h_thetax_recotruep_test", "h_thetax_recotruep_test", 50, -3.14,3.14);  		
            h_mom_nonrecop = new TH1D("h_mom_nonrecop", "h_mom_nonrecop", 50, 0.0, 2.0);       								

            h_mom_selectedtruep = new TH1D("h_mom_selectedtruep", "h_mom_selectedtruep", 50, 0.0, 2.0);
            h_mom_trkscoretruep = new TH1D("h_mom_trkscoretruep", "h_mom_trkscoretruep", 50, 0.0, 2.0);
            h_mom_tmdqdxtruep = new TH1D("h_mom_tmdqdxtruep", "h_mom_tmdqdxtruep", 50, 0.0, 2.0);
            h_mom_chi2truep = new TH1D("h_mom_chi2truep", "h_mom_chi2truep", 50, 0.0, 2.0);

            h_mom_selectedtruepipm = new TH1D("h_mom_selectedtruepipm", "h_mom_selectedtruepipm", 50, 0.0, 1.2);
            h_mom_trkscoretruepipm = new TH1D("h_mom_trkscoretruepipm", "h_mom_trkscoretruepipm", 50, 0.0, 1.2);
            h_mom_tmdqdxtruepipm = new TH1D("h_mom_tmdqdxtruepipm", "h_mom_tmdqdxtruepipm", 50, 0.0, 1.2);
            h_mom_chi2truepipm = new TH1D("h_mom_chi2truepipm", "h_mom_chi2truepipm", 50, 0.0, 1.2);





	    h_mom_gentruepionpm = new TH1D("h_mom_gentruepionpm", "h_mom_gentruepionpm", 50, 0.0,1.2);  		
	    h_mom_recotruepionpm = new TH1D("h_mom_recotruepionpm", "h_mom_recotruepionpm", 50, 0.0, 1.2);  		
	    h_mom_gentrueother = new TH1D("h_mom_gentrueother", "h_mom_gentrueother", 50, 0.0,1.2);  		
	    h_mom_recotrueother = new TH1D("h_mom_recotrueother", "h_mom_recotrueother", 50, 0.0, 1.2);  		
		 

	    h_tmdqdx=new TH1D("h_tmdqdx", "h_tmdqdx", 100, 0.0, 10.);
	    h_tmdqdx_proton=new TH1D("h_tmdqdx_proton", "h_tmdqdx_proton", 100, 0.0, 10.);
	    h_tmdqdx_pionpm=new TH1D("h_tmdqdx_pionpm", "h_tmdqdx_pionpm", 100, 0.0, 10.);
	    

            h_tmdqdxvsrange_proton=new TH2D("h_tmdqdxvsrange_proton", "h_tmdqdxvsrange_proton", 100, 0.0, 10., 100, 0.0, 300.0);
            h_tmdqdxvsrange_pionpm=new TH2D("h_tmdqdxvsrange_pionpm", "h_tmdqdxvsrange_pionpm", 100, 0.0, 10., 100, 0.0, 300.0);

            /*h_tmdqdx=new TH1D("h_tmdqdx", "h_tmdqdx", 100, 0.0, 1200);
	    h_tmdqdx_proton=new TH1D("h_tmdqdx_proton", "h_tmdqdx_proton", 100, 0.0, 1200);
	    h_tmdqdx_pionpm=new TH1D("h_tmdqdx_pionpm", "h_tmdqdx_pionpm", 100, 0.0, 1200);
	    

            h_tmdqdxvsrange_proton=new TH2D("h_tmdqdxvsrange_proton", "h_tmdqdxvsrange_proton", 100, 0.0, 1200, 100, 0.0, 300.0);
            h_tmdqdxvsrange_pionpm=new TH2D("h_tmdqdxvsrange_pionpm", "h_tmdqdxvsrange_pionpm", 100, 0.0, 1200, 100, 0.0, 300.0);
            */

	    h_nhits_shwlike_photon=new TH1D("h_nhits_shwlike_photon", "h_nhits_shwlike_photon", 600, -0.5, 599.5);
	    h_nhits_shwlike_electron=new TH1D("h_nhits_shwlike_electron", "h_nhits_shwlike_electron", 600, -0.5, 599.5);
	    h_nhits_shwlike_nonphoton=new TH1D("h_nhits_shwlike_nonphoton", "h_nhits_shwlike_nonphoton", 600, -0.5, 599.5);
	    h_nhits_shwlike_all=new TH1D("h_nhits_shwlike_all", "h_nhits_shwlike_all", 600, -0.5, 599.5);
 
	    h_shweng_photon=new TH1D("h_shweng_photon", "h_shweng_photon", 60, -0., 600);
	    h_shweng_electron=new TH1D("h_shweng_electron", "h_shweng_electron", 60, -0., 600);
	    h_shweng_nonphoton=new TH1D("h_shweng_nonphoton", "h_shweng_nonphoton", 60, -0., 600);
	    h_shweng_all=new TH1D("h_shweng_all", "h_shweng_all", 60, -0.0, 600);
 
	    h_shwang_photon=new TH1D("h_shwang_photon", "h_shwang_photon", 60, -0., 3.14);
	    h_shwang_electron=new TH1D("h_shwang_electron", "h_shwang_electron", 60, -0., 3.14);
	    h_shwang_nonphoton=new TH1D("h_shwang_nonphoton", "h_shwang_nonphoton", 60, -0., 3.14);
	    h_shwang_all=new TH1D("h_shwang_all", "h_shwang_all", 60, -0.0, 3.14);
 
	    h_shwdis_photon=new TH1D("h_shwdis_photon", "h_shwdis_photon", 60, -0., 150);
	    h_shwdis_electron=new TH1D("h_shwdis_electron", "h_shwdis_electron", 60, -0., 150);
	    h_shwdis_nonphoton=new TH1D("h_shwdis_nonphoton", "h_shwdis_nonphoton", 60, -0., 150);
	    h_shwdis_all=new TH1D("h_shwdis_all", "h_shwdis_all", 60, -0.0, 150);
 

	    h_trackscore_shwlike_photon=new TH1D("h_trackscore_shwlike_photon", "h_trackscore_shwlike_photon", 100, -0.0, 0.3);
	    h_trackscore_shwlike_electron=new TH1D("h_trackscore_shwlike_electron", "h_trackscore_shwlike_electron", 100, -0.0, 0.3);
	    h_trackscore_shwlike_nonphoton=new TH1D("h_trackscore_shwlike_nonphoton", "h_trackscore_shwlike_nonphoton", 100, -0.0, 0.3);
	    h_trackscore_shwlike_all=new TH1D("h_trackscore_shwlike_all", "h_trackscore_shwlike_all", 100, -0.0, 0.3);
 


	    h_trackscore_all=new TH1D("h_trackscore_all", "h_trackscore_all", 100, -0.0, 1.0);
	    h_trackscore_proton=new TH1D("h_trackscore_proton", "h_trackscore_proton", 100, -0.0, 1.0);
	    h_trackscore_pionpm=new TH1D("h_trackscore_pionpm", "h_trackscore_pionpm", 100, -0.0, 1.0);
	    h_trackscore_electron=new TH1D("h_trackscore_electron", "h_trackscore_electron", 100, -0.0, 1.0);
	    h_trackscore_pion0=new TH1D("h_trackscore_pion0", "h_trackscore_pion0", 100, -0.0, 1.0);
	    h_trackscore_muon=new TH1D("h_trackscore_muon", "h_trackscore_muon", 100, -0.0, 1.0);
	    h_trackscore_other=new TH1D("h_trackscore_other", "h_trackscore_other", 100, -0.0, 1.0);
	    h_trackscore_photon=new TH1D("h_trackscore_photon", "h_trackscore_photon", 100, -0.0, 1.0);

            h_nonrecop_momvsnhits = new TH2D("h_nonrecop_momvsnhits", "h_nonrecop_momvsnhits", 100, 0.0, 1.5, 100, 0.0, 100);
            h_nonrecop_momvsnhits -> GetXaxis()->SetTitle("True Momentum (proton)[GeV]");
            h_nonrecop_momvsnhits -> GetYaxis()->SetTitle("True number of hits(proton)");
            h_nonrecop_momvsnhits -> SetTitle("");
            //========================================================================================
            h_lownhitsp_thetax = new TH1D("h_lownhitsp_thetax", "h_lownhitsp_thetax", 100, -1.0, 1.);       								
            h_lownhitsp_thetay = new TH1D("h_lownhitsp_thetay", "h_lownhitsp_thetay", 100, -1.0, 1.);       								
            h_lownhitsp_thetaz = new TH1D("h_lownhitsp_thetaz", "h_lownhitsp_thetaz", 100, -1.0, 1.);       								

            h_lownhitsp_startX = new TH1D("h_lownhitsp_startX", "h_lownhitsp_startX", 100, -400.0, 400.);
            h_lownhitsp_startY = new TH1D("h_lownhitsp_startY", "h_lownhitsp_startY", 100, -.0, 700.); 
            h_lownhitsp_startZ = new TH1D("h_lownhitsp_startZ", "h_lownhitsp_startZ", 100, -.0, 800.);       							
            h_lownhitsp_endX = new TH1D("h_lownhitsp_endX", "h_lownhitsp_endX", 100, -400.0, 400.); 
            h_lownhitsp_endY = new TH1D("h_lownhitsp_endY", "h_lownhitsp_endY", 100, -.0, 700.);   
            h_lownhitsp_endZ = new TH1D("h_lownhitsp_endZ", "h_lownhitsp_endZ", 100, -.0, 800.);       								

            h_chi2cutp_thetax = new TH1D("h_chi2cutp_thetax", "h_chi2cutp_thetax", 100, -1.0, 1.);  
            h_chi2cutp_thetay = new TH1D("h_chi2cutp_thetay", "h_chi2cutp_thetay", 100, -1.0, 1.);     
            h_chi2cutp_thetaz = new TH1D("h_chi2cutp_thetaz", "h_chi2cutp_thetaz", 100, -1.0, 1.);       								
            

            h_mom_gentruepion0 = new TH1D("h_mom_gentruepion0", "h_mom_gentruepion0", 50, 0.0, 1.2);
            h_mom_recotruepion0 = new TH1D("h_mom_recotruepion0", "h_mom_recotruepion0", 50, 0.0, 1.2);
            h_mom_recotruenonpi0 = new TH1D("h_mom_recotruenonpi0", "h_mom_recotruenonpi0", 50, 0.0, 1.2);



	    h_mom_gentruepionpm = new TH1D("h_mom_gentruepionpm", "h_mom_gentruepionpm", 50, 0.0, 1.2);
            h_mom_recotruepionpm = new TH1D("h_mom_recotruepionpm", "h_mom_recotruepionpm", 50, 0.0, 1.2);

            h_PiAbs_sel_pmom = new TH1D ("h_PiAbs_sel_pmom", "h_PiAbs_sel_pmom", 50, 0.0, 1.2);
            h_PiAbs_sig_pmom = new TH1D ("h_PiAbs_sig_pmom", "h_PiAbs_sig_pmom", 50, 0.0, 1.2);
            h_PiAbs_chxbac_pmom = new TH1D ("h_PiAbs_chxbac_pmom", "h_PiAbs_chxbac_pmom", 50, 0.0, 1.2);
            h_PiAbs_reabac_pmom = new TH1D ("h_PiAbs_reabac_pmom", "h_PiAbs_reabac_pmom", 50, 0.0, 1.2);
            h_PiAbs_other_pmom = new TH1D ("h_PiAbs_other_pmom", "h_PiAbs_other_pmom", 50, 0.0, 1.2);

            h_PiAbs_sel_pcostheta = new TH1D ("h_PiAbs_sel_pcostheta", "h_PiAbs_sel_pcostheta", 50, -1., 1.0);
            h_PiAbs_sig_pcostheta = new TH1D ("h_PiAbs_sig_pcostheta", "h_PiAbs_sig_pcostheta", 50, -1., 1.0);
            h_PiAbs_chxbac_pcostheta = new TH1D ("h_PiAbs_chxbac_pcostheta", "h_PiAbs_chxbac_pcostheta", 50, -1., 1.0);
            h_PiAbs_reabac_pcostheta = new TH1D ("h_PiAbs_reabac_pcostheta", "h_PiAbs_reabac_pcostheta", 50, -1., 1.0);
            h_PiAbs_other_pcostheta = new TH1D ("h_PiAbs_other_pcostheta", "h_PiAbs_other_pcostheta", 50, -1., 1.0);

            h_PiAbs_sel_ptheta = new TH1D ("h_PiAbs_sel_ptheta", "h_PiAbs_sel_ptheta", 50, -0., 3.14);
            h_PiAbs_sig_ptheta = new TH1D ("h_PiAbs_sig_ptheta", "h_PiAbs_sig_ptheta", 50, -0., 3.14);
            h_PiAbs_chxbac_ptheta = new TH1D ("h_PiAbs_chxbac_ptheta", "h_PiAbs_chxbac_ptheta", 50, -0., 3.14);
            h_PiAbs_reabac_ptheta = new TH1D ("h_PiAbs_reabac_ptheta", "h_PiAbs_reabac_ptheta", 50, -0., 3.14);
            h_PiAbs_other_ptheta = new TH1D ("h_PiAbs_other_ptheta", "h_PiAbs_other_ptheta", 50, -0., 3.14);


            h_PiAbs_sel_pphi = new TH1D ("h_PiAbs_sel_pphi", "h_PiAbs_sel_pphi", 50, -3.14, 3.14);
            h_PiAbs_sig_pphi = new TH1D ("h_PiAbs_sig_pphi", "h_PiAbs_sig_pphi", 50, -3.14, 3.14);
            h_PiAbs_chxbac_pphi = new TH1D ("h_PiAbs_chxbac_pphi", "h_PiAbs_chxbac_pphi", 50, -3.14, 3.14);
            h_PiAbs_reabac_pphi = new TH1D ("h_PiAbs_reabac_pphi", "h_PiAbs_reabac_pphi", 50, -3.14, 3.14);
            h_PiAbs_other_pphi = new TH1D ("h_PiAbs_other_pphi", "h_PiAbs_other_pphi", 50, -3.14, 3.14);

            h_PiAbs_sel_pcosthetax= new TH1D("h_PiAbs_sel_pcosthetax", "h_PiAbs_sel_pcosthetax", 20, -1., 1.);

  

            h_PiAbs_sel_pmult = new TH1D ("h_PiAbs_sel_pmult", "h_PiAbs_sel_pmult", 5, -0.5, 4.5);
            h_PiAbs_sig_pmult = new TH1D ("h_PiAbs_sig_pmult", "h_PiAbs_sig_pmult", 5, -0.5, 4.5);
            h_PiAbs_chxbac_pmult = new TH1D ("h_PiAbs_chxbac_pmult", "h_PiAbs_chxbac_pmult", 5, -0.5, 4.5);
            h_PiAbs_reabac_pmult = new TH1D ("h_PiAbs_reabac_pmult", "h_PiAbs_reabac_pmult", 5, -0.5, 4.5);

            h_PiAbs_sig_energeticproton_reco_mom = new TH1D("h_PiAbs_sig_energeticproton_reco_mom", "h_PiAbs_sig_energeticproton_reco_mom", 15, 0.0, 1.2);
            h_PiAbs_sig_energeticproton_reco_costheta = new TH1D("h_PiAbs_sig_energeticproton_reco_costheta", "h_PiAbs_sig_energeticproton_reco_costheta", 15, -1.0, 1.0);
            h_PiAbs_sig_energeticproton_reco_phi = new TH1D("h_PiAbs_sig_energeticproton_reco_phi", "h_PiAbs_sig_energeticproton_reco_phi", 15, -3.14, 3.14);

            h_PiAbs_sig_energeticproton_truevsreco_mom = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_mom", "h_PiAbs_sig_energeticproton_truevsreco_mom", 30, 0.0, 1.2, 30, 0.0, 1.2);
            h_PiAbs_sig_energeticproton_truevsreco_costheta = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_costheta", "h_PiAbs_sig_energeticproton_truevsreco_costheta", 30, -1.0, 1.0, 30, -1.0, 1.0);
            h_PiAbs_sig_energeticproton_truevsreco_phi = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_phi", "h_PiAbs_sig_energeticproton_truevsreco_phi", 30, -3.14, 3.14, 30, -3.14, 3.14);


            h_PiAbs_gen_sig_energeticproton_mom = new TH1D("h_PiAbs_gen_sig_energeticproton_mom", "h_PiAbs_gen_sig_energeticproton_mom", 30, 0.0, 1.2);
            h_PiAbs_gen_sig_energeticproton_costheta = new TH1D("h_PiAbs_gen_sig_energeticproton_costheta", "h_PiAbs_gen_sig_energeticproton_costheta", 30, -1.0, 1.0);
            h_PiAbs_gen_sig_energeticproton_phi = new TH1D("h_PiAbs_gen_sig_energeticproton_phi", "h_PiAbs_gen_sig_energeticproton_phi", 30, -3.14, 3.14);







            h_PiAbs_sig_nmult = new TH1D ("h_PiAbs_sig_nmult", "h_PiAbs_sig_nmult", 15, -0.5, 14.5);
            h_PiAbs_sig_nmult ->GetXaxis()->SetTitle("Neutron Multiplicity[True]");
            h_PiAbs_sig_nmult ->GetYaxis()->SetTitle("No. of Events");

            h_PiAbs_sig_truevsreco_pmult = new TH2D ("h_PiAbs_sig_truevsreco_pmult", "h_PiAbs_sig_truevsreco_pmult", 10, -0.5, 9.5, 10, -0.5, 9.5);
            h_PiAbs_sig_truevsreco_pmult ->GetXaxis()->SetTitle("Proton Multiplicity[True]");
            h_PiAbs_sig_truevsreco_pmult ->GetYaxis()->SetTitle("Proton_Multiplicity[Reco]");

            h_sig_pvsnmult = new TH2D("h_sig_pvsnmult", "h_sig_pvsnmult", 15, -0.5, 14.5, 15, -0.5, 14.5);
            h_sig_pvsnmult->GetXaxis()->SetTitle("No. of Protons");
            h_sig_pvsnmult->GetYaxis()->SetTitle("No. of Neutrons");

            h_PiAbs_sig_true_totalKE_protonvsneutron = new TH2D("h_PiAbs_sig_true_totalKE_protonvsneutron", "h_PiAbs_sig_true_totalKE_protonvsneutron", 50, 0., 1., 50, 0., 1.);
            h_PiAbs_sig_true_totalKE_protonvsneutron->GetXaxis()->SetTitle("Total KE:Protons[GeV]");
            h_PiAbs_sig_true_totalKE_protonvsneutron->GetYaxis()->SetTitle("Total KE:Neutrons[GeV]");


            h_PiAbs_sig_proton_mom = new TH1D("h_PiAbs_sig_proton_mom", "h_PiAbs_sig_proton_mom", 30, 0.0, 1.0);

            h_PiAbs_sig_proton_mom->GetXaxis()->SetTitle("True Momentum of Proton [GeV]");
            h_PiAbs_sig_proton_mom->GetYaxis()->SetTitle("No. of Tracks");

          

            h_PiAbs_sig_neutron_mom = new TH1D("h_PiAbs_sig_neutron_mom", "h_PiAbs_sig_neutron_mom", 30, 0.0, 1.0);

            h_PiAbs_sig_neutron_mom->GetXaxis()->SetTitle("True Momentum of Neutron [GeV]");
            h_PiAbs_sig_neutron_mom->GetYaxis()->SetTitle("No. of Tracks");

            h_PiAbs_sig_proton_costheta = new TH1D("h_PiAbs_sig_proton_costheta", "h_PiAbs_sig_proton_costheta", 30, -1., 1.0);

            h_PiAbs_sig_proton_costheta->GetXaxis()->SetTitle("True cos#theta of Proton [GeV]");
            h_PiAbs_sig_proton_costheta->GetYaxis()->SetTitle("No. of Tracks");

          

            h_PiAbs_sig_neutron_costheta = new TH1D("h_PiAbs_sig_neutron_costheta", "h_PiAbs_sig_neutron_costheta", 30, -1., 1.0);

            h_PiAbs_sig_neutron_costheta->GetXaxis()->SetTitle("True cos#theta of Neutron [GeV]");
            h_PiAbs_sig_neutron_costheta->GetYaxis()->SetTitle("No. of Tracks");




            h_PiAbs_sig_true_Ptmissing = new TH1D("h_PiAbs_sig_true_Ptmissing","h_PiAbs_sig_true_Ptmissing", 30, 0.0, 0.5);
            h_PiAbs_sig_true_Ptmissing->GetXaxis()->SetTitle("True missing-P_{T}[GeV]");
            h_PiAbs_sig_true_Ptmissing->GetYaxis()->SetTitle("No. of Tracks");

            h_PiAbs_sig_true_Ptmissing_onlyproton = new TH1D("h_PiAbs_sig_true_Ptmissing_onlyproton","h_PiAbs_sig_true_Ptmissing_onlyproton", 30, 0.0, 1.0);
            h_PiAbs_sig_true_Ptmissing_onlyproton->GetXaxis()->SetTitle("True missing_onlyproton-P_{T}[GeV]");
            h_PiAbs_sig_true_Ptmissing_onlyproton->GetYaxis()->SetTitle("No. of Tracks");


            h_sel_Ptmissing = new TH1D ("h_sel_Ptmissing", "h_sel_Ptmissing", 30, 0.0, 1.2);
            h_sel_Pmissing = new TH1D ("h_sel_Pmissing", "h_sel_Pmissing", 30, 0.0, 1.6);
            h_sel_Emissing = new TH1D ("h_sel_Emissing", "h_sel_Emissing", 30, 0.0, 1.2);
            h_sel_Plongit = new TH1D ("h_sel_Plongit", "h_sel_Plongit", 30, 0.0, 1.2);
 

            h_sig_Ptmissing = new TH1D ("h_sig_Ptmissing", "h_sig_Ptmissing", 30, 0.0, 1.2);
            h_sig_Pmissing = new TH1D ("h_sig_Pmissing", "h_sig_Pmissing", 30, 0.0, 1.6);
            h_sig_Emissing = new TH1D ("h_sig_Emissing", "h_sig_Emissing", 30, 0.0, 1.2);
            h_sig_Plongit = new TH1D ("h_sig_Plongit", "h_sig_Plongit", 30, 0.0, 1.2);

            h_sig_Ptmissing_vs_pcand = new TH2D("h_sig_Ptmissing_vs_pcand", "h_sig_Ptmissing_vs_pcand", 30, 0.0, 1.2, 5, -0.5, 4.5);
            
            h_chxbac_Ptmissing = new TH1D ("h_chxbac_Ptmissing", "h_chxbac_Ptmissing", 30, 0.0, 1.2);
            h_chxbac_Pmissing = new TH1D ("h_chxbac_Pmissing", "h_chxbac_Pmissing", 30, 0.0, 1.6);
            h_chxbac_Emissing = new TH1D ("h_chxbac_Emissing", "h_chxbac_Emissing", 30, 0.0, 1.2);
            h_chxbac_Plongit = new TH1D ("h_chxbac_Plongit", "h_chxbac_Plongit", 30, 0.0, 1.2);
            
            h_reabac_Ptmissing = new TH1D ("h_reabac_Ptmissing", "h_reabac_Ptmissing", 30, 0.0, 1.2);
            h_reabac_Pmissing = new TH1D ("h_reabac_Pmissing", "h_reabac_Pmissing", 30, 0.0, 1.6);
            h_reabac_Emissing = new TH1D ("h_reabac_Emissing", "h_reabac_Emissing", 30, 0.0, 1.2);
            h_reabac_Plongit = new TH1D ("h_reabac_Plongit", "h_reabac_Plongit", 30, 0.0, 1.2);
            
            h_other_Ptmissing = new TH1D ("h_other_Ptmissing", "h_other_Ptmissing", 30, 0.0, 1.2);
            h_other_Pmissing = new TH1D ("h_other_Pmissing", "h_other_Pmissing", 30, 0.0, 1.6);
            h_other_Emissing = new TH1D ("h_other_Emissing", "h_other_Emissing", 30, 0.0, 1.2);
            h_other_Plongit = new TH1D ("h_other_Plongit", "h_other_Plongit", 30, 0.0, 1.2);


            h_reco_pionpm_rea = new TH1D("h_reco_pionpm_rea", "h_reco_pionpm_rea ", 10, -0.5, 9.5);
            h_reco_pionpm_rea -> GetXaxis() ->SetTitle("No. of #pi^{+}/#pi^{-}");
            h_reco_pionpm_rea -> GetYaxis() ->SetTitle("No. of Events");

            h_reco_photon_rea = new TH1D("h_reco_photon_rea", "h_reco_photon_rea ", 10, -0.5, 9.5);
            h_reco_photon_rea -> GetXaxis() ->SetTitle("No. of #gamma");
            h_reco_photon_rea -> GetYaxis() ->SetTitle("No. of Events");
 
            h_reco_photon_chx = new TH1D("h_reco_photon_chx", "h_reco_photon_chx ", 10, -0.5, 9.5);
            h_reco_photon_chx -> GetXaxis() ->SetTitle("No. of #gamma");
            h_reco_photon_chx -> GetYaxis() ->SetTitle("No. of Events");

            h_daughter_beam_dist = new TH1D("h_daughter_beam_dist", "h_daughter_beam_dist", 30, 0.0, 30.0); 

            h_reco_grand_daughter_beam_dist = new TH1D("h_reco_grand_daughter_beam_dist", "h_reco_grand_daughter_beam_dist", 30, 0.0, 30.0); 
            h_reco_daughter_beam_dist = new TH1D("h_reco_daughter_beam_dist", "h_reco_daughter_beam_dist", 30, 0.0, 30.0); 

            h_reco_grand_daughter_angle = new TH1D("h_reco_grand_daughter_angle", "h_reco_grand_daughter_angle", 30, 0.0, 3.14); 
            h_reco_daughter_angle = new TH1D("h_reco_daughter_angle", "h_reco_daughter_angle", 30, 0.0, 3.14); 

            h_reco_daughter_deltax = new TH1D("h_reco_daughter_deltax", "h_reco_daughter_deltax", 30, -20., 40); 
            h_reco_daughter_deltay = new TH1D("h_reco_daughter_deltay", "h_reco_daughter_deltay", 30, -20., 40); 
            h_reco_daughter_deltaz = new TH1D("h_reco_daughter_deltaz", "h_reco_daughter_deltaz", 30, -20., 40); 






            h_reco_daughter_distvsangle=new TH2D("h_reco_daughter_distvsangle", "h_reco_daughter_distvsangle", 30, 0.0, 30.0, 30, 0.0, 3.14);


            h_gdfromproton= new TH1D("h_gdfromproton", "h_gdfromproton", 30, 0.0, 0.8);
            h_gdfromother= new TH1D("h_gdfromother", "h_gdfromother", 30, 0.0, 0.8);

            h_selproton_momreso = new TH2D("h_selproton_momreso", "h_selproton_momreso", 50, 0.0, 1.0, 50, 0.0, 1.0);
            h_selproton_momreso_new = new TH1D("h_selproton_momreso_new", "h_selproton_momreso_new", 60, -0.3, 0.3);
            h_selproton_momreso_new->GetXaxis()->SetTitle("(reco P - true P)/true P");
            h_selproton_momreso_new->GetYaxis()->SetTitle("No. of Tracks");

            h_selproton_momPxreso_new = new TH1D("h_selproton_momPxreso_new", "h_selproton_momPxreso_new", 60, -0.3, 0.3);
            h_selproton_momPxreso_new->GetXaxis()->SetTitle("(recoPx - truePx)/truePx");
            h_selproton_momPxreso_new->GetYaxis()->SetTitle("No. of Tracks");

            h_selproton_momPyreso_new = new TH1D("h_selproton_momPyreso_new", "h_selproton_momPyreso_new", 60, -0.3, 0.3);
            h_selproton_momPyreso_new->GetXaxis()->SetTitle("(recoPy - truePy)/truePy");
            h_selproton_momPyreso_new->GetYaxis()->SetTitle("No. of Tracks");

            h_selproton_momPzreso_new = new TH1D("h_selproton_momPzreso_new", "h_selproton_momPzreso_new", 60, -0.3, 0.3);
            h_selproton_momPzreso_new->GetXaxis()->SetTitle("(recoPz - truePz)/truePz");
            h_selproton_momPzreso_new->GetYaxis()->SetTitle("No. of Tracks");

            h_selproton_costhetareso_new = new TH1D("h_selproton_costhetareso_new", "h_selproton_costhetareso_new", 60, -0.3, 0.3);
            h_selproton_costhetareso_new->GetXaxis()->SetTitle("(reco cos#theta - true cos#theta)/true cos#theta");
            h_selproton_costhetareso_new->GetYaxis()->SetTitle("No. of Tracks");





            h_selproton_costhetareso = new TH2D("h_selproton_costhetareso", "h_selproton_costhetareso", 50, -1.0, 1.0, 50, -1.0, 1.0);
            h_ngamma_frompi0=new TH1D("h_ngamma_frompi0", "h_ngamma_frompi0", 2.0, 0.5, 2.5);
            h_pgamma_frompi0=new TH1D("h_pgamma_frompi0", "h_pgamma_frompi0", 20, 0.0, 1.5);

            //selected_signal_percut = new TH1D("selected_signal_percut", "selected_percut", 4, 0, 4);
            //selected_percut = new TH1D("selected_percut", "selected_percut", 4, 0, 4);

            h_true_sig_trkmult = new TH1D("h_true_sig_trkmult", "h_true_sig_trkmult", 6, -0.5, 5.5);
            h_true_chxbac_trkmult = new TH1D("h_true_chxbac_trkmult", "h_true_chxbac_trkmult", 6, -0.5, 5.5);
            h_true_reabac_trkmult = new TH1D("h_true_reabac_trkmult", "h_true_reabac_trkmult", 6, -0.5, 5.5);
            h_selproton_dedxRR = new TH2D ("h_selproton_dedxRR", "h_selproton_dedxRR", 50, 0, 20, 50, 0, 30);


            h_shwcut_mult=new TH1D("h_shwcut_mult", "h_shwcut_mult", 10, -0.5, 9.5);
            h_tmdqdxcut_mult=new TH1D("h_tmdqdxcut_mult", "h_tmdqdxcut_mult", 10, -0.5, 9.5);
            h_chi2cut_mult=new TH1D("h_chi2cut_mult", "h_chi2cut_mult", 10, -0.5, 9.5);

            hsig_shwcut_mult=new TH1D("hsig_shwcut_mult", "hsig_shwcut_mult", 10, -0.5, 9.5);
            hsig_tmdqdxcut_mult=new TH1D("hsig_tmdqdxcut_mult", "hsig_tmdqdxcut_mult", 10, -0.5, 9.5);
            hsig_chi2cut_mult=new TH1D("hsig_chi2cut_mult", "hsig_chi2cut_mult", 10, -0.5, 9.5);

            hchxbac_shwcut_mult=new TH1D("hchxbac_shwcut_mult", "hchxbac_shwcut_mult", 10, -0.5, 9.5);
            hchxbac_tmdqdxcut_mult=new TH1D("hchxbac_tmdqdxcut_mult", "hchxbac_tmdqdxcut_mult", 10, -0.5, 9.5);
            hchxbac_chi2cut_mult=new TH1D("hchxbac_chi2cut_mult", "hchxbac_chi2cut_mult", 10, -0.5, 9.5);

            hreabac_shwcut_mult=new TH1D("hreabac_shwcut_mult", "hreabac_shwcut_mult", 10, -0.5, 9.5);
            hreabac_tmdqdxcut_mult=new TH1D("hreabac_tmdqdxcut_mult", "hreabac_tmdqdxcut_mult", 10, -0.5, 9.5);
            hreabac_chi2cut_mult=new TH1D("hreabac_chi2cut_mult", "hreabac_chi2cut_mult", 10, -0.5, 9.5);

            hsig_nrecogamma = new TH1D("hsig_nrecogamma", "hsig_nrecogamma", 5, -0.5, 4.5);
            hchxbac_nrecogamma = new TH1D("hchxbac_nrecogamma", "hchxbac_nrecogamma", 5, -0.5, 4.5);
            hreabac_nrecogamma = new TH1D("hreabac_nrecogamma", "hreabac_nrecogamma", 5, -0.5, 4.5);
            hother_nrecogamma = new TH1D("hother_nrecogamma", "hother_nrecogamma", 5, -0.5, 4.5);


            h_totalp_mom=new TH1D("h_totalp_mom", "h_totalp_mom", 60, 0.0, 1.2);
            h_reintp_mom=new TH1D("h_reintp_mom", "h_reintp_mom", 60, 0.0, 1.2);

          
            h_trklen_reso = new TH1D("h_trklen_reso", "h_trklen_reso", 60, -0.3, 0.3); 
            h_trklen_reso->GetXaxis()->SetTitle("(Reco L-True L)/True L");
            h_trklen_reso->GetYaxis()->SetTitle("No. of Tracks");

            h_seltrk_ptheta_cosmic = new TH1D("h_seltrk_ptheta_cosmic", "h_seltrk_ptheta_cosmic", 50, 0.0, 3.14);
            h_seltrk_ptheta_beam_daughter = new TH1D("h_seltrk_ptheta_beam_daughter", "h_seltrk_ptheta_beam_daughter", 50, 0.0, 3.14);
            h_seltrk_ptheta_beam_granddaughter = new TH1D("h_seltrk_ptheta_beam_granddaughter", "h_seltrk_ptheta_beam_granddaughter", 50, 0.0, 3.14);
            h_seltrk_pphi_cosmic = new TH1D("h_seltrk_pphi_cosmic", "h_seltrk_pphi_cosmic", 50, -3.14, 3.14);
            h_seltrk_pphi_beam_daughter = new TH1D("h_seltrk_pphi_beam_daughter", "h_seltrk_pphi_beam_daughter", 50, -3.14, 3.14);
            h_seltrk_pphi_beam_granddaughter = new TH1D("h_seltrk_pphi_beam_granddaughter", "h_seltrk_pphi_beam_granddaughter", 50, -3.14, 3.14);

            h_seltrk_distvtx_daughter = new TH1D("h_seltrk_distvtx_daughter","h_seltrk_distvtx_daughter",50, 0.0, 30);
            h_seltrk_distvtx_granddaughter = new TH1D("h_seltrk_distvtx_granddaughter","h_seltrk_distvtx_granddaughter",50, 0.0, 30);

            h_seltrk_angle_daughter = new TH1D("h_seltrk_angle_daughter","h_seltrk_angle_daughter",50, -0.0, 3.14);
            h_seltrk_angle_granddaughter = new TH1D("h_seltrk_angle_granddaughter","h_seltrk_angle_granddaughter",50, -0.0, 3.14);


            h_true_beam_endE_den = new TH1D("h_true_beam_endE_den", "h_true_beam_endE_den", nbins_beamE, bins_beamE);
            h_true_beam_endE_num = new TH1D("h_true_beam_endE_num", "h_true_beam_endE_num", nbins_beamE, bins_beamE);

            h_sel_gdvsd=new TH2D("h_sel_gdvsd", "h_sel_gdvsd", 6, -0.5, 5.5, 6, -0.5, 5.5); 

    }

  }

#endif
