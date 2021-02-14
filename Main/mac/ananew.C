#define ananew_cxx
#include "ananew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <vector>

//For MC from Owen Goodwins studies
double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

double cutAPA3_Z = 226.;

bool manual_beamPos_mc(double beam_startX, double beam_startY,
                       double beam_startZ, double beam_dirX,
                       double beam_dirY,   double beam_dirZ, 
                       double true_dirX,   double true_dirY,
                       double true_dirZ,   double true_startX,
                       double true_startY, double true_startZ) {
  double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;
  
  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;
  
  if ( (beam_startY - projectY) < ylow )
    return false;
  
  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;
  
  return true;
  
};

bool endAPA3(double reco_beam_endZ){
  return(reco_beam_endZ < cutAPA3_Z);
}

void ananew::LoadHist(){

  TFile *infile = TFile::Open("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");

  //Load in files
  TH3F* hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
  TH3F* hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
  TH3F* hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
  TH3F* hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
  TH3F* hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
  TH3F* hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
  
  TH3F* hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
  TH3F* hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
  TH3F* hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
  TH3F* hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
  TH3F* hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
  TH3F* hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
  
  TH3F* hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
  TH3F* hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
  TH3F* hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
  TH3F* hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
  TH3F* hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
  TH3F* hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");
  
  TH3F* hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
  TH3F* hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
  TH3F* hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
  TH3F* hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
  TH3F* hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
  TH3F* hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");
  
  hDx_sim_pos->SetDirectory(0);
  hDy_sim_pos->SetDirectory(0);
  hDz_sim_pos->SetDirectory(0);
  hEx_sim_pos->SetDirectory(0);
  hEy_sim_pos->SetDirectory(0);
  hEz_sim_pos->SetDirectory(0);
  
  hDx_sim_neg->SetDirectory(0);
  hDy_sim_neg->SetDirectory(0);
  hDz_sim_neg->SetDirectory(0);
  hEx_sim_neg->SetDirectory(0);
  hEy_sim_neg->SetDirectory(0);
  hEz_sim_neg->SetDirectory(0);
  
  for(int y = 1; y <= 31; y++){
    for(int z = 1; z <= 37; z++){
      spline_dx_fwd_neg[y-1][z-1] = MakeSpline(hDx_sim_neg,1,y,z,1,1);
      spline_dx_fwd_pos[y-1][z-1] = MakeSpline(hDx_sim_pos,1,y,z,1,2);
      spline_dEx_neg[y-1][z-1] = MakeSpline(hEx_sim_neg,1,y,z,3,1);
      spline_dEx_pos[y-1][z-1] = MakeSpline(hEx_sim_pos,1,y,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int z = 1; z <= 37; z++){
      spline_dy_fwd_neg[x-1][z-1] = MakeSpline(hDy_sim_neg,2,x,z,1,1);
      spline_dy_fwd_pos[x-1][z-1] = MakeSpline(hDy_sim_pos,2,x,z,1,2);
      spline_dEy_neg[x-1][z-1] = MakeSpline(hEy_sim_neg,2,x,z,3,1);
      spline_dEy_pos[x-1][z-1] = MakeSpline(hEy_sim_pos,2,x,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int y = 1; y <= 31; y++){
      spline_dz_fwd_neg[x-1][y-1] = MakeSpline(hDz_sim_neg,3,x,y,1,1);
      spline_dz_fwd_pos[x-1][y-1] = MakeSpline(hDz_sim_pos,3,x,y,1,2);
      spline_dEz_neg[x-1][y-1] = MakeSpline(hEz_sim_neg,3,x,y,3,1);
      spline_dEz_pos[x-1][y-1] = MakeSpline(hEz_sim_pos,3,x,y,3,2);
    }
  }
}

TSpline3* ananew::MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const
{
  TSpline3 *spline = 0;
  
  if(dim1 == 1)
  {
    double a[19];
    double b[19];
    for(int x = 1; x <= 19; x++)
    {
      a[x-1] = spline_hist->GetXaxis()->GetBinCenter(x);
      b[x-1] = spline_hist->GetBinContent(x,dim2_bin,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 2)
  {
    double a[31];
    double b[31];
    for(int y = 1; y <= 31; y++)
    {
      a[y-1] = spline_hist->GetYaxis()->GetBinCenter(y);
      b[y-1] = spline_hist->GetBinContent(dim2_bin,y,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 3)
  {
    double a[37];
    double b[37];
    for(int z = 1; z <= 37; z++)
    {
      a[z-1] = spline_hist->GetZaxis()->GetBinCenter(z);
      b[z-1] = spline_hist->GetBinContent(dim2_bin,dim3_bin,z);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }

  return spline;
}

double ananew::InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const
{
  int bin_x = interp_hist->GetXaxis()->FindBin(xVal);
  int bin_y = interp_hist->GetYaxis()->FindBin(yVal);
  int bin_z = interp_hist->GetZaxis()->FindBin(zVal);

  int bincenter_x = interp_hist->GetXaxis()->GetBinCenter(bin_x);
  int bincenter_y = interp_hist->GetYaxis()->GetBinCenter(bin_y);
  int bincenter_z = interp_hist->GetZaxis()->GetBinCenter(bin_z);

  int max_x = interp_hist->GetNbinsX();
  int max_y = interp_hist->GetNbinsY();
  int max_z = interp_hist->GetNbinsZ();
  
  int low_x;
  int high_x;
  if(bin_x <= 1)
  {
    low_x = 1;
    high_x = 2;
  }
  else if(bin_x >= max_x)
  {
    low_x = max_x-1;
    high_x = max_x;
  }
  else if(xVal > bincenter_x)
  {
    low_x = bin_x;
    high_x = bin_x+1;
  }
  else
  {
    low_x = bin_x-1;
    high_x = bin_x;
  }

  int low_y;
  int high_y;
  if(bin_y <= 1)
  {
    low_y = 1;
    high_y = 2;
  }
  else if(bin_y >= max_y)
  {
    low_y = max_y-1;
    high_y = max_y;
  }
  else if(yVal > bincenter_y)
  {
    low_y = bin_y;
    high_y = bin_y+1;
  }
  else
  {
    low_y = bin_y-1;
    high_y = bin_y;
  }

  int low_z;
  int high_z;
  if(bin_z <= 1)
  {
    low_z = 1;
    high_z = 2;
  }
  else if(bin_z >= max_z)
  {
    low_z = max_z-1;
    high_z = max_z;
  }
  else if(zVal > bincenter_z)
  {
    low_z = bin_z;
    high_z = bin_z+1;
  }
  else
  {
    low_z = bin_z-1;
    high_z = bin_z;
  }

  double interp_val = 0.0;
  
  if(dim == 1)
  {
    double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_neg[high_y-1][high_z-1]->Eval(xVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_pos[high_y-1][high_z-1]->Eval(xVal);
      }
    }

    interp_val = (f_11*(a_2-yVal)*(b_2-zVal) + f_21*(yVal-a_1)*(b_2-zVal) + f_12*(a_2-yVal)*(zVal-b_1) + f_22*(yVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 2)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_neg[high_x-1][high_z-1]->Eval(yVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_pos[high_x-1][high_z-1]->Eval(yVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-zVal) + f_21*(xVal-a_1)*(b_2-zVal) + f_12*(a_2-xVal)*(zVal-b_1) + f_22*(xVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 3)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_neg[high_x-1][high_y-1]->Eval(zVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_pos[high_x-1][high_y-1]->Eval(zVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-yVal) + f_21*(xVal-a_1)*(b_2-yVal) + f_12*(a_2-xVal)*(yVal-b_1) + f_22*(xVal-a_1)*(yVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }

  return interp_val;
}

bool ananew::IsInsideBoundaries(TVector3 const& point) const
{
  return !(
           (TMath::Abs(point.X()) <= 0.0) || (TMath::Abs(point.X()) >= 360.0)
           || (point.Y()             <= 5.2) || (point.Y()             >= 604.0)
           || (point.Z()             <= -0.5) || (point.Z()             >= 695.3)
           );
} 
  
bool ananew::IsTooFarFromBoundaries(TVector3 const& point) const
{
  return (
          (TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X())  >= 360.0)
          || (point.Y()             < -14.8) || (point.Y()              >  624.0)
          || (point.Z()             < -20.5) || (point.Z()              >  715.3)
          );
}

TVector3 ananew::PretendAtBoundary(TVector3 const& point) const
{
  double x = point.X(), y = point.Y(), z = point.Z();
  
  
  if      (TMath::Abs(point.X()) ==    0.0    ) x =                           -0.00001;
  else if (TMath::Abs(point.X()) <	 0.00001) x =   TMath::Sign(point.X(),1)*0.00001; 
  else if (TMath::Abs(point.X()) >=    360.0  ) x = TMath::Sign(point.X(),1)*359.99999;
  
  if      (point.Y() <=   5.2) y =   5.20001;
  else if (point.Y() >= 604.0) y = 603.99999;
  
  if      (point.Z() <=   -0.5) z =   -0.49999;
  else if (point.Z() >= 695.3) z = 695.29999;
  
  return {x, y, z};
}


void ananew::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ananew.C
//      root> ananew t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   LoadHist();

   TFile f2("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
   TH3F *RecoFwd_Displacement_Z_Neg = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Neg");
   TH3F *RecoFwd_Displacement_Z_Pos = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Pos");

   const int nwires_in_slice = 20;
   const int nslices = 480/nwires_in_slice;
   Int_t nbinse=12; 
   Int_t nbinsthickness = 100;

   double NA=6.02214076e23;
   double MAr=39.95; //gmol
   double Density = 1.39; // g/cm^3

   TFile f("histnew.root","recreate");
   double interaction[nslices];
   double true_interaction[nslices];
   double true_abs[nslices];
   double true_cex[nslices];
   double signal[nslices];
   double incident[nslices];
   double true_incident[nslices];
   TH1D *incE[nslices];
   TH1D *pitch[nslices];
   TH1D *dEdx[nslices];
   TH2D *wirenum = new TH2D("wirenum","wirenum",480,0,480,480,0,480);
   for (int i = 0; i<nslices; ++i){
     interaction[i] = 0;
     true_interaction[i] = 0;
     true_abs[i] = 0;
     true_cex[i] = 0;
     signal[i] = 0;
     incident[i] = 0;
     true_incident[i] = 0;
     incE[i] = new TH1D(Form("incE_%d",i),Form("Incident energy, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     pitch[i] = new TH1D(Form("pitch_%d",i),Form("Slice thickness, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
     dEdx[i] = new TH1D(Form("dEdx_%d",i),Form("dE/dx, %d<wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
   }

   TH1D *dslcID = new TH1D("dslcID","reco slice ID - true slice ID",20,-10,10);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     if(abs(true_beam_PDG) != 211) continue;
     //if(abs(true_beam_PDG) != 2212) continue;
     
     // truth studes
     
     int true_sliceID = -1;
     
     double true_endz = true_beam_endZ;
     double offset_z = 0;
     TVector3 point(true_beam_endX, true_beam_endY, true_beam_endZ);
     if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
     
     if (point.X()>0){
       //true_endz += RecoFwd_Displacement_Z_Pos->GetBinContent(RecoFwd_Displacement_Z_Pos->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
       offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Pos, point.X(), point.Y(), point.Z(), 3, 1, 2);
     }
     else{
       //true_endz += RecoFwd_Displacement_Z_Neg->GetBinContent(RecoFwd_Displacement_Z_Neg->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
       offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Neg, point.X(), point.Y(), point.Z(), 3, 1, 1);
     }
     //std::cout<<true_beam_endX<<" "<<true_beam_endY<<" "<<true_beam_endZ<<" "<<offset_z<<std::endl;
     true_endz += offset_z;
     true_sliceID = int((true_endz-0.5603500-0.479/2)/0.479/nwires_in_slice); //z=0.56 is z coordinate for wire 0
     wirenum->Fill(int(true_endz-0.5603500-0.479/2)/0.479, reco_beam_calo_wire->back());
     if (true_sliceID <0) true_sliceID = 0;

     for (int i = 0; i<=true_sliceID; ++i){
       if (i<nslices) ++true_incident[i];
     }

//     if (true_sliceID == 5){
//       if ((*true_beam_endProcess) != "pi+Inelastic"){
//         std::cout<<(*true_beam_endProcess)<<" "<<run<<" "<<subrun<<" "<<event<<" "<<true_beam_endZ<<" "<<true_endz<<" "<<int(true_endz-0.5603500-0.479/2)/0.479<<std::endl;
//       }
//     }
     //std::cout<<(*true_beam_endProcess)<<std::endl;
     if ((*true_beam_endProcess) == "pi+Inelastic"){
     //if ((*true_beam_endProcess) == "protonInelastic"){
       if (true_sliceID < nslices){
         ++true_interaction[true_sliceID];
         if (true_daughter_nPi0 == 0 && true_daughter_nPiPlus == 0 && true_daughter_nPiMinus == 0){//true absorption
           ++true_abs[true_sliceID];
         }
         if (true_daughter_nPi0 == 1 && true_daughter_nPiPlus == 0  && true_daughter_nPiMinus == 0){//true absorption
           ++true_cex[true_sliceID];
         }
       }
     }
     
     if(!manual_beamPos_mc(reco_beam_startX, reco_beam_startY, reco_beam_startZ, 
                           reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ,
                           true_beam_startDirX, true_beam_startDirY, true_beam_startDirZ,
                           true_beam_startX, true_beam_startY, true_beam_startZ)) continue;
     if(!endAPA3(reco_beam_endZ) )continue;
     //std::cout<<*true_beam_endProcess<<" "<<true_beam_endX<<" "<<true_beam_endY<<" "<<true_beam_endZ<<std::endl;
     
     
     //std::cout<<(*true_beam_endProcess)<<std::endl;
     //if ((*true_beam_endProcess) == "pi+Inelastic" && abs(true_beam_PDG) == 211){  //std::cout<<"signal"<<std::endl;
     
     //interaction slice ID based on the track end wire number
     int sliceID = reco_beam_calo_wire->back()/nwires_in_slice;

//     if (sliceID == 5 && true_sliceID!=5){
//       std::cout<<sliceID<<" "<<true_sliceID<<" "<<(*true_beam_endProcess)<<" "<<run<<" "<<subrun<<" "<<event<<" "<<true_beam_endZ<<" "<<true_endz<<" "<<int(true_endz-0.5603500-0.479/2)/0.479<<std::endl;
//     }
     //Add code to determine if it is absorption
     bool isAbs = false;
     
     if (sliceID>=nslices) continue;
     
     if (true_sliceID!=-1){
       dslcID->Fill(sliceID - true_sliceID);
       if (sliceID == true_sliceID){
         ++signal[sliceID];
       }
     }
     
     //increment interaction counter
     ++interaction[sliceID];
     //increment incident counter
     for (int i = 0; i<=sliceID; ++i){
       ++incident[i];
     }
     
     std::vector<std::vector<double>> vpitch(nslices);
     std::vector<std::vector<double>> vincE(nslices);
     std::vector<std::vector<double>> vdEdx(nslices);
     for (size_t i = 0; i<reco_beam_calo_wire->size(); ++i){
       int this_wire = (*reco_beam_calo_wire)[i];
       int this_sliceID = this_wire/nwires_in_slice;
       //ignore the last slice for pitch and incident energy calculations
       if (this_sliceID>=sliceID) continue;
       
       double this_incE = (*reco_beam_incidentEnergies)[i];
       double this_pitch = (*reco_beam_TrkPitch_SCE)[i];
       double this_dEdx = (*reco_beam_dEdX_SCE)[i];
       //std::cout<<this_pitch<<std::endl;
       vpitch[this_sliceID].push_back(this_pitch);
       vincE[this_sliceID].push_back(this_incE);
       vdEdx[this_sliceID].push_back(this_dEdx);
     }
     for (size_t i = 0; i<vpitch.size(); ++i){
       if (!vpitch[i].empty()){
         double sum_pitch = 0;
         for (size_t j = 0; j<vpitch[i].size(); ++j){
           //std::cout<<vpitch[i][j]<<std::endl;
           sum_pitch += vpitch[i][j];
         }
         //std::cout<<sum_pitch<<" "<<vpitch[i].size()<<std::endl;
         pitch[i]->Fill(sum_pitch/vpitch[i].size()*nwires_in_slice);
       }
     }
     for (size_t i = 0; i<vincE.size(); ++i){
       if (!vincE[i].empty()){
         double sum_incE = 0;
         for (size_t j = 0; j<vincE[i].size(); ++j){
           sum_incE += vincE[i][j];
         }
         incE[i]->Fill(sum_incE/vincE[i].size());
       }
     }
     for (size_t i = 0; i<vdEdx.size(); ++i){
       if (!vdEdx[i].empty()){
         double sum_dEdx = 0;
         for (size_t j = 0; j<vdEdx[i].size(); ++j){
           sum_dEdx += vdEdx[i][j];
         }
         dEdx[i]->Fill(sum_dEdx/vdEdx[i].size());
       }
     }
   }
   
   double slcid[nslices];
   double avg_incE[nslices];
   double avg_pitch[nslices];
   double avg_dEdx[nslices];
   double err_dEdx[nslices];
   double eff_int[nslices];
   double eff_inc[nslices];
   double pur[nslices];
   double xs[nslices];
   double xs_abs[nslices];
   double xs_cex[nslices];
   double exs[nslices];
   double exs_abs[nslices];
   double exs_cex[nslices];
   for (int i = 0; i<nslices; ++i){
     slcid[i] = i;
     avg_incE[i] = incE[i]->GetMean();
     avg_pitch[i] = pitch[i]->GetMean();
     //avg_pitch[i] = 11;
     avg_dEdx[i] = dEdx[i]->GetMean();
     err_dEdx[i] = dEdx[i]->GetRMS()/sqrt(dEdx[i]->GetEntries());
     if (interaction[i]){
       eff_int[i] = signal[i]/true_interaction[i];
       eff_inc[i] = true_incident[i]/incident[i];
       pur[i] = signal[i]/interaction[i];
     }
     if (avg_pitch[i]&&true_incident[i]){
//       xs[i] = MAr/(Density*NA*avg_pitch[i])*true_interaction[i]/true_incident[i]*1e27;
//       xs_abs[i] = MAr/(Density*NA*avg_pitch[i])*true_abs[i]/true_incident[i]*1e27;
//       xs_cex[i] = MAr/(Density*NA*avg_pitch[i])*true_cex[i]/true_incident[i]*1e27;

       xs[i] = MAr/(Density*NA*avg_pitch[i])*log(true_incident[i]/(true_incident[i]-true_interaction[i]))*1e27;
       xs_abs[i] = MAr/(Density*NA*avg_pitch[i])*log(true_incident[i]/(true_incident[i]-true_abs[i]))*1e27;
       xs_cex[i] = MAr/(Density*NA*avg_pitch[i])*log(true_incident[i]/(true_incident[i]-true_cex[i]))*1e27;

       exs[i] = MAr/(Density*NA*avg_pitch[i])*1e27*sqrt(true_interaction[i]+pow(true_interaction[i],2)/true_incident[i])/true_incident[i];
       exs_abs[i] = MAr/(Density*NA*avg_pitch[i])*1e27*sqrt(true_abs[i]+pow(true_abs[i],2)/true_incident[i])/true_incident[i];
       exs_cex[i] = MAr/(Density*NA*avg_pitch[i])*1e27*sqrt(true_cex[i]+pow(true_cex[i],2)/true_incident[i])/true_incident[i];
     }
     else xs[i] = 0;
     //std::cout<<MAr<<" "<<Density<<" "<<NA<<" "<<avg_pitch[i]<<" "<<signal[i]<<" "<<incident[i]<<" "<<xs[i]<<std::endl;
   }

   TGraph *gr_int_slc = new TGraph(nslices-5, &(slcid[4]), &(interaction[4]));
   TGraph *gr_inc_slc = new TGraph(nslices-5, &(slcid[4]), &(incident[4]));
   TGraph *gr_incE_slc = new TGraph(nslices-5, &(slcid[4]), &(avg_incE[4]));
   TGraph *gr_pitch_slc = new TGraph(nslices-5, &(slcid[4]), &(avg_pitch[4]));
   TGraphErrors *gr_dEdx_slc = new TGraphErrors(nslices-5, &(slcid[4]), &(avg_dEdx[4]), 0, &(err_dEdx[4]));
   TGraph *gr_effint_slc = new TGraph(nslices-5, &(slcid[4]), &(eff_int[4]));
   TGraph *gr_effinc_slc = new TGraph(nslices-5, &(slcid[4]), &(eff_inc[4]));
   TGraph *gr_pur_slc = new TGraph(nslices-5, &(slcid[4]), &(pur[4]));
   TGraphErrors *gr_xs_incE = new TGraphErrors(nslices-5, &(avg_incE[4]), &(xs[4]), 0, &exs[4]);
   TGraphErrors *gr_xsabs_incE = new TGraphErrors(nslices-5, &(avg_incE[4]), &(xs_abs[4]), 0, &exs_abs[4]);
   TGraphErrors *gr_xscex_incE = new TGraphErrors(nslices-5, &(avg_incE[4]), &(xs_cex[4]), 0, &exs_cex[4]);
   TGraph *gr_trueint_slc = new TGraph(nslices-5, &(slcid[4]), &(true_interaction[4]));
   TGraph *gr_trueinc_slc = new TGraph(nslices-5, &(slcid[4]), &(true_incident[4]));

   f.Write();
   gr_int_slc->Write("gr_int_slc");
   gr_trueint_slc->Write("gr_trueint_slc");
   gr_inc_slc->Write("gr_inc_slc");
   gr_trueinc_slc->Write("gr_trueinc_slc");
   gr_incE_slc->Write("gr_incE_slc");
   gr_pitch_slc->Write("gr_pitch_slc");
   gr_dEdx_slc->Write("gr_dEdx_slc");
   gr_effint_slc->Write("gr_effint_slc");
   gr_effinc_slc->Write("gr_effinc_slc");
   gr_pur_slc->Write("gr_pur_slc");
   gr_xs_incE->Write("gr_xs_incE");
   gr_xsabs_incE->Write("gr_xsabs_incE");
   gr_xscex_incE->Write("gr_xscex_incE");
   f.Close();

}
