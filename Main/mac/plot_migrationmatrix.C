TMatrix CalculateMigrationMatrix(TH2D* _h_true_reco_mom){

TMatrix _S;
_S.Clear();
int _n=_h_true_reco_mom->GetNbinsY();
int _m=_h_true_reco_mom->GetNbinsX();

std::cout<<" _m bins of X axis is "<<_m<<std::endl;
std::cout<<" _n bins of Y axis is "<<_n<<std::endl;
_S.ResizeTo(_n + 1, _m + 1);

for (int j = 0; j < _m + 2; j++) {    // True bin
   int true_idx = j-1;
      if (j == 0)    true_idx = _m;
      if (j == _m+1) true_idx = _m;
      std::vector<double> p_v;
      p_v.resize(_n + 2);

      double sum = 0;

      for (int i = 0; i < _n + 2; i++) {      // Reco bin
       p_v.at(i) = _h_true_reco_mom->GetBinContent(j, i);
        sum += p_v.at(i);
      } //reco bin
     double tot_prob = 0;

      for (int i = 1; i < _n + 1; i++) {

        if (sum == 0 || std::isnan(sum))
          p_v.at(i) = 0;
        else
          p_v.at(i) /= sum;

        tot_prob += p_v.at(i);

        _S[i - 1][true_idx] += p_v.at(i);
      }
      // Add over/under-flow
      _S[_n][true_idx] = p_v.at(0) / sum + p_v.at(_n + 1) / sum;
} // true bin
     return _S;

}
void plot_migrationmatrix(){
gStyle->SetOptStat(0);
TFile *inputfile = new TFile("output_test.root");

    int n_bins_mumom = 50;
    int n_bins_mucostheta = 50;
    int n_bins_muphi=50;

    double bins_mumom[51]; 
    double bins_mucostheta[51]; 
    double bins_muphi[51];

    for(int i=0; i<n_bins_mumom+1; i++){
         bins_mumom[i]=0.0+1.2/n_bins_mumom*i;
    }
    for(int i=0; i<n_bins_mucostheta+1; i++){
         bins_mucostheta[i]=-1.0+2.0/n_bins_mucostheta*i;
    }
    for(int i=0; i<n_bins_muphi+1; i++){
         bins_muphi[i]=-TMath::Pi()+2.0*TMath::Pi()/n_bins_muphi*i;
    }



TLatex *prelim_Right = new TLatex(0.9, 0.93, "protoDUNE Simulation");
 
prelim_Right->SetTextFont(62);
prelim_Right->SetTextColor(kGray+1); 
prelim_Right->SetNDC(); 
prelim_Right->SetTextSize(0.05); 
prelim_Right->SetTextAlign(32); 
//prelim->Draw();

TLatex* prelim_Left = new TLatex(0.4,0.8, "protoDUNE Simulation");
prelim_Left->SetTextFont(62);
prelim_Left->SetTextColor(kGray+1); 
prelim_Left->SetNDC(); 
prelim_Left->SetTextSize(0.05); 
prelim_Left->SetTextAlign(32); 
 

TH2D *h_true_reco_mom;
h_true_reco_mom=(TH2D*)inputfile->Get("h_true_reco_mom");

TCanvas *c_mom=new TCanvas("c_mom", "c_mom", 900, 900);
h_true_reco_mom->Draw("colz");

TMatrix _S_mom = CalculateMigrationMatrix(h_true_reco_mom);

int bins_x = h_true_reco_mom->GetNbinsX()+1;
int bins_y = h_true_reco_mom->GetNbinsY()+1;
int _m = bins_x - 1;
int _n = bins_y - 1;
TH2D * smearing_matrix_mom= new TH2D("smearing_matrix_mom", "", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
smearing_matrix_mom->GetXaxis()->SetTitle("p_{#mu} (Truth) [GeV/c]");
smearing_matrix_mom->GetYaxis()->SetTitle("p_{#mu} (Reco) [GeV/c]");
smearing_matrix_mom->GetYaxis()->SetTitleOffset(1.4);
 for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_mom->SetBinContent(j+1, i+1, _S_mom[i][j]);

      }
    }
//smearing_matrix_mom->GetXaxis()->SetRangeUser(0.0,2.5);
//smearing_matrix_mom->GetYaxis()->SetRangeUser(0.0,2.5);
smearing_matrix_mom->Draw("colz");
prelim_Right->Draw("same");
c_mom->SaveAs("Fractional_MigrationMatrix_trkmom.png");




TH2D *h_true_reco_costheta;
h_true_reco_costheta=(TH2D*)inputfile->Get("h_true_reco_costheta");
TCanvas *c_costheta=new TCanvas("c_costheta", "c_costheta", 900, 900);
h_true_reco_costheta->Draw("colz");

TMatrix _S_costheta = CalculateMigrationMatrix(h_true_reco_costheta);

bins_x = h_true_reco_costheta->GetNbinsX()+1;
bins_y = h_true_reco_costheta->GetNbinsY()+1;
_m = bins_x - 1;
_n = bins_y - 1;
TH2D * smearing_matrix_costheta= new TH2D("smearing_matrix_costheta", "", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
smearing_matrix_costheta->GetXaxis()->SetTitle("cos#theta_{#mu} (Truth)");
smearing_matrix_costheta->GetYaxis()->SetTitle("cos#theta_{#mu} (Reco)");
smearing_matrix_costheta->GetYaxis()->SetTitleOffset(1.4);
 for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_costheta->SetBinContent(j+1, i+1, _S_costheta[i][j]);

      }
    }
smearing_matrix_costheta->Draw("colz");
prelim_Right->Draw("same");
c_costheta->SaveAs("Fractional_MigrationMatrix_trkcostheta.png");

TH2D *h_true_reco_phi;
h_true_reco_phi=(TH2D*)inputfile->Get("h_true_reco_phi");
TCanvas *c_phi=new TCanvas("c_phi", "c_phi", 900, 900);
h_true_reco_phi->Draw("colz");

TMatrix _S_phi = CalculateMigrationMatrix(h_true_reco_phi);

bins_x = h_true_reco_phi->GetNbinsX()+1;
bins_y = h_true_reco_phi->GetNbinsY()+1;
_m = bins_x - 1;
_n = bins_y - 1;
TH2D * smearing_matrix_phi= new TH2D("smearing_matrix_phi", "", n_bins_muphi, bins_muphi, n_bins_muphi, bins_muphi);
smearing_matrix_phi->GetXaxis()->SetTitle("cos#theta_{#mu} (Truth)");
smearing_matrix_phi->GetYaxis()->SetTitle("cos#theta_{#mu} (Reco)");
smearing_matrix_phi->GetYaxis()->SetTitleOffset(1.4);
 for (int i = 0; i < _n; i++) {
      for (int j = 0; j < _m; j++) {

        smearing_matrix_phi->SetBinContent(j+1, i+1, _S_phi[i][j]);

      }
    }
smearing_matrix_phi->Draw("colz");
prelim_Right->Draw("same");
c_phi->SaveAs("Fractional_MigrationMatrix_trkphi.png");




}
