#include "mu_run1_3_ly_ratio.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
using namespace std;

void mu_run1_3_ly_ratio() {
    
    //Open run1 and run 3 selection output root files
    TFile *run1_output = TFile::Open("run1_analyz_outputs.root");
    TFile *run3_output = TFile::Open("run3_analyz_outputs.root");
    
    //Get the 2 muon light yield distribution in x for run1 and run3
    TProfile* h_mu_ly_reco_x_run1 = (TProfile*)(run1_output->Get("h_mu_ly_reco_x"));
    TProfile* h_mu_ly_reco_x_run3 = (TProfile*)(run3_output->Get("h_mu_ly_reco_x"));
    
    //divide the run3 by run1 ly to get our muon ly ratio plot

    TH1D *p1 = h_mu_ly_reco_x_run3->ProjectionX();
    TH1D *p2 = h_mu_ly_reco_x_run1->ProjectionX();
    //p1->Divide(p2);
    TH1D* ratio = new TH1D("ratio","Reconstructed Cosmic Muon Run3/Run1 Light Yield Ratio",10,0,256.4);

    //get the bins for the ratio TProfile
    int bins = ratio->GetNbinsX();
    vector <double> ratio_vals;
    vector <double> ratio_err;

    for(int n=1; n<=bins; n++) {
        //get the values of the ratio for each bin
        double ratio_vals_single = p1->GetBinContent(n)/p2->GetBinContent(n);
        ratio_vals.push_back(ratio_vals_single);
        //now the error propogation to calculateing the error
        double ratio_err_single = TMath::Sqrt(pow((p2->GetBinError(n)/p2->GetBinContent(n)),2) + pow((p1->GetBinError(n)/p1->GetBinContent(n)),2))*ratio_vals_single;
        ratio_err.push_back(ratio_err_single);
        
        ratio->SetBinContent(n,(p1->GetBinContent(n)/p2->GetBinContent(n)));
        ratio -> SetBinError(n,ratio_err_single);
    }
    TF1 *func = new TF1("func", "[0]*x+[1]", 25.64, 230.76);
    ratio->Fit(func, "R");

    
    //Calculate the values and the error of % change of PE/cm as calculated in Patrick's analysis
    //first get the value of ly for cathode and anode piercing muons in run 1 and 3
    double mu_ly_run1_anode = p2 -> GetBinContent(1);
    double mu_ly_run3_anode = p1 -> GetBinContent(1);
    double mu_ly_run1_cathode = p2 -> GetBinContent(10);
    double mu_ly_run3_cathode = p1 -> GetBinContent(10);
    
    
    //now calculate the % change of PE/cm
    //double temp1 = mu_ly_run1_cathode - mu_ly_run3_cathode;
    double perc_change_PE_per_cm_cathode = 1- (mu_ly_run3_cathode / mu_ly_run1_cathode);
        
    //double temp2 = mu_ly_run1_anode - mu_ly_run3_anode;
    double perc_change_PE_per_cm_anode = 1- (mu_ly_run3_anode / mu_ly_run1_anode);
        
    //now error propagation
    //double cathode_ly_difference_err = TMath::Sqrt( pow(p2->GetBinError(10),2) + pow(p1->GetBinError(10),2) );
    double cathode_ly_perc_change_err = perc_change_PE_per_cm_cathode * TMath::Sqrt( pow( (p1->GetBinError(10)/mu_ly_run3_cathode) ,2) + pow( (p2->GetBinError(10)/mu_ly_run1_cathode) ,2) );
        
    //double anode_ly_difference_err = TMath::Sqrt( pow(p2->GetBinError(1),2) + pow(p1->GetBinError(1),2) );
    double anode_ly_perc_change_err = perc_change_PE_per_cm_anode * TMath::Sqrt( pow( (p1->GetBinError(1)/mu_ly_run3_anode),2) + pow( (p2->GetBinError(1)/mu_ly_run1_anode),2) );
        
    cout << "Cathode % change of PE/cm: " << perc_change_PE_per_cm_cathode << " Error for Cathode % change of PE/cm: " << cathode_ly_perc_change_err << endl;
    cout << "Anode % change of PE/cm: " << perc_change_PE_per_cm_anode << " Error for Anode & change of PE/cm: " << anode_ly_perc_change_err << endl;
    
    //create a new plot for my points and try to overlap it on patrick's plot
    double x[2] = {4,4.03};
    double y[2] = {-perc_change_PE_per_cm_cathode*100,-perc_change_PE_per_cm_anode*100};
    double ex[2] = {0,0};
    double ey[2] = {cathode_ly_perc_change_err*100,anode_ly_perc_change_err*100};
    auto g = new TGraphErrors(2,x,y,ex,ey);

    TCanvas* c2 = new TCanvas("c2");
    c2 -> cd();
    ratio -> Draw();
    ratio -> GetXaxis() -> SetTitle("Track Mid-point (x)(cm)");
    ratio -> GetYaxis() -> SetTitle("mu_LY_run3/mu_LY_run1");
    func->Draw("SAME");
    gStyle -> SetOptStat(10);
    
    TCanvas* c3 = new TCanvas("c3");
    c3 -> cd();
    g -> GetXaxis() -> SetLimits(0.,9.);
    g -> GetHistogram() -> SetMaximum(10.);
    g -> GetHistogram() -> SetMinimum(-60.);
    g -> Draw();
    
    //Overplotting different histograms
    
    
}
