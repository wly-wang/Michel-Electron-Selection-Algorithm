#define Myclass_cxx
#include "Myclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
using namespace std;

//function that calculates the track length
float calc_len(float stpx,float stpy,float stpz,float enpx,float enpy,float enpz)
{
    return (sqrt(pow(enpx-stpx, 2) + pow(enpy-stpy, 2) + pow(enpz-stpz, 2)));
}
            
float calc_k(float y1, float y2, float x1, float x2)
{
    return ((y2-y1)/(x2-x1));
}

float calc_comp_coord_4_intersec(float x, float y, float k)
{
    return ((115 - (y + k*x))/k);
}

float calc_flash_evt_dist(float flash_ycent, float flash_zcent, float starty, float endy,float startz, float endz)
{
    return (sqrt( pow(flash_zcent - ((startz+endz)/2.),2) + pow(flash_ycent - ((starty+endy)/2.),2) ));
}

float calc_devi_flash_trk(float flash_zcent, float flash_zwid, float startz, float endz)
{
    return (abs(flash_zcent - ((startz+endz)/2.))/flash_zwid);
}

void Myclass::Loop()
{
    //   In a ROOT session, you can do:
    //      root> .L Myclass.C
    //      root> Myclass t
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
    
    
    //defining new histograms
    TH1F * startposy = new TH1F("startposy","y-Start Position",44,-116.5,116.5);
    TH1F * startposy_FV = new TH1F("startposy_FV","startposy_FV",44,-116.5,116.5);
    TH1F * startposy_stp_trk = new TH1F("startposy_stp_trk","startposy_FV",44,-116.5,116.5);
    TH1F * startposy_FV2 = new TH1F("startposy_FV2","startposy_FV2",44,-116.5,116.5);
    TH1F * startposy_ME = new TH1F("startposy_ME","Michel Electron y-Start Position",5,-116.5,116.5);
    startposy_FV->SetLineColor(kRed);
    startposy_FV2->SetLineColor(kRed);
    startposy_stp_trk->SetLineColor(kMagenta);
    
    TH1F * startposx = new TH1F("startposx","x-Start Position",44,0,256.4);
    TH1F * startposx_FV = new TH1F("startposx_FV","startposx_FV",44,0,256.4);
    TH1F * startposx_stp_trk = new TH1F("startposx_stp_trk","startposx_stp_trk",44,0,256.4);
    TH1F * startposx_FV2 = new TH1F("startposx_FV2","startposx_FV2",44,0,256.4);
    TH1F * startposx_ME = new
    TH1F("startposx_ME","Michel Electron x-Start Position",5,0,256.4);
    startposx_FV->SetLineColor(kRed);
    startposx_FV2->SetLineColor(kRed);
    startposx_stp_trk->SetLineColor(kMagenta);
    
    TH1F * startposz = new TH1F("startposz","z-Start Position",44,0,1036.8);
    TH1F * startposz_FV = new TH1F("startposz_FV","startposz_FV",44,0,1036.8);
    TH1F * startposz_stp_trk = new TH1F("startposz_stp_trk","startposz_stp_trk",44,0,1036.8);
    TH1F * startposz_FV2 = new TH1F("startposz_FV2","startposz_FV2",44,0,1036.8);
    TH1F * startposz_ME = new TH1F("startposz_ME","Michel Electron z-Start Position",5,0,1036.8);
    startposz_FV->SetLineColor(kRed);
    startposz_FV2->SetLineColor(kRed);
    startposz_stp_trk->SetLineColor(kMagenta);
    
    TH1F * endposx = new TH1F("endposx","x-End Position",44,0,256.4);
    TH1F * endposx_FV = new TH1F("endposx_FV","endposx_FV",44,0,256.4);
    TH1F * endposx_stp_trk = new TH1F("endposx_stp_trk","endposx_stp_trk",44,0,256.4);
    TH1F * endposx_muon_selection = new TH1F("endposx_muon_selection","x-End Position of Selected Muons",13,0,256.4);
    TH1F * endposx_FV2 = new TH1F("endposx_FV2","endposx_FV2",44,0,256.4);
    endposx_FV->SetLineColor(kRed);
    endposx_FV2->SetLineColor(kRed);
    endposx_stp_trk->SetLineColor(kMagenta);
    endposx_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * endposz = new TH1F("endposz","z-End Position",44,0,1036.8);
    TH1F * endposz_FV = new TH1F("endposz_FV","endposz_FV",44,0,1036.8);
    TH1F * endposz_stp_trk = new TH1F("endposz_stp_trk","endposz_stp_trk",44,0,1036.8);
    TH1F * endposz_muon_selection = new TH1F("endposz_muon_selection","z-End Position of Selected Muons",52,0,1036.8);
    TH1F * endposz_FV2 = new TH1F("endposz_FV2","endposz_FV2",44,0,1036.8);
    endposz_FV->SetLineColor(kRed);
    endposz_FV2->SetLineColor(kRed);
    endposz_stp_trk->SetLineColor(kMagenta);
    endposz_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * endposy = new TH1F("endposy","y-End Position",44,-116.5,116.5);
    TH1F * endposy_FV = new TH1F("endposy_FV","endposy_FV",44,-116.5,116.5);
    TH1F * endposy_stp_trk = new TH1F("endposy_stp_trk","endposy_stp_trk",44,-116.5,116.5);
    TH1F * endposy_muon_selection = new TH1F("endposy_muon_selection","y-End Position of Selected Muons",12,-116.5,116.5);
    TH1F * endposy_FV2 = new TH1F("endposy_FV2","endposy_FV2",44,-116.5,116.5);
    endposy_FV->SetLineColor(kRed);
    endposy_FV2->SetLineColor(kRed);
    endposy_stp_trk->SetLineColor(kMagenta);
    endposy_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * track_length = new
    TH1F ("track_length", "Track Length", 80, 0, 450);
    TH1F * track_length_FV = new
    TH1F ("track_length_FV", "track_length_FV", 80, 0, 450);
    TH1F * track_length_stp_trk = new
    TH1F ("track_length_stp_trk", "Track Length", 80, 0, 450);
    TH1F * track_length_score_cuts = new
    TH1F ("track_length_score_cuts", "Track Length w/ Score Cuts", 80, 0, 450);
    TH1F * track_length_FV2 = new
    TH1F ("track_length_FV2", "track_length_FV2", 80, 0, 450);
    TH1F * track_length_ME = new
    TH1F ("track_length_ME", "track_length_ME", 80, 0, 450);
    TH1F * track_length_ME_no_n_showers = new
    TH1F ("track_length_ME_no_n_showers", "track_length_ME_no_n_showers", 80, 0, 450);
    track_length_FV->SetLineColor(kRed);
    track_length_FV2->SetLineColor(kRed);
    track_length_stp_trk->SetLineColor(kMagenta);
    track_length_score_cuts->SetLineColor(kViolet+3);
    track_length_ME->SetLineColor(kViolet+3);
    track_length_ME_no_n_showers->SetLineColor(kViolet+3);
    
    
    TH1F * track_score = new
    TH1F ("track_score", "Track Score", 50, 0, 1);
    TH1F * track_score_cut = new TH1F ("track_score_cut", "track_score_cut", 50, 0, 1);
    TH1F * track_score_ME = new
    TH1F ("track_score_ME", "Track Score of Selected Michel Electron Candidates", 50, 0, 1);
    track_score_cut->SetLineColor(kRed);
    track_score_ME->SetLineColor(kRed);
    
    TH1F * track_PID = new
    TH1F ("track_PID", "Track PID Score", 50, -1, 1);
    TH1F * track_PID_cut = new TH1F ("track_PID_cut", "track_PID_cut", 50, -1, 1);
    TH1F * track_PID_ME = new TH1F ("track_PID_ME", "track_PID_ME", 50, -1, 1);
    track_PID_cut->SetLineColor(kRed);
    track_PID_ME->SetLineColor(kRed);
    
    TH1F * MuonME_gap = new
    TH1F ("MuonME_gap", "Muon/Michel Electron Candidates Gap Distance", 100, 0, 150);
    TH1F * MuonME_gap_zoom = new
    TH1F ("MuonME_gap_zoom", "Zoomed in Muon/Michel Electron Candidates Gap Distance", 100, 0, 20);
    //TH1F * MuonME_start_gap = new
    //TH1F ("MuonME_start_gap", "MuonME_start_gap", 80, 0, 450);
    MuonME_gap->SetLineColor(kRed);
    MuonME_gap_zoom->SetLineColor(kMagenta);
    //MuonME_start_gap->SetLineColor(kViolet+3);
    
    TH1F * Track_Energy = new
    TH1F ("Track_Energy", "Track Energy", 100, 0, 4000);
    TH1F * Track_Energy_Muon = new
    TH1F ("Track_Energy_Muon", "Muon Track Energy", 100, 0, 4000);
    TH1F * Track_Energy_ME = new
    TH1F ("Track_Energy_ME", "Michel Electron Track Energy", 100, 0, 4000);
    TH1F * Track_Energy_u = new
    TH1F ("Track_Energy_u", "Track Energy (u plane)", 100, 0, 4000);
    TH1F * Track_Energy_Muon_u = new
    TH1F ("Track_Energy_Muon_u", "Muon Track Energy (u Plane)", 100, 0, 4000);
    TH1F * Track_Energy_ME_u = new
    TH1F ("Track_Energy_ME_u", "Michel Electron Track Energy (u Plane)", 100, 0, 4000);
    TH1F * Track_Energy_v = new
    TH1F ("Track_Energy_v", "Track Energy (v plane)", 100, 0, 4000);
    TH1F * Track_Energy_Muon_v = new
    TH1F ("Track_Energy_Muon_v", "Muon Track Energy (v Plane)", 100, 0, 4000);
    TH1F * Track_Energy_ME_v = new
    TH1F ("Track_Energy_ME_v", "Michel Electron Track Energy (u Plane)", 100, 0, 4000);
    Track_Energy_Muon->SetLineColor(kRed);
    Track_Energy_ME->SetLineColor(kMagenta);
    Track_Energy_Muon_u->SetLineColor(kRed);
    Track_Energy_ME_u->SetLineColor(kMagenta);
    Track_Energy_Muon_v->SetLineColor(kRed);
    Track_Energy_ME_v->SetLineColor(kMagenta);
    

    TH1F * mc_end_x = new TH1F("mc_end_x","Monte Carlo Truth End Position (x)",44,0,256.4);
    TH1F * mc_end_x_reco = new TH1F("mc_end_x_reco","Monte Carlo True End Position of Reconstructed Tracks (x)",44,0,256.4);
    TH1F * ME_mc_end_x = new TH1F("ME_mc_end_x","x-End Position for Monte Carlo Michel Electrons",44,0,256.4);
    mc_end_x_reco->SetLineColor(kRed);
    
    TH1F * mc_end_y = new TH1F("mc_end_y","Monte Carlo Truth End Position (y)",44,-116.5,116.5);
    TH1F * mc_end_y_reco = new TH1F("mc_end_y_reco","Monte Carlo True End Position of Reconstructed Tracks (y)",44,-116.5,116.5);
    TH1F * mc_end_y_st_deep= new TH1F("mc_end_y_st_deep","Monte Carlo True End Position of Reconstructed Tracks for Events Starting Deep in the Detector (y)",44,-116.5,116.5);
    TH1F * mc_end_y_end_shallow= new TH1F("mc_end_y_end_shallow","Monte Carlo True End Position of Reconstructed Tracks for Events Ending Shallow in the Detector (y)",44,-116.5,116.5);
    TH1F * ME_mc_end_y = new TH1F("ME_mc_end_y","y-End Position for Monte Carlo Michel Electrons ",44,-116.5,116.5);
    mc_end_y_reco->SetLineColor(kRed);
    mc_end_y_st_deep->SetLineColor(kRed);
    mc_end_y_end_shallow -> SetLineColor(kMagenta);
    
    TH1F * mc_end_z = new TH1F("mc_end_z","Monte Carlo Truth End Position (z)",44,0,1036.8);
    TH1F * mc_end_z_reco = new TH1F("mc_end_z_reco","Monte Carlo True End Position of Reconstructed Tracks (z)",44,0,1036.8);
    TH1F * ME_mc_end_z = new TH1F("ME_mc_end_z","z-End Position for Monte Carlo Michel Electrons",44,0,1036.8);
    mc_end_z_reco->SetLineColor(kRed);
    
    TH1F * mc_v_x = new TH1F("mc_v_x","Monte Carlo Truth Start Position (x)",44,0,256.4);
    TH1F * mc_v_x_st_deep = new TH1F("mc_v_x_st_deep","Monte Carlo Truth Start Position for Events Starting Deep in the Detector (x)",44,0,256.4);
    TH1F * mc_v_x_end_shallow = new TH1F("mc_v_x_end_shallow","Monte Carlo Truth Start Position for Events Ending Shallow in the Detector (x)",44,0,256.4);
    mc_v_x_st_deep->SetLineColor(kRed);
    mc_v_x_end_shallow -> SetLineColor(kMagenta);
    
    TH1F * mc_v_y = new TH1F("mc_v_y","Monte Carlo Truth Start Position (y)",44,-116.5,120);

    TH1F * mc_v_z = new TH1F("mc_v_z","Monte Carlo Truth Start Position (z)",44,0,1036.8);
    TH1F * mc_v_z_st_deep = new TH1F("mc_v_z_st_deep","Monte Carlo Truth Start Position for Events Starting Deep in the Detector (z)",44,0,1036.8);
    TH1F * mc_v_z_end_shallow = new TH1F("mc_v_z_end_shallow","Monte Carlo Truth Start Position for Events Ending Shallow in the Detector (z)",44,0,1036.8);
    mc_v_z_st_deep->SetLineColor(kRed);
    mc_v_z_end_shallow -> SetLineColor(kMagenta);
    
    TH1F * mc_energy = new TH1F("mc_energy","Muon Monte Carlo True Energy",100,0.2, 0.6);
    TH1F * mc_energy_st_deep = new TH1F("mc_energy_st_deep","Monte Carlo True Energy for Muon Events Starting Deep in the Detector",100,0.2,0.6);
    TH1F * mc_energy_end_shallow = new TH1F("mc_energy_end_shallow","Monte Carlo True Energy for Muon Events Ending Shallow in the Detector",100,0.2, 0.6);
    mc_energy_st_deep->SetLineColor(kRed);
    mc_energy_end_shallow -> SetLineColor(kMagenta);
    
    TH1F * orig_me_mc_energy = new TH1F("me_mc_energy", "Michel Electron Monte Carlo True Energy", 100, 0, 60);
    TH1F * me_mc_energy = new TH1F("me_mc_energy", "Michel Electron Monte Carlo True Energy", 100, 0, 60);
    me_mc_energy -> SetLineColor(kRed);
    

    TH2F * mc_endy_vs_y_end_reco_st_deep = new TH2F("mc_endy_vs_y_end_reco_st_deep", "Monte Carlo True End Position (y) vs. y End Position for Reconstructed Tracks for Events Starting Deep in the Detector", 44, -116.5, 116.5, 44, -116.5, 116.5);
    TH2F * mc_endy_vs_y_end_reco_end_shallow = new TH2F("mc_endy_vs_y_end_reco_end_shallow", "Monte Carlo True End Position (y) vs. y End Position for Reconstructed Tracks for Events Ending Shallow in the Detector", 44, -116.5, 116.5, 44, -116.5, 116.5);
    
    TH2F * mc_E_vs_reco_E = new TH2F("mc_E_vs_reco_E", "Muon Monte Carlo True Energy vs. Track Calorimetric Energy for Reconstructed Tracks", 30, 200, 600, 30, 0, 600);
    TH2F * me_mc_E_v_reco_E = new TH2F("me_mc_E_v_reco_E", "Michel Electron Monte Carlo True Energy vs. Track Calorimetric Energy for Reconstructed Tracks",20, 0, 100, 20, 0, 500);
    
    TH2F* end_y_v_z_st_deep = new TH2F("end_y_v_z_st_deep", "The Reconstructed yz Plane End Position for Events Starting Deep in the Detector", 100, 0, 1036.8, 100, -116.5, 116.5);
    TH2F* end_y_v_z_end_shallow= new TH2F("end_y_v_z_end_shallow", "The Reconstructed yz Plane End Position for Events Ending Shallow in the Detector", 100, 0, 1036.8, 100, -116.5, 116.5);
    
    TH1F * mc_pdg_code = new TH1F("mc_pdg_code", "Monte Carlo PDG Code", 100, -15, 15);
    
    TH1F * h_me_dist_flsh_trk = new TH1F("h_me_dist_flsh_trk", "Distance Between Flash and Michel Electron Event Tracks", 100, 0, 500);
    TH1F * h_mu_dist_flsh_trk = new TH1F("h_mu_dist_flsh_trk", "Distance Between Flash and Muon Event Tracks", 100, 0, 500);
    TH1F * h_mc_mu_dist_flsh_trk = new TH1F("h_mc_mu_dist_flsh_trk", "Distance Between Flash and MC Muon Information Event Tracks", 100, 0, 500);
    
    TProfile* h_tot_ly_v_reco_x_recME = new TProfile("h_tot_ly_v_reco_x_recME", "Total Light Yield Measurent Using Reconstructed Michel Electrons", 10, 0, 256.4, "");
    TH1F* h_reco_x_counter_recME = new TH1F("h_reco_x_counter_recME", "Event Counter for reconstructed x position", 10, 0, 256.4);
    TH1F* h_tot_ly_recME = new TH1F("h_tot_ly_recME", "Total Light Yield of Reconstructed Michel Electrons", 10, 0, 256.4);
    
    TProfile* h_tot_ly_v_reco_x_trueMEe = new  TProfile("h_tot_ly_v_reco_x_trueMEe", "Total Light Yield Measurent Using Reconstructed Michel Electrons Using True Energy", 10, 0, 256.4, "");
    TH1F* h_tot_ly_trueMEe = new TH1F("h_tot_ly_trueMEe", "Total Light Yield of Reconstructed Michel Electrons Using True Energy", 10, 0, 256.4);
    
    
    TProfile* h_tot_ly_v_reco_x_recmu = new  TProfile("h_tot_ly_v_reco_x_recmu", "Total Light Yield Measurent Using Reconstructed Muons", 10, 0, 256.4, "");
    TH1F* h_reco_x_counter_recmu = new TH1F("h_reco_x_counter_recmu", "Event Counter for reconstructed x position", 10, 0, 256.4);
    TH1F* h_tot_ly_recmu = new TH1F("h_tot_ly_recmu", "Total Light Yield of Reconstructed Muons", 10, 0, 256.4);

    TProfile* h_tot_ly_v_reco_x_truemue = new  TProfile("h_tot_ly_v_reco_x_truemue", "Total Light Yield Measurent Using Reconstructed muons Using True Energy", 10, 0, 256.4, "");
    TH1F* h_tot_ly_truemue = new TH1F("h_tot_ly_truemue", "Total Light Yield of Reconstructed muons Using True Energy", 10, 0, 256.4);
    
    //open particle_selection txt file (for event display)
    fstream file;
    file.open("selected_muon_info.txt", ios::app | ios::out);
    
    //open me_info_len_score_pid (printing out to see if values are sensible
    fstream file1;
    file1.open("selected_me_info.txt", ios::app | ios::out);
    
    //open txt file to look at ME with energy >50 MeV
    fstream file2;
    file2.open("ME_E_more_50mev.txt", ios::app | ios::out);
    
    //open txt file for storing info about the events that end at the top of the detector
    fstream file3;
    file3.open("evt_end_top.txt", ios::app | ios::out);
    
    //open txt file for storing info about the events that start deep within the detector
    fstream file4;
    file4.open("evt_start_deep.txt", ios::app | ios::out);
    
    //create a root file that stores all the histogram/results
    std::unique_ptr<TFile> myFile(TFile::Open("me_selection.root", "RECREATE"));
    cout<< "ROOT file for storing all results is created" <<endl;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        //for the original histogram stored in TTree, fill on spot
        for(unsigned int i=0;
               i<(*mc_endx).size(); i++)
        {
            //helper variables
            float mcendx, mcendy, mcendz, mcvx, mcvy, mcvz, mcpdg, mcE, me_mc_E;
            mcendx = (*mc_endx)[i];
            mcendy =(*mc_endy)[i];
            mcendz =(*mc_endz)[i];
            mcvx = (*mc_vx)[i];
            mcvy = (*mc_vy)[i];
            mcvz = (*mc_vz)[i];
            mcpdg = (*mc_pdg)[i];
            mcE = (*mc_E)[i];
            me_mc_E = endmuonmichel;
            
            mc_end_x -> Fill(mcendx);
            mc_end_y -> Fill(mcendy);
            mc_end_z -> Fill(mcendz);
            mc_v_x -> Fill(mcvx);
            mc_v_y -> Fill(mcvy);
            mc_v_z -> Fill(mcvz);
            mc_energy -> Fill(mcE);
            orig_me_mc_energy -> Fill(me_mc_E);
            mc_pdg_code -> Fill(mcpdg);
            
            
        }
        //Defining new vectors for storing muon information
        vector<float> startposx_1;
        vector<float> startposy_1;
        vector<float> startposz_1;
        vector<float> endposx_1;
        vector<float> endposy_1;
        vector<float> endposz_1;
        vector<float> reco_me_energy;
        vector<float> trk_length;
        vector<float> tpc_mc_vx;
        vector<float> tpc_mc_vz;
        
        bool GoodMuonExists=false;
        //bool GoodMichelExists=false;
        vector<int> goodmuonindex;
        
        
        //cut out unphysical regions in our position variables
        for(unsigned int itrack=0;
            itrack<(*trk_sce_start_x_v).size(); itrack++)
        {
            
            
            //helper variabls
            float stposx, stposy, stposz, enposx, enposy, enposz, trklen, trk_score, trk_pid, trk_E_y, trk_E_u, trk_E_v, mcendx, mcendy, mcendz, mcvx, mcvy, mcvz, mcpdg, kx, kz, x_coord_st_intersec, z_coord_st_intersec, mc_trklen, mcE;
            stposx=(*trk_sce_start_x_v)[itrack];
            stposy=(*trk_sce_start_y_v)[itrack];
            stposz=(*trk_sce_start_z_v)[itrack];
            enposx=(*trk_sce_end_x_v)[itrack];
            enposy=(*trk_sce_end_y_v)[itrack];
            enposz=(*trk_sce_end_z_v)[itrack];
            trk_score=(*trk_score_v)[itrack];
            trk_pid=(*trk_llr_pid_score_v)[itrack];
            trk_E_y=(*trk_calo_energy_y_v)[itrack];
            trk_E_u=(*trk_calo_energy_u_v)[itrack];
            trk_E_v=(*trk_calo_energy_v_v)[itrack];
            mcendx = (*mc_endx)[itrack];
            mcendy =(*mc_endy)[itrack];
            mcendz =(*mc_endz)[itrack];
            mcvx = (*mc_vx)[itrack];
            mcvy = (*mc_vy)[itrack];
            mcvz = (*mc_vz)[itrack];
            mcpdg = (*mc_pdg)[itrack];
            mcE = (*mc_E)[itrack];
            
            //cout << mcendy << endl;
            
            //fill out mc end xyz since we want to check the original first
           
            /*
            if (itrack < 100)
            {
                cout << "mc_end_x: " << mcendx << endl;
                cout << "mc_end_y: " << mcendy << endl;
                cout << "mc_end_z: " << mcendz << endl;
                cout << "mc_pdg: " << mcpdg << endl;
            }
             */
             
            /* if (run == 6959 && sub == 201 && evt == 10056)
             {
             cout << "x start: " << stposx << " y start: " << stposy << " z start: " << stposz << " x end: " << enposx << " y end: " << enposy << " z end: " << enposz << endl;
             }
             else continue; */
            
            // use only sensible vectors
            if(stposx > 0 && stposy > -115 && stposz >0 && enposx >0 && enposy > -115 && enposz > 0 && trk_E_y > -1 )
            {
                trklen = calc_len(stposx,stposy,stposz,enposx,enposy,enposz); //calculate the track length of the event tracks
                
                //Now, for the simulated muons, calculate the x, z coordinates of the start point of their tracks first entering the detector.
                //start with the x-y plane
                //calculate the gradient of the tracks using the mc_vx, y and mc_endx, y
                //kx = calc_k(mcvy, mcendy, mcvx, mcendx);
                //calculate the x coordinate of the start point
                //x_coord_st_intersec = calc_comp_coord_4_intersec(mcendx, mcendy, kx);
                //pushback this value into the vector created for this everytime looping over each track
                //tpc_mc_vx.push_back(x_coord_st_intersec);
                
                //now repeat same step for z-y plane
                //calculate the gradient of the tracks using the mc_vz, y and mc_endz, y
                // kz = calc_k(mcvy, mcendy, mcvz, mcendz);
                //calculate the z coordinate of the start point
                //z_coord_st_intersec = calc_comp_coord_4_intersec(mcendz, mcendy, kz);
                //pushback this value into the vector created for this everytime looping over each track
                //tpc_mc_vz.push_back(z_coord_st_intersec);
                
                //Now, calculate the track length of the MC muon tracks within the detector
                mc_trklen = calc_len(mcvx, mcvy, mcvz, mcendx, mcendy, mcendz);
                //if (mc_trklen > 10)
                {
                    
                    if(enposy > stposy) // if end point is higher - flip.
                    {
                        float tempx,tempy,tempz;
                        tempy=enposy;enposy=stposy;stposy=tempy;
                        tempz=enposz;enposz=stposz;stposz=tempz;
                        tempx=enposx;enposx=stposx;stposx=tempx;
                    } // now end point is proper endpoint.
                    
                    //here can fill vectors or histograms (with no cuts/only non crazy values).
                    startposx -> Fill(stposx);
                    startposy -> Fill(stposy);
                    startposz -> Fill(stposz);
                    endposx -> Fill(enposx);
                    endposy -> Fill(enposy);
                    endposz -> Fill(enposz);
                    track_length -> Fill(trklen);
                    Track_Energy -> Fill(trk_E_y);
                    Track_Energy_u -> Fill(trk_E_u);
                    Track_Energy_v -> Fill(trk_E_v);
                    //request non-zero, non-crazy PDG code for Monte Carlo sample, 0 pdg means it's a real data particle and anything above 250 will not appear in this analysis
                    if (mcpdg != 0 && mcpdg < 100)
                    {
                        mc_end_x_reco -> Fill(mcendx);
                        mc_end_y_reco -> Fill(mcendy);
                        mc_end_z_reco -> Fill(mcendz);
                    }
                    
                    // add fiducial volume cuts
                    if(stposx > 10 && stposx < 246 && enposx > 10 && enposx < 246
                       && stposz >10 && stposz < 1026 && enposz > 10 && enposz < 1026 && trk_E_y != 0)
                    {
                        // fill second tier of histograms (start/end/length)
                        
                        startposx_FV -> Fill(stposx);
                        startposy_FV -> Fill(stposy);
                        startposz_FV -> Fill(stposz);
                        endposx_FV -> Fill(enposx);
                        endposy_FV -> Fill(enposy);
                        endposz_FV -> Fill(enposz);
                        track_length_FV -> Fill(trklen);
                
                        
                        if (mcpdg != 0 && mcpdg < 100 && stposy < 110)
                        {
                            mc_end_y_st_deep -> Fill(mcendy);
                            mc_v_x_st_deep -> Fill(mcvx);
                            mc_v_z_st_deep -> Fill(mcvz);
                            mc_energy_st_deep -> Fill(mcE);
                            mc_endy_vs_y_end_reco_st_deep -> Fill(enposy, mcendy);
                            end_y_v_z_st_deep -> Fill(enposz, enposy);
                            
                            file4 << left << setw(10) << run << setw(10) << sub << setw(10) << evt << setw(15) << stposx << setw(10) << stposy << setw(10) << stposz << setw(10) << enposx << setw(10) << enposy << setw(10)<< enposz << setw(10) << trklen << setw(15) << trk_score << setw(15) << trk_pid << setw(15) << mcvx << setw(15) << mcvy << setw(15) << mcvz << setw(15) << mcE << setw(15) << mcendx << setw(15) << mcendy << setw(15) << mcendz << endl;
                        }
                        
                        //select only tracks coming in from the top and stopping
                        // this could the if statement to select stopping track.
                        if(stposy > 110 && stposy < 115 && enposy > -110 )
                        {
                            //fill third tier of histograms,...
                            
                            startposx_stp_trk -> Fill(stposx);
                            startposy_stp_trk -> Fill(stposy);
                            startposz_stp_trk -> Fill(stposz);
                            endposx_stp_trk -> Fill(enposx);
                            endposy_stp_trk -> Fill(enposy);
                            endposz_stp_trk -> Fill(enposz);
                            track_length_stp_trk -> Fill(trklen);
                            track_score -> Fill(trk_score);
                            track_PID -> Fill(trk_pid);
                            
                            if (enposy > 105 && mcpdg != 0 && mcpdg < 100 && mcE > 0)
                            {
                                mc_end_y_end_shallow -> Fill(mcendy);
                                mc_v_x_end_shallow -> Fill(mcvx);
                                mc_v_z_end_shallow -> Fill(mcvz);
                                end_y_v_z_end_shallow -> Fill(enposz, enposy);
                                mc_energy_end_shallow -> Fill(mcE);
                                mc_endy_vs_y_end_reco_end_shallow -> Fill(enposy, mcendy);
                                
                                file3 << left << setw(10) << run << setw(10) << sub << setw(10) << evt << setw(15) << stposx << setw(10) << stposy << setw(10) << stposz << setw(10) << enposx << setw(10) << enposy << setw(10)<< enposz << setw(10) << trklen << setw(15) << trk_score << setw(15) << trk_pid << setw(15) << mcvx << setw(15) << mcvy << setw(15) << mcvz << setw(15) << mcE << setw(15) << mcendx << setw(15) << mcendy << setw(15) << mcendz << endl;
                            }
                            
                            // when all cuts satisfied, set GoodMuonExists=true;
                            
                            if(trk_score > 0.5 && trk_pid > 0 && trklen > 10 && mcpdg != 0 && mcpdg < 100 && mcE > 0)
                            {
                                track_score_cut -> Fill(trk_score);
                                track_PID_cut -> Fill(trk_pid);
                                track_length_score_cuts -> Fill(trklen);
                                endposx_muon_selection -> Fill(enposx);
                                endposy_muon_selection -> Fill(enposy);
                                endposz_muon_selection -> Fill(enposz);
                                Track_Energy_Muon -> Fill(trk_E_y);
                                Track_Energy_Muon_u -> Fill(trk_E_u);
                                Track_Energy_Muon_v -> Fill(trk_E_v);
                                if(mcpdg != 0 && mcpdg < 100 && mcE > 0)
                                {
                                    mc_E_vs_reco_E -> Fill(mcE*1000., trk_E_y);
                                    
                                    
                                }
                                
                                //calculate the distance between flash and track (mu)
                                float mu_distance_flash_trk = calc_flash_evt_dist(flash_ycenter, flash_zcenter, stposy, enposy, stposz, enposz);
                                h_mu_dist_flsh_trk -> Fill(mu_distance_flash_trk);
                                //flash matching condition
                                if (mu_distance_flash_trk < 160. && calc_devi_flash_trk(flash_zcenter, flash_zwidth, stposz, enposz) <= 1.)
                                {
                                    //create a variable with value 0 for the total photoelectrons (PEs)
                                    
                                    float tot_pe_per_trk = 0;
                                    //looping over all the pmt for the total PE per event
                                    for (unsigned int PMT_no = 0; PMT_no != 32; ++PMT_no)
                                    {
                                        tot_pe_per_trk += flash_pe_v->at(PMT_no)*(20./gain_ampl_v->at(PMT_no));
                                    }
                                    //now calculate the total LY for 4 cases
                                    //total LY for true Michel electron energy
                                    if (tot_pe_per_trk >0 && trk_E_y > 0)
                                    {
                                        float tot_ly_per_trk = tot_pe_per_trk/trk_E_y;
                                        float tot_ly_per_trk_mcE = tot_pe_per_trk/(mcE*1000.);
                                        
                                        float trk_mid_x = (stposx+enposx)/2.;
                                        
                                        //first fill rough ly plot (different no. of tracks in each bin
                                        h_tot_ly_v_reco_x_recmu -> Fill(trk_mid_x, tot_ly_per_trk);
                                        h_tot_ly_v_reco_x_truemue-> Fill(trk_mid_x, tot_ly_per_trk_mcE);
                                        
                                        /*
                                         for (unsigned int bin_no = 1; bin_no != 100; ++bin_no)
                                         {
                                         cout << h_tot_ly_v_reco_x_recME->GetBinContent(bin_no) << endl;
                                         }
                                         */
                                        //int bin = h_tot_ly_v_reco_x_recME->FindBin(trk_mid_x);
                                        //h_tot_ly_v_reco_x_recME->AddBinContent(bin, tot_ly_per_trk);
                                        // Next fill a reco_x midpoint counter to count no of evt in each bin
                        
                                            h_reco_x_counter_recmu -> Fill(trk_mid_x);
                                            
                                    
                                            //next get bin content of the two histo, divide the first by the second to get a normalized and final light yield plot
                                        
                                            for (unsigned int bin_no = 1; bin_no != 100; ++bin_no)
                                            {
                                                if(h_reco_x_counter_recmu->GetBinContent(bin_no) != 0)
                                                {
                                                    h_tot_ly_recmu->SetBinContent(bin_no,(h_tot_ly_v_reco_x_recmu->GetBinContent(bin_no)/h_reco_x_counter_recmu->GetBinContent(bin_no)));
                                                    float err = sqrt(h_tot_ly_v_reco_x_recmu->GetBinContent(bin_no))/h_reco_x_counter_recmu->GetBinContent(bin_no);
                                                    h_tot_ly_recmu->SetBinError(bin_no, err);
                                                    
                                                    h_tot_ly_truemue->SetBinContent(bin_no,(h_tot_ly_v_reco_x_truemue->GetBinContent(bin_no)/h_reco_x_counter_recmu->GetBinContent(bin_no)));
                                                    float err_tru = sqrt(h_tot_ly_v_reco_x_truemue->GetBinContent(bin_no))/h_reco_x_counter_recmu->GetBinContent(bin_no);
                                                    h_tot_ly_truemue->SetBinError(bin_no, err_tru);
                                                }
                                         
                                            }//end of filling recomu LY plot loop
                                         
                                    
                                    }
                                }//end of flash matching condition
                                
                                GoodMuonExists=true;
                                startposx_1.push_back(stposx);
                                startposy_1.push_back(stposy);
                                startposz_1.push_back(stposz);
                                endposx_1.push_back(enposx);
                                endposy_1.push_back(enposy);
                                endposz_1.push_back(enposz);
                                
                                goodmuonindex.push_back(itrack);
                                
                                file << run << " "<< sub << " "<< evt << " " << stposx << " " << stposy << " " << stposz << " " << enposx << " " << enposy << " " << enposz << " " << trklen << " " << trk_score << " " << trk_pid << " " << n_tracks << " " << n_showers << " " << endl;
                                
                                
                            }
                        }
                        
                        /*else if(stposy < 110 && enposy > -110)
                         {
                         // fill different third tier of histos
                         
                         // next cut go here (e.g. length/shower/
                         // when all satsified, set GoodMichelExists=true;
                         
                         
                         }*/
                        
                        
                        
                    }
                    
                    
                }
            }
        }
        
        ///////second loop to find Michel candidates:
        //cut out unphysical regions in our position variables
        
        
        
        
        if(GoodMuonExists)
        {
            //std::cout << "size of index vector " << goodmuonindex.size() << std::endl;
            
            for(unsigned int itrack=0;
                itrack<(*trk_sce_start_x_v).size(); itrack++)
            {
                
                if(goodmuonindex[0]==itrack)
                    continue;
                
                //helper variables/variable definitions
                float stposx,stposy,stposz,enposx,enposy,enposz,trklen,mu_enposx,mu_enposy,mu_enposz,mu_me_trklen,mu_mest_trklen,mu_meen_trklen,trk_score,trk_pid,trk_E_y,trk_E_u,trk_E_v, mcendx, mcendy, mcendz, mcpdg, me_mcE;
                stposx=(*trk_sce_start_x_v)[itrack];
                stposy=(*trk_sce_start_y_v)[itrack];
                stposz=(*trk_sce_start_z_v)[itrack];
                enposx=(*trk_sce_end_x_v)[itrack];
                enposy=(*trk_sce_end_y_v)[itrack];
                enposz=(*trk_sce_end_z_v)[itrack];
                trk_score=(*trk_score_v)[itrack];
                trk_pid=(*trk_llr_pid_score_v)[itrack];
                trk_E_y=(*trk_calo_energy_y_v)[itrack];
                trk_E_u=(*trk_calo_energy_u_v)[itrack];
                trk_E_v=(*trk_calo_energy_v_v)[itrack];
                mu_enposx=endposx_1[goodmuonindex[0]];
                mu_enposy=endposy_1[goodmuonindex[0]];
                mu_enposz=endposz_1[goodmuonindex[0]];
                mcendx = (*mc_endx)[itrack];
                mcendy =(*mc_endy)[itrack];
                mcendz =(*mc_endz)[itrack];
                mcpdg = (*mc_pdg)[itrack];
                me_mcE = endmuonmichel;

            
                me_mc_energy -> Fill(me_mcE*1000.);
                //trk_score=(*trk_score_v)[itrack];
                //trk_pid=(*trk_llr_pid_score_v)[itrack];
                
                
                
                // use only sensible vectors
                if(stposx > 0 && stposy > -115 && stposz >0 && enposx >0 && enposy > -115 && enposz > 0)
                {
                    trklen = calc_len(stposx,stposy,stposz,enposx,enposy,enposz);
                    
                    //here can fill vectors or histograms (with no cuts/only non crazy values).
                    /* startposx -> Fill(stposx);
                     startposy -> Fill(stposy);
                     startposz -> Fill(stposz);
                     endposx -> Fill(enposx);
                     endposy -> Fill(enposy);
                     endposz -> Fill(enposz);
                     track_length -> Fill(trklen);
                     */
                    
                    // add fiducial volume cuts
                    if(stposx > 10 && stposx < 246 && enposx > 10 && enposx < 246
                       && stposz >10 && stposz < 1026 && enposz > 10 && enposz < 1026 && stposy<105 && stposy>-105 && enposy > -105 && enposy < 105)
                    {
                        // fill second tier of histograms
                        startposx_FV2->Fill(stposx);
                        startposy_FV2->Fill(stposy);
                        startposz_FV2->Fill(stposz);
                        endposx_FV2->Fill(enposx);
                        endposy_FV2->Fill(enposy);
                        endposz_FV2->Fill(enposz);
                        track_length_FV2->Fill(trklen);
                        
                        // here calculate distance between stpos/enpos and
                        // goldmuonend points and plot histogram of smmaller of the two
                        mu_mest_trklen = calc_len(mu_enposx,mu_enposy,mu_enposz,stposx,stposy,stposz);
                        
                        /* if (run == 6959 && sub == 201 && evt == 10056)
                         {
                         std::cout << "GoodMuon " << GoodMuonExists << " " << goodmuonindex[0] << " size: " << endposx_1.size() << std::endl;
                         std::cout << " start points: " << mu_enposx << " " << stposx << " y " << mu_enposy << " " << stposy << " z " << mu_enposz << " " << stposz << " " << mu_mest_trklen << std::endl;
                         } */
                        
                        mu_meen_trklen = calc_len(mu_enposx,mu_enposy,mu_enposz,enposx,enposy,enposz);
                        /* if (run == 6959 && sub == 201 && evt == 10056)
                         {
                         std::cout << " start points: " << mu_enposx << " " << enposx << " y " << mu_enposy << " " << enposy << " z " << mu_enposz << " " << enposz << " " << mu_meen_trklen << std::endl;
                         } */
                        
                        
                        
                        //compare the two distance, shorter distance is the muon end point & ME start point, swap the two if mixed
                        
                        float shorter_trklen = (mu_meen_trklen < mu_mest_trklen) ? mu_meen_trklen : mu_mest_trklen;
                        
                        /*if (run == 6959 && sub == 201 && evt == 10056)
                         {
                         cout << " Gap distance: " << shorter_trklen << endl;
                         }*/
                        
                        
                        MuonME_gap_zoom->Fill(shorter_trklen);
                        
                        //now set a maximum gap limit to select our ME
                        if(shorter_trklen < 10)
                        {
                            track_length_ME_no_n_showers->Fill(trklen);
                            
                            if(trk_score < 0.5)
                            {
                                track_length_ME->Fill(trklen);
                                track_score_ME->Fill(trk_score);
                                track_PID_ME->Fill(trk_pid);
                                Track_Energy_ME->Fill(trk_E_y);
                                Track_Energy_ME_u->Fill(trk_E_u);
                                Track_Energy_ME_v->Fill(trk_E_v);
                                startposx_ME->Fill(stposx);
                                startposy_ME->Fill(stposy);
                                startposz_ME->Fill(stposz);
                                MuonME_gap->Fill(shorter_trklen);
                                
                                reco_me_energy.push_back(trk_E_y);
                                me_mc_E_v_reco_E -> Fill(me_mcE*1000., trk_E_y);
                                if(mcpdg != 0 && mcpdg < 100)
                                {
                                    ME_mc_end_x->Fill(mcendx);
                                    ME_mc_end_y->Fill(mcendy);
                                    ME_mc_end_z->Fill(mcendz);
                                    
                                }
                            
                                file1 << run << " "<< sub << " "<< evt << " " << stposx << " " << stposy << " " << stposz << " " << enposx << " " << enposy << " " << enposz << " " << trklen << " " << trk_score << " " << trk_pid << " " << shorter_trklen << " " << n_tracks << " " << n_showers << " " << " " << trk_E_y << endl;
                                
                                if(trk_E_y > 50)
                                {
                                    file2 << run << " " << sub << " "<< evt << " " << stposx << " " << stposy << " " << stposz << " " << enposx << " " << enposy << " " << enposz << " " << trklen << " " << trk_score << " " << trk_pid << " " << trk_E_y << " " << endl;
                                }
                                
                                //calculate the distance between flash and track (michel)
                                float me_distance_flash_trk = calc_flash_evt_dist(flash_ycenter, flash_zcenter, stposy, enposy, stposz, enposz);
                                h_me_dist_flsh_trk -> Fill(me_distance_flash_trk);
                                
                                //flash/trk matching selection
                                if (me_distance_flash_trk < 160. && calc_devi_flash_trk(flash_zcenter, flash_zwidth, stposz, enposz) <= 1. && me_mcE > 0)
                                {
                                    //create a variable with value 0 for the total photoelectrons (PEs)
                                    
                                    float tot_pe_per_trk = 0;
                                    //looping over all the pmt for the total PE per event
                                    for (unsigned int PMT_no = 0; PMT_no != 32; ++PMT_no)
                                    {
                                        tot_pe_per_trk += flash_pe_v->at(PMT_no)*(20./gain_ampl_v->at(PMT_no));
                                    }
                                    //now calculate the total LY for 4 cases
                                    //total LY for true Michel electron energy
                                    if (tot_pe_per_trk >0 && trk_E_y > 0)
                                    {
                                        float tot_ly_per_trk = tot_pe_per_trk/trk_E_y;
                                        
                                        float tot_ly_per_trk_mcE = tot_pe_per_trk/(me_mcE*1000.);
                                        
                                        float trk_mid_x = (stposx+enposx)/2.;
                                        
                                        //first fill rough ly plot (different no. of tracks in each bin
                                        h_tot_ly_v_reco_x_recME -> Fill(trk_mid_x, tot_ly_per_trk);
                                        h_tot_ly_v_reco_x_trueMEe -> Fill(trk_mid_x, tot_ly_per_trk_mcE);
                                        
                                        /*
                                         for (unsigned int bin_no = 1; bin_no != 100; ++bin_no)
                                         {
                                         cout << h_tot_ly_v_reco_x_recME->GetBinContent(bin_no) << endl;
                                         }
                                         */
                                        //int bin = h_tot_ly_v_reco_x_recME->FindBin(trk_mid_x);
                                        //h_tot_ly_v_reco_x_recME->AddBinContent(bin, tot_ly_per_trk);
                                        // Next fill a reco_x midpoint counter to count no of evt in each bin
                        
                                            h_reco_x_counter_recME -> Fill(trk_mid_x);
                                            
                                    
                                            //next get bin content of the two histo, divide the first by the second to get a normalized and final light yield plot
                                            for (unsigned int bin_no = 1; bin_no != 100; ++bin_no)
                                            {
                                                if(h_reco_x_counter_recME->GetBinContent(bin_no) != 0)
                                                {
                                                    if(h_tot_ly_v_reco_x_recME->GetBinContent(bin_no) < 700)
                                                    {
                                                        h_tot_ly_recME->SetBinContent(bin_no,(h_tot_ly_v_reco_x_recME->GetBinContent(bin_no)/h_reco_x_counter_recME->GetBinContent(bin_no)));
                                                        float err = sqrt(h_tot_ly_v_reco_x_recME->GetBinContent(bin_no))/h_reco_x_counter_recME->GetBinContent(bin_no);
                                                        h_tot_ly_recME->SetBinError(bin_no, err);
                                                    }
                                                    if(h_tot_ly_v_reco_x_trueMEe->GetBinContent(bin_no) < 100000)
                                                    {
                                                        h_tot_ly_trueMEe->SetBinContent(bin_no,(h_tot_ly_v_reco_x_trueMEe->GetBinContent(bin_no)/h_reco_x_counter_recME->GetBinContent(bin_no)));
                                                        float err_tru = sqrt(h_tot_ly_v_reco_x_trueMEe->GetBinContent(bin_no))/h_reco_x_counter_recME->GetBinContent(bin_no);
                                                        h_tot_ly_trueMEe->SetBinError(bin_no, err_tru);
                                                    }
                                                }
                                            }//end of filling recoME LY plot loop
                                    
                                    }
                                }//end of flash matching condition
                                
                            } //michel electron selection
                        }//end of max gap distance condition
                        
                    }
                }
                //define another helper variable for the selected reco ME energy
                float reco_meE;
                reco_meE = reco_me_energy[itrack];
                
                if (mcpdg != 0 && mcpdg < 100) {
                    
                }
                
                
            }  // end of for loop on tracks
        } // end of "if GoodMuonExists is true"
        
        
        // if (Cut(ientry) < 0) continue;
    }
    mc_E_vs_reco_E->Fit("pol1","","",200,600);
    
    
    //closing particle information file
    file.close();
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    cout<< "Selected muon info file for evt displayed is filled." <<endl;
    cout<< "Selected Michel electron info file filled." <<endl;
    cout<< "Selected Michel electrons info w/ E>50MeV file filled." <<endl;
    cout<< ".txt file filled with info for events that end very shallow in the detector." <<endl;
    cout<< ".txt file filled with info for events that start deep in the detector." <<endl;
    
    //Visualization Drawing histo part
    //defining new canvases
    TCanvas *canvas1 = new TCanvas("canvas1");
    TCanvas *canvas2 = new TCanvas("canvas2");
    TCanvas *canvas3 = new TCanvas("canvas3");
    TCanvas *canvas4 = new TCanvas("canvas4");
    TCanvas *canvas5 = new TCanvas("canvas5");
    TCanvas *canvas6 = new TCanvas("canvas6");
    TCanvas *canvas7 = new TCanvas("canvas7");
    TCanvas *canvas8 = new TCanvas("canvas8");
    TCanvas *canvas9 = new TCanvas("canvas9");
    TCanvas *canvas10 = new TCanvas("canvas10");
    TCanvas *canvas11 = new TCanvas("canvas11");
    TCanvas *canvas12 = new TCanvas("canvas12");
    TCanvas *canvas13 = new TCanvas("canvas13");
    TCanvas *canvas14 = new TCanvas("canvas14");
    TCanvas *canvas15 = new TCanvas("canvas15");
    TCanvas *canvas16 = new TCanvas("canvas16");
    TCanvas *canvas17 = new TCanvas("canvas17");
    TCanvas *canvas18 = new TCanvas("canvas18");
    TCanvas *canvas19 = new TCanvas("canvas19");
    TCanvas *canvas20 = new TCanvas("canvas20");
    TCanvas *canvas21 = new TCanvas("canvas21");
    TCanvas *canvas22 = new TCanvas("canvas22");
    TCanvas *canvas23 = new TCanvas("canvas23");
    TCanvas *canvas24 = new TCanvas("canvas24");
    TCanvas *canvas25 = new TCanvas("canvas25");
    TCanvas *canvas26 = new TCanvas("canvas26");
    TCanvas *canvas27 = new TCanvas("canvas27");
    TCanvas *canvas28 = new TCanvas("canvas28");
    TCanvas *canvas29 = new TCanvas("canvas29");
    TCanvas *canvas30 = new TCanvas("canvas30");
    TCanvas *canvas31 = new TCanvas("canvas31");
    TCanvas *canvas32 = new TCanvas("canvas32");
    TCanvas *canvas33 = new TCanvas("canvas33");
    TCanvas *canvas34 = new TCanvas("canvas34");
    TCanvas *canvas35 = new TCanvas("canvas35");
    TCanvas *canvas36 = new TCanvas("canvas36");
    TCanvas *canvas37 = new TCanvas("canvas37");
    TCanvas *canvas38 = new TCanvas("canvas38");
    TCanvas *canvas39 = new TCanvas("canvas39");
    TCanvas *canvas40 = new TCanvas("canvas40");
    TCanvas *canvas41 = new TCanvas("canvas41");
    TCanvas *canvas42 = new TCanvas("canvas42");
    TCanvas *canvas43 = new TCanvas("canvas43");
    TCanvas *canvas44 = new TCanvas("canvas44");
    TCanvas *canvas45 = new TCanvas("canvas45");
    TCanvas *canvas46 = new TCanvas("canvas46");
    TCanvas *canvas47 = new TCanvas("canvas47");
    TCanvas *canvas48 = new TCanvas("canvas48");
    TCanvas *canvas49 = new TCanvas("canvas49");
    TCanvas *canvas50 = new TCanvas("canvas50");
    TCanvas *canvas51 = new TCanvas("canvas51");
    TCanvas *canvas52 = new TCanvas("canvas52");
    TCanvas *canvas53 = new TCanvas("canvas53");
    
    //actual drawing
    canvas47->cd();
    h_me_dist_flsh_trk->Draw();
    h_me_dist_flsh_trk->GetYaxis()->SetTitle("Number of Events");
    h_me_dist_flsh_trk->GetXaxis()->SetTitle("Distance (cm)");
    gStyle->SetOptStat(10);
    canvas47->Print("plots.ps(");
    myFile->WriteObject(canvas47, "Distance between Flash and ME Event Track");
    
    canvas48->cd();
    h_mu_dist_flsh_trk->Draw();
    h_mu_dist_flsh_trk->GetYaxis()->SetTitle("Number of Events");
    h_mu_dist_flsh_trk->GetXaxis()->SetTitle("Distance (cm)");
    gStyle->SetOptStat(10);
    canvas48->Print("plots.ps");
    myFile->WriteObject(canvas48, "Distance between Flash and Muon Event Track");

    canvas49->cd();
    h_tot_ly_v_reco_x_recME->Draw("E1");
    h_tot_ly_v_reco_x_recME->GetYaxis()->SetTitle("Number of PE/MeV");
    h_tot_ly_v_reco_x_recME->GetXaxis()->SetTitle("Track mid-x position (cm)");
    gStyle->SetOptStat(10);
    canvas49->Print("plots.ps");
    myFile->WriteObject(canvas49, "Total Light Yield of Reconstructed Michel Electrons (Reco E)");
    
    canvas50->cd();
    h_tot_ly_v_reco_x_trueMEe->Draw("E1");
    h_tot_ly_v_reco_x_trueMEe->GetYaxis()->SetTitle("Number of PE/MeV");
    h_tot_ly_v_reco_x_trueMEe->GetXaxis()->SetTitle("Track mid-x position (cm)");
    gStyle->SetOptStat(10);
    canvas50->Print("plots.ps");
    myFile->WriteObject(canvas50, "Total Light Yield of Reconstructed Michel Electrons (True E)");
    
    canvas51->cd();
    h_tot_ly_v_reco_x_recmu->Draw("E1");
    h_tot_ly_v_reco_x_recmu->GetYaxis()->SetTitle("Number of PE/MeV");
    h_tot_ly_v_reco_x_recmu->GetXaxis()->SetTitle("Track mid-x position (cm)");
    gStyle->SetOptStat(10);
    canvas51->Print("plots.ps");
    myFile->WriteObject(canvas51, "Total Light Yield of Reconstructed Muons (Reco E)");
    
    canvas52->cd();
    h_tot_ly_v_reco_x_truemue->Draw("E1");
    h_tot_ly_v_reco_x_truemue->GetYaxis()->SetTitle("Number of PE/MeV");
    h_tot_ly_v_reco_x_truemue->GetXaxis()->SetTitle("Track mid-x position (cm)");
    gStyle->SetOptStat(10);
    canvas52->Print("plots.ps");
    myFile->WriteObject(canvas52, "Total Light Yield of Reconstructed Michel Electrons (True E)");

    /*
    canvas49->cd();
    h_tot_ly_recME->Draw();
    h_tot_ly_recME->GetYaxis()->SetTitle("Number of PE/MeV");
    h_tot_ly_recME->GetXaxis()->SetTitle("Track mid-x position");
    gStyle->SetOptStat(10);
    canvas49->Print("plots.ps");
    myFile->WriteObject(canvas49, "Total Light Yield of Reconstructed Michel Electrons");
    */
     
    
    canvas29->cd();
    mc_end_x->Draw();
    mc_end_x->GetXaxis()->SetTitle("x-Position (cm)");
    mc_end_x->GetYaxis()->SetTitle("Number of Events");
    mc_end_x_reco->Draw("same");
    mc_end_x->GetYaxis()->SetRangeUser(-10, mc_end_x->GetMaximum()*1.1);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg20 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg20->SetBorderSize(0);
    leg20->SetFillStyle(0);
    leg20->AddEntry(mc_end_x, "Original Monte Carlo true end position (x)", "lep");
    leg20->AddEntry(mc_end_x_reco,"Monte Carlo true end position for reconstructed tracks (x)", "lep");
    leg20->SetTextFont(62);
    leg20->Draw();
    canvas29->Print("plots.ps");
    myFile->WriteObject(canvas29, "Monte Carlo Truth End Position (x)");
    

    
    canvas30->cd();
    mc_end_y->Draw();
    mc_end_y->GetXaxis()->SetTitle("y-Position (cm)");
    mc_end_y->GetYaxis()->SetTitle("Number of Events");
    mc_end_y_reco->Draw("same");
    mc_end_y->GetYaxis()->SetRangeUser(-10, mc_end_y->GetMaximum()*1.1);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg21 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg21->SetBorderSize(0);
    leg21->SetFillStyle(0);
    leg21->AddEntry(mc_end_y, "Original Monte Carlo true end position (y)", "lep");
    leg21->AddEntry(mc_end_y_reco,"Monte Carlo true end position for reconstructed tracks (y)", "lep");
    leg21->SetTextFont(62);
    leg21->Draw();
    canvas30->Print("plots.pdf)","Monte Carlo true end position (y)");
    canvas30->Print("plots.ps");
    myFile->WriteObject(canvas30, "Monte Carlo Truth End Position (y)");
    
    
    canvas31->cd();
    mc_end_z->Draw();
    mc_end_z->GetXaxis()->SetTitle("z-Position (cm)");
    mc_end_z->GetYaxis()->SetTitle("Number of Events");
    mc_end_z_reco->Draw("same");
    mc_end_z->GetYaxis()->SetRangeUser(-10, mc_end_z->GetMaximum()*1.1);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg22 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg22->SetBorderSize(0);
    leg22->SetFillStyle(0);
    leg22->AddEntry(mc_end_z, "Original Monte Carlo true end position (z)", "lep");
    leg22->AddEntry(mc_end_z_reco,"Monte Carlo true end position for reconstructed tracks (z)", "lep");
    leg22->SetTextFont(62);
    leg22->Draw();
    canvas31->Print("plots.ps");
    myFile->WriteObject(canvas31, "Monte Carlo Truth End Position (z)");
    
    canvas32->cd();
    mc_v_x->Draw();
    mc_v_x->GetXaxis()->SetTitle("x-Position (cm)");
    mc_v_x->GetYaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas32->Print("plots.ps");
    myFile->WriteObject(canvas32, "Monte Carlo True Start Position (x)");
    
    canvas33->cd();
    mc_v_y->Draw();
    mc_v_y->GetXaxis()->SetTitle("y-Position (cm)");
    mc_v_y->GetYaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas33->Print("plots.ps");
    myFile->WriteObject(canvas33, "Monte Carlo True Start Position (y)");
    
    canvas34->cd();
    mc_v_z->Draw();
    mc_v_z->GetXaxis()->SetTitle("z-Position (cm)");
    mc_v_z->GetYaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas34->Print("plots.ps");
    myFile->WriteObject(canvas34, "Monte Carlo True Start Position (z)");
    
    canvas35->cd();
    mc_v_x->Draw();
    mc_v_x_st_deep->Draw("same");
    mc_v_x_end_shallow->Draw("same");
    mc_v_x->GetXaxis()->SetTitle("x-Position (cm)");
    mc_v_x->GetYaxis()->SetTitle("Number of Events");
    mc_v_x->GetYaxis()->SetRangeUser(0., 1000.);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg23 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg23->SetBorderSize(0);
    leg23->SetFillStyle(0);
    leg23->AddEntry(mc_v_x, "Original Monte Carlo true start position (x)", "lep");
    leg23->AddEntry(mc_v_x_st_deep,"Monte Carlo true start position for events that start deep in the detector (x)", "lep");
    leg23->AddEntry(mc_v_x_end_shallow,"Monte Carlo true start position for events that end shallow in the detector (x)", "lep");
    leg23->SetTextFont(62);
    leg23->Draw();
    canvas35->Print("plots.ps");
    myFile->WriteObject(canvas35, "Monte Carlo True Start Position for Events in Question (x)");

    canvas36->cd();
    mc_v_z->Draw();
    mc_v_z_st_deep->Draw("same");
    mc_v_z_end_shallow->Draw("same");
    mc_v_z->GetXaxis()->SetTitle("z-Position (cm)");
    mc_v_z->GetYaxis()->SetTitle("Number of Events");
    mc_v_z->GetYaxis()->SetRangeUser(0., 1000.);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg24 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg24->SetBorderSize(0);
    leg24->SetFillStyle(0);
    leg24->AddEntry(mc_v_z, "Original Monte Carlo true start position (z)", "lep");
    leg24->AddEntry(mc_v_z_st_deep,"Monte Carlo true start position for events that start deep in the detector (z)", "lep");
    leg24->AddEntry(mc_v_z_end_shallow,"Monte Carlo true start position for events that end shallow in the detector (z)", "lep");
    leg24->SetTextFont(62);
    leg24->Draw();
    canvas36->Print("plots.ps");
    myFile->WriteObject(canvas36, "Monte Carlo True Start Position for Events in Question (z)");
    
    canvas37->cd();
    mc_end_y->Draw();
    mc_end_y->GetXaxis()->SetTitle("y-Position (cm)");
    mc_end_y->GetYaxis()->SetTitle("Number of Events");
    mc_end_y_st_deep->Draw("same");
    mc_end_y_end_shallow->Draw("same");
    mc_end_y->GetYaxis()->SetRangeUser(-10, mc_end_y->GetMaximum()*1.1);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg25 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg25->SetBorderSize(0);
    leg25->SetFillStyle(0);
    leg25->AddEntry(mc_end_y, "Original Monte Carlo true end position (y)", "lep");
    leg25->AddEntry(mc_end_y_st_deep,"Monte Carlo true end position for events that start deep in the detector (y)", "lep");
    leg25->AddEntry(mc_end_y_end_shallow,"Monte Carlo true end position for events that end shallow in the detector (y)", "lep");
    leg25->SetTextFont(62);
    leg25->Draw();
    canvas37->Print("plots.ps");
    myFile->WriteObject(canvas37, "Monte Carlo True End Position for Events in Question (y)");
    
    canvas38->cd();
    mc_energy->Draw();
    mc_energy->GetXaxis()->SetTitle("Energy (MeV)");
    mc_energy->GetYaxis()->SetTitle("Number of Events");
    mc_energy_st_deep->Draw("same");
    mc_energy_end_shallow->Draw("same");
    mc_energy->GetYaxis()->SetRangeUser(0., 500.);
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg26 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg26->SetBorderSize(0);
    leg26->SetFillStyle(0);
    leg26->AddEntry(mc_energy, "Original Monte Carlo true energy", "lep");
    leg26->AddEntry(mc_energy_st_deep,"Monte Carlo true energy for events that start deep in the detector", "lep");
    leg26->AddEntry(mc_energy_end_shallow,"Monte Carlo true energy for events that end shallow in the detector", "lep");
    leg26->SetTextFont(62);
    leg26->Draw();
    canvas38->Print("plots.ps");
    myFile->WriteObject(canvas38, "Monte Carlo True Energy for Events in Question");
    
    canvas39->cd();
    mc_endy_vs_y_end_reco_st_deep -> Draw("colz");
    mc_endy_vs_y_end_reco_st_deep->GetYaxis()->SetTitle("MC y-Position (cm)");
    mc_endy_vs_y_end_reco_st_deep->GetXaxis()->SetTitle("Reco y-Position (cm)");
    mc_endy_vs_y_end_reco_st_deep->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas39->Print("plots.ps");
    myFile->WriteObject(canvas39, "MC End (y) vs Reconstructed End Position (y) Comparison Plot for Events that Start Deep in the Detector");
    
    canvas40->cd();
    mc_endy_vs_y_end_reco_end_shallow->Draw("colz");
    mc_endy_vs_y_end_reco_end_shallow->GetYaxis()->SetTitle("MC y-Position (cm)");
    mc_endy_vs_y_end_reco_end_shallow->GetXaxis()->SetTitle("Reco y-Position (cm)");
    mc_endy_vs_y_end_reco_end_shallow->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas40->Print("plots.ps");
    myFile->WriteObject(canvas40, "MC End (y) vs Reconstructed End Position (y) Comparison Plot for Events that End Shallow in the Detector");
    
    canvas41->cd();
    mc_E_vs_reco_E->Draw("colz");
    mc_E_vs_reco_E->GetYaxis()->SetTitle("Reconstructed Track Energy (MeV)");
    mc_E_vs_reco_E->GetXaxis()->SetTitle("MC Energy (MeV)");
    mc_E_vs_reco_E->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas41->Print("plots.ps");
    myFile->WriteObject(canvas41, "Muon Monte Carlo True Energy vs. Track Calorimetric Energy for Reconstructed Tracks");
    
    canvas42->cd();
    orig_me_mc_energy->Draw();
    me_mc_energy->Draw("same");
    orig_me_mc_energy->GetXaxis()->SetTitle("Energy (MeV)");
    orig_me_mc_energy->GetYaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    //legend
    TLegend *leg27= new TLegend(0.55, 0.8, 0.9, 0.9);
    leg27->SetBorderSize(0);
    leg27->SetFillStyle(0);
    leg27->AddEntry(orig_me_mc_energy, "Original Monte Carlo Michel electron energy", "lep");
    leg27->AddEntry(me_mc_energy,"Monte Carlo Michel electron energy after selection", "lep");
    leg27->SetTextFont(62);
    leg27->Draw();
    canvas42->Print("plots.ps");
    myFile->WriteObject(canvas42, "Monte Carlo True Energy of Michel Electrons");
    
    canvas43->cd();
    me_mc_E_v_reco_E->Draw("colz");
    me_mc_E_v_reco_E->GetYaxis()->SetTitle("Reconstructed Track Energy (MeV)");
    me_mc_E_v_reco_E->GetXaxis()->SetTitle("MC Energy (MeV)");
    me_mc_E_v_reco_E->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas43->Print("plots.ps");
    myFile->WriteObject(canvas43, "Michel Electron Monte Carlo True Energy vs. Track Calorimetric Energy for Reconstructed Tracks");
     
    canvas44->cd();
    mc_pdg_code->Draw();
    mc_pdg_code->GetYaxis()->SetTitle("Number of Events");
    mc_pdg_code->GetXaxis()->SetTitle("PDG Code (AU)");
    gStyle->SetOptStat(10);
    canvas44->Print("plots.ps");
    myFile->WriteObject(canvas44, "Monte Carlo PDG Code");
    
    canvas45->cd();
    end_y_v_z_st_deep->Draw("colz");
    end_y_v_z_st_deep->GetYaxis()->SetTitle("Reconstructed y End Position (cm)");
    end_y_v_z_st_deep->GetXaxis()->SetTitle("Reconstructed z End Position (cm)");
    end_y_v_z_st_deep->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas45->Print("plots.ps");
    myFile->WriteObject(canvas45, "The Reconstructed yz Plane End Position for Events Starting Deep in the Detector");
                                       
    canvas46->cd();
    end_y_v_z_end_shallow->Draw("colz");
    end_y_v_z_end_shallow->GetYaxis()->SetTitle("Reconstructed y End Position (cm)");
    end_y_v_z_end_shallow->GetXaxis()->SetTitle("Reconstructed z End Position (cm)");
    end_y_v_z_end_shallow->GetZaxis()->SetTitle("Number of Events");
    gStyle->SetOptStat(10);
    canvas46->Print("plots.ps");
    myFile->WriteObject(canvas46, "The Reconstructed yz Plane End Position for Events Ending Shallow in the Detector");
    
    
    //fChain->Get(mc_endx)->GetXaxis()->SetTitle("x-position (cm)");
    //mc_endx->GetYaxis()->SetRangeUser(-10, startposx->GetMaximum()*1.1);
    //mc_endx->GetYaxis()->SetTitle("Number of Events");
    //orig_mc_endx->Draw("HISTO");
    //canvas32->Draw();
    
    /*
    canvas29->cd();
    mc_end_x->GetXaxis()->SetTitle("x-position (cm)");
    mc_end_x->GetYaxis()->SetTitle("Number of Events");
    mc_end_x->Draw();
    ME_mc_end_x->Draw("same");
    //legend
    TLegend *leg20 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg20->SetBorderSize(0);
    leg20->SetFillStyle(0);
    leg20->AddEntry(mc_end_x, "x-position of the Monte Carlo muon sample", "lep");
    leg20->AddEntry(ME_mc_end_x,"x-position of the selected Monte Carlo Michel electron sample", "lep");
    leg20->SetTextFont(62);
    leg20->Draw();
    
    
    canvas30->cd();
    mc_end_y->GetXaxis()->SetTitle("y-position (cm)");
    mc_end_y->GetYaxis()->SetTitle("Number of Events");
    mc_end_y->Draw();
    ME_mc_end_y->Draw("same");
    //legend
    TLegend *leg21 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg21->SetBorderSize(0);
    leg21->SetFillStyle(0);
    leg21->AddEntry(mc_end_y, "y-position of the Monte Carlo muon sample", "lep");
    leg21->AddEntry(ME_mc_end_y,"y-position of the selected Monte Carlo Michel electron sample", "lep");
    leg21->SetTextFont(62);
    leg21->Draw();
    
    canvas31->cd();
    mc_end_z->GetXaxis()->SetTitle("z-position (cm)");
    mc_end_z->GetYaxis()->SetTitle("Number of Events");
    mc_end_z->Draw();
    ME_mc_end_z->Draw("same");
    //legend
    TLegend *leg22 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg22->SetBorderSize(0);
    leg22->SetFillStyle(0);
    leg22->AddEntry(mc_end_z, "z-position of the Monte Carlo muon sample", "lep");
    leg22->AddEntry(ME_mc_end_z,"z-position of the selected Monte Carlo Michel electron sample", "lep");
    leg22->SetTextFont(62);
    leg22->Draw();
    */
    
    
    canvas1->cd();
    startposx->GetXaxis()->SetTitle("x-position (cm)");
    startposx->GetYaxis()->SetRangeUser(-10, startposx->GetMaximum()*1.1);
    startposx->GetYaxis()->SetTitle("Number of Events");
    startposx->Draw();
    startposx_FV->Draw("same");
    startposx_stp_trk->Draw("same");
    //legend
    TLegend *leg1 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(startposx, "Original", "lep");
    leg1->AddEntry(startposx_FV,"After fiducial volume", "lep");
    leg1->AddEntry(startposx_stp_trk,"F.V + stopping trks selection", "lep");
    leg1->SetTextFont(62);
    leg1->Draw();
    gStyle->SetOptStat(10);
    canvas1->Print("plots.ps");
    myFile->WriteObject(canvas1, "Muon Selection Start Position of the Reconstructed Event Tracks (x)");
    
    //TLegend *leg1 = InitializeLegend(canvas1);
    //DrawSingleCanvas(canvas1,startposx,"x-position (cm)", "x position, no cuts" ,leg1);
    //DrawSecondHistOnCanvas(canvas1,startposx_FV,"fiducial volume",kRed, leg1);
    
    canvas2->cd();
    startposy->GetXaxis()->SetTitle("y-position (cm)");
    startposy->GetYaxis()->SetRangeUser(-10, startposy->GetMaximum()*1.1);
    startposy->GetYaxis()->SetTitle("Number of Events");
    startposy->Draw();
    startposy_FV->Draw("same");
    startposy_stp_trk->Draw("same");
    //legend
    TLegend *leg2 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(startposy, "Original", "lep");
    leg2->AddEntry(startposy_FV,"After fiducial volume", "lep");
    leg2->AddEntry(startposy_stp_trk,"F.V + stopping trks selection", "lep");
    leg2->SetTextFont(62);
    leg2->Draw();
    gStyle->SetOptStat(10);
    canvas2->Print("plots.ps");
    myFile->WriteObject(canvas2, "Muon Selection Start Position of the Reconstructed Event Tracks (y)");
    
    canvas3->cd();
    startposz->GetXaxis()->SetTitle("z-position (cm)");
    startposz->GetYaxis()->SetRangeUser(-10, startposz->GetMaximum()*1.1);
    startposz->GetYaxis()->SetTitle("Number of Events");
    startposz->Draw();
    startposz_FV->Draw("same");
    startposz_stp_trk->Draw("same");
    //legend
    TLegend *leg3 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(startposz, "Original", "lep");
    leg3->AddEntry(startposz_FV,"After fiducial volume", "lep");
    leg3->AddEntry(startposz_stp_trk,"F.V + stopping trks selection", "lep");
    leg3->SetTextFont(62);
    leg3->Draw();
    gStyle->SetOptStat(10);
    canvas3->Print("plots.ps");
    myFile->WriteObject(canvas3, "Muon Selection Start Position of the Reconstructed Event Tracks (z)");
    
    canvas4->cd();
    endposx->GetXaxis()->SetTitle("x-position (cm)");
    endposx->GetYaxis()->SetRangeUser(-10, endposx->GetMaximum()*1.1);
    endposx->GetYaxis()->SetTitle("Number of Events");
    endposx->Draw();
    endposx_FV->Draw("same");
    endposx_stp_trk->Draw("same");
    //legend
    TLegend *leg4 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->AddEntry(endposx, "Original", "lep");
    leg4->AddEntry(endposx_FV,"After fiducial volume", "lep");
    leg4->AddEntry(endposx_stp_trk,"F.V + stopping trks selection", "lep");
    leg4->SetTextFont(62);
    leg4->Draw();
    gStyle->SetOptStat(10);
    canvas4->Print("plots.ps");
    myFile->WriteObject(canvas4, "Muon Selection End Position of the Reconstructed Event Tracks (x)");
    
    canvas5->cd();
    endposy->GetXaxis()->SetTitle("y-position (cm)");
    endposy->GetYaxis()->SetTitle("Number of Events");
    endposy->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposy->Draw();
    endposy_FV->Draw("same");
    endposy_stp_trk->Draw("same");
    //legend
    TLegend *leg5 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);
    leg5->AddEntry(endposy, "Original", "lep");
    leg5->AddEntry(endposy_FV,"After fiducial volume", "lep");
    leg5->AddEntry(endposy_stp_trk,"F.V + stopping trks selection", "lep");
    leg5->SetTextFont(62);
    leg5->Draw();
    gStyle->SetOptStat(10);
    canvas5->Print("plots.ps");
    myFile->WriteObject(canvas5, "Muon Selection End Position of the Reconstructed Event Tracks (y)");
    
    canvas6->cd();
    endposz->GetXaxis()->SetTitle("z-position (cm)");
    endposz->GetYaxis()->SetRangeUser(-10, endposz->GetMaximum()*1.1);
    endposz->GetYaxis()->SetTitle("Number of Events");
    endposz->Draw();
    endposz_FV->Draw("same");
    endposz_stp_trk->Draw("same");
    //legend
    TLegend *leg6 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg6->SetBorderSize(0);
    leg6->SetFillStyle(0);
    leg6->AddEntry(endposz, "Original", "lep");
    leg6->AddEntry(endposz_FV,"After fiducial volume", "lep");
    leg6->AddEntry(endposz_stp_trk,"F.V + stopping trks selection", "lep");
    leg6->SetTextFont(62);
    leg6->Draw();
    gStyle->SetOptStat(10);
    canvas6->Print("plots.ps");
    myFile->WriteObject(canvas6, "Muon Selection End Position of the Reconstructed Event Tracks (z)");
    
    canvas7->cd();
    track_length->GetXaxis()->SetTitle("track length (cm)");
    track_length->GetYaxis()->SetTitle("Number of Events");
    //track_length->GetYaxis()->SetRangeUser(0, endposz->GetMaximum()*1.1);
    track_length->Draw();
    track_length_FV->Draw("same");
    track_length_stp_trk->Draw("same");
    track_length_score_cuts->Draw("same");
    //legend
    TLegend *leg7 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg7->SetBorderSize(0);
    leg7->SetFillStyle(0);
    leg7->AddEntry(track_length, "Original", "lep");
    leg7->AddEntry(track_length_FV,"After fiducial volume", "lep");
    leg7->AddEntry(track_length_stp_trk,"F.V + stopping trks selection", "lep");
    leg7->AddEntry(track_length_score_cuts,"+ score/PID cuts", "lep");
    leg7->SetTextFont(62);
    leg7->Draw();
    gStyle->SetOptStat(10);
    canvas7->Print("plots.ps");
    myFile->WriteObject(canvas7, "Muon Selection Track Length of the Reconstructed Event Tracks");
    
    canvas8->cd();
    track_score->GetXaxis()->SetTitle("track score (AU)");
    track_score->GetYaxis()->SetTitle("Number of Events");
    track_score->Draw();
    track_score_cut->Draw("same");
    //legend
    TLegend *leg8 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg8->SetBorderSize(0);
    leg8->SetFillStyle(0);
    leg8->AddEntry(track_score, "F.V + stopping trks selection", "lep");
    leg8->AddEntry(track_score_cut,"+ score/PID + track length cuts", "lep");
    leg8->SetTextFont(62);
    leg8->Draw();
    gStyle->SetOptStat(10);
    canvas8->Print("plots.ps");
    myFile->WriteObject(canvas8, "Muon Selection Track Score of the Reconstructed Event Tracks");
    
    canvas9->cd();
    track_PID->GetXaxis()->SetTitle("track PID score (AU)");
    track_PID->GetYaxis()->SetTitle("Number of Events");
    track_PID->Draw();
    track_PID_cut->Draw("same");
    //legend
    TLegend *leg9 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->AddEntry(track_PID, "F.V + stopping trks selection", "lep");
    leg9->AddEntry(track_PID_cut,"+ score/PID + track length cuts", "lep");
    leg9->SetTextFont(62);
    leg9->Draw();
    gStyle->SetOptStat(10);
    canvas9->Print("plots.ps");
    myFile->WriteObject(canvas9, "Muon Selection Track PID of the Reconstructed Event Tracks");
    
    canvas10->cd();
    endposx_muon_selection->GetXaxis()->SetTitle("x-position(cm)");
    endposx_muon_selection->GetYaxis()->SetRangeUser(0, endposx_muon_selection->GetMaximum()*1.1);
    endposx_muon_selection->GetYaxis()->SetTitle("Number of Events");
    endposx_muon_selection->Draw();
    gStyle->SetOptStat(10);
    canvas10->Print("plots.ps");
    myFile->WriteObject(canvas10, "End Position of the Selected Muon Candidates (x)");
    
    canvas11->cd();
    endposy_muon_selection->GetXaxis()->SetTitle("y-position(cm)");
    endposy_muon_selection->GetYaxis()->SetRangeUser(0, endposy_muon_selection->GetMaximum()*1.1);
    endposy_muon_selection->GetYaxis()->SetTitle("Number of Events");
    endposy_muon_selection->Draw();
    gStyle->SetOptStat(10);
    canvas11->Print("plots.ps");
    myFile->WriteObject(canvas11, "End Position of the Selected Muon Candidates (y)");
    
    canvas12->cd();
    endposz_muon_selection->GetXaxis()->SetTitle("z-position(cm)");
    endposz_muon_selection->GetYaxis()->SetRangeUser(0, endposz_muon_selection->GetMaximum()*1.1);
    endposz_muon_selection->GetYaxis()->SetTitle("Number of Events");
    endposz_muon_selection->Draw();
    gStyle->SetOptStat(10);
    canvas12->Print("plots.ps");
    myFile->WriteObject(canvas12, "End Position of the Selected Muon Candidates (z)");
    
    canvas13->cd();
    startposx->GetXaxis()->SetTitle("x-position (cm)");
    startposx->GetYaxis()->SetRangeUser(-10, startposx->GetMaximum()*1.1);
    startposx->GetYaxis()->SetTitle("Number of Events");
    startposx->Draw();
    startposx_FV2->Draw("same");
    //legend
    TLegend *leg10 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(startposx, "Original", "lep");
    leg10->AddEntry(startposx_FV2,"F.V + Michel electron selection", "lep");
    leg10->SetTextFont(62);
    leg10->Draw();
    gStyle->SetOptStat(10);
    canvas13->Print("plots.ps");
    myFile->WriteObject(canvas13, "Michel Electron Selection Start Position of the Reconstructed Event Tracks (x)");
    
    
    canvas14->cd();
    startposy->GetXaxis()->SetTitle("y-position (cm)");
    startposy->GetYaxis()->SetRangeUser(-10, startposy->GetMaximum()*1.1);
    startposy->GetYaxis()->SetTitle("Number of Events");
    startposy->Draw();
    startposy_FV2->Draw("same");
    //legend
    TLegend *leg11 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg11->SetBorderSize(0);
    leg11->SetFillStyle(0);
    leg11->AddEntry(startposy, "Original", "lep");
    leg11->AddEntry(startposy_FV2,"F.V + Michel electron selection", "lep");
    leg11->SetTextFont(62);
    leg11->Draw();
    gStyle->SetOptStat(10);
    canvas14->Print("plots.ps");
    myFile->WriteObject(canvas14, "Michel Electron Selection Start Position of the Reconstructed Event Tracks (y)");
    
    canvas15->cd();
    startposz->GetXaxis()->SetTitle("z-position (cm)");
    startposz->GetYaxis()->SetRangeUser(-10, startposz->GetMaximum()*1.1);
    startposz->GetYaxis()->SetTitle("Number of Events");
    startposz->Draw();
    startposz_FV2->Draw("same");
    //legend
    TLegend *leg12 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg12->SetBorderSize(0);
    leg12->SetFillStyle(0);
    leg12->AddEntry(startposz, "Original", "lep");
    leg12->AddEntry(startposz_FV2,"F.V + Michel electron selection", "lep");
    leg12->SetTextFont(62);
    leg12->Draw();
    gStyle->SetOptStat(10);
    canvas15->Print("plots.ps");
    myFile->WriteObject(canvas15, "Michel Electron Selection Start Position of the Reconstructed Event Tracks (z)");
    
    canvas16->cd();
    endposx->GetXaxis()->SetTitle("x-position (cm)");
    endposx->GetYaxis()->SetRangeUser(-10, endposx->GetMaximum()*1.1);
    endposx->GetYaxis()->SetTitle("Number of Events");
    endposx->Draw();
    endposx_FV2->Draw("same");
    //legend
    TLegend *leg13 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg13->SetBorderSize(0);
    leg13->SetFillStyle(0);
    leg13->AddEntry(startposz, "Original", "lep");
    leg13->AddEntry(startposz_FV2,"F.V + Michel electron selection", "lep");
    leg13->SetTextFont(62);
    leg13->Draw();
    gStyle->SetOptStat(10);
    canvas16->Print("plots.ps");
    myFile->WriteObject(canvas16, "Michel Electron Selection End Position of the Reconstructed Event Tracks (x)");
    
    canvas17->cd();
    endposy->GetXaxis()->SetTitle("y-position (cm)");
    endposy->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposy->GetYaxis()->SetTitle("Number of Events");
    endposy->Draw();
    endposy_FV2->Draw("same");
    //legend
    TLegend *leg14 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg14->SetBorderSize(0);
    leg14->SetFillStyle(0);
    leg14->AddEntry(startposz, "Original", "lep");
    leg14->AddEntry(startposz_FV2,"F.V + Michel electron selection", "lep");
    leg14->SetTextFont(62);
    leg14->Draw();
    gStyle->SetOptStat(10);
    canvas17->Print("plots.ps");
    myFile->WriteObject(canvas17, "Michel Electron Selection End Position of the Reconstructed Event Tracks (y)");
    
    canvas18->cd();
    endposz->GetXaxis()->SetTitle("z-position (cm)");
    endposz->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposz->GetYaxis()->SetTitle("Number of Events");
    endposz->Draw();
    endposz_FV2->Draw("same");
    //legend
    TLegend *leg15 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg15->SetBorderSize(0);
    leg15->SetFillStyle(0);
    leg15->AddEntry(startposz, "Original", "lep");
    leg15->AddEntry(startposz_FV2,"F.V + Michel electron selection", "lep");
    leg15->SetTextFont(62);
    leg15->Draw();
    gStyle->SetOptStat(10);
    canvas18->Print("plots.ps");
    myFile->WriteObject(canvas18, "Michel Electron Selection End Position of the Reconstructed Event Tracks (z)");
    
    canvas19->cd();
    track_length->GetXaxis()->SetTitle("track length (cm)");
    //track_length->GetYaxis()->SetRangeUser(0, endposz->GetMaximum()*1.1);
    track_length->GetYaxis()->SetTitle("Number of Events");
    track_length->Draw();
    track_length_FV2->Draw("same");
    track_length_ME->Draw("same");
    track_length_ME_no_n_showers->Draw("same");
    //legend
    TLegend *leg16 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg16->SetBorderSize(0);
    leg16->SetFillStyle(0);
    leg16->AddEntry(track_length, "Original", "lep");
    leg16->AddEntry(track_length_FV,"F.V + Michel electron selection", "lep");
    leg16->AddEntry(track_length_ME,"Track length of Michel electron candidates", "lep");
    leg16->AddEntry(track_length_ME_no_n_showers, "ME candidates track length without the n_shower selection", "lep");
    leg16->SetTextFont(62);
    leg16->Draw();
    gStyle->SetOptStat(10);
    canvas19->Print("plots.ps");
    myFile->WriteObject(canvas19, "Michel Electron Selection Track Length of the Reconstructed Event Tracks");
    
    canvas20->cd();
    MuonME_gap->GetXaxis()->SetTitle("Distance (cm)");
    MuonME_gap->GetYaxis()->SetTitle("Number of Events");
    MuonME_gap->Draw();
    gStyle->SetOptStat(10);
    canvas20->Print("plots.ps");
    myFile->WriteObject(canvas20, "Muon/Michel Electron Gap Distance");
    
    canvas21->cd();
    MuonME_gap_zoom->GetXaxis()->SetTitle("Distance (cm)");
    MuonME_gap_zoom->GetYaxis()->SetTitle("Number of Events");
    MuonME_gap_zoom->Draw();
    gStyle->SetOptStat(10);
    canvas21->Print("plots.ps");
    myFile->WriteObject(canvas21, "Muon/Michel Electron Gap Distance (Zoomed in)");
    
    canvas22->cd();
    track_score_ME->GetXaxis()->SetTitle("Track Score (AU)");
    track_score_ME->GetYaxis()->SetTitle("Number of Events");
    track_score_ME->Draw();
    gStyle->SetOptStat(10);
    canvas22->Print("plots.ps");
    myFile->WriteObject(canvas22, "Track Score of the Michel Electron Candidates");
    
    canvas23->cd();
    Track_Energy->GetXaxis()->SetTitle("Energy (MeV)");
    Track_Energy->GetYaxis()->SetTitle("Number of Events");
    Track_Energy->Draw();
    Track_Energy_Muon->Draw("same");
    Track_Energy_ME->Draw("same");
    TLegend *leg17 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg17->SetBorderSize(0);
    leg17->SetFillStyle(0);
    leg17->AddEntry(Track_Energy, "track energy", "lep");
    leg17->AddEntry(Track_Energy_Muon,"Muon track energy", "lep");
    leg17->AddEntry(Track_Energy_ME,"Michel electron track energy", "lep");
    leg17->SetTextFont(62);
    leg17->Draw();
    gStyle->SetOptStat(10);
    canvas23->Print("plots.ps");
    myFile->WriteObject(canvas23, "Track Energy");
    
    canvas24->cd();
    Track_Energy_u->GetXaxis()->SetTitle("Energy (MeV)");
    Track_Energy_u->GetYaxis()->SetTitle("Number of Events");
    Track_Energy_u->Draw();
    Track_Energy_Muon_u->Draw("same");
    Track_Energy_ME_u->Draw("same");
    TLegend *leg18 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg18->SetBorderSize(0);
    leg18->SetFillStyle(0);
    leg18->AddEntry(Track_Energy_u, "track energy (u)", "lep");
    leg18->AddEntry(Track_Energy_Muon_u,"Muon track energy (u)", "lep");
    leg18->AddEntry(Track_Energy_ME_u,"Michel electron track energy (u)", "lep");
    leg18->SetTextFont(62);
    leg18->Draw();
    gStyle->SetOptStat(10);
    canvas24->Print("plots.ps");
    myFile->WriteObject(canvas24, "Track Energy (u Plane)");
    
    canvas25->cd();
    Track_Energy_v->GetXaxis()->SetTitle("Energy (MeV)");
    Track_Energy_v->GetYaxis()->SetTitle("Number of Events");
    Track_Energy_v->Draw();
    Track_Energy_Muon_v->Draw("same");
    Track_Energy_ME_v->Draw("same");
    TLegend *leg19 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg19->SetBorderSize(0);
    leg19->SetFillStyle(0);
    leg19->AddEntry(Track_Energy_v, "track energy (v)", "lep");
    leg19->AddEntry(Track_Energy_Muon_v,"Muon track energy (v)", "lep");
    leg19->AddEntry(Track_Energy_ME_v,"Michel electron track energy (v)", "lep");
    leg19->SetTextFont(62);
    leg19->Draw();
    gStyle->SetOptStat(10);
    canvas25->Print("plots.ps");
    myFile->WriteObject(canvas25, "Track Energy (v Plane)");
    
    canvas26->cd();
    startposx_ME->GetXaxis()->SetTitle("x-position(cm)");
    startposx_ME->GetYaxis()->SetTitle("Number of Events");
    startposx_ME->GetYaxis()->SetRangeUser(0, startposx_ME->GetMaximum()*1.1);
    startposx_ME->Draw();
    gStyle->SetOptStat(10);
    canvas26->Print("plots.ps");
    myFile->WriteObject(canvas26, "Start Position of the Selected Michel Electron Candidates (x)");
    
    canvas27->cd();
    startposy_ME->GetXaxis()->SetTitle("y-position(cm)");
    startposy_ME->GetYaxis()->SetTitle("Number of Events");
    startposy_ME->GetYaxis()->SetRangeUser(0, startposy_ME->GetMaximum()*1.1);
    startposy_ME->Draw();
    gStyle->SetOptStat(10);
    canvas27->Print("plots.ps");
    myFile->WriteObject(canvas27, "Start Position of the Selected Michel Electron Candidates (y)");
    
    canvas28->cd();
    startposz_ME->GetXaxis()->SetTitle("z-position(cm)");
    startposz_ME->GetYaxis()->SetTitle("Number of Events");
    startposz_ME->GetYaxis()->SetRangeUser(0, startposz_ME->GetMaximum()*1.1);
    startposz_ME->Draw();
    gStyle->SetOptStat(10);
    canvas28->Print("plots.ps)");
    myFile->WriteObject(canvas28, "Start Position of the Selected Michel Electron Candidates (z)");
    
    
}

/*
    
void Myclass::DrawSingleCanvas(TCanvas *canvas,TH1F *histo, TString axis_title ,TString label, TLegend *leg1)
    {
        canvas->cd();
        histo->GetXaxis()->SetTitle(axis_title);
        histo->GetYaxis()->SetRangeUser(-10, histo->GetMaximum()*1.1);
        histo->GetYaxis()->SetTitle("Number of Events");
        histo->Draw();
        // startposx_FV->Draw("same");
        //  startposx_stp_trk->Draw("same");
        //legend
        //leg1->AddEntry(histo, label, "lep");
        //std::cout << "first: " << "histo "<< histo << "label" << label << std::endl;
        //leg1->AddEntry(startposx_FV,"x position, Fid Vol", "lep");
        // leg1->AddEntry(startposx_stp_trk,"x position, Stopping Trks", "lep");
       // leg1->Draw();
        
    }
    
    
void Myclass::DrawSecondHistOnCanvas(TCanvas *canvas,TH1F *histo,TString label, Color_t color, TLegend *leg1)
    {
        canvas->cd();
        histo->SetLineColor(color);
        histo->Draw("same");
        //leg1->AddEntry(histo, label, "lep");
        //std::cout << "second: " << "histo "<< histo << "label" << label << std::endl;
        //leg1->Draw();
        //canvas->GetLegend->AddEntry(histo, label, "lep");
    }
    
void Myclass::InitializeLegend(TCanvas *canvas, TH1F *histo, TString label)
    {
        canvas -> cd();
        TLegend *leg1 = new TLegend(0.55, 0.8, 0.9, 0.9);
        leg1 -> SetBorderSize(0);
        leg1 -> SetFillStyle(0);
        leg1->SetTextFont(58);
    //a while loop here? but what should i be importing in the parameter?
        leg1->AddEntry(histo, label, "lep");
    
    }
*/
