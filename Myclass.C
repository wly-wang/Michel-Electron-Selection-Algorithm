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
    TH1F * startposy = new TH1F("startposy","startposy",44,-116,116);
    TH1F * startposy_FV = new TH1F("startposy_FV","startposy_cut",44,-116,116);
    TH1F * startposy_stp_trk = new TH1F("startposy_stp_trk","startposy_stp_trk",44,-116,116);
    TH1F * startposy_FV2 = new TH1F("startposy_FV2","startposy_FV2",44,-116,116);
    TH1F * startposy_ME = new TH1F("startposy_ME","Michel Electron y-Start Position",5,-116,116);
    startposy_FV->SetLineColor(kRed);
    startposy_FV2->SetLineColor(kRed);
    startposy_stp_trk->SetLineColor(kMagenta);
    
    TH1F * startposx = new TH1F("startposx","startposx",44,0,256);
    TH1F * startposx_FV = new TH1F("startposx_FV","startposx_FV",44,0,256);
    TH1F * startposx_stp_trk = new TH1F("startposx_stp_trk","startposx_stp_trk",44,0,256);
    TH1F * startposx_FV2 = new TH1F("startposx_FV2","startposx_FV2",44,0,256);
    TH1F * startposx_ME = new
    TH1F("startposx_ME","Michel Electron x-Start Position",5,0,256);
    startposx_FV->SetLineColor(kRed);
    startposx_FV2->SetLineColor(kRed);
    startposx_stp_trk->SetLineColor(kMagenta);
    
    TH1F * startposz = new TH1F("startposz","startposz",44,0,1036);
    TH1F * startposz_FV = new TH1F("startposz_FV","startposz_FV",44,0,1036);
    TH1F * startposz_stp_trk = new TH1F("startposz_stp_trk","startposz_stp_trk",44,0,1036);
    TH1F * startposz_FV2 = new TH1F("startposz_FV2","startposz_FV2",44,0,1036);
    TH1F * startposz_ME = new TH1F("startposz_ME","Michel Electron z-Start Position",5,0,1036);
    startposz_FV->SetLineColor(kRed);
    startposz_FV2->SetLineColor(kRed);
    startposz_stp_trk->SetLineColor(kMagenta);
    
    TH1F * endposx = new TH1F("endposx","endposx",44,0,256);
    TH1F * endposx_FV = new TH1F("endposx_FV","endposx_FV",44,0,256);
    TH1F * endposx_stp_trk = new TH1F("endposx_stp_trk","endposx_stp_trk",44,0,256);
    TH1F * endposx_muon_selection = new TH1F("endposx_muon_selection","endposx_muon_selection",13,0,256);
    TH1F * endposx_FV2 = new TH1F("endposx_FV2","endposx_FV2",44,0,256);
    endposx_FV->SetLineColor(kRed);
    endposx_FV2->SetLineColor(kRed);
    endposx_stp_trk->SetLineColor(kMagenta);
    endposx_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * endposz = new TH1F("endposz","endposz",44,0,1036);
    TH1F * endposz_FV = new TH1F("endposz_FV","endposz_FV",44,0,1036);
    TH1F * endposz_stp_trk = new TH1F("endposz_stp_trk","endposz_stp_trk",44,0,1036);
    TH1F * endposz_muon_selection = new TH1F("endposz_muon_selection","endposz_muon_selection",52,0,1036);
    TH1F * endposz_FV2 = new TH1F("endposz_FV2","endposz_FV2",44,0,1036);
    endposz_FV->SetLineColor(kRed);
    endposz_FV2->SetLineColor(kRed);
    endposz_stp_trk->SetLineColor(kMagenta);
    endposz_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * endposy = new
    TH1F("endposy","endposy",44,-116,116);
    TH1F * endposy_FV = new TH1F("endposy_FV","endposy_FV",44,-116,116);
    TH1F * endposy_stp_trk = new TH1F("endposy_stp_trk","endposy_stp_trk",44,-116,116);
    TH1F * endposy_muon_selection = new TH1F("endposy_muon_selection","endposy_muon_selection",12,-116,116);
    TH1F * endposy_FV2 = new TH1F("endposy_FV2","endposy_FV2",44,-116,116);
    endposy_FV->SetLineColor(kRed);
    endposy_FV2->SetLineColor(kRed);
    endposy_stp_trk->SetLineColor(kMagenta);
    endposy_muon_selection->SetLineColor(kViolet+3);
    
    TH1F * track_length = new
    TH1F ("track_length", "track_length", 80, 0, 450);
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
    TH1F ("track_score", "track_score", 50, 0, 1);
    TH1F * track_score_cut = new TH1F ("track_score_cut", "track_score_cut", 50, 0, 1);
    TH1F * track_score_ME = new
    TH1F ("track_score_ME", "track_score_ME", 50, 0, 1);
    track_score_cut->SetLineColor(kRed);
    track_score_ME->SetLineColor(kRed);
    
    TH1F * track_PID = new
    TH1F ("track_PID", "track_PID", 50, -1, 1);
    TH1F * track_PID_cut = new TH1F ("track_PID_cut", "track_PID_cut", 50, -1, 1);
    TH1F * track_PID_ME = new TH1F ("track_PID_ME", "track_PID_ME", 50, -1, 1);
    track_PID_cut->SetLineColor(kRed);
    track_PID_ME->SetLineColor(kRed);
    
    TH1F * MuonME_gap = new
    TH1F ("MuonME_gap", "Muon/Michel Electron Candidates Gap Distance", 100, 0, 150);
    TH1F * MuonME_gap_zoom = new
    TH1F ("MuonME_gap_zoom", "MuonME_gap_zoom", 100, 0, 20);
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
    
    //open particle_selection txt file (for event display)
    fstream file;
    file.open("selected_muon_info.txt", ios::app | ios::out);
    
    //open me_info_len_score_pid (printing out to see if values are sensible
    fstream file1;
    file1.open("selected_me_info.txt", ios::app | ios::out);
    
    //open txt file to look at ME with energy >50 MeV
    fstream file2;
    file2.open("ME_E_more_50mev.txt", ios::app | ios::out);
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        //Defining new vectors for storing muon information
        vector<float> startposx_1;
        vector<float> startposy_1;
        vector<float> startposz_1;
        vector<float> endposx_1;
        vector<float> endposy_1;
        vector<float> endposz_1;
        vector<float> trk_length;
        
        bool GoodMuonExists=false;
        //bool GoodMichelExists=false;
        vector<int> goodmuonindex;
        
        //cut out unphysical regions in our position variables
        for(unsigned int itrack=0;
            itrack<(*trk_sce_start_x_v).size(); itrack++)
        {
            
            
            //helper variabls
            float stposx,stposy,stposz,enposx,enposy,enposz,trklen, trk_score,trk_pid,trk_E_y,trk_E_u,trk_E_v;
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
            
            /* if (run == 6959 && sub == 201 && evt == 10056)
             {
             cout << "x start: " << stposx << " y start: " << stposy << " z start: " << stposz << " x end: " << enposx << " y end: " << enposy << " z end: " << enposz << endl;
             }
             else continue; */
            
            // use only sensible vectors
            if(stposx > 0 && stposy > -115 && stposz >0 && enposx >0 && enposy > -115 && enposz > 0 && trk_E_y > -1)
            {
                trklen = calc_len(stposx,stposy,stposz,enposx,enposy,enposz);
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
                
                
                // add fiducial volume cuts
                if(stposx > 10 && stposx < 246 && enposx > 10 && enposx < 246
                   && stposz >10 && stposz < 1026 && enposz > 10 && enposz < 1026)
                {
                    // fill second tier of histograms (start/end/length)
                    
                    startposx_FV -> Fill(stposx);
                    startposy_FV -> Fill(stposy);
                    startposz_FV -> Fill(stposz);
                    endposx_FV -> Fill(enposx);
                    endposy_FV -> Fill(enposy);
                    endposz_FV -> Fill(enposz);
                    track_length_FV -> Fill(trklen);
                    
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
                        
                        // when all cuts satisfied, set GoodMuonExists=true;
                        if(trk_score > 0.5 && trk_pid > 0 && trklen > 10)
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
                float stposx,stposy,stposz,enposx,enposy,enposz,trklen,mu_enposx,mu_enposy,mu_enposz,mu_me_trklen,mu_mest_trklen,mu_meen_trklen,trk_score,trk_pid,trk_E_y,trk_E_u,trk_E_v;
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
                        
                        MuonME_gap->Fill(shorter_trklen);
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
                                
                                
                                file1 << run << " "<< sub << " "<< evt << " " << stposx << " " << stposy << " " << stposz << " " << enposx << " " << enposy << " " << enposz << " " << trklen << " " << trk_score << " " << trk_pid << " " << shorter_trklen << " " << n_tracks << " " << n_showers << " " << " " << trk_E_y << endl;
                                
                                if(trk_E_y > 50)
                                {
                                    file2 << run << " " << sub << " "<< evt << " " << stposx << " " << stposy << " " << stposz << " " << enposx << " " << enposy << " " << enposz << " " << trklen << " " << trk_score << " " << trk_pid << " " << trk_E_y << " " << endl;
                                }
                                
                                
                            }
                        }
                        
                    }
                }
            }  // end of for loop on tracks
        } // end of "if GoodMuonExists is true"
        
        // if (Cut(ientry) < 0) continue;
    }
    //closing particle information file
    file.close();
    file1.close();
    file2.close();
    cout<< "Selected muon info file for evt displayed is filled" <<endl;
    cout<< "Selected Michel electron info file filled" <<endl;
    cout<< "Selected Michel electrons info w/ E>50MeV file filled" <<endl;
    
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
    
    //actual drawing
    
    //canvas1->cd();
    //startposx->GetXaxis()->SetTitle("x-position (cm)");
    //startposx->GetYaxis()->SetRangeUser(-10, startposx->GetMaximum()*1.1);
    //startposx->Draw();
    //startposx_FV->Draw("same");
    //startposx_stp_trk->Draw("same");
    //legend
    //TLegend *leg1 = new TLegend(0.55, 0.8, 0.9, 0.9);
    //leg1->SetBorderSize(0);
    //leg1->SetFillStyle(0);
    //leg1->AddEntry(startposx, "x position, no cuts", "lep");
    //leg1->AddEntry(startposx_FV,"x position, Fid Vol", "lep");
    //leg1->AddEntry(startposx_stp_trk,"x position, Stopping Trks", "lep");
    //leg1->SetTextFont(62);
    //leg1->Draw();
    
    TLegend *leg1 = InitializeLegend(canvas1);
    DrawSingleCanvas(canvas1,startposx,"x-position (cm)", leg1);
    DrawSecondHistOnCanvas(canvas1,startposx_FV,"fiducial volume",kRed, leg1);
    

    
    return;
    
    canvas2->cd();
    startposy->GetXaxis()->SetTitle("y-position (cm)");
    startposy->GetYaxis()->SetRangeUser(-10, startposy->GetMaximum()*1.1);
    startposy->Draw();
    startposy_FV->Draw("same");
    startposy_stp_trk->Draw("same");
    //legend
    TLegend *leg2 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(startposy, "y position, no cuts", "lep");
    leg2->AddEntry(startposy_FV,"y position, Fid Volume", "lep");
    leg2->AddEntry(startposy_stp_trk,"y position, Stopping Trks", "lep");
    leg2->SetTextFont(62);
    leg2->Draw();
    
    canvas3->cd();
    startposz->GetXaxis()->SetTitle("z-position (cm)");
    startposz->GetYaxis()->SetRangeUser(-10, startposz->GetMaximum()*1.1);
    startposz->Draw();
    startposz_FV->Draw("same");
    startposz_stp_trk->Draw("same");
    //legend
    TLegend *leg3 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(startposz, "z position, no cuts", "lep");
    leg3->AddEntry(startposz_FV,"z position, Fid Volume", "lep");
    leg3->AddEntry(startposz_stp_trk,"z position, Stopping Trks", "lep");
    leg3->SetTextFont(62);
    leg3->Draw();
    
    canvas4->cd();
    endposx->GetXaxis()->SetTitle("x-position (cm)");
    endposx->GetYaxis()->SetRangeUser(-10, endposx->GetMaximum()*1.1);
    endposx->Draw();
    endposx_FV->Draw("same");
    endposx_stp_trk->Draw("same");
    //legend
    TLegend *leg4 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->AddEntry(endposx, "x position, no cuts", "lep");
    leg4->AddEntry(endposx_FV,"x position, Fid Vol", "lep");
    leg4->AddEntry(endposx_stp_trk,"x position, Stopping Trks", "lep");
    leg4->SetTextFont(62);
    leg4->Draw();
    
    canvas5->cd();
    endposy->GetXaxis()->SetTitle("y-position (cm)");
    endposy->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposy->Draw();
    endposy_FV->Draw("same");
    endposy_stp_trk->Draw("same");
    //legend
    TLegend *leg5 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);
    leg5->AddEntry(endposy, "y position, no cuts", "lep");
    leg5->AddEntry(endposy_FV,"y position, Fid Vol", "lep");
    leg5->AddEntry(endposy_stp_trk,"y position, Stopping Trks", "lep");
    leg5->SetTextFont(62);
    leg5->Draw();
    
    canvas6->cd();
    endposz->GetXaxis()->SetTitle("z-position (cm)");
    endposz->GetYaxis()->SetRangeUser(-10, endposz->GetMaximum()*1.1);
    endposz->Draw();
    endposz_FV->Draw("same");
    endposz_stp_trk->Draw("same");
    //legend
    TLegend *leg6 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg6->SetBorderSize(0);
    leg6->SetFillStyle(0);
    leg6->AddEntry(endposz, "z position, no cuts", "lep");
    leg6->AddEntry(endposz_FV,"z position, Fid Vol", "lep");
    leg6->AddEntry(endposz_stp_trk,"z position, Stopping Trks", "lep");
    leg6->SetTextFont(62);
    leg6->Draw();
    
    canvas7->cd();
    track_length->GetXaxis()->SetTitle("track length (cm)");
    //track_length->GetYaxis()->SetRangeUser(0, endposz->GetMaximum()*1.1);
    track_length->Draw();
    track_length_FV->Draw("same");
    track_length_stp_trk->Draw("same");
    track_length_score_cuts->Draw("same");
    //legend
    TLegend *leg7 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg7->SetBorderSize(0);
    leg7->SetFillStyle(0);
    leg7->AddEntry(track_length, "track length, no cuts", "lep");
    leg7->AddEntry(track_length_FV,"track length, Fid Vol", "lep");
    leg7->AddEntry(track_length_stp_trk,"track length, Stopping Trks", "lep");
    leg7->AddEntry(track_length_score_cuts,"track length, score/length cuts", "lep");
    leg7->SetTextFont(62);
    leg7->Draw();
    
    canvas8->cd();
    track_score->GetXaxis()->SetTitle("track score (AU)");
    track_score->Draw();
    track_score_cut->Draw("same");
    //legend
    TLegend *leg8 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg8->SetBorderSize(0);
    leg8->SetFillStyle(0);
    leg8->AddEntry(track_score, "track score, FV & stopping trks", "lep");
    leg8->AddEntry(track_score_cut,"track score > 0.5 & track PID > 0.5 & track length > 10cm", "lep");
    leg8->SetTextFont(62);
    leg8->Draw();
    
    canvas9->cd();
    track_PID->GetXaxis()->SetTitle("track PID score (AU)");
    track_PID->Draw();
    track_PID_cut->Draw("same");
    //legend
    TLegend *leg9 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->AddEntry(track_PID, "track PID score, FV & stopping trks", "lep");
    leg9->AddEntry(track_PID_cut,"track score > 0.5 & track PID > 0.5 & track length > 10cm", "lep");
    leg9->SetTextFont(62);
    leg9->Draw();
    
    canvas10->cd();
    endposx_muon_selection->GetXaxis()->SetTitle("x-position(cm)");
    endposx_muon_selection->GetYaxis()->SetRangeUser(0, endposx_muon_selection->GetMaximum()*1.1);
    endposx_muon_selection->Draw();
    
    canvas11->cd();
    endposy_muon_selection->GetXaxis()->SetTitle("y-position(cm)");
    endposy_muon_selection->GetYaxis()->SetRangeUser(0, endposy_muon_selection->GetMaximum()*1.1);
    endposy_muon_selection->Draw();
    
    canvas12->cd();
    endposz_muon_selection->GetXaxis()->SetTitle("z-position(cm)");
    endposz_muon_selection->GetYaxis()->SetRangeUser(0, endposz_muon_selection->GetMaximum()*1.1);
    endposz_muon_selection->Draw();
    
    canvas13->cd();
    startposx->GetXaxis()->SetTitle("x-position (cm)");
    startposx->GetYaxis()->SetRangeUser(-10, startposx->GetMaximum()*1.1);
    startposx->Draw();
    startposx_FV2->Draw("same");
    //legend
    TLegend *leg10 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(startposx, "startpos x, no cuts", "lep");
    leg10->AddEntry(startposx_FV2,"startpos x, fiducial volume w/ ME Selection", "lep");
    leg10->SetTextFont(62);
    leg10->Draw();
    
    canvas14->cd();
    startposy->GetXaxis()->SetTitle("y-position (cm)");
    startposy->GetYaxis()->SetRangeUser(-10, startposy->GetMaximum()*1.1);
    startposy->Draw();
    startposy_FV2->Draw("same");
    //legend
    TLegend *leg11 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg11->SetBorderSize(0);
    leg11->SetFillStyle(0);
    leg11->AddEntry(startposy, "startpos y, no cuts", "lep");
    leg11->AddEntry(startposy_FV2,"startpos y, fiducial volume w/ ME Selection", "lep");
    leg11->SetTextFont(62);
    leg11->Draw();
    
    canvas15->cd();
    startposz->GetXaxis()->SetTitle("z-position (cm)");
    startposz->GetYaxis()->SetRangeUser(-10, startposz->GetMaximum()*1.1);
    startposz->Draw();
    startposz_FV2->Draw("same");
    //legend
    TLegend *leg12 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg12->SetBorderSize(0);
    leg12->SetFillStyle(0);
    leg12->AddEntry(startposz, "startpos z, no cuts", "lep");
    leg12->AddEntry(startposz_FV2,"startpos z, fiducial volume w/ ME Selection", "lep");
    leg12->SetTextFont(62);
    leg12->Draw();
    
    canvas16->cd();
    endposx->GetXaxis()->SetTitle("x-position (cm)");
    endposx->GetYaxis()->SetRangeUser(-10, endposx->GetMaximum()*1.1);
    endposx->Draw();
    endposx_FV2->Draw("same");
    //legend
    TLegend *leg13 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg13->SetBorderSize(0);
    leg13->SetFillStyle(0);
    leg13->AddEntry(startposz, "endpos x, no cuts", "lep");
    leg13->AddEntry(startposz_FV2,"endpos x, fiducial volume w/ ME Selection", "lep");
    leg13->SetTextFont(62);
    leg13->Draw();
    
    canvas17->cd();
    endposy->GetXaxis()->SetTitle("y-position (cm)");
    endposy->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposy->Draw();
    endposy_FV2->Draw("same");
    //legend
    TLegend *leg14 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg14->SetBorderSize(0);
    leg14->SetFillStyle(0);
    leg14->AddEntry(startposz, "endpos y, no cuts", "lep");
    leg14->AddEntry(startposz_FV2,"endpos y, fiducial volume w/ ME Selection", "lep");
    leg14->SetTextFont(62);
    leg14->Draw();
    
    canvas18->cd();
    endposz->GetXaxis()->SetTitle("z-position (cm)");
    endposz->GetYaxis()->SetRangeUser(-10, endposy->GetMaximum()*1.1);
    endposz->Draw();
    endposz_FV2->Draw("same");
    //legend
    TLegend *leg15 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg15->SetBorderSize(0);
    leg15->SetFillStyle(0);
    leg15->AddEntry(startposz, "endpos z, no cuts", "lep");
    leg15->AddEntry(startposz_FV2,"endpos z, fiducial volume w/ ME Selection", "lep");
    leg15->SetTextFont(62);
    leg15->Draw();
    
    canvas19->cd();
    track_length->GetXaxis()->SetTitle("track length (cm)");
    //track_length->GetYaxis()->SetRangeUser(0, endposz->GetMaximum()*1.1);
    track_length->Draw();
    track_length_FV2->Draw("same");
    track_length_ME->Draw("same");
    track_length_ME_no_n_showers->Draw("same");
    //legend
    TLegend *leg16 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg16->SetBorderSize(0);
    leg16->SetFillStyle(0);
    leg16->AddEntry(track_length, "track length, no cuts", "lep");
    leg16->AddEntry(track_length_FV,"track length, Fid Vol w/ ME selection", "lep");
    leg16->AddEntry(track_length_ME,"track length of ME candidates", "lep");
    leg16->AddEntry(track_length_ME_no_n_showers, "ME candidates track length w/o n_shower selection", "lep");
    leg16->SetTextFont(62);
    leg16->Draw();
    
    canvas20->cd();
    MuonME_gap->GetXaxis()->SetTitle("Distance (cm)");
    MuonME_gap->Draw();
    
    canvas21->cd();
    MuonME_gap_zoom->GetXaxis()->SetTitle("Distance (cm)");
    MuonME_gap_zoom->Draw();
    
    canvas22->cd();
    track_score_ME->GetXaxis()->SetTitle("Track Score (AU)");
    track_score_ME->Draw();
    
    canvas23->cd();
    Track_Energy->GetXaxis()->SetTitle("Energy (MeV)");
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
    
    canvas24->cd();
    Track_Energy_u->GetXaxis()->SetTitle("Energy (MeV)");
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
    
    canvas25->cd();
    Track_Energy_v->GetXaxis()->SetTitle("Energy (MeV)");
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
    
    canvas26->cd();
    startposx_ME->GetXaxis()->SetTitle("x-position(cm)");
    startposx_ME->GetYaxis()->SetTitle("Number of Events");
    startposx_ME->GetYaxis()->SetRangeUser(0, startposx_ME->GetMaximum()*1.1);
    startposx_ME->Draw();
    
    canvas27->cd();
    startposy_ME->GetXaxis()->SetTitle("y-position(cm)");
    startposy_ME->GetYaxis()->SetTitle("Number of Events");
    startposy_ME->GetYaxis()->SetRangeUser(0, startposy_ME->GetMaximum()*1.1);
    startposy_ME->Draw();
    
    canvas28->cd();
    startposz_ME->GetXaxis()->SetTitle("z-position(cm)");
    startposz_ME->GetYaxis()->SetTitle("Number of Events");
    startposz_ME->GetYaxis()->SetRangeUser(0, startposz_ME->GetMaximum()*1.1);
    startposz_ME->Draw();
}
    
    
void Myclass::DrawSingleCanvas(TCanvas *canvas1,TH1F *histo,TString label, TLegend *leg1)
    {
        canvas1->cd();
        histo->GetXaxis()->SetTitle(label);
        histo->GetYaxis()->SetRangeUser(-10, histo->GetMaximum()*1.1);
        histo->Draw();
        // startposx_FV->Draw("same");
        //  startposx_stp_trk->Draw("same");
        //legend
        leg1->AddEntry(histo, label, "lep");
        //leg1->AddEntry(startposx_FV,"x position, Fid Vol", "lep");
        // leg1->AddEntry(startposx_stp_trk,"x position, Stopping Trks", "lep");
        leg1->Draw();
        
    }
    
    
void Myclass::DrawSecondHistOnCanvas(TCanvas *canvas1,TH1F *histo,TString label, Color_t color, TLegend *leg1)
    {
        canvas1->cd();
        histo->SetLineColor(color);
        histo->Draw("same");
        leg1->AddEntry(histo, label, "lep");
        leg1->Draw();
        //canvas->GetLegend->AddEntry(histo, label, "lep");
    }
    
TLegend* Myclass::InitializeLegend(TCanvas *canvas1)
    {
        canvas1 -> cd();
        TLegend *leg1 = new TLegend(0.55, 0.8, 0.9, 0.9);
        leg1 -> SetBorderSize(0);
        leg1 -> SetFillStyle(0);
        leg1->SetTextFont(58);
        return leg1;
    }

