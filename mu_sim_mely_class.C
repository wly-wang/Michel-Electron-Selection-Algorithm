#define mu_sim_mely_class_cxx
#include "mu_sim_mely_class.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
using namespace std;

/*———————————————————————————————————————————————————————————————
 here contains functions written that would be called and used
 within the main algorithm
 ————————————————————————————————————————————————————————————————*/
//Track length function (general distance formula)
float calc_len(float stpx,float stpy,float stpz,float enpx,float enpy,float enpz)
{
    return (sqrt(pow(enpx-stpx, 2) + pow(enpy-stpy, 2) + pow(enpz-stpz, 2)));
}

//function that calculates the distance between flash center and event track mid-pt
float calc_flash_evt_dist(float flash_ycent, float flash_zcent, float starty, float endy,float startz, float endz)
{
    return (sqrt( pow(flash_zcent - ((startz+endz)/2.),2) + pow(flash_ycent - ((starty+endy)/2.),2) ));
}

//function that calculates the z-axis proportion of deviation between flash center and event track mid-point
float calc_devi_flash_trk(float flash_zcent, float flash_zwid, float startz, float endz)
{
    return (abs(flash_zcent - ((startz+endz)/2.))/flash_zwid);
}

/*——————————————————————————————————————————————————————————————*/

void mu_sim_mely_class::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mu_sim_mely_class.C
//      root> mu_sim_mely_class <choice-of-analyzer-name>
//      root> <choice-of-analyzer-name>.Loop(); //loop on all entries
//      root> <choice-of-analyzer-name>.GetEntry(12); // Fill analyzer data members with entry number 12
//      root> <choice-of-analyzer-name>.Show();       // Show values of entry 12
//      root> <choice-of-analyzer-name>.Show(16);     // Read and show values of entry 16
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
    
    //define new empty plots for filling later in the analyzer.
    
    //Empty histo for storing Truth information for the simulated muon tracks + truth energy for Michel
    TH1F* h_mc_vx = new TH1F("h_mc_vx","MC Start Position (x)",100,-1,257);
    TH1F* h_mc_vy = new TH1F("h_mc_vy","MC Start Position (y)",100,-117,117);
    TH1F* h_mc_vz = new TH1F("h_mc_vz","MC Start Position (z)",100,-1,1037);
    TH1F* h_mc_endx = new TH1F("h_mc_endx","MC End Position (x)",100,-1,257);
    TH1F* h_mc_endy = new TH1F("h_mc_endy","MC End Position (y)",100,-117,117);
    TH1F* h_mc_endz = new TH1F("h_mc_endz","MC End Position (z)",100,-1,1037);
    TH1F* h_mc_pdg = new TH1F("h_mc_pdg","MC PDG Code",100,-15,15);
    TH1F* h_mc_E = new TH1F("h_mc_E","MC Energy",100,190,610);
    TH1F* h_mc_michel_E = new TH1F("h_mc_michel_E","MC Michel Electron Energy",100,-1,70);
    
    //Empty histo for storing original reco start/end position
    TH1F* h_startposx = new TH1F("h_startposx","Reconstructed Track Start Position (x)",100,-1,257);
    TH1F* h_startposy = new TH1F("h_startposy","Reconstructed Track Start Position (y)",100,-117,117);
    TH1F* h_startposz = new TH1F("h_startposz","Reconstructed Track Start Position (z)",100,-1,1037);
    TH1F* h_endposx = new TH1F("h_endposx","Reconstructed Track End Position (x)",100,-1,257);
    TH1F* h_endposy = new TH1F("h_endposy","Reconstructed Track End Position (y)",100,-117,117);
    TH1F* h_endposz = new TH1F("h_endposz","Reconstructed Track End Position (z)",100,-1,1037);
    
    //Empty histo for non-selected reconstructed track length and reconstructed track (calorimetric) energy
    TH1F* h_reco_trk_len = new TH1F("h_reco_trk_len","Reconstructed Track Length",100,-1,500);
    TH1F* h_reco_trk_E = new TH1F("h_reco_trk_E","Reconstructed Track Energy",100,-1,1000);

    //Empty histo for storing unselected track score and track PID information
    TH1F* h_trk_score = new TH1F("h_trk_score","Track Score",100,0,1);
    TH1F* h_trk_pid = new TH1F("h_trk_pid","Track PID",100,-1,1);
    
    //Empty histo for start/end position and track length after fiducial volume selection
    TH1F* h_startposx_fv = new TH1F("h_startposx_fv","Fiducial Volume Reconstructed Start Position (x)",100,-1,257);
    TH1F* h_startposy_fv = new TH1F("h_startposy_fv","Fiducial Volume Reconstructed Start Position (y)",100,-117,117);
    TH1F* h_startposz_fv = new TH1F("h_startposz_fv","Fiducial Volume Reconstructed Start Position (z)",100,-1,1037);
    TH1F* h_endposx_fv = new TH1F("h_endposx_fv","Fiducial Volume Reconstructed End Position (x)",100,-1,257);
    TH1F* h_endposy_fv = new TH1F("h_endposy_fv","Fiducial Volume Reconstructed End Position (y)",100,-117,117);
    TH1F* h_endposz_fv = new TH1F("h_endposz_fv","Fiducial Volume Reconstructed End Position (z)",100,-1,1037);
    TH1F* h_reco_trk_len_fv = new TH1F("h_reco_trk_len_fv","Fiducial Volume Reconstructed Track Length",100,0,500);
    h_startposx_fv -> SetLineColor(kRed);
    h_startposy_fv -> SetLineColor(kRed);
    h_startposz_fv -> SetLineColor(kRed);
    h_endposx_fv -> SetLineColor(kRed);
    h_endposy_fv -> SetLineColor(kRed);
    h_endposz_fv -> SetLineColor(kRed);
    h_reco_trk_len_fv -> SetLineColor(kRed);
    
    //empty histo for storing information after selecting tracks coming in from the top of detector and stopping inside the AV
    TH1F* h_startposx_st = new TH1F("h_startposx_st","Reconstructed Start Position (x) After Stopping Tracks Selection",100,-1,257);
    TH1F* h_startposy_st = new TH1F("h_startposy_st","Reconstructed Start Position (y) After Stopping Tracks Selection",100,-117,117);
    TH1F* h_startposz_st = new TH1F("h_startposz_st","Reconstructed Start Position (z) After Stopping Tracks Selection",100,-1,1037);
    TH1F* h_endposx_st = new TH1F("h_endposx_st","Reconstructed End Position (x) After Stopping Tracks Selection",100,-1,257);
    TH1F* h_endposy_st = new TH1F("h_endposy_st","Reconstructed End Position (y) After Stopping Tracks Selection",100,-117,117);
    TH1F* h_endposz_st = new TH1F("h_endposz_st","Reconstructed End Position (z) After Stopping Tracks Selection",100,-1,1037);
    TH1F* h_reco_trklen_st = new TH1F("h_reco_trklen_st","Reconstructed Track Length After Stopping Tracks Selection",100,0,500);
    h_startposx_st -> SetLineColor(kMagenta);
    h_startposy_st -> SetLineColor(kMagenta);
    h_startposz_st -> SetLineColor(kMagenta);
    h_endposx_st -> SetLineColor(kMagenta);
    h_endposy_st -> SetLineColor(kMagenta);
    h_endposz_st -> SetLineColor(kMagenta);
    h_reco_trklen_st -> SetLineColor(kMagenta);
    
    //empty histo for storing muon selection result
    TH1F* h_startposx_mu = new TH1F("h_startposx_mu","Reconstructed Start Position (x) for Selected Muons",100,-1,257);
    TH1F* h_startposy_mu = new TH1F("h_startposy_mu","Reconstructed Start Position (y) for Selected Muons",100,-117,117);
    TH1F* h_startposz_mu = new TH1F("h_startposz_mu","Reconstructed Start Position (z) for Selected Muons",100,-1,1037);
    TH1F* h_endposx_mu = new TH1F("h_endposx_mu","Reconstructed End Position (x) for Selected Muons",100,-1,257);
    TH1F* h_endposy_mu = new TH1F("h_endposy_mu","Reconstructed End Position (y) for Selected Muons",100,-117,117);
    TH1F* h_endposz_mu = new TH1F("h_endposz_mu","Reconstructed End Position (z) for Selected Muons",100,-1,1037);
    TH1F* h_reco_trklen_mu = new TH1F("h_reco_trklen_mu","Reconstructed Track Length for Selected Muons",100,0,500);
    TH1F* h_trk_score_mu = new TH1F("h_trk_store_mu","Muon Track Score",100,0,1);
    TH1F* h_trk_pid_mu = new TH1F("h_trk_pid_mu","Muon Track PID",100,-1,1);
    TH1F* h_mu_trk_E = new TH1F("h_mu_trk_E","Muon Track Energy (y)",100,0,1000);
    TH2F* h2d_mu_mc_v_reco_E = new TH2F("h2d_mu_mc_v_reco_E","Muon MC-Reco Energy Correlation",30,200,600,30,20,540);
    h_startposx_mu -> SetLineColor(kAzure+3);
    h_startposy_mu -> SetLineColor(kAzure+3);
    h_startposz_mu -> SetLineColor(kAzure+3);
    h_endposx_mu -> SetLineColor(kAzure+3);
    h_endposy_mu -> SetLineColor(kAzure+3);
    h_endposz_mu -> SetLineColor(kAzure+3);
    h_reco_trklen_mu -> SetLineColor(kAzure+3);
    h_trk_score_mu -> SetLineColor(kRed);
    h_trk_pid_mu -> SetLineColor(kRed);
    h_mu_trk_E -> SetLineColor(kRed);
    
    //empty histos for storing related histograms for the muon light yield measurement
    TH1F* h_mu_flash_trk_dist = new TH1F("h_mu_flash_trk_dist","Distance Between the Muon Tracks and Flashes",100,0,350);
    TProfile* h_mu_ly_reco_x = new TProfile("h_mu_ly_reco_x","Total Light Yield Measurement Using Reconstructed MC Muons",10,-1,257,"");
    TProfile* h_mu_ly_mc_x = new TProfile("h_mu_ly_mc_x","Total Light Yield Measurement Using MC Muon (True)",10,-1,257,"");
    
    //empty histo for Muon-Michel gap distance histo
    TH1F* h_mu_michel_gap = new TH1F("h_mu_michel_gap","Muon-Michel electron Gap Distance",100,0,1000);
    
    //empty histo for storing selected Michel information
    TH1F* h_startposx_michel = new TH1F("h_startposx_michel","Michel Electron Start Position (x)",100,-1,257);
    TH1F* h_startposy_michel = new TH1F("h_startposy_michel","Michel Electron Start Position (y)",100,-117,117);
    TH1F* h_startposz_michel = new TH1F("h_startposz_michel","Michel Electron Start Position (z)",100,-1,1037);
    TH1F* h_reco_trklen_michel = new TH1F("h_reco_trklen_michel","Reconstructed Track Length for Selected Michel Electrons",100,0,100);
    TH1F* h_trk_score_michel = new TH1F("h_trk_store_michel","Michel Electron Track Score",100,0,1);
    TH1F* h_trk_pid_michel = new TH1F("h_trk_pid_michel","Michel Electron Track PID",100,-1,1);
    
    TH1F* h_michel_trk_E = new TH1F("h_michel_trk_E","Michel Electron Reco Track Energy (y)",100,0,300);
    h_michel_trk_E -> SetLineColor(kMagenta);
    
    TH2F* h2d_michel_mc_v_reco_E = new TH2F("h2d_michel_mc_v_reco_E","Michel Electron MC-Reco Energy Correlation",100,0,150,100,0,150);
    
    //empty histos for storing related histograms for the Michel light yield measurement
    TH1F* h_michel_flash_trk_dist = new TH1F("h_michel_flash_trk_dist","Distance Between the Michel Electron Tracks and Flashes",100,0,350);
    TProfile* h_michel_ly_reco_x = new TProfile("h_michel_ly_reco_x","Total Light Yield Measurement Using Reconstructed MC Michel Electrons",10,0,256.4,"");
    TProfile* h_michel_ly_mc_x = new TProfile("h_michel_ly_mc_x","Total Light Yield Measurement Using MC Michel Electrons (True)",10,0,256.4,"");
    
    //empty histo for shr_dedx vs the reconstructed track energy of the Michel
    TH2F* h2d_michel_shrdedx_v_reco_E = new TH2F("h2d_michel_shrdedx_v_reco_E","Michel Electron Shr_dedx vs Reconstructed Energy",100,0,150,0,15);
    TH2F* h2d_michel_shrdedxcali_v_reco_E = new TH2F("h2d_michel_shrdedxcali_v_reco_E","Michel Electron Shr_dedx_cali vs Reconstructed Energy",100,0,150,0,15);
    
    //Create a root file that stores all output plots
    unique_ptr<TFile> myFile(TFile::Open("sim_analyz_outputs.root", "RECREATE"));
    cout << "An empty .root file created for storing output plots." << endl;
    
    //now define a set of geometric variables that'll be used for the selection
    float low_edge_x = 0;
    float high_edge_x = 256.4;
    float bottom_edge_y = -116.5;
    float top_edge_y = 116.5;
    float low_edge_z = 0;
    float high_edge_z = 1036.8;
    
    float FV_low_x = 10;
    float FV_high_x = 246.4;
    float FV_low_z = 10;
    float FV_high_z = 1026.8;
    
    float mu_st_from_top_low_y = 110;
    float mu_st_from_top_high_y = 116.5;
    float mu_st_end_in_FV_y = -110;
    
    float michel_FV_low_y = -106.5;
    float michel_FV_high_y = 106.5;
    
    //create and open txt file for storing Michel electron events in the first tall bin in reco LY
    fstream file1;
    file1.open("1st_bin_reco_LY_Michel_info.txt", ios::app | ios::out);
    
    //create and open txt file for storing information of the Michel electron events that has an higher-than-expected reconstructed energy (>50 MeV)
    fstream greater_than_50;
    greater_than_50.open("michel_greater_than_50mev.txt", ios::app | ios::out);
    
    
   //getting all the event entries from the TTree
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   //the loop that loops through all of the event entries in the TTree
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
       for(unsigned int i=0;
              i<(*mc_endx).size(); i++)
       {
           //local helper variables
           //Truth information for the ~40k simulated muon events, and truth energy for Michel electrons.
           float mcvx, mcvy, mcvz, mcendx, mcendy, mcendz, mcpdg, mcE, me_mc_E;
           mcvx = (*mc_vx)[i];
           mcvy = (*mc_vy)[i];
           mcvz = (*mc_vz)[i];
           mcendx = (*mc_endx)[i];
           mcendy =(*mc_endy)[i];
           mcendz =(*mc_endz)[i];
           mcpdg = (*mc_pdg)[i];
           mcE = (*mc_E)[i] * 1000.; //convert to MeV
           me_mc_E = endmuonmichel * 1000.; //convert to MeV
           
           //fill the above values into the respective histograms.
           h_mc_vx -> Fill(mcvx);
           h_mc_vy -> Fill(mcvy);
           h_mc_vz -> Fill(mcvz);
           h_mc_endx -> Fill(mcendx);
           h_mc_endy -> Fill(mcendy);
           h_mc_endz -> Fill(mcendz);
           h_mc_pdg -> Fill(mcpdg);
           h_mc_E -> Fill(mcE);
           h_mc_michel_E -> Fill(me_mc_E);
       }
           
       //Define new vectors for storing muon position information later after the muon selection track loop
       vector<float> mu_endposx;
       vector<float> mu_endposy;
       vector<float> mu_endposz;
       
       //at this point before any selection cuts, we set a boolean parameter indicating whether good muons exist or not
       bool GoodMuExists = false;
       //create a vector containing index which refers to the selected muons
       vector<int> goodmuonindex;
       
       //loop through all the RECO tracks that's contained within each entry. This RECO track loop is also our muon selection track loop
       for(unsigned int itrack=0;
           itrack<(*trk_sce_start_x_v).size(); itrack++)
       {
           //define local variables inside the loop for computing and calculation usage.
           float stposx, stposy, stposz, endposx, endposy, endposz, trklen, trk_score, trk_pid, trk_recoE_y, mcvx, mcvy, mcvz, mcendx, mcendy, mcendz, mcE, mcpdg;
           stposx = (*trk_sce_start_x_v)[itrack];
           stposy = (*trk_sce_start_y_v)[itrack];
           stposz = (*trk_sce_start_z_v)[itrack];
           endposx = (*trk_sce_end_x_v)[itrack];
           endposy = (*trk_sce_end_y_v)[itrack];
           endposz = (*trk_sce_end_z_v)[itrack];
           trk_score = (*trk_score_v)[itrack];
           trk_pid = (*trk_llr_pid_score_v)[itrack];
           trk_recoE_y = (*trk_calo_energy_y_v)[itrack];
           mcvx = (*mc_vx)[itrack];
           mcvy = (*mc_vy)[itrack];
           mcvz = (*mc_vz)[itrack];
           mcendx = (*mc_endx)[itrack];
           mcendy = (*mc_endy)[itrack];
           mcendz = (*mc_endz)[itrack];
           mcpdg = (*mc_pdg)[itrack];
           mcE = (*mc_E)[itrack] * 1000.; //convert to MeV
           
           //Selection on the reconstruction variables (i.e. start/end position, recoE in this case for now) to make sure all our position values are reasonable and within the detector volume
           if(stposx>low_edge_x && stposy>bottom_edge_y && stposz>low_edge_z && endposx>low_edge_x && endposy>bottom_edge_y && endposz>low_edge_z && trk_recoE_y>-1)
           {
               //track length of each tracks can be calculated with default distance formula (made into a function)
               trklen = calc_len(stposx,stposy,stposz,endposx,endposy,endposz);
               
               //store in the track score and PID value before any proper selection
               h_trk_score -> Fill(trk_score);
               h_trk_pid -> Fill(trk_pid);
               
               //identify the proper endpoint (start/end indistinguishable at this stage even there are start/end variable)
               if(endposy>stposy) //if endpoint is higher than startpoint, flip them
               {
                   float tempx, tempy, tempz; //define temporary vars
                   tempx = endposx; endposx = stposx; stposx = tempx;
                   tempy = endposy; endposy = stposy; stposy = tempy;
                   tempz = endposz; endposz = stposz; stposz = tempz;
               } //now endpos is the proper endpoint
               
               //here fill non-selected reco start/end position with the respective reco track length and reco track energy
               h_startposx -> Fill(stposx);
               h_startposy -> Fill(stposy);
               h_startposz -> Fill(stposz);
               h_endposx -> Fill(endposx);
               h_endposy -> Fill(endposy);
               h_endposz -> Fill(endposz);
               h_reco_trk_len -> Fill(trklen);
               h_reco_trk_E -> Fill(trk_recoE_y);
               
               //fiducial volume selection (AV) (drift & beam dir)
               if(stposx>FV_low_x && stposx<FV_high_x && endposx>FV_low_x && endposx<FV_high_x && stposz>FV_low_z && stposz<FV_high_z && endposz>FV_low_z && endposz<FV_high_z && trk_recoE_y!=0)
               {
                   h_startposx_fv -> Fill(stposx);
                   h_startposy_fv -> Fill(stposy);
                   h_startposz_fv -> Fill(stposz);
                   h_endposx_fv -> Fill(endposx);
                   h_endposy_fv -> Fill(endposy);
                   h_endposz_fv -> Fill(endposz);
                   h_reco_trk_len_fv -> Fill(trklen);
               
                   //we would like to select only tracks coming into the detector from the top and stopping within the AV
                   if(stposy>mu_st_from_top_low_y && stposy<mu_st_from_top_high_y && endposy>mu_st_end_in_FV_y)
                   {
                       h_startposx_st -> Fill(stposx);
                       h_startposy_st -> Fill(stposy);
                       h_startposz_st -> Fill(stposz);
                       h_endposx_st -> Fill(endposx);
                       h_endposy_st -> Fill(endposy);
                       h_endposz_st -> Fill(endposz);
                       h_reco_trklen_st -> Fill(trklen);
                   
                       //Now, from these stopping tracks in the detector AV, identify the muon tracks
                       if(trk_score>0.5 && trk_pid>0 && trklen>10 && mcpdg!=0 && mcpdg<100 && mcE>0)
                       {
                           h_trk_score_mu -> Fill(trk_score);
                           h_trk_pid_mu -> Fill(trk_pid);
                           h_startposx_mu -> Fill(stposx);
                           h_startposy_mu -> Fill(stposy);
                           h_startposz_mu -> Fill(stposz);
                           h_endposx_mu -> Fill(endposx);
                           h_endposy_mu -> Fill(endposy);
                           h_endposz_mu -> Fill(endposz);
                           h_reco_trklen_mu -> Fill(trklen);
                           h_mu_trk_E -> Fill(trk_recoE_y);
                           h2d_mu_mc_v_reco_E -> Fill(mcE,trk_recoE_y);
                           
                           
            /*————————————————————————————————————————————————
              here starts the light yield measurement for the
              muons after we have selected the muon candidates,
              the michel electron light yield measurement will be
              using the exact same method.
             ————————————————————————————————————————————————*/
                           //here we first calculate the distance between the optically reconstructed flash and the muon track (midpt)
                           float mu_flash_trk_dist = calc_flash_evt_dist(flash_ycenter,flash_zcenter,stposy,endposy,stposz,endposz);
                           h_mu_flash_trk_dist -> Fill(mu_flash_trk_dist);
                           
                           //now, match the flash with its corresponding muon track, and the light yield measurement will be done within this condition
                           if (mu_flash_trk_dist<160. && calc_devi_flash_trk(flash_zcenter,flash_zwidth,stposz,endposz)<=1.)
                           {
                               //define a var for the total number of PEs
                               float tot_pe_per_trk = 0;
                               //now calculate this value by looping through all 32 pmt and add the PEs together
                               for (unsigned int PMT_no=0; PMT_no!=32; ++PMT_no)
                               {
                                   tot_pe_per_trk += flash_pe_v->at(PMT_no)*(20./gain_ampl_v->at(PMT_no));
                               }//end of all pmt loop for calculating the tot no. of PEs
                               
                               //now calculate the total LY with reconstructed and MC energy
                               if(tot_pe_per_trk>0 && trk_recoE_y>0 )
                               {
                                   float mu_reco_ly = tot_pe_per_trk/trk_recoE_y;
                                   float mu_mc_ly = tot_pe_per_trk/mcE;
                                   float trk_mid_x = (stposx+endposx)/2.;
                                   
                                   //fill LY plot for the muons
                                   h_mu_ly_reco_x -> Fill(trk_mid_x,mu_reco_ly);
                                   h_mu_ly_mc_x -> Fill(trk_mid_x,mu_mc_ly);
                                   
                               }//end of ensuring no division by 0
                            
                           }//end of muon flash-matching condition
                           
                           //here when all selections are satisfied, set goodmuexists to true and store the muon endpoints into pre-defined vector for Michel electron selection usage
                           GoodMuExists = true;
                           mu_endposx.push_back(endposx);
                           mu_endposy.push_back(endposy);
                           mu_endposz.push_back(endposz);
                           
                           goodmuonindex.push_back(itrack);
                           
                       }//end of muon identification condition
                       
                   } //end of stopping track selection condition
                   
               }//end of fiducial volume selection
               
           } //end of sensible value condition (reasonable position/energy values)
       
       } //end of muon RECO track loop
       /*——————————————————————————————————————————————————————————
         here starts the Michel electron selection loop (track loop)
        —————————————————————————————————————————————————————————*/
       
       if(GoodMuExists)
       {
           //Another loop which loops through all the reconstructed tracks for the Michel selection
           for(unsigned int itrack=0;
               itrack<(*trk_sce_start_x_v).size(); itrack++)
           {
               if(goodmuonindex[0]==itrack)
                   continue;
               
               //define local variables inside the loop for computing and calculation usage.
               float stposx, stposy, stposz, endposx, endposy, endposz, trklen, trk_score, trk_pid, trk_recoE_u, trk_recoE_v, trk_recoE_y, mcvx, mcvy, mcvz, mcendx, mcendy, mcendz, mcE, mc_michel_E, mcpdg, val_mu_endposx, val_mu_endposy, val_mu_endposz, shr_stx, shr_sty, shr_stz, shr_dedx, shr_dedx_cali;
               stposx = (*trk_sce_start_x_v)[itrack];
               stposy = (*trk_sce_start_y_v)[itrack];
               stposz = (*trk_sce_start_z_v)[itrack];
               endposx = (*trk_sce_end_x_v)[itrack];
               endposy = (*trk_sce_end_y_v)[itrack];
               endposz = (*trk_sce_end_z_v)[itrack];
               trk_score = (*trk_score_v)[itrack];
               trk_pid = (*trk_llr_pid_score_v)[itrack];
               trk_recoE_u = (*trk_calo_energy_u_v)[itrack];
               trk_recoE_v = (*trk_calo_energy_v_v)[itrack];
               trk_recoE_y = (*trk_calo_energy_y_v)[itrack];
               mcvx = (*mc_vx)[itrack];
               mcvy = (*mc_vy)[itrack];
               mcvz = (*mc_vz)[itrack];
               mcendx = (*mc_endx)[itrack];
               mcendy = (*mc_endy)[itrack];
               mcendz = (*mc_endz)[itrack];
               mcpdg = (*mc_pdg)[itrack];
               mcE = (*mc_E)[itrack] * 1000.; //convert to MeV
               mc_michel_E = endmuonmichel*1000.; //convert to MeV
               val_mu_endposx = mu_endposx[goodmuonindex[0]];
               val_mu_endposy = mu_endposy[goodmuonindex[0]];
               val_mu_endposz = mu_endposz[goodmuonindex[0]];
               shr_stx = shr_start_x;
               shr_sty = shr_start_y;
               shr_stz = shr_start_z;
               shr_dedx = shr_dedx_Y;
               shr_dedx_cali = shr_dedx_Y_cali;
               
               
               //Selection to make sure no glitching non-realistic values
               if(stposx>low_edge_x && stposy>bottom_edge_y && stposz>low_edge_z && endposx>low_edge_x && endposy>bottom_edge_y && endposz>low_edge_z)
               {
                   //track length of each tracks can be calculated with default distance formula
                   trklen = calc_len(stposx,stposy,stposz,endposx,endposy,endposz);
                   
                   
                   
                   //the fiducial volume selection, make sure the y-value are all within the detector (since ME are expected to be within the detector)
                   if(stposx>FV_low_x && stposx<FV_high_x && endposx>FV_low_x && endposx<FV_high_x && stposz>FV_low_z && stposz<FV_high_z && endposz>FV_low_z && endposz<FV_high_z && stposy<michel_FV_high_y && stposy>michel_FV_low_y && endposy>michel_FV_low_y && endposy<michel_FV_high_y)
                   {
                       
                       //create two var for storing the distance value between the muon endpoint we stored and the start/end positions. The shorter distance var will be the gap between reco muon track and reco Michel shower
                       float mu_michel_st_dist = calc_len(val_mu_endposx,val_mu_endposy,val_mu_endposz,stposx,stposy,stposz);
                       
                       
                       float mu_michel_end_dist = calc_len(val_mu_endposx,val_mu_endposy,val_mu_endposz,endposx,endposy,endposz);
                       
                       
                       
                       //compare the two dist and determine the correct gap between muon and michel
                       float mu_michel_gap = (mu_michel_end_dist < mu_michel_st_dist) ? mu_michel_end_dist : mu_michel_st_dist;
                       
                       
                       h_mu_michel_gap -> Fill(mu_michel_gap);
                       
                       //set a maximum gap limit as part of michel selection
                       if(mu_michel_gap<10 && trk_score<0.5)
                       {
                           h_startposx_michel -> Fill(stposx);
                           h_startposy_michel -> Fill(stposy);
                           h_startposz_michel -> Fill(stposz);
                           h_reco_trklen_michel -> Fill(trklen);
                           h_trk_score_michel -> Fill(trk_score);
                           h_trk_pid_michel -> Fill(trk_pid);
                           h_michel_trk_E -> Fill(trk_recoE_y);
                           if(mc_michel_E!=0)
                           {
                               h2d_michel_mc_v_reco_E -> Fill(mc_michel_E,trk_recoE_y);
                           }
                           
                           /*——————————————————————————————————————
                             here starts the light yield measurement for the
                             Michel electrons after we have selected the Michel candidates
                            —————————————————————————————————————*/
                           //similarly, first calculate the dist b/w flash and Michel track (midpt)
                           float michel_flash_trk_dist = calc_flash_evt_dist(flash_ycenter,flash_zcenter,stposy,endposy,stposz,endposz);
                           h_michel_flash_trk_dist -> Fill(michel_flash_trk_dist);
                           
                           //now, match the flash with its corresponding muon track, and the light yield measurement will be done within this condition
                           if (michel_flash_trk_dist<160. && calc_devi_flash_trk(flash_zcenter,flash_zwidth,stposz,endposz)<=1.)
                           {
                               //define a var for the total number of PEs of this track
                               float tot_pe_per_trk = 0;
                               //calculate and sum up all the PEs from all 32 PMTs
                               for (unsigned int PMT_no=0; PMT_no!=32; ++PMT_no)
                               {
                                   tot_pe_per_trk += flash_pe_v->at(PMT_no)*(20./gain_ampl_v->at(PMT_no));
                               }//end of all pmt loop for calculating the tot no. of PEs in the current track
                               
                               //now calculate the total LY with reconstructed and MC energy
                               if(tot_pe_per_trk>0 && trk_recoE_y>0 && mc_michel_E>0)
                               {
                                   float michel_reco_ly = tot_pe_per_trk/trk_recoE_y;
                                   float michel_mc_ly = tot_pe_per_trk/mc_michel_E;
                                   float trk_mid_x = (stposx+endposx)/2.;
                                                                      
                                   //fill LY plot for the michel electrons
                                   h_michel_ly_reco_x -> Fill(trk_mid_x,michel_reco_ly);
                                   h_michel_ly_mc_x -> Fill(trk_mid_x,michel_mc_ly);
                                   
                                   //here investigate the events in the first bin of h_michel_ly_reco_x
                                   if(trk_mid_x<25)
                                   {
                                       file1 << left << setw(10) << run << setw(10) << sub << setw(10) << evt << setw(15) << stposx << setw(15) << stposy << setw(15) << stposz << setw(15) << endposx << setw(15) << endposy << setw(15) << endposz << setw(15) << shr_stx << setw(15) << shr_sty << setw(15) << shr_stz << setw(15) << shr_dedx_cali << setw(15) << trk_score << setw(15) << mu_michel_gap << setw(15) << trk_recoE_u << setw(15) << trk_recoE_v << setw(15) << trk_recoE_y << endl;
                                   }
                                   
                                   //here fill the txt file with info of Michel events with energy >50 MeV
                                   if(trk_recoE_y>50)
                                   {
                                       greater_than_50 << left << setw(10) << run << setw(10) << sub << setw(10) << evt << setw(15) << stposx << setw(15) << stposy << setw(15) << stposz << setw(15) << endposx << setw(15) << endposy << setw(15) << endposz << setw(15) << shr_stx << setw(15) << shr_sty << setw(15) << shr_stz << setw(15) << shr_dedx_cali << setw(15) << trk_score << setw(15) << mu_michel_gap << setw(10) << trk_recoE_u << setw(10) << trk_recoE_v << setw(10) << trk_recoE_y << endl;
                                   }
                                   
                                   h2d_michel_shrdedx_v_reco_E -> Fill(trk_recoE_y,shr_dedx);
                                   h2d_michel_shrdedxcali_v_reco_E -> Fill(trk_recoE_y,shr_dedx_cali);
                                
                               }//end of ensuring no division by 0
                               
                           }//end of flash-matching condition
                           
                       }//end of max gap limit condition
                       
                   }//End of Michel Fiducial Volume selection
                   
               }//end of sensible value condition (reasonable position/energy values)
               
           }//end of Michel electron RECO track loop
           
       }//end of goodmuexist condition
       
   } //end of the entry loop that loops through all event entries
    
    //here, linearly fit the 2d muon mc-reco energy histogram
    h2d_mu_mc_v_reco_E -> Fit("pol1","","",200,600);
    
    //closing particle information .txt files
        file1.close();
        greater_than_50.close();
    
    //Drawing and formatting plots/histograms starts here
    TCanvas* c1 = new TCanvas("c1");
    c1 -> cd();
    h_mc_vx -> Draw();
    h_mc_vx -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_mc_vx -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c1 -> Print("plots.ps(");
    myFile -> WriteObject(c1, "Monte Carlo True Start Position (x)");
    
    TCanvas* c2 = new TCanvas("c2");
    c2 -> cd();
    h_mc_vy -> Draw();
    h_mc_vy -> GetXaxis() -> SetTitle("y-Position (cm)");
    h_mc_vy -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c2 -> Print("plots.ps");
    myFile -> WriteObject(c2, "Monte Carlo True Start Position (y)");
    
    TCanvas* c3 = new TCanvas("c3");
    c3 -> cd();
    h_mc_vz -> Draw();
    h_mc_vz -> GetXaxis() -> SetTitle("z-Position (cm)");
    h_mc_vz -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c3 -> Print("plots.ps");
    myFile -> WriteObject(c3, "Monte Carlo True Start Position (z)");
    
    TCanvas* c4 = new TCanvas("c4");
    c4 -> cd();
    h_mc_endx -> Draw();
    h_mc_endx -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_mc_endx -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c4 -> Print("plots.ps");
    myFile -> WriteObject(c4, "Monte Carlo True End Position (x)");
    
    TCanvas* c5 = new TCanvas("c5");
    c5 -> cd();
    h_mc_endy -> Draw();
    h_mc_endy -> GetXaxis() -> SetTitle("y-Position (cm)");
    h_mc_endy -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c5 -> Print("plots.ps");
    myFile -> WriteObject(c5, "Monte Carlo True End Position (y)");
    
    TCanvas* c6 = new TCanvas("c6");
    c6 -> cd();
    h_mc_endz -> Draw();
    h_mc_endz -> GetXaxis() -> SetTitle("z-Position (cm)");
    h_mc_endz -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c6 -> Print("plots.ps");
    myFile -> WriteObject(c6, "Monte Carlo True End Position (z)");
    
    TCanvas* c7 = new TCanvas("c7");
    c7 -> cd();
    h_mc_pdg -> Draw();
    h_mc_pdg -> GetXaxis()-> SetTitle("MC PDG Code (AU)");
    h_mc_pdg -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c7 -> Print("plots.ps");
    myFile -> WriteObject(c7, "Monte Carlo PDG Code");
    
    TCanvas* c8 = new TCanvas("c8");
    c8 -> cd();
    h_mc_E -> Draw();
    h_mc_E -> GetXaxis()-> SetTitle("Energy (MeV)");
    h_mc_E -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c8 -> Print("plots.ps");
    myFile -> WriteObject(c8, "Monte Carlo Energy");
    
    TCanvas* c9 = new TCanvas("c9");
    c9 -> cd();
    h_mc_michel_E -> Draw();
    h_mc_michel_E -> GetXaxis()-> SetTitle("Energy (MeV)");
    h_mc_michel_E -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c9 -> Print("plots.ps");
    myFile -> WriteObject(c9, "Monte Carlo Michel Electron Energy");
    
    TCanvas* c10 = new TCanvas("c10");
    c10 -> cd();
    h_startposx -> Draw();
    h_startposx_fv -> Draw("same");
    h_startposx_st -> Draw("same");
    h_startposx_mu -> Draw("same");
    h_startposx -> GetXaxis()-> SetTitle("x-Position (cm)");
    h_startposx -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L1 = new TLegend(0.55,0.8,0.9,0.9);
    L1 -> SetBorderSize(0);
    L1 -> SetFillStyle(0);
    L1 -> AddEntry(h_startposx,"No Selection Applied","lep");
    L1 -> AddEntry(h_startposx_fv,"With Fiducial Volume Selection","lep");
    L1 -> AddEntry(h_startposx_st,"With Stopping Track Selection","lep");
    L1 -> AddEntry(h_startposx_mu,"For the Selected Muons","lep");
    L1 -> SetTextFont(62);
    L1 -> Draw();
    //print to .root file and pdf
    c10 -> Print("plots.ps");
    myFile -> WriteObject(c10, "Reco Track Start Position (x)");
    
    TCanvas* c11 = new TCanvas("c11");
    c11 -> cd();
    h_startposy -> Draw();
    h_startposy_fv -> Draw("same");
    h_startposy_st -> Draw("same");
    h_startposy_mu -> Draw("same");
    h_startposy -> GetXaxis()-> SetTitle("y-Position (cm)");
    h_startposy -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L2 = new TLegend(0.55,0.8,0.9,0.9);
    L2 -> SetBorderSize(0);
    L2 -> SetFillStyle(0);
    L2 -> AddEntry(h_startposy,"No Selection Applied","lep");
    L2 -> AddEntry(h_startposy_fv,"With Fiducial Volume Selection","lep");
    L2 -> AddEntry(h_startposy_st,"With Stopping Track Selection","lep");
    L2 -> AddEntry(h_startposy_mu,"For the Selected Muons","lep");
    L2 -> SetTextFont(62);
    L2 -> Draw();
    //print to .root file and pdf
    c11 -> Print("plots.ps");
    myFile -> WriteObject(c11, "Reco Track Start Position (y)");
    
    TCanvas* c12 = new TCanvas("c12");
    c12 -> cd();
    h_startposz -> Draw();
    h_startposz_fv -> Draw("same");
    h_startposz_st -> Draw("same");
    h_startposz_mu -> Draw("same");
    gStyle -> SetOptStat(10);
    h_startposz -> GetXaxis()-> SetTitle("z-Position (cm)");
    h_startposz -> GetYaxis() -> SetTitle("Number of Events");
    //legend
    TLegend* L3 = new TLegend(0.55,0.8,0.9,0.9);
    L3 -> SetBorderSize(0);
    L3 -> SetFillStyle(0);
    L3 -> AddEntry(h_startposz,"No Selection Applied","lep");
    L3 -> AddEntry(h_startposz_fv,"With Fiducial Volume Selection","lep");
    L3 -> AddEntry(h_startposz_st,"With Stopping Track Selection","lep");
    L3 -> AddEntry(h_startposz_mu,"For the Selected Muons","lep");
    L3 -> SetTextFont(62);
    L3 -> Draw();
    //print to .root file and pdf
    c12 -> Print("plots.ps");
    myFile -> WriteObject(c12, "Reco Track Start Position (z)");
    
    TCanvas* c13 = new TCanvas("c13");
    c13 -> cd();
    h_endposx -> Draw();
    h_endposx_fv -> Draw("same");
    h_endposx_st -> Draw("same");
    h_endposx_mu -> Draw("same");
    h_endposx -> GetXaxis()-> SetTitle("x-Position (cm)");
    h_endposx -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L4 = new TLegend(0.55,0.8,0.9,0.9);
    L4 -> SetBorderSize(0);
    L4 -> SetFillStyle(0);
    L4 -> AddEntry(h_endposx,"No Selection Applied","lep");
    L4 -> AddEntry(h_endposx_fv,"With Fiducial Volume Selection","lep");
    L4 -> AddEntry(h_endposx_st,"With Stopping Track Selection","lep");
    L4 -> AddEntry(h_endposx_mu,"For the Selected Muons","lep");
    L4 -> SetTextFont(62);
    L4 -> Draw();
    //print to .root file and pdf
    c13 -> Print("plots.ps");
    myFile -> WriteObject(c13, "Reco Track End Position (x)");
    
    TCanvas* c14 = new TCanvas("c14");
    c14 -> cd();
    h_endposy -> Draw();
    h_endposy_fv -> Draw("same");
    h_endposy_st -> Draw("same");
    h_endposy_mu -> Draw("same");
    h_endposy -> GetXaxis()-> SetTitle("y-Position (cm)");
    h_endposy -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L5 = new TLegend(0.55,0.8,0.9,0.9);
    L5 -> SetBorderSize(0);
    L5 -> SetFillStyle(0);
    L5 -> AddEntry(h_endposy,"No Selection Applied","lep");
    L5 -> AddEntry(h_endposy_fv,"With Fiducial Volume Selection","lep");
    L5 -> AddEntry(h_endposy_st,"With Stopping Track Selection","lep");
    L5 -> AddEntry(h_endposy_mu,"For the Selected Muons","lep");
    L5 -> SetTextFont(62);
    L5 -> Draw();
    //print to .root file and pdf
    c14 -> Print("plots.ps");
    myFile -> WriteObject(c14, "Reco Track End Position (y)");
    
    TCanvas* c15 = new TCanvas("c15");
    c15 -> cd();
    h_endposz -> Draw();
    h_endposz_fv -> Draw("same");
    h_endposz_st -> Draw("same");
    h_endposz_mu -> Draw("same");
    h_endposz -> GetXaxis()-> SetTitle("z-Position (cm)");
    h_endposz -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L6 = new TLegend(0.55,0.8,0.9,0.9);
    L6 -> SetBorderSize(0);
    L6 -> SetFillStyle(0);
    L6 -> AddEntry(h_endposz,"No Selection Applied","lep");
    L6 -> AddEntry(h_endposz_fv,"With Fiducial Volume Selection","lep");
    L6 -> AddEntry(h_endposz_st,"With Stopping Track Selection","lep");
    L6 -> AddEntry(h_endposz_mu,"For the Selected Muons","lep");
    L6 -> SetTextFont(62);
    L6 -> Draw();
    //print to .root file and pdf
    c15 -> Print("plots.ps");
    myFile -> WriteObject(c15, "Reco Track End Position (z)");
    
    TCanvas* c16 = new TCanvas("c16");
    c16 -> cd();
    h_reco_trk_len -> Draw();
    h_reco_trk_len_fv -> Draw("same");
    h_reco_trklen_st -> Draw("same");
    h_reco_trklen_mu -> Draw("same");
    h_reco_trk_len -> GetXaxis()-> SetTitle("Track Length (cm)");
    h_reco_trk_len -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L7 = new TLegend(0.55,0.8,0.9,0.9);
    L7 -> SetBorderSize(0);
    L7 -> SetFillStyle(0);
    L7 -> AddEntry(h_reco_trk_len,"No Selection Applied","lep");
    L7 -> AddEntry(h_reco_trk_len_fv,"With Fiducial Volume Selection","lep");
    L7 -> AddEntry(h_reco_trklen_st,"With Stopping Track Selection","lep");
    L7 -> AddEntry(h_reco_trklen_mu,"For the Selected Muons","lep");
    L7 -> SetTextFont(62);
    L7 -> Draw();
    //print to .root file and pdf
    c16 -> Print("plots.ps");
    myFile -> WriteObject(c16, "Reco Track Length");
    
    TCanvas* c17 = new TCanvas("c17");
    c17 -> cd();
    h_reco_trk_E -> Draw();
    h_mu_trk_E -> Draw("same");
    h_michel_trk_E -> Draw("same");
    h_reco_trk_E -> GetXaxis()-> SetTitle("Energy (MeV)");
    h_reco_trk_E -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L8 = new TLegend(0.55,0.8,0.9,0.9);
    L8 -> SetBorderSize(0);
    L8 -> SetFillStyle(0);
    L8 -> AddEntry(h_reco_trk_E,"No Selection Applied","lep");
    L8 -> AddEntry(h_mu_trk_E,"For Selected Muons","lep");
    L8 -> AddEntry(h_michel_trk_E,"For Selected Michel Electrons","lep");
    L8 -> SetTextFont(62);
    L8 -> Draw();
    //print to .root file and pdf
    c17 -> Print("plots.ps");
    myFile -> WriteObject(c17, "Reco Track Energy");
    
    TCanvas* c18 = new TCanvas("c18");
    c18 -> cd();
    h_trk_score -> Draw();
    h_trk_score_mu -> Draw("same");
    h_trk_score_michel ->Draw("same");
    h_trk_score -> GetXaxis() -> SetTitle("Track Score (AU)");
    h_trk_score -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L9 = new TLegend(0.55,0.8,0.9,0.9);
    L9 -> SetBorderSize(0);
    L9 -> SetFillStyle(0);
    L9 -> AddEntry(h_trk_score,"No Selection Applied","lep");
    L9 -> AddEntry(h_trk_score_mu,"For Selected Muons","lep");
    L9 -> AddEntry(h_trk_score_michel,"For Selected Michel Electrons","lep");
    L9 -> SetTextFont(62);
    L9 -> Draw();
    //print to .root file and pdf
    c18 -> Print("plots.ps");
    myFile -> WriteObject(c18, "Track Score");
    
    TCanvas* c19 = new TCanvas("c19");
    c19 -> cd();
    h_trk_pid -> Draw();
    h_trk_pid_mu -> Draw("same");
    h_trk_pid_michel ->Draw("same");
    h_trk_pid -> GetXaxis() -> SetTitle("Track PID (AU)");
    h_trk_pid -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //legend
    TLegend* L10 = new TLegend(0.55,0.8,0.9,0.9);
    L10 -> SetBorderSize(0);
    L10 -> SetFillStyle(0);
    L10 -> AddEntry(h_trk_pid,"No Selection Applied","lep");
    L10 -> AddEntry(h_trk_pid_mu,"For Selected Muons","lep");
    L10 -> AddEntry(h_trk_pid_michel,"For Selected Michel Electrons","lep");
    L10 -> SetTextFont(62);
    L10 -> Draw();
    //print to .root file and pdf
    c19 -> Print("plots.ps");
    myFile -> WriteObject(c19, "Track Score");
    
    TCanvas* c20 = new TCanvas("c20");
    c20 -> cd();
    h2d_mu_mc_v_reco_E -> Draw("colz");
    h2d_mu_mc_v_reco_E -> GetXaxis() -> SetTitle("MC Energy (MeV)");
    h2d_mu_mc_v_reco_E -> GetYaxis() -> SetTitle("Reco Track Energy (MeV)");
    h2d_mu_mc_v_reco_E -> GetZaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c20 -> Print("plots.ps");
    myFile -> WriteObject(c20, "Muon MC vs. Reco Track Energy");
    
    TCanvas* c21 = new TCanvas("c21");
    c21 -> cd();
    h2d_michel_mc_v_reco_E -> Draw("colz");
    h2d_michel_mc_v_reco_E -> GetXaxis() -> SetTitle("MC Energy (MeV)");
    h2d_michel_mc_v_reco_E -> GetYaxis() -> SetTitle("Reco Track Energy (MeV)");
    h2d_michel_mc_v_reco_E -> GetZaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c21 -> Print("plots.ps");
    myFile -> WriteObject(c21, "Michel Electron MC vs. Reco Track Energy");
    
    TCanvas* c22 = new TCanvas("c22");
    c22 -> cd();
    h_mu_michel_gap -> Draw();
    h_mu_michel_gap -> GetXaxis() -> SetTitle("Gap Distance (cm)");
    h_mu_michel_gap -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c22 -> Print("plots.ps");
    myFile -> WriteObject(c22, "Muon-Michel Gap Distance");
    
    TCanvas* c23 = new TCanvas("c23");
    c23 -> cd();
    h_mu_flash_trk_dist -> Draw();
    h_mu_flash_trk_dist -> GetXaxis() -> SetTitle("Distance (cm)");
    h_mu_flash_trk_dist -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c23 -> Print("plots.ps");
    myFile -> WriteObject(c23, "Muon Flash-Track Distance");
    
    TCanvas* c24 = new TCanvas("c24");
    c24 -> cd();
    h_michel_flash_trk_dist -> Draw();
    h_michel_flash_trk_dist -> GetXaxis() -> SetTitle("Distance (cm)");
    h_michel_flash_trk_dist -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c24 -> Print("plots.ps");
    myFile -> WriteObject(c24, "Michel Flash-Track Distance");
    
    TCanvas* c25 = new TCanvas("c25");
    c25 -> cd();
    h_mu_ly_reco_x -> Draw("E1");
    h_mu_ly_reco_x -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_mu_ly_reco_x -> GetYaxis() -> SetTitle("Number of PE/MeV");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c25 -> Print("plots.ps");
    myFile -> WriteObject(c25, "Muon Reconstructed Light Yield Measurement (x)");
    
    TCanvas* c26 = new TCanvas("c26");
    c26 -> cd();
    h_mu_ly_mc_x -> Draw("E1");
    h_mu_ly_mc_x -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_mu_ly_mc_x -> GetYaxis() -> SetTitle("Number of PE/MeV");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c26 -> Print("plots.ps");
    myFile -> WriteObject(c26, "Muon MC Light Yield Measurement (x)");
    
    TCanvas* c27 = new TCanvas("c27");
    c27 -> cd();
    h_michel_ly_reco_x -> Draw("E1");
    h_michel_ly_reco_x -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_michel_ly_reco_x -> GetYaxis() -> SetTitle("Number of PE/MeV");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c27 -> Print("plots.ps");
    myFile -> WriteObject(c27, "Michel Reconstructed Light Yield Measurement (x)");
    
    TCanvas* c28 = new TCanvas("c28");
    c28 -> cd();
    h_michel_ly_mc_x -> Draw("E1");
    h_michel_ly_mc_x -> GetXaxis() -> SetTitle("x-Position (cm)");
    h_michel_ly_mc_x -> GetYaxis() -> SetTitle("Number of PE/MeV");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c28-> Print("plots.ps");
    myFile -> WriteObject(c28, "Michel MC Light Yield Measurement (x)");
    
    TCanvas* c29 = new TCanvas("c29");
    c29 -> cd();
    h_michel_trk_E -> Draw();
    h_michel_trk_E-> GetXaxis()-> SetTitle("Energy (MeV)");
    h_michel_trk_E -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c29-> Print("plots.ps");
    myFile -> WriteObject(c29, "Michel Reco Track Energy");
    
    TCanvas* c30 = new TCanvas("c30");
    c30 -> cd();
    h_reco_trklen_michel -> Draw();
    h_reco_trklen_michel-> GetXaxis()-> SetTitle("Track Length (cm)");
    h_reco_trklen_michel -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c30-> Print("plots.ps");
    myFile -> WriteObject(c30, "Michel Reco Track Length");
    
    TCanvas* c31 = new TCanvas("c31");
    c31 -> cd();
    h_mu_trk_E -> Draw();
    h_mu_trk_E-> GetXaxis()-> SetTitle("Energy (MeV)");
    h_mu_trk_E -> GetYaxis() -> SetTitle("Number of Events");
    gStyle -> SetOptStat(10);
    //print to .root file and pdf
    c31 -> Print("plots.ps");
    myFile -> WriteObject(c31, "Muon Reco Track Energy");
    
    TCanvas* c32 = new TCanvas("c32");
    c32 -> cd();
    h2d_michel_shrdedx_v_reco_E -> Draw("colz");
    
    TCanvas* c33 = new TCanvas("c33");
    c33 -> cd();
    h2d_michel_shrdedxcali_v_reco_E -> Draw("colz");
    
} //end of the selection analyzer class
