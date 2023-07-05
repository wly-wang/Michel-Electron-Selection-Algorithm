//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 24 14:38:39 2023 by ROOT version 6.26/06
// from TTree NeutrinoSelectionFilter/Neutrino Selection TTree
// found on file: ub_mu_sim_40k.root
//////////////////////////////////////////////////////////

#ifndef mu_sim_mely_class_h
#define mu_sim_mely_class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "string"
#include "vector"
#include "vector"

class mu_sim_mely_class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           selected;
   Int_t           run;
   Int_t           sub;
   Int_t           evt;
   UInt_t          trk_id;
   UInt_t          shr_id;
   UInt_t          trk2_id;
   UInt_t          shr2_id;
   Float_t         shr_energy_tot;
   Float_t         shr_energy;
   Float_t         shr_energy_tot_cali;
   Float_t         shr_energy_cali;
   Float_t         shr_theta;
   Float_t         shr_phi;
   Float_t         shr_pca_0;
   Float_t         shr_pca_1;
   Float_t         shr_pca_2;
   Float_t         shr_px;
   Float_t         shr_py;
   Float_t         shr_pz;
   Float_t         shr_openangle;
   Float_t         shr_tkfit_start_x;
   Float_t         shr_tkfit_start_y;
   Float_t         shr_tkfit_start_z;
   Float_t         shr_tkfit_theta;
   Float_t         shr_tkfit_phi;
   Float_t         shr_start_x;
   Float_t         shr_start_y;
   Float_t         shr_start_z;
   Float_t         shr_dedx_Y;
   Float_t         shr_dedx_V;
   Float_t         shr_dedx_U;
   Float_t         shr_dedx_Y_cali;
   Float_t         shr_dedx_V_cali;
   Float_t         shr_dedx_U_cali;
   Float_t         shr_tkfit_dedx_Y;
   Float_t         shr_tkfit_dedx_V;
   Float_t         shr_tkfit_dedx_U;
   Float_t         shr_tkfit_dedx_max;
   UInt_t          shr_tkfit_nhits_Y;
   UInt_t          shr_tkfit_nhits_V;
   UInt_t          shr_tkfit_nhits_U;
   Float_t         shr_llrpid_dedx_Y;
   Float_t         shr_llrpid_dedx_V;
   Float_t         shr_llrpid_dedx_U;
   Float_t         shr_llrpid_dedx;
   Float_t         shr_tkfit_dedx_Y_alt;
   Float_t         shr_tkfit_dedx_V_alt;
   Float_t         shr_tkfit_dedx_U_alt;
   UInt_t          shr_tkfit_nhits_Y_alt;
   UInt_t          shr_tkfit_nhits_V_alt;
   UInt_t          shr_tkfit_nhits_U_alt;
   Float_t         trkfit;
   UInt_t          shr_tkfit_npoints;
   UInt_t          shr_tkfit_npointsvalid;
   Float_t         shr_trkfitmedangle;
   Float_t         shrmoliereavg;
   Float_t         shrmoliererms;
   Float_t         shr1shr2moliereavg;
   Float_t         shr1shr2moliererms;
   Float_t         shr1trk1moliereavg;
   Float_t         shr1trk1moliererms;
   Float_t         shr1trk2moliereavg;
   Float_t         shr1trk2moliererms;
   UChar_t         ismerged;
   Float_t         merge_bestdot;
   Float_t         merge_bestdist;
   Float_t         merge_vtx_x;
   Float_t         merge_vtx_y;
   Float_t         merge_vtx_z;
   UInt_t          merge_tk_ipfp;
   Float_t         shr_tkfit_2cm_dedx_Y;
   Float_t         shr_tkfit_2cm_dedx_V;
   Float_t         shr_tkfit_2cm_dedx_U;
   UInt_t          shr_tkfit_2cm_nhits_Y;
   UInt_t          shr_tkfit_2cm_nhits_V;
   UInt_t          shr_tkfit_2cm_nhits_U;
   Float_t         shr_tkfit_gap05_dedx_Y;
   Float_t         shr_tkfit_gap05_dedx_V;
   Float_t         shr_tkfit_gap05_dedx_U;
   UInt_t          shr_tkfit_gap05_nhits_Y;
   UInt_t          shr_tkfit_gap05_nhits_V;
   UInt_t          shr_tkfit_gap05_nhits_U;
   Float_t         shr_tkfit_gap10_dedx_Y;
   Float_t         shr_tkfit_gap10_dedx_V;
   Float_t         shr_tkfit_gap10_dedx_U;
   UInt_t          shr_tkfit_gap10_nhits_Y;
   UInt_t          shr_tkfit_gap10_nhits_V;
   UInt_t          shr_tkfit_gap10_nhits_U;
   Float_t         shr_chipr;
   Float_t         shr_chimu;
   Float_t         shr_bragg_p;
   Float_t         shr_bragg_mu;
   Float_t         shr_bragg_mip;
   Float_t         shr_bragg_kaon;
   Float_t         shr_bragg_pion;
   Float_t         tksh_distance;
   Float_t         tksh_angle;
   Float_t         shr_distance;
   Float_t         shr_score;
   Int_t           shr_bkt_pdg;
   Float_t         shr_bkt_purity;
   Float_t         shr_bkt_completeness;
   Float_t         shr_bkt_E;
   Float_t         trk_len;
   Float_t         trk_theta;
   Float_t         trk_phi;
   Float_t         trk_energy;
   Float_t         trk_energy_muon;
   Float_t         trk_energy_muon_mcs;
   Float_t         trk_energy_tot;
   Float_t         trk_energy_muon_tot;
   Float_t         trk_distance;
   Float_t         trk_score;
   Int_t           trk_bkt_pdg;
   Float_t         trk_bkt_purity;
   Float_t         trk_bkt_completeness;
   Float_t         trk_bkt_E;
   Float_t         trk_chipr_best;
   Float_t         trk_chipr_worst;
   Float_t         trk_chimu_best;
   Float_t         trk_chimu_worst;
   Float_t         trk_chipr;
   Float_t         trk_chimu;
   Float_t         trk_pida;
   Float_t         trk_bragg_p;
   Float_t         trk_bragg_mu;
   Float_t         trk_bragg_mip;
   Float_t         trk_bragg_kaon;
   Float_t         trk_bragg_pion;
   UInt_t          trk_hits_max;
   UInt_t          shr_hits_max;
   vector<int>     *all_shr_hits;
   vector<int>     *all_trk_hits;
   vector<float>   *all_shr_energies;
   vector<float>   *all_trk_energies;
   //UInt_t          shr_hits_max;
   UInt_t          trk_hits_2nd;
   UInt_t          shr_hits_2nd;
   Float_t         trkshrhitdist0;
   Float_t         trkshrhitdist1;
   Float_t         trkshrhitdist2;
   Float_t         trk2shrhitdist0;
   Float_t         trk2shrhitdist1;
   Float_t         trk2shrhitdist2;
   Float_t         trk1trk2hitdist0;
   Float_t         trk1trk2hitdist1;
   Float_t         trk1trk2hitdist2;
   UInt_t          total_hits_y;
   Float_t         extra_energy_y;
   Float_t         trk_energy_hits_tot;
   UInt_t          subcluster;
   UInt_t          shrsubclusters0;
   UInt_t          shrsubclusters1;
   UInt_t          shrsubclusters2;
   Float_t         shrclusfrac0;
   Float_t         shrclusfrac1;
   Float_t         shrclusfrac2;
   Float_t         shrclusdir0;
   Float_t         shrclusdir1;
   Float_t         shrclusdir2;
   UInt_t          shr_hits_tot;
   UInt_t          shr_hits_y_tot;
   UInt_t          shr_hits_u_tot;
   UInt_t          shr_hits_v_tot;
   UInt_t          trk_hits_tot;
   UInt_t          trk_hits_y_tot;
   UInt_t          trk_hits_u_tot;
   UInt_t          trk_hits_v_tot;
   Float_t         _elecclusters_U_charge;
   Float_t         _elecclusters_V_charge;
   Float_t         _elecclusters_Y_charge;
   Int_t           _elecclusters_U_N;
   Int_t           _elecclusters_V_N;
   Int_t           _elecclusters_Y_N;
   UInt_t          n_tracks_contained;
   UInt_t          n_showers_contained;
   Float_t         matched_E;
   Float_t         hits_ratio;
   Float_t         contained_fraction;
   Float_t         sps_contained_fraction;
   Float_t         pt;
   Float_t         p;
   Float_t         pt_assume_muon;
   Float_t         p_assume_muon;
   Float_t         reco_e;
   Float_t         dvtx;
   Float_t         dtrk;
   Float_t         contained_sps_ratio;
   vector<double>  *dtrk_x_boundary;
   vector<double>  *dtrk_y_boundary;
   vector<double>  *dtrk_z_boundary;
   vector<double>  *dshr_x_boundary;
   vector<double>  *dshr_y_boundary;
   vector<double>  *dshr_z_boundary;
   vector<double>  *dvtx_x_boundary;
   vector<double>  *dvtx_y_boundary;
   vector<double>  *dvtx_z_boundary;
   vector<vector<double> > *dtrk_boundary;
   vector<vector<double> > *dvtx_boundary;
   vector<vector<double> > *dshr_boundary;
   vector<vector<double> > *dmc_boundary;
   Float_t         leeweight;
   Float_t         true_pt;
   Float_t         true_pt_visible;
   Float_t         true_p;
   Float_t         true_p_visible;
   Float_t         true_e_visible;
   Float_t         _opfilter_pe_beam;
   Float_t         _opfilter_pe_veto;
   Int_t           nu_pdg;
   Int_t           ccnc;
   Int_t           nu_parent_pdg;
   Int_t           nu_hadron_pdg;
   Int_t           nu_decay_mode;
   Int_t           interaction;
   Float_t         nu_e;
   Float_t         nu_l;
   Float_t         nu_pt;
   Float_t         theta;
   Bool_t          isVtxInFiducial;
   Bool_t          truthFiducial;
   Float_t         true_nu_vtx_t;
   Float_t         true_nu_vtx_x;
   Float_t         true_nu_vtx_y;
   Float_t         true_nu_vtx_z;
   Float_t         true_nu_vtx_sce_x;
   Float_t         true_nu_vtx_sce_y;
   Float_t         true_nu_vtx_sce_z;
   Float_t         reco_nu_vtx_x;
   Float_t         reco_nu_vtx_y;
   Float_t         reco_nu_vtx_z;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           nmuon;
   Float_t         muon_e;
   Float_t         muon_c;
   Float_t         muon_p;
   Int_t           nelec;
   Float_t         elec_e;
   Float_t         elec_c;
   Float_t         elec_p;
   Float_t         elec_vx;
   Float_t         elec_vy;
   Float_t         elec_vz;
   Float_t         elec_px;
   Float_t         elec_py;
   Float_t         elec_pz;
   Int_t           npi0;
   Float_t         pi0_e;
   Float_t         pi0_c;
   Float_t         pi0_p;
   Int_t           nneutron;
   Int_t           nproton;
   Float_t         proton_e;
   Float_t         proton_c;
   Float_t         proton_p;
   Int_t           npion;
   Float_t         pion_e;
   Float_t         pion_c;
   Float_t         pion_p;
   Int_t           neta;
   Float_t         eta_e;
   Int_t           nslice;
   Int_t           crtveto;
   Float_t         crthitpe;
   vector<int>     *pfp_slice_idx;
   Int_t           category;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<int>     *backtracked_tid;
   vector<float>   *backtracked_purity;
   vector<float>   *backtracked_completeness;
   vector<float>   *backtracked_overlay_purity;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *backtracked_start_x;
   vector<float>   *backtracked_start_y;
   vector<float>   *backtracked_start_z;
   vector<float>   *backtracked_start_t;
   vector<float>   *backtracked_start_U;
   vector<float>   *backtracked_start_V;
   vector<float>   *backtracked_start_Y;
   vector<float>   *backtracked_sce_start_x;
   vector<float>   *backtracked_sce_start_y;
   vector<float>   *backtracked_sce_start_z;
   vector<float>   *backtracked_sce_start_U;
   vector<float>   *backtracked_sce_start_V;
   vector<float>   *backtracked_sce_start_Y;
   Float_t         lep_e;
   Int_t           pass;
   Int_t           swtrig;
   Int_t           evnhits;
   Int_t           slpdg;
   Int_t           slnhits;
   Int_t           n_pfps;
   Int_t           n_tracks;
   Int_t           n_showers;
   vector<unsigned int> *pfp_generation_v;
   vector<unsigned int> *pfp_trk_daughters_v;
   vector<unsigned int> *pfp_shr_daughters_v;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   vector<int>     *pfpplanesubclusters_U;
   vector<int>     *pfpplanesubclusters_V;
   vector<int>     *pfpplanesubclusters_Y;
   vector<float>   *pfpplanesubhitfracmax_U;
   vector<float>   *pfpplanesubhitfracmax_V;
   vector<float>   *pfpplanesubhitfracmax_Y;
   UInt_t          hits_u;
   UInt_t          hits_v;
   UInt_t          hits_y;
   Float_t         topological_score;
   Float_t         slclustfrac;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<float>   *mc_vx;
   vector<float>   *mc_vy;
   vector<float>   *mc_vz;
   vector<float>   *mc_endx;
   vector<float>   *mc_endy;
   vector<float>   *mc_endz;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_completeness;
   vector<float>   *mc_purity;
   string          *endmuonprocess;
   Float_t         endmuonmichel;
   Float_t         flash_pe;
   vector<float>   *flash_pe_v;
   vector<float>   *slice_pe_v;
   Float_t         flash_time;
   Float_t         flash_y;
   Float_t         flash_z;
   Float_t         flash_timewidth;
   Float_t         flash_ywidth;
   Float_t         flash_zwidth;
   Float_t         nu_flashmatch_score;
   Float_t         nu_centerX;
   Float_t         nu_centerY;
   Float_t         nu_centerZ;
   Float_t         nu_totalCharge;
   Float_t         best_cosmic_flashmatch_score;
   Float_t         best_obviouscosmic_flashmatch_score;
   vector<float>   *cosmic_flashmatch_score_v;
   vector<float>   *cosmic_topological_score_v;
   vector<float>   *cosmic_centerX_v;
   vector<float>   *cosmic_centerY_v;
   vector<float>   *cosmic_centerZ_v;
   vector<float>   *cosmic_totalCharge_v;
   vector<int>     *cosmic_nhits_v;
   vector<int>     *cosmic_nunhits_v;
   vector<int>     *cosmic_isclear_v;
   //Float_t         flash_pe;
   Float_t         flash_pe_calib;
   //Float_t         flash_time;
   Float_t         flash_zcenter;
   Float_t         flash_ycenter;
   //Float_t         flash_zwidth;
   //Float_t         flash_ywidth;
   //vector<float>   *flash_pe_v;
   vector<float>   *flash_pe_calib_v;
   vector<float>   *gain_area_v;
   vector<float>   *gain_ampl_v;
   vector<short>   *waveform_00;
   vector<short>   *waveform_01;
   vector<short>   *waveform_02;
   vector<short>   *waveform_03;
   vector<short>   *waveform_04;
   vector<short>   *waveform_05;
   vector<short>   *waveform_06;
   vector<short>   *waveform_07;
   vector<short>   *waveform_08;
   vector<short>   *waveform_09;
   vector<short>   *waveform_10;
   vector<short>   *waveform_11;
   vector<short>   *waveform_12;
   vector<short>   *waveform_13;
   vector<short>   *waveform_14;
   vector<short>   *waveform_15;
   vector<short>   *waveform_16;
   vector<short>   *waveform_17;
   vector<short>   *waveform_18;
   vector<short>   *waveform_19;
   vector<short>   *waveform_20;
   vector<short>   *waveform_21;
   vector<short>   *waveform_22;
   vector<short>   *waveform_23;
   vector<short>   *waveform_24;
   vector<short>   *waveform_25;
   vector<short>   *waveform_26;
   vector<short>   *waveform_27;
   vector<short>   *waveform_28;
   vector<short>   *waveform_29;
   vector<short>   *waveform_30;
   vector<short>   *waveform_31;
   Float_t         secondshower_U_charge;
   Int_t           secondshower_U_nhit;
   Float_t         secondshower_U_vtxdist;
   Float_t         secondshower_U_eigenratio;
   Float_t         secondshower_U_dot;
   Float_t         secondshower_U_dir;
   Float_t         secondshower_V_charge;
   Int_t           secondshower_V_nhit;
   Float_t         secondshower_V_vtxdist;
   Float_t         secondshower_V_eigenratio;
   Float_t         secondshower_V_dot;
   Float_t         secondshower_V_dir;
   Float_t         secondshower_Y_charge;
   Int_t           secondshower_Y_nhit;
   Float_t         secondshower_Y_vtxdist;
   Float_t         secondshower_Y_eigenratio;
   Float_t         secondshower_Y_dot;
   Float_t         secondshower_Y_dir;
   vector<float>   *simphoton_number_v;
   vector<float>   *simphoton_tmin_v;
   vector<float>   *simphoton_tmax_v;
   vector<float>   *trk_bragg_p_v;
   vector<float>   *trk_bragg_mu_v;
   vector<float>   *trk_bragg_mip_v;
   vector<float>   *trk_pida_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_pid_chipi_v;
   vector<float>   *trk_pid_chika_v;
   vector<float>   *trk_pid_chimu_v;
   vector<float>   *trk_bragg_p_u_v;
   vector<float>   *trk_bragg_mu_u_v;
   vector<float>   *trk_bragg_mip_u_v;
   vector<float>   *trk_pida_u_v;
   vector<float>   *trk_pid_chipr_u_v;
   vector<float>   *trk_pid_chipi_u_v;
   vector<float>   *trk_pid_chika_u_v;
   vector<float>   *trk_pid_chimu_u_v;
   vector<float>   *trk_bragg_p_v_v;
   vector<float>   *trk_bragg_mu_v_v;
   vector<float>   *trk_bragg_mip_v_v;
   vector<float>   *trk_pida_v_v;
   vector<float>   *trk_pid_chipr_v_v;
   vector<float>   *trk_pid_chipi_v_v;
   vector<float>   *trk_pid_chika_v_v;
   vector<float>   *trk_pid_chimu_v_v;
   vector<unsigned long> *trk_pfp_id_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_start_x_v;
   vector<float>   *trk_start_y_v;
   vector<float>   *trk_start_z_v;
   vector<float>   *trk_sce_start_x_v;
   vector<float>   *trk_sce_start_y_v;
   vector<float>   *trk_sce_start_z_v;
   vector<float>   *trk_end_x_v;
   vector<float>   *trk_end_y_v;
   vector<float>   *trk_end_z_v;
   vector<float>   *trk_sce_end_x_v;
   vector<float>   *trk_sce_end_y_v;
   vector<float>   *trk_sce_end_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_theta_v;
   vector<float>   *trk_phi_v;
   vector<float>   *trk_len_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_range_muon_mom_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_energy_muon_v;
   vector<float>   *trk_calo_energy_u_v;
   vector<float>   *trk_calo_energy_v_v;
   vector<float>   *trk_calo_energy_y_v;
   vector<float>   *trk_llr_pid_u_v;
   vector<float>   *trk_llr_pid_v_v;
   vector<float>   *trk_llr_pid_y_v;
   vector<float>   *trk_llr_pid_v;
   vector<float>   *trk_llr_pid_score_v;

   // List of branches
   TBranch        *b_selected;   //!
   TBranch        *b_run;   //!
   TBranch        *b_sub;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_trk_pfp_id;   //!
   TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_trk2_pfp_id;   //!
   TBranch        *b_shr2_pfp_id;   //!
   TBranch        *b_shr_energy_tot;   //!
   TBranch        *b_shr_energy;   //!
   TBranch        *b_shr_energy_tot_cali;   //!
   TBranch        *b_shr_energy_cali;   //!
   TBranch        *b_shr_theta;   //!
   TBranch        *b_shr_phi;   //!
   TBranch        *b_shr_pca_0;   //!
   TBranch        *b_shr_pca_1;   //!
   TBranch        *b_shr_pca_2;   //!
   TBranch        *b_shr_px;   //!
   TBranch        *b_shr_py;   //!
   TBranch        *b_shr_pz;   //!
   TBranch        *b_shr_openangle;   //!
   TBranch        *b_shr_tkfit_start_x;   //!
   TBranch        *b_shr_tkfit_start_y;   //!
   TBranch        *b_shr_tkfit_start_z;   //!
   TBranch        *b_shr_tkfit_theta;   //!
   TBranch        *b_shr_tkfit_phi;   //!
   TBranch        *b_shr_start_x;   //!
   TBranch        *b_shr_start_y;   //!
   TBranch        *b_shr_start_z;   //!
   TBranch        *b_shr_dedx_Y;   //!
   TBranch        *b_shr_dedx_V;   //!
   TBranch        *b_shr_dedx_U;   //!
   TBranch        *b_shr_dedx_Y_cali;   //!
   TBranch        *b_shr_dedx_V_cali;   //!
   TBranch        *b_shr_dedx_U_cali;   //!
   TBranch        *b_shr_tkfit_dedx_Y;   //!
   TBranch        *b_shr_tkfit_dedx_V;   //!
   TBranch        *b_shr_tkfit_dedx_U;   //!
   TBranch        *b_shr_tkfit_dedx_max;   //!
   TBranch        *b_shr_tkfit_nhits_Y;   //!
   TBranch        *b_shr_tkfit_nhits_V;   //!
   TBranch        *b_shr_tkfit_nhits_U;   //!
   TBranch        *b_shr_llrpid_dedx_Y;   //!
   TBranch        *b_shr_llrpid_dedx_V;   //!
   TBranch        *b_shr_llrpid_dedx_U;   //!
   TBranch        *b_shr_llrpid_dedx;   //!
   TBranch        *b_shr_tkfit_dedx_Y_alt;   //!
   TBranch        *b_shr_tkfit_dedx_V_alt;   //!
   TBranch        *b_shr_tkfit_dedx_U_alt;   //!
   TBranch        *b_shr_tkfit_nhits_Y_alt;   //!
   TBranch        *b_shr_tkfit_nhits_V_alt;   //!
   TBranch        *b_shr_tkfit_nhits_U_alt;   //!
   TBranch        *b__trkfit;   //!
   TBranch        *b_shr_tkfit_npoints;   //!
   TBranch        *b_shr_tkfit_npointsvalid;   //!
   TBranch        *b_shr_trkfitmedangle;   //!
   TBranch        *b_shrmoliereavg;   //!
   TBranch        *b_shrmoliererms;   //!
   TBranch        *b_shr1shr2moliereavg;   //!
   TBranch        *b_shr1shr2moliererms;   //!
   TBranch        *b_shr1trk1moliereavg;   //!
   TBranch        *b_shr1trk1moliererms;   //!
   TBranch        *b_shr1trk2moliereavg;   //!
   TBranch        *b_shr1trk2moliererms;   //!
   TBranch        *b_ismerged;   //!
   TBranch        *b_merge_bestdot;   //!
   TBranch        *b_merge_bestdist;   //!
   TBranch        *b_merge_vtx_x;   //!
   TBranch        *b_merge_vtx_y;   //!
   TBranch        *b_merge_vtx_z;   //!
   TBranch        *b_merge_tk_ipfp;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_Y;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_V;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_U;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_Y;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_V;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_U;   //!
   TBranch        *b_shr_chipr;   //!
   TBranch        *b_shr_chimu;   //!
   TBranch        *b_shr_bragg_p;   //!
   TBranch        *b_shr_bragg_mu;   //!
   TBranch        *b_shr_bragg_mip;   //!
   TBranch        *b_shr_bragg_kaon;   //!
   TBranch        *b_shr_bragg_pion;   //!
   TBranch        *b_tksh_distance;   //!
   TBranch        *b_tksh_angle;   //!
   TBranch        *b_shr_distance;   //!
   TBranch        *b_shr_score;   //!
   TBranch        *b_shr_bkt_pdg;   //!
   TBranch        *b_shr_bkt_purity;   //!
   TBranch        *b_shr_bkt_completeness;   //!
   TBranch        *b_shr_bkt_E;   //!
   TBranch        *b_trk_len;   //!
   TBranch        *b_trk_theta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_energy;   //!
   TBranch        *b_trk_energy_muon;   //!
   TBranch        *b_trk_energy_muon_mcs;   //!
   TBranch        *b_trk_energy_tot;   //!
   TBranch        *b_trk_energy_muon_tot;   //!
   TBranch        *b_trk_distance;   //!
   TBranch        *b_trk_score;   //!
   TBranch        *b_trk_bkt_pdg;   //!
   TBranch        *b_trk_bkt_purity;   //!
   TBranch        *b_trk_bkt_completeness;   //!
   TBranch        *b_trk_bkt_E;   //!
   TBranch        *b_trk_chipr_best;   //!
   TBranch        *b_trk_chipr_worst;   //!
   TBranch        *b_trk_chimu_best;   //!
   TBranch        *b_trk_chimu_worst;   //!
   TBranch        *b_trk_chipr;   //!
   TBranch        *b_trk_chimu;   //!
   TBranch        *b_trk_pida;   //!
   TBranch        *b_trk_bragg_p;   //!
   TBranch        *b_trk_bragg_mu;   //!
   TBranch        *b_trk_bragg_mip;   //!
   TBranch        *b_trk_bragg_kaon;   //!
   TBranch        *b_trk_bragg_pion;   //!
   TBranch        *b_trk_hits_max;   //!
   TBranch        *b_shr_hits_max;   //!
   TBranch        *b_all_shr_hits;   //!
   TBranch        *b_all_trk_hits;   //!
   TBranch        *b_all_shr_energies;   //!
   TBranch        *b_all_trk_energies;   //!
   //TBranch        *b_shr_hits_max;   //!
   TBranch        *b_trk_hits_2nd;   //!
   TBranch        *b_shr_hits_2nd;   //!
   TBranch        *b_trkshrhitdist0;   //!
   TBranch        *b_trkshrhitdist1;   //!
   TBranch        *b_trkshrhitdist2;   //!
   TBranch        *b_trk2shrhitdist0;   //!
   TBranch        *b_trk2shrhitdist1;   //!
   TBranch        *b_trk2shrhitdist2;   //!
   TBranch        *b_trk1trk2hitdist0;   //!
   TBranch        *b_trk1trk2hitdist1;   //!
   TBranch        *b_trk1trk2hitdist2;   //!
   TBranch        *b_total_hits_y;   //!
   TBranch        *b_extra_energy_y;   //!
   TBranch        *b_trk_energy_hits_tot;   //!
   TBranch        *b_subcluster;   //!
   TBranch        *b_shrsubclusters0;   //!
   TBranch        *b_shrsubclusters1;   //!
   TBranch        *b_shrsubclusters2;   //!
   TBranch        *b_shrclusfrac0;   //!
   TBranch        *b_shrclusfrac1;   //!
   TBranch        *b_shrclusfrac2;   //!
   TBranch        *b_shrclusdir0;   //!
   TBranch        *b_shrclusdir1;   //!
   TBranch        *b_shrclusdir2;   //!
   TBranch        *b_shr_hits_tot;   //!
   TBranch        *b_shr_hits_y_tot;   //!
   TBranch        *b_shr_hits_u_tot;   //!
   TBranch        *b_shr_hits_v_tot;   //!
   TBranch        *b_trk_hits_tot;   //!
   TBranch        *b_trk_hits_y_tot;   //!
   TBranch        *b_trk_hits_u_tot;   //!
   TBranch        *b_trk_hits_v_tot;   //!
   TBranch        *b_elecclusters_U_charge;   //!
   TBranch        *b_elecclusters_V_charge;   //!
   TBranch        *b_elecclusters_Y_charge;   //!
   TBranch        *b_elecclusters_U_N;   //!
   TBranch        *b_elecclusters_V_N;   //!
   TBranch        *b_elecclusters_Y_N;   //!
   TBranch        *b_n_tracks_contained;   //!
   TBranch        *b_n_showers_contained;   //!
   TBranch        *b_matched_E;   //!
   TBranch        *b_hits_ratio;   //!
   TBranch        *b_contained_fraction;   //!
   TBranch        *b_sps_contained_fraction;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt_assume_muon;   //!
   TBranch        *b_p_assume_muon;   //!
   TBranch        *b_reco_e;   //!
   TBranch        *b_dvtx;   //!
   TBranch        *b_dtrk;   //!
   TBranch        *b_contained_sps_ratio;   //!
   TBranch        *b_dtrk_x_boundary;   //!
   TBranch        *b_dtrk_y_boundary;   //!
   TBranch        *b_dtrk_z_boundary;   //!
   TBranch        *b_dshr_x_boundary;   //!
   TBranch        *b_dshr_y_boundary;   //!
   TBranch        *b_dshr_z_boundary;   //!
   TBranch        *b_dvtx_x_boundary;   //!
   TBranch        *b_dvtx_y_boundary;   //!
   TBranch        *b_dvtx_z_boundary;   //!
   TBranch        *b_dtrk_boundary;   //!
   TBranch        *b_dvtx_boundary;   //!
   TBranch        *b_dshr_boundary;   //!
   TBranch        *b_dmc_boundary;   //!
   TBranch        *b_leeweight;   //!
   TBranch        *b_true_pt;   //!
   TBranch        *b_true_pt_visible;   //!
   TBranch        *b_true_p;   //!
   TBranch        *b_true_p_visible;   //!
   TBranch        *b_true_e_visible;   //!
   TBranch        *b_opfilter_pe_beam;   //!
   TBranch        *b_opfilter_pe_veto;   //!
   TBranch        *b_nu_pdg;   //!
   TBranch        *b_ccnc;   //!
   TBranch        *b_nu_parent_pdg;   //!
   TBranch        *b_nu_hadron_pdg;   //!
   TBranch        *b_nu_decay_mode;   //!
   TBranch        *b_interaction;   //!
   TBranch        *b_nu_e;   //!
   TBranch        *b_nu_l;   //!
   TBranch        *b_nu_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_isVtxInFiducial;   //!
   TBranch        *b_truthFiducial;   //!
   TBranch        *b_true_nu_vtx_t;   //!
   TBranch        *b_true_nu_vtx_x;   //!
   TBranch        *b_true_nu_vtx_y;   //!
   TBranch        *b_true_nu_vtx_z;   //!
   TBranch        *b_true_nu_vtx_sce_x;   //!
   TBranch        *b_true_nu_vtx_sce_y;   //!
   TBranch        *b_true_nu_vtx_sce_z;   //!
   TBranch        *b_reco_nu_vtx_x;   //!
   TBranch        *b_reco_nu_vtx_y;   //!
   TBranch        *b_reco_nu_vtx_z;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_e;   //!
   TBranch        *b_muon_c;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_nelec;   //!
   TBranch        *b_elec_e;   //!
   TBranch        *b_elec_c;   //!
   TBranch        *b_elec_p;   //!
   TBranch        *b_elec_vx;   //!
   TBranch        *b_elec_vy;   //!
   TBranch        *b_elec_vz;   //!
   TBranch        *b_elec_px;   //!
   TBranch        *b_elec_py;   //!
   TBranch        *b_elec_pz;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_c;   //!
   TBranch        *b_pi0_p;   //!
   TBranch        *b_nneutron;   //!
   TBranch        *b_nproton;   //!
   TBranch        *b_proton_e;   //!
   TBranch        *b_proton_c;   //!
   TBranch        *b_proton_p;   //!
   TBranch        *b_npion;   //!
   TBranch        *b_pion_e;   //!
   TBranch        *b_pion_c;   //!
   TBranch        *b_pion_p;   //!
   TBranch        *b_neta;   //!
   TBranch        *b_eta_e;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_crtveto;   //!
   TBranch        *b_crthitpe;   //!
   TBranch        *b_pfp_slice_idx;   //!
   TBranch        *b_category;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_tid;   //!
   TBranch        *b_backtracked_purity;   //!
   TBranch        *b_backtracked_completeness;   //!
   TBranch        *b_backtracked_overlay_purity;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_backtracked_start_x;   //!
   TBranch        *b_backtracked_start_y;   //!
   TBranch        *b_backtracked_start_z;   //!
   TBranch        *b_backtracked_start_t;   //!
   TBranch        *b_backtracked_start_U;   //!
   TBranch        *b_backtracked_start_V;   //!
   TBranch        *b_backtracked_start_Y;   //!
   TBranch        *b_backtracked_sce_start_x;   //!
   TBranch        *b_backtracked_sce_start_y;   //!
   TBranch        *b_backtracked_sce_start_z;   //!
   TBranch        *b_backtracked_sce_start_U;   //!
   TBranch        *b_backtracked_sce_start_V;   //!
   TBranch        *b_backtracked_sce_start_Y;   //!
   TBranch        *b_lep_e;   //!
   TBranch        *b_pass;   //!
   TBranch        *b_swtrig;   //!
   TBranch        *b_evnhits;   //!
   TBranch        *b_slpdg;   //!
   TBranch        *b_slnhits;   //!
   TBranch        *b_n_pfps;   //!
   TBranch        *b_n_tracks;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_pfp_generation_v;   //!
   TBranch        *b_pfp_trk_daughters_v;   //!
   TBranch        *b_pfp_shr_daughters_v;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_pfpplanesubclusters_U;   //!
   TBranch        *b_pfpplanesubclusters_V;   //!
   TBranch        *b_pfpplanesubclusters_Y;   //!
   TBranch        *b_pfpplanesubhitfracmax_U;   //!
   TBranch        *b_pfpplanesubhitfracmax_V;   //!
   TBranch        *b_pfpplanesubhitfracmax_Y;   //!
   TBranch        *b_hits_u;   //!
   TBranch        *b_hits_v;   //!
   TBranch        *b_hits_y;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_slclustfrac;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_vx;   //!
   TBranch        *b_mc_vy;   //!
   TBranch        *b_mc_vz;   //!
   TBranch        *b_mc_endx;   //!
   TBranch        *b_mc_endy;   //!
   TBranch        *b_mc_endz;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_completeness;   //!
   TBranch        *b_mc_purity;   //!
   TBranch        *b_endmuonprocess;   //!
   TBranch        *b_endmuonmichel;   //!
   TBranch        *b_flash_pe;   //!
   TBranch        *b_flash_pe_v;   //!
   TBranch        *b_slice_pe_v;   //!
   TBranch        *b_flash_time;   //!
   TBranch        *b_flash_y;   //!
   TBranch        *b_flash_z;   //!
   TBranch        *b_flash_timewidth;   //!
   TBranch        *b_flash_ywidth;   //!
   TBranch        *b_flash_zwidth;   //!
   TBranch        *b_nu_flashmatch_score;   //!
   TBranch        *b_nu_centerX;   //!
   TBranch        *b_nu_centerY;   //!
   TBranch        *b_nu_centerZ;   //!
   TBranch        *b_nu_totalCharge;   //!
   TBranch        *b_best_cosmic_flashmatch_score;   //!
   TBranch        *b_best_obviouscosmic_flashmatch_score;   //!
   TBranch        *b_cosmic_flashmatch_score_v;   //!
   TBranch        *b_cosmic_topological_score_v;   //!
   TBranch        *b_cosmic_centerX_v;   //!
   TBranch        *b_cosmic_centerY_v;   //!
   TBranch        *b_cosmic_centerZ_v;   //!
   TBranch        *b_cosmic_totalCharge_v;   //!
   TBranch        *b_cosmic_nhits_v;   //!
   TBranch        *b_cosmic_nunhits_v;   //!
   TBranch        *b_cosmic_isclear_v;   //!
   //TBranch        *b_flash_pe;   //!
   TBranch        *b_flash_pe_calib;   //!
   //TBranch        *b_flash_time;   //!
   TBranch        *b_flash_zcenter;   //!
   TBranch        *b_flash_ycenter;   //!
   //TBranch        *b_flash_zwidth;   //!
   //TBranch        *b_flash_ywidth;   //!
   //TBranch        *b_flash_pe_v;   //!
   TBranch        *b_flash_pe_calib_v;   //!
   TBranch        *b_gain_area_v;   //!
   TBranch        *b_gain_ampl_v;   //!
   TBranch        *b_waveform_00;   //!
   TBranch        *b_waveform_01;   //!
   TBranch        *b_waveform_02;   //!
   TBranch        *b_waveform_03;   //!
   TBranch        *b_waveform_04;   //!
   TBranch        *b_waveform_05;   //!
   TBranch        *b_waveform_06;   //!
   TBranch        *b_waveform_07;   //!
   TBranch        *b_waveform_08;   //!
   TBranch        *b_waveform_09;   //!
   TBranch        *b_waveform_10;   //!
   TBranch        *b_waveform_11;   //!
   TBranch        *b_waveform_12;   //!
   TBranch        *b_waveform_13;   //!
   TBranch        *b_waveform_14;   //!
   TBranch        *b_waveform_15;   //!
   TBranch        *b_waveform_16;   //!
   TBranch        *b_waveform_17;   //!
   TBranch        *b_waveform_18;   //!
   TBranch        *b_waveform_19;   //!
   TBranch        *b_waveform_20;   //!
   TBranch        *b_waveform_21;   //!
   TBranch        *b_waveform_22;   //!
   TBranch        *b_waveform_23;   //!
   TBranch        *b_waveform_24;   //!
   TBranch        *b_waveform_25;   //!
   TBranch        *b_waveform_26;   //!
   TBranch        *b_waveform_27;   //!
   TBranch        *b_waveform_28;   //!
   TBranch        *b_waveform_29;   //!
   TBranch        *b_waveform_30;   //!
   TBranch        *b_waveform_31;   //!
   TBranch        *b_secondshower_U_charge;   //!
   TBranch        *b_secondshower_U_nhit;   //!
   TBranch        *b_secondshower_U_vtxdist;   //!
   TBranch        *b_secondshower_U_eigenratio;   //!
   TBranch        *b_secondshower_U_dot;   //!
   TBranch        *b_secondshower_U_dir;   //!
   TBranch        *b_secondshower_V_charge;   //!
   TBranch        *b_secondshower_V_nhit;   //!
   TBranch        *b_secondshower_V_vtxdist;   //!
   TBranch        *b_secondshower_V_eigenratio;   //!
   TBranch        *b_secondshower_V_dot;   //!
   TBranch        *b_secondshower_V_dir;   //!
   TBranch        *b_secondshower_Y_charge;   //!
   TBranch        *b_secondshower_Y_nhit;   //!
   TBranch        *b_secondshower_Y_vtxdist;   //!
   TBranch        *b_secondshower_Y_eigenratio;   //!
   TBranch        *b_secondshower_Y_dot;   //!
   TBranch        *b_secondshower_Y_dir;   //!
   TBranch        *b_simphoton_number_v;   //!
   TBranch        *b_simphoton_tmin_v;   //!
   TBranch        *b_simphoton_tmax_v;   //!
   TBranch        *b_trk_bragg_p_v;   //!
   TBranch        *b_trk_bragg_mu_v;   //!
   TBranch        *b_trk_bragg_mip_v;   //!
   TBranch        *b_trk_pida_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_pid_chipi_v;   //!
   TBranch        *b_trk_pid_chika_v;   //!
   TBranch        *b_trk_pid_chimu_v;   //!
   TBranch        *b_trk_bragg_p_u_v;   //!
   TBranch        *b_trk_bragg_mu_u_v;   //!
   TBranch        *b_trk_bragg_mip_u_v;   //!
   TBranch        *b_trk_pida_u_v;   //!
   TBranch        *b_trk_pid_chipr_u_v;   //!
   TBranch        *b_trk_pid_chipi_u_v;   //!
   TBranch        *b_trk_pid_chika_u_v;   //!
   TBranch        *b_trk_pid_chimu_u_v;   //!
   TBranch        *b_trk_bragg_p_v_v;   //!
   TBranch        *b_trk_bragg_mu_v_v;   //!
   TBranch        *b_trk_bragg_mip_v_v;   //!
   TBranch        *b_trk_pida_v_v;   //!
   TBranch        *b_trk_pid_chipr_v_v;   //!
   TBranch        *b_trk_pid_chipi_v_v;   //!
   TBranch        *b_trk_pid_chika_v_v;   //!
   TBranch        *b_trk_pid_chimu_v_v;   //!
   TBranch        *b_trk_pfp_id_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_start_x_v;   //!
   TBranch        *b_trk_start_y_v;   //!
   TBranch        *b_trk_start_z_v;   //!
   TBranch        *b_trk_sce_start_x_v;   //!
   TBranch        *b_trk_sce_start_y_v;   //!
   TBranch        *b_trk_sce_start_z_v;   //!
   TBranch        *b_trk_end_x_v;   //!
   TBranch        *b_trk_end_y_v;   //!
   TBranch        *b_trk_end_z_v;   //!
   TBranch        *b_trk_sce_end_x_v;   //!
   TBranch        *b_trk_sce_end_y_v;   //!
   TBranch        *b_trk_sce_end_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_theta_v;   //!
   TBranch        *b_trk_phi_v;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_range_muon_mom_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_energy_muon_v;   //!
   TBranch        *b_trk_calo_energy_u_v;   //!
   TBranch        *b_trk_calo_energy_v_v;   //!
   TBranch        *b_trk_calo_energy_y_v;   //!
   TBranch        *b_trk_llr_pid_u_v;   //!
   TBranch        *b_trk_llr_pid_v_v;   //!
   TBranch        *b_trk_llr_pid_y_v;   //!
   TBranch        *b_trk_llr_pid_v;   //!
   TBranch        *b_trk_llr_pid_score_v;   //!

   mu_sim_mely_class(TTree *tree=0);
   virtual ~mu_sim_mely_class();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mu_sim_mely_class_cxx
mu_sim_mely_class::mu_sim_mely_class(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ub_mu_sim_40k.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ub_mu_sim_40k.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ub_mu_sim_40k.root:/nuselection");
      dir->GetObject("NeutrinoSelectionFilter",tree);

   }
   Init(tree);
}

mu_sim_mely_class::~mu_sim_mely_class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mu_sim_mely_class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mu_sim_mely_class::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mu_sim_mely_class::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   all_shr_hits = 0;
   all_trk_hits = 0;
   all_shr_energies = 0;
   all_trk_energies = 0;
   dtrk_x_boundary = 0;
   dtrk_y_boundary = 0;
   dtrk_z_boundary = 0;
   dshr_x_boundary = 0;
   dshr_y_boundary = 0;
   dshr_z_boundary = 0;
   dvtx_x_boundary = 0;
   dvtx_y_boundary = 0;
   dvtx_z_boundary = 0;
   dtrk_boundary = 0;
   dvtx_boundary = 0;
   dshr_boundary = 0;
   dmc_boundary = 0;
   pfp_slice_idx = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_tid = 0;
   backtracked_purity = 0;
   backtracked_completeness = 0;
   backtracked_overlay_purity = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   backtracked_start_x = 0;
   backtracked_start_y = 0;
   backtracked_start_z = 0;
   backtracked_start_t = 0;
   backtracked_start_U = 0;
   backtracked_start_V = 0;
   backtracked_start_Y = 0;
   backtracked_sce_start_x = 0;
   backtracked_sce_start_y = 0;
   backtracked_sce_start_z = 0;
   backtracked_sce_start_U = 0;
   backtracked_sce_start_V = 0;
   backtracked_sce_start_Y = 0;
   pfp_generation_v = 0;
   pfp_trk_daughters_v = 0;
   pfp_shr_daughters_v = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   pfpplanesubclusters_U = 0;
   pfpplanesubclusters_V = 0;
   pfpplanesubclusters_Y = 0;
   pfpplanesubhitfracmax_U = 0;
   pfpplanesubhitfracmax_V = 0;
   pfpplanesubhitfracmax_Y = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_vx = 0;
   mc_vy = 0;
   mc_vz = 0;
   mc_endx = 0;
   mc_endy = 0;
   mc_endz = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_completeness = 0;
   mc_purity = 0;
   endmuonprocess = 0;
   flash_pe_v = 0;
   slice_pe_v = 0;
   cosmic_flashmatch_score_v = 0;
   cosmic_topological_score_v = 0;
   cosmic_centerX_v = 0;
   cosmic_centerY_v = 0;
   cosmic_centerZ_v = 0;
   cosmic_totalCharge_v = 0;
   cosmic_nhits_v = 0;
   cosmic_nunhits_v = 0;
   cosmic_isclear_v = 0;
   //flash_pe_v = 0;
   flash_pe_calib_v = 0;
   gain_area_v = 0;
   gain_ampl_v = 0;
   waveform_00 = 0;
   waveform_01 = 0;
   waveform_02 = 0;
   waveform_03 = 0;
   waveform_04 = 0;
   waveform_05 = 0;
   waveform_06 = 0;
   waveform_07 = 0;
   waveform_08 = 0;
   waveform_09 = 0;
   waveform_10 = 0;
   waveform_11 = 0;
   waveform_12 = 0;
   waveform_13 = 0;
   waveform_14 = 0;
   waveform_15 = 0;
   waveform_16 = 0;
   waveform_17 = 0;
   waveform_18 = 0;
   waveform_19 = 0;
   waveform_20 = 0;
   waveform_21 = 0;
   waveform_22 = 0;
   waveform_23 = 0;
   waveform_24 = 0;
   waveform_25 = 0;
   waveform_26 = 0;
   waveform_27 = 0;
   waveform_28 = 0;
   waveform_29 = 0;
   waveform_30 = 0;
   waveform_31 = 0;
   simphoton_number_v = 0;
   simphoton_tmin_v = 0;
   simphoton_tmax_v = 0;
   trk_bragg_p_v = 0;
   trk_bragg_mu_v = 0;
   trk_bragg_mip_v = 0;
   trk_pida_v = 0;
   trk_pid_chipr_v = 0;
   trk_pid_chipi_v = 0;
   trk_pid_chika_v = 0;
   trk_pid_chimu_v = 0;
   trk_bragg_p_u_v = 0;
   trk_bragg_mu_u_v = 0;
   trk_bragg_mip_u_v = 0;
   trk_pida_u_v = 0;
   trk_pid_chipr_u_v = 0;
   trk_pid_chipi_u_v = 0;
   trk_pid_chika_u_v = 0;
   trk_pid_chimu_u_v = 0;
   trk_bragg_p_v_v = 0;
   trk_bragg_mu_v_v = 0;
   trk_bragg_mip_v_v = 0;
   trk_pida_v_v = 0;
   trk_pid_chipr_v_v = 0;
   trk_pid_chipi_v_v = 0;
   trk_pid_chika_v_v = 0;
   trk_pid_chimu_v_v = 0;
   trk_pfp_id_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_start_x_v = 0;
   trk_start_y_v = 0;
   trk_start_z_v = 0;
   trk_sce_start_x_v = 0;
   trk_sce_start_y_v = 0;
   trk_sce_start_z_v = 0;
   trk_end_x_v = 0;
   trk_end_y_v = 0;
   trk_end_z_v = 0;
   trk_sce_end_x_v = 0;
   trk_sce_end_y_v = 0;
   trk_sce_end_z_v = 0;
   trk_distance_v = 0;
   trk_theta_v = 0;
   trk_phi_v = 0;
   trk_len_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_range_muon_mom_v = 0;
   trk_energy_proton_v = 0;
   trk_energy_muon_v = 0;
   trk_calo_energy_u_v = 0;
   trk_calo_energy_v_v = 0;
   trk_calo_energy_y_v = 0;
   trk_llr_pid_u_v = 0;
   trk_llr_pid_v_v = 0;
   trk_llr_pid_y_v = 0;
   trk_llr_pid_v = 0;
   trk_llr_pid_score_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("sub", &sub, &b_sub);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("trk_id", &trk_id, &b_trk_pfp_id);
   fChain->SetBranchAddress("shr_id", &shr_id, &b_shr_pfp_id);
   fChain->SetBranchAddress("trk2_id", &trk2_id, &b_trk2_pfp_id);
   fChain->SetBranchAddress("shr2_id", &shr2_id, &b_shr2_pfp_id);
   fChain->SetBranchAddress("shr_energy_tot", &shr_energy_tot, &b_shr_energy_tot);
   fChain->SetBranchAddress("shr_energy", &shr_energy, &b_shr_energy);
   fChain->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali, &b_shr_energy_tot_cali);
   fChain->SetBranchAddress("shr_energy_cali", &shr_energy_cali, &b_shr_energy_cali);
   fChain->SetBranchAddress("shr_theta", &shr_theta, &b_shr_theta);
   fChain->SetBranchAddress("shr_phi", &shr_phi, &b_shr_phi);
   fChain->SetBranchAddress("shr_pca_0", &shr_pca_0, &b_shr_pca_0);
   fChain->SetBranchAddress("shr_pca_1", &shr_pca_1, &b_shr_pca_1);
   fChain->SetBranchAddress("shr_pca_2", &shr_pca_2, &b_shr_pca_2);
   fChain->SetBranchAddress("shr_px", &shr_px, &b_shr_px);
   fChain->SetBranchAddress("shr_py", &shr_py, &b_shr_py);
   fChain->SetBranchAddress("shr_pz", &shr_pz, &b_shr_pz);
   fChain->SetBranchAddress("shr_openangle", &shr_openangle, &b_shr_openangle);
   fChain->SetBranchAddress("shr_tkfit_start_x", &shr_tkfit_start_x, &b_shr_tkfit_start_x);
   fChain->SetBranchAddress("shr_tkfit_start_y", &shr_tkfit_start_y, &b_shr_tkfit_start_y);
   fChain->SetBranchAddress("shr_tkfit_start_z", &shr_tkfit_start_z, &b_shr_tkfit_start_z);
   fChain->SetBranchAddress("shr_tkfit_theta", &shr_tkfit_theta, &b_shr_tkfit_theta);
   fChain->SetBranchAddress("shr_tkfit_phi", &shr_tkfit_phi, &b_shr_tkfit_phi);
   fChain->SetBranchAddress("shr_start_x", &shr_start_x, &b_shr_start_x);
   fChain->SetBranchAddress("shr_start_y", &shr_start_y, &b_shr_start_y);
   fChain->SetBranchAddress("shr_start_z", &shr_start_z, &b_shr_start_z);
   fChain->SetBranchAddress("shr_dedx_Y", &shr_dedx_Y, &b_shr_dedx_Y);
   fChain->SetBranchAddress("shr_dedx_V", &shr_dedx_V, &b_shr_dedx_V);
   fChain->SetBranchAddress("shr_dedx_U", &shr_dedx_U, &b_shr_dedx_U);
   fChain->SetBranchAddress("shr_dedx_Y_cali", &shr_dedx_Y_cali, &b_shr_dedx_Y_cali);
   fChain->SetBranchAddress("shr_dedx_V_cali", &shr_dedx_V_cali, &b_shr_dedx_V_cali);
   fChain->SetBranchAddress("shr_dedx_U_cali", &shr_dedx_U_cali, &b_shr_dedx_U_cali);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, &b_shr_tkfit_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, &b_shr_tkfit_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, &b_shr_tkfit_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_dedx_max", &shr_tkfit_dedx_max, &b_shr_tkfit_dedx_max);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y, &b_shr_tkfit_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V, &b_shr_tkfit_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U, &b_shr_tkfit_nhits_U);
   fChain->SetBranchAddress("shr_llrpid_dedx_Y", &shr_llrpid_dedx_Y, &b_shr_llrpid_dedx_Y);
   fChain->SetBranchAddress("shr_llrpid_dedx_V", &shr_llrpid_dedx_V, &b_shr_llrpid_dedx_V);
   fChain->SetBranchAddress("shr_llrpid_dedx_U", &shr_llrpid_dedx_U, &b_shr_llrpid_dedx_U);
   fChain->SetBranchAddress("shr_llrpid_dedx", &shr_llrpid_dedx, &b_shr_llrpid_dedx);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y_alt", &shr_tkfit_dedx_Y_alt, &b_shr_tkfit_dedx_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_V_alt", &shr_tkfit_dedx_V_alt, &b_shr_tkfit_dedx_V_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_U_alt", &shr_tkfit_dedx_U_alt, &b_shr_tkfit_dedx_U_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y_alt", &shr_tkfit_nhits_Y_alt, &b_shr_tkfit_nhits_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_V_alt", &shr_tkfit_nhits_V_alt, &b_shr_tkfit_nhits_V_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_U_alt", &shr_tkfit_nhits_U_alt, &b_shr_tkfit_nhits_U_alt);
   fChain->SetBranchAddress("trkfit", &trkfit, &b__trkfit);
   fChain->SetBranchAddress("shr_tkfit_npoints", &shr_tkfit_npoints, &b_shr_tkfit_npoints);
   fChain->SetBranchAddress("shr_tkfit_npointsvalid", &shr_tkfit_npointsvalid, &b_shr_tkfit_npointsvalid);
   fChain->SetBranchAddress("shr_trkfitmedangle", &shr_trkfitmedangle, &b_shr_trkfitmedangle);
   fChain->SetBranchAddress("shrmoliereavg", &shrmoliereavg, &b_shrmoliereavg);
   fChain->SetBranchAddress("shrmoliererms", &shrmoliererms, &b_shrmoliererms);
   fChain->SetBranchAddress("shr1shr2moliereavg", &shr1shr2moliereavg, &b_shr1shr2moliereavg);
   fChain->SetBranchAddress("shr1shr2moliererms", &shr1shr2moliererms, &b_shr1shr2moliererms);
   fChain->SetBranchAddress("shr1trk1moliereavg", &shr1trk1moliereavg, &b_shr1trk1moliereavg);
   fChain->SetBranchAddress("shr1trk1moliererms", &shr1trk1moliererms, &b_shr1trk1moliererms);
   fChain->SetBranchAddress("shr1trk2moliereavg", &shr1trk2moliereavg, &b_shr1trk2moliereavg);
   fChain->SetBranchAddress("shr1trk2moliererms", &shr1trk2moliererms, &b_shr1trk2moliererms);
   fChain->SetBranchAddress("ismerged", &ismerged, &b_ismerged);
   fChain->SetBranchAddress("merge_bestdot", &merge_bestdot, &b_merge_bestdot);
   fChain->SetBranchAddress("merge_bestdist", &merge_bestdist, &b_merge_bestdist);
   fChain->SetBranchAddress("merge_vtx_x", &merge_vtx_x, &b_merge_vtx_x);
   fChain->SetBranchAddress("merge_vtx_y", &merge_vtx_y, &b_merge_vtx_y);
   fChain->SetBranchAddress("merge_vtx_z", &merge_vtx_z, &b_merge_vtx_z);
   fChain->SetBranchAddress("merge_tk_ipfp", &merge_tk_ipfp, &b_merge_tk_ipfp);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y, &b_shr_tkfit_2cm_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_V", &shr_tkfit_2cm_dedx_V, &b_shr_tkfit_2cm_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_U", &shr_tkfit_2cm_dedx_U, &b_shr_tkfit_2cm_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_Y", &shr_tkfit_2cm_nhits_Y, &b_shr_tkfit_2cm_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_V", &shr_tkfit_2cm_nhits_V, &b_shr_tkfit_2cm_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_U", &shr_tkfit_2cm_nhits_U, &b_shr_tkfit_2cm_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_Y", &shr_tkfit_gap05_dedx_Y, &b_shr_tkfit_gap05_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_V", &shr_tkfit_gap05_dedx_V, &b_shr_tkfit_gap05_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_U", &shr_tkfit_gap05_dedx_U, &b_shr_tkfit_gap05_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_Y", &shr_tkfit_gap05_nhits_Y, &b_shr_tkfit_gap05_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_V", &shr_tkfit_gap05_nhits_V, &b_shr_tkfit_gap05_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_U", &shr_tkfit_gap05_nhits_U, &b_shr_tkfit_gap05_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y, &b_shr_tkfit_gap10_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_V", &shr_tkfit_gap10_dedx_V, &b_shr_tkfit_gap10_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_U", &shr_tkfit_gap10_dedx_U, &b_shr_tkfit_gap10_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_Y", &shr_tkfit_gap10_nhits_Y, &b_shr_tkfit_gap10_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_V", &shr_tkfit_gap10_nhits_V, &b_shr_tkfit_gap10_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_U", &shr_tkfit_gap10_nhits_U, &b_shr_tkfit_gap10_nhits_U);
   fChain->SetBranchAddress("shr_chipr", &shr_chipr, &b_shr_chipr);
   fChain->SetBranchAddress("shr_chimu", &shr_chimu, &b_shr_chimu);
   fChain->SetBranchAddress("shr_bragg_p", &shr_bragg_p, &b_shr_bragg_p);
   fChain->SetBranchAddress("shr_bragg_mu", &shr_bragg_mu, &b_shr_bragg_mu);
   fChain->SetBranchAddress("shr_bragg_mip", &shr_bragg_mip, &b_shr_bragg_mip);
   fChain->SetBranchAddress("shr_bragg_kaon", &shr_bragg_kaon, &b_shr_bragg_kaon);
   fChain->SetBranchAddress("shr_bragg_pion", &shr_bragg_pion, &b_shr_bragg_pion);
   fChain->SetBranchAddress("tksh_distance", &tksh_distance, &b_tksh_distance);
   fChain->SetBranchAddress("tksh_angle", &tksh_angle, &b_tksh_angle);
   fChain->SetBranchAddress("shr_distance", &shr_distance, &b_shr_distance);
   fChain->SetBranchAddress("shr_score", &shr_score, &b_shr_score);
   fChain->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg, &b_shr_bkt_pdg);
   fChain->SetBranchAddress("shr_bkt_purity", &shr_bkt_purity, &b_shr_bkt_purity);
   fChain->SetBranchAddress("shr_bkt_completeness", &shr_bkt_completeness, &b_shr_bkt_completeness);
   fChain->SetBranchAddress("shr_bkt_E", &shr_bkt_E, &b_shr_bkt_E);
   fChain->SetBranchAddress("trk_len", &trk_len, &b_trk_len);
   fChain->SetBranchAddress("trk_theta", &trk_theta, &b_trk_theta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_energy", &trk_energy, &b_trk_energy);
   fChain->SetBranchAddress("trk_energy_muon", &trk_energy_muon, &b_trk_energy_muon);
   fChain->SetBranchAddress("trk_energy_muon_mcs", &trk_energy_muon_mcs, &b_trk_energy_muon_mcs);
   fChain->SetBranchAddress("trk_energy_tot", &trk_energy_tot, &b_trk_energy_tot);
   fChain->SetBranchAddress("trk_energy_muon_tot", &trk_energy_muon_tot, &b_trk_energy_muon_tot);
   fChain->SetBranchAddress("trk_distance", &trk_distance, &b_trk_distance);
   fChain->SetBranchAddress("trk_score", &trk_score, &b_trk_score);
   fChain->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg, &b_trk_bkt_pdg);
   fChain->SetBranchAddress("trk_bkt_purity", &trk_bkt_purity, &b_trk_bkt_purity);
   fChain->SetBranchAddress("trk_bkt_completeness", &trk_bkt_completeness, &b_trk_bkt_completeness);
   fChain->SetBranchAddress("trk_bkt_E", &trk_bkt_E, &b_trk_bkt_E);
   fChain->SetBranchAddress("trk_chipr_best", &trk_chipr_best, &b_trk_chipr_best);
   fChain->SetBranchAddress("trk_chipr_worst", &trk_chipr_worst, &b_trk_chipr_worst);
   fChain->SetBranchAddress("trk_chimu_best", &trk_chimu_best, &b_trk_chimu_best);
   fChain->SetBranchAddress("trk_chimu_worst", &trk_chimu_worst, &b_trk_chimu_worst);
   fChain->SetBranchAddress("trk_chipr", &trk_chipr, &b_trk_chipr);
   fChain->SetBranchAddress("trk_chimu", &trk_chimu, &b_trk_chimu);
   fChain->SetBranchAddress("trk_pida", &trk_pida, &b_trk_pida);
   fChain->SetBranchAddress("trk_bragg_p", &trk_bragg_p, &b_trk_bragg_p);
   fChain->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu, &b_trk_bragg_mu);
   fChain->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip, &b_trk_bragg_mip);
   fChain->SetBranchAddress("trk_bragg_kaon", &trk_bragg_kaon, &b_trk_bragg_kaon);
   fChain->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion, &b_trk_bragg_pion);
   fChain->SetBranchAddress("trk_hits_max", &trk_hits_max, &b_trk_hits_max);
   fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("all_shr_hits", &all_shr_hits, &b_all_shr_hits);
   fChain->SetBranchAddress("all_trk_hits", &all_trk_hits, &b_all_trk_hits);
   fChain->SetBranchAddress("all_shr_energies", &all_shr_energies, &b_all_shr_energies);
   fChain->SetBranchAddress("all_trk_energies", &all_trk_energies, &b_all_trk_energies);
//    fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("trk_hits_2nd", &trk_hits_2nd, &b_trk_hits_2nd);
   fChain->SetBranchAddress("shr_hits_2nd", &shr_hits_2nd, &b_shr_hits_2nd);
   fChain->SetBranchAddress("trkshrhitdist0", &trkshrhitdist0, &b_trkshrhitdist0);
   fChain->SetBranchAddress("trkshrhitdist1", &trkshrhitdist1, &b_trkshrhitdist1);
   fChain->SetBranchAddress("trkshrhitdist2", &trkshrhitdist2, &b_trkshrhitdist2);
   fChain->SetBranchAddress("trk2shrhitdist0", &trk2shrhitdist0, &b_trk2shrhitdist0);
   fChain->SetBranchAddress("trk2shrhitdist1", &trk2shrhitdist1, &b_trk2shrhitdist1);
   fChain->SetBranchAddress("trk2shrhitdist2", &trk2shrhitdist2, &b_trk2shrhitdist2);
   fChain->SetBranchAddress("trk1trk2hitdist0", &trk1trk2hitdist0, &b_trk1trk2hitdist0);
   fChain->SetBranchAddress("trk1trk2hitdist1", &trk1trk2hitdist1, &b_trk1trk2hitdist1);
   fChain->SetBranchAddress("trk1trk2hitdist2", &trk1trk2hitdist2, &b_trk1trk2hitdist2);
   fChain->SetBranchAddress("total_hits_y", &total_hits_y, &b_total_hits_y);
   fChain->SetBranchAddress("extra_energy_y", &extra_energy_y, &b_extra_energy_y);
   fChain->SetBranchAddress("trk_energy_hits_tot", &trk_energy_hits_tot, &b_trk_energy_hits_tot);
   fChain->SetBranchAddress("subcluster", &subcluster, &b_subcluster);
   fChain->SetBranchAddress("shrsubclusters0", &shrsubclusters0, &b_shrsubclusters0);
   fChain->SetBranchAddress("shrsubclusters1", &shrsubclusters1, &b_shrsubclusters1);
   fChain->SetBranchAddress("shrsubclusters2", &shrsubclusters2, &b_shrsubclusters2);
   fChain->SetBranchAddress("shrclusfrac0", &shrclusfrac0, &b_shrclusfrac0);
   fChain->SetBranchAddress("shrclusfrac1", &shrclusfrac1, &b_shrclusfrac1);
   fChain->SetBranchAddress("shrclusfrac2", &shrclusfrac2, &b_shrclusfrac2);
   fChain->SetBranchAddress("shrclusdir0", &shrclusdir0, &b_shrclusdir0);
   fChain->SetBranchAddress("shrclusdir1", &shrclusdir1, &b_shrclusdir1);
   fChain->SetBranchAddress("shrclusdir2", &shrclusdir2, &b_shrclusdir2);
   fChain->SetBranchAddress("shr_hits_tot", &shr_hits_tot, &b_shr_hits_tot);
   fChain->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot, &b_shr_hits_y_tot);
   fChain->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot, &b_shr_hits_u_tot);
   fChain->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot, &b_shr_hits_v_tot);
   fChain->SetBranchAddress("trk_hits_tot", &trk_hits_tot, &b_trk_hits_tot);
   fChain->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot, &b_trk_hits_y_tot);
   fChain->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot, &b_trk_hits_u_tot);
   fChain->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot, &b_trk_hits_v_tot);
   fChain->SetBranchAddress("_elecclusters_U_charge", &_elecclusters_U_charge, &b_elecclusters_U_charge);
   fChain->SetBranchAddress("_elecclusters_V_charge", &_elecclusters_V_charge, &b_elecclusters_V_charge);
   fChain->SetBranchAddress("_elecclusters_Y_charge", &_elecclusters_Y_charge, &b_elecclusters_Y_charge);
   fChain->SetBranchAddress("_elecclusters_U_N", &_elecclusters_U_N, &b_elecclusters_U_N);
   fChain->SetBranchAddress("_elecclusters_V_N", &_elecclusters_V_N, &b_elecclusters_V_N);
   fChain->SetBranchAddress("_elecclusters_Y_N", &_elecclusters_Y_N, &b_elecclusters_Y_N);
   fChain->SetBranchAddress("n_tracks_contained", &n_tracks_contained, &b_n_tracks_contained);
   fChain->SetBranchAddress("n_showers_contained", &n_showers_contained, &b_n_showers_contained);
   fChain->SetBranchAddress("matched_E", &matched_E, &b_matched_E);
   fChain->SetBranchAddress("hits_ratio", &hits_ratio, &b_hits_ratio);
   fChain->SetBranchAddress("contained_fraction", &contained_fraction, &b_contained_fraction);
   fChain->SetBranchAddress("sps_contained_fraction", &sps_contained_fraction, &b_sps_contained_fraction);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pt_assume_muon", &pt_assume_muon, &b_pt_assume_muon);
   fChain->SetBranchAddress("p_assume_muon", &p_assume_muon, &b_p_assume_muon);
   fChain->SetBranchAddress("reco_e", &reco_e, &b_reco_e);
   fChain->SetBranchAddress("dvtx", &dvtx, &b_dvtx);
   fChain->SetBranchAddress("dtrk", &dtrk, &b_dtrk);
   fChain->SetBranchAddress("contained_sps_ratio", &contained_sps_ratio, &b_contained_sps_ratio);
   fChain->SetBranchAddress("dtrk_x_boundary", &dtrk_x_boundary, &b_dtrk_x_boundary);
   fChain->SetBranchAddress("dtrk_y_boundary", &dtrk_y_boundary, &b_dtrk_y_boundary);
   fChain->SetBranchAddress("dtrk_z_boundary", &dtrk_z_boundary, &b_dtrk_z_boundary);
   fChain->SetBranchAddress("dshr_x_boundary", &dshr_x_boundary, &b_dshr_x_boundary);
   fChain->SetBranchAddress("dshr_y_boundary", &dshr_y_boundary, &b_dshr_y_boundary);
   fChain->SetBranchAddress("dshr_z_boundary", &dshr_z_boundary, &b_dshr_z_boundary);
   fChain->SetBranchAddress("dvtx_x_boundary", &dvtx_x_boundary, &b_dvtx_x_boundary);
   fChain->SetBranchAddress("dvtx_y_boundary", &dvtx_y_boundary, &b_dvtx_y_boundary);
   fChain->SetBranchAddress("dvtx_z_boundary", &dvtx_z_boundary, &b_dvtx_z_boundary);
   fChain->SetBranchAddress("dtrk_boundary", &dtrk_boundary, &b_dtrk_boundary);
   fChain->SetBranchAddress("dvtx_boundary", &dvtx_boundary, &b_dvtx_boundary);
   fChain->SetBranchAddress("dshr_boundary", &dshr_boundary, &b_dshr_boundary);
   fChain->SetBranchAddress("dmc_boundary", &dmc_boundary, &b_dmc_boundary);
   fChain->SetBranchAddress("leeweight", &leeweight, &b_leeweight);
   fChain->SetBranchAddress("true_pt", &true_pt, &b_true_pt);
   fChain->SetBranchAddress("true_pt_visible", &true_pt_visible, &b_true_pt_visible);
   fChain->SetBranchAddress("true_p", &true_p, &b_true_p);
   fChain->SetBranchAddress("true_p_visible", &true_p_visible, &b_true_p_visible);
   fChain->SetBranchAddress("true_e_visible", &true_e_visible, &b_true_e_visible);
   fChain->SetBranchAddress("_opfilter_pe_beam", &_opfilter_pe_beam, &b_opfilter_pe_beam);
   fChain->SetBranchAddress("_opfilter_pe_veto", &_opfilter_pe_veto, &b_opfilter_pe_veto);
   fChain->SetBranchAddress("nu_pdg", &nu_pdg, &b_nu_pdg);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   fChain->SetBranchAddress("nu_parent_pdg", &nu_parent_pdg, &b_nu_parent_pdg);
   fChain->SetBranchAddress("nu_hadron_pdg", &nu_hadron_pdg, &b_nu_hadron_pdg);
   fChain->SetBranchAddress("nu_decay_mode", &nu_decay_mode, &b_nu_decay_mode);
   fChain->SetBranchAddress("interaction", &interaction, &b_interaction);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_nu_e);
   fChain->SetBranchAddress("nu_l", &nu_l, &b_nu_l);
   fChain->SetBranchAddress("nu_pt", &nu_pt, &b_nu_pt);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("isVtxInFiducial", &isVtxInFiducial, &b_isVtxInFiducial);
   fChain->SetBranchAddress("truthFiducial", &truthFiducial, &b_truthFiducial);
   fChain->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t, &b_true_nu_vtx_t);
   fChain->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x, &b_true_nu_vtx_x);
   fChain->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y, &b_true_nu_vtx_y);
   fChain->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z, &b_true_nu_vtx_z);
   fChain->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x, &b_true_nu_vtx_sce_x);
   fChain->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y, &b_true_nu_vtx_sce_y);
   fChain->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z, &b_true_nu_vtx_sce_z);
   fChain->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x, &b_reco_nu_vtx_x);
   fChain->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y, &b_reco_nu_vtx_y);
   fChain->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z, &b_reco_nu_vtx_z);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_e", &muon_e, &b_muon_e);
   fChain->SetBranchAddress("muon_c", &muon_c, &b_muon_c);
   fChain->SetBranchAddress("muon_p", &muon_p, &b_muon_p);
   fChain->SetBranchAddress("nelec", &nelec, &b_nelec);
   fChain->SetBranchAddress("elec_e", &elec_e, &b_elec_e);
   fChain->SetBranchAddress("elec_c", &elec_c, &b_elec_c);
   fChain->SetBranchAddress("elec_p", &elec_p, &b_elec_p);
   fChain->SetBranchAddress("elec_vx", &elec_vx, &b_elec_vx);
   fChain->SetBranchAddress("elec_vy", &elec_vy, &b_elec_vy);
   fChain->SetBranchAddress("elec_vz", &elec_vz, &b_elec_vz);
   fChain->SetBranchAddress("elec_px", &elec_px, &b_elec_px);
   fChain->SetBranchAddress("elec_py", &elec_py, &b_elec_py);
   fChain->SetBranchAddress("elec_pz", &elec_pz, &b_elec_pz);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("pi0_e", &pi0_e, &b_pi0_e);
   fChain->SetBranchAddress("pi0_c", &pi0_c, &b_pi0_c);
   fChain->SetBranchAddress("pi0_p", &pi0_p, &b_pi0_p);
   fChain->SetBranchAddress("nneutron", &nneutron, &b_nneutron);
   fChain->SetBranchAddress("nproton", &nproton, &b_nproton);
   fChain->SetBranchAddress("proton_e", &proton_e, &b_proton_e);
   fChain->SetBranchAddress("proton_c", &proton_c, &b_proton_c);
   fChain->SetBranchAddress("proton_p", &proton_p, &b_proton_p);
   fChain->SetBranchAddress("npion", &npion, &b_npion);
   fChain->SetBranchAddress("pion_e", &pion_e, &b_pion_e);
   fChain->SetBranchAddress("pion_c", &pion_c, &b_pion_c);
   fChain->SetBranchAddress("pion_p", &pion_p, &b_pion_p);
   fChain->SetBranchAddress("neta", &neta, &b_neta);
   fChain->SetBranchAddress("eta_e", &eta_e, &b_eta_e);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("crtveto", &crtveto, &b_crtveto);
   fChain->SetBranchAddress("crthitpe", &crthitpe, &b_crthitpe);
   fChain->SetBranchAddress("pfp_slice_idx", &pfp_slice_idx, &b_pfp_slice_idx);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_tid", &backtracked_tid, &b_backtracked_tid);
   fChain->SetBranchAddress("backtracked_purity", &backtracked_purity, &b_backtracked_purity);
   fChain->SetBranchAddress("backtracked_completeness", &backtracked_completeness, &b_backtracked_completeness);
   fChain->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity, &b_backtracked_overlay_purity);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("backtracked_start_x", &backtracked_start_x, &b_backtracked_start_x);
   fChain->SetBranchAddress("backtracked_start_y", &backtracked_start_y, &b_backtracked_start_y);
   fChain->SetBranchAddress("backtracked_start_z", &backtracked_start_z, &b_backtracked_start_z);
   fChain->SetBranchAddress("backtracked_start_t", &backtracked_start_t, &b_backtracked_start_t);
   fChain->SetBranchAddress("backtracked_start_U", &backtracked_start_U, &b_backtracked_start_U);
   fChain->SetBranchAddress("backtracked_start_V", &backtracked_start_V, &b_backtracked_start_V);
   fChain->SetBranchAddress("backtracked_start_Y", &backtracked_start_Y, &b_backtracked_start_Y);
   fChain->SetBranchAddress("backtracked_sce_start_x", &backtracked_sce_start_x, &b_backtracked_sce_start_x);
   fChain->SetBranchAddress("backtracked_sce_start_y", &backtracked_sce_start_y, &b_backtracked_sce_start_y);
   fChain->SetBranchAddress("backtracked_sce_start_z", &backtracked_sce_start_z, &b_backtracked_sce_start_z);
   fChain->SetBranchAddress("backtracked_sce_start_U", &backtracked_sce_start_U, &b_backtracked_sce_start_U);
   fChain->SetBranchAddress("backtracked_sce_start_V", &backtracked_sce_start_V, &b_backtracked_sce_start_V);
   fChain->SetBranchAddress("backtracked_sce_start_Y", &backtracked_sce_start_Y, &b_backtracked_sce_start_Y);
   fChain->SetBranchAddress("lep_e", &lep_e, &b_lep_e);
   fChain->SetBranchAddress("pass", &pass, &b_pass);
   fChain->SetBranchAddress("swtrig", &swtrig, &b_swtrig);
   fChain->SetBranchAddress("evnhits", &evnhits, &b_evnhits);
   fChain->SetBranchAddress("slpdg", &slpdg, &b_slpdg);
   fChain->SetBranchAddress("slnhits", &slnhits, &b_slnhits);
   fChain->SetBranchAddress("n_pfps", &n_pfps, &b_n_pfps);
   fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("pfp_generation_v", &pfp_generation_v, &b_pfp_generation_v);
   fChain->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v, &b_pfp_trk_daughters_v);
   fChain->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v, &b_pfp_shr_daughters_v);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("pfpplanesubclusters_U", &pfpplanesubclusters_U, &b_pfpplanesubclusters_U);
   fChain->SetBranchAddress("pfpplanesubclusters_V", &pfpplanesubclusters_V, &b_pfpplanesubclusters_V);
   fChain->SetBranchAddress("pfpplanesubclusters_Y", &pfpplanesubclusters_Y, &b_pfpplanesubclusters_Y);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_U", &pfpplanesubhitfracmax_U, &b_pfpplanesubhitfracmax_U);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_V", &pfpplanesubhitfracmax_V, &b_pfpplanesubhitfracmax_V);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_Y", &pfpplanesubhitfracmax_Y, &b_pfpplanesubhitfracmax_Y);
   fChain->SetBranchAddress("hits_u", &hits_u, &b_hits_u);
   fChain->SetBranchAddress("hits_v", &hits_v, &b_hits_v);
   fChain->SetBranchAddress("hits_y", &hits_y, &b_hits_y);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("slclustfrac", &slclustfrac, &b_slclustfrac);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_vx", &mc_vx, &b_mc_vx);
   fChain->SetBranchAddress("mc_vy", &mc_vy, &b_mc_vy);
   fChain->SetBranchAddress("mc_vz", &mc_vz, &b_mc_vz);
   fChain->SetBranchAddress("mc_endx", &mc_endx, &b_mc_endx);
   fChain->SetBranchAddress("mc_endy", &mc_endy, &b_mc_endy);
   fChain->SetBranchAddress("mc_endz", &mc_endz, &b_mc_endz);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_completeness", &mc_completeness, &b_mc_completeness);
   fChain->SetBranchAddress("mc_purity", &mc_purity, &b_mc_purity);
   fChain->SetBranchAddress("endmuonprocess", &endmuonprocess, &b_endmuonprocess);
   fChain->SetBranchAddress("endmuonmichel", &endmuonmichel, &b_endmuonmichel);
   fChain->SetBranchAddress("flash_pe", &flash_pe, &b_flash_pe);
   fChain->SetBranchAddress("flash_pe_v", &flash_pe_v, &b_flash_pe_v);
   fChain->SetBranchAddress("slice_pe_v", &slice_pe_v, &b_slice_pe_v);
   fChain->SetBranchAddress("flash_time", &flash_time, &b_flash_time);
   fChain->SetBranchAddress("flash_y", &flash_y, &b_flash_y);
   fChain->SetBranchAddress("flash_z", &flash_z, &b_flash_z);
   fChain->SetBranchAddress("flash_timewidth", &flash_timewidth, &b_flash_timewidth);
   fChain->SetBranchAddress("flash_ywidth", &flash_ywidth, &b_flash_ywidth);
   fChain->SetBranchAddress("flash_zwidth", &flash_zwidth, &b_flash_zwidth);
   fChain->SetBranchAddress("nu_flashmatch_score", &nu_flashmatch_score, &b_nu_flashmatch_score);
   fChain->SetBranchAddress("nu_centerX", &nu_centerX, &b_nu_centerX);
   fChain->SetBranchAddress("nu_centerY", &nu_centerY, &b_nu_centerY);
   fChain->SetBranchAddress("nu_centerZ", &nu_centerZ, &b_nu_centerZ);
   fChain->SetBranchAddress("nu_totalCharge", &nu_totalCharge, &b_nu_totalCharge);
   fChain->SetBranchAddress("best_cosmic_flashmatch_score", &best_cosmic_flashmatch_score, &b_best_cosmic_flashmatch_score);
   fChain->SetBranchAddress("best_obviouscosmic_flashmatch_score", &best_obviouscosmic_flashmatch_score, &b_best_obviouscosmic_flashmatch_score);
   fChain->SetBranchAddress("cosmic_flashmatch_score_v", &cosmic_flashmatch_score_v, &b_cosmic_flashmatch_score_v);
   fChain->SetBranchAddress("cosmic_topological_score_v", &cosmic_topological_score_v, &b_cosmic_topological_score_v);
   fChain->SetBranchAddress("cosmic_centerX_v", &cosmic_centerX_v, &b_cosmic_centerX_v);
   fChain->SetBranchAddress("cosmic_centerY_v", &cosmic_centerY_v, &b_cosmic_centerY_v);
   fChain->SetBranchAddress("cosmic_centerZ_v", &cosmic_centerZ_v, &b_cosmic_centerZ_v);
   fChain->SetBranchAddress("cosmic_totalCharge_v", &cosmic_totalCharge_v, &b_cosmic_totalCharge_v);
   fChain->SetBranchAddress("cosmic_nhits_v", &cosmic_nhits_v, &b_cosmic_nhits_v);
   fChain->SetBranchAddress("cosmic_nunhits_v", &cosmic_nunhits_v, &b_cosmic_nunhits_v);
   fChain->SetBranchAddress("cosmic_isclear_v", &cosmic_isclear_v, &b_cosmic_isclear_v);
//    fChain->SetBranchAddress("flash_pe", &flash_pe, &b_flash_pe);
   fChain->SetBranchAddress("flash_pe_calib", &flash_pe_calib, &b_flash_pe_calib);
//    fChain->SetBranchAddress("flash_time", &flash_time, &b_flash_time);
   fChain->SetBranchAddress("flash_zcenter", &flash_zcenter, &b_flash_zcenter);
   fChain->SetBranchAddress("flash_ycenter", &flash_ycenter, &b_flash_ycenter);
//    fChain->SetBranchAddress("flash_zwidth", &flash_zwidth, &b_flash_zwidth);
//    fChain->SetBranchAddress("flash_ywidth", &flash_ywidth, &b_flash_ywidth);
//    fChain->SetBranchAddress("flash_pe_v", &flash_pe_v, &b_flash_pe_v);
   fChain->SetBranchAddress("flash_pe_calib_v", &flash_pe_calib_v, &b_flash_pe_calib_v);
   fChain->SetBranchAddress("gain_area_v", &gain_area_v, &b_gain_area_v);
   fChain->SetBranchAddress("gain_ampl_v", &gain_ampl_v, &b_gain_ampl_v);
   fChain->SetBranchAddress("waveform_00", &waveform_00, &b_waveform_00);
   fChain->SetBranchAddress("waveform_01", &waveform_01, &b_waveform_01);
   fChain->SetBranchAddress("waveform_02", &waveform_02, &b_waveform_02);
   fChain->SetBranchAddress("waveform_03", &waveform_03, &b_waveform_03);
   fChain->SetBranchAddress("waveform_04", &waveform_04, &b_waveform_04);
   fChain->SetBranchAddress("waveform_05", &waveform_05, &b_waveform_05);
   fChain->SetBranchAddress("waveform_06", &waveform_06, &b_waveform_06);
   fChain->SetBranchAddress("waveform_07", &waveform_07, &b_waveform_07);
   fChain->SetBranchAddress("waveform_08", &waveform_08, &b_waveform_08);
   fChain->SetBranchAddress("waveform_09", &waveform_09, &b_waveform_09);
   fChain->SetBranchAddress("waveform_10", &waveform_10, &b_waveform_10);
   fChain->SetBranchAddress("waveform_11", &waveform_11, &b_waveform_11);
   fChain->SetBranchAddress("waveform_12", &waveform_12, &b_waveform_12);
   fChain->SetBranchAddress("waveform_13", &waveform_13, &b_waveform_13);
   fChain->SetBranchAddress("waveform_14", &waveform_14, &b_waveform_14);
   fChain->SetBranchAddress("waveform_15", &waveform_15, &b_waveform_15);
   fChain->SetBranchAddress("waveform_16", &waveform_16, &b_waveform_16);
   fChain->SetBranchAddress("waveform_17", &waveform_17, &b_waveform_17);
   fChain->SetBranchAddress("waveform_18", &waveform_18, &b_waveform_18);
   fChain->SetBranchAddress("waveform_19", &waveform_19, &b_waveform_19);
   fChain->SetBranchAddress("waveform_20", &waveform_20, &b_waveform_20);
   fChain->SetBranchAddress("waveform_21", &waveform_21, &b_waveform_21);
   fChain->SetBranchAddress("waveform_22", &waveform_22, &b_waveform_22);
   fChain->SetBranchAddress("waveform_23", &waveform_23, &b_waveform_23);
   fChain->SetBranchAddress("waveform_24", &waveform_24, &b_waveform_24);
   fChain->SetBranchAddress("waveform_25", &waveform_25, &b_waveform_25);
   fChain->SetBranchAddress("waveform_26", &waveform_26, &b_waveform_26);
   fChain->SetBranchAddress("waveform_27", &waveform_27, &b_waveform_27);
   fChain->SetBranchAddress("waveform_28", &waveform_28, &b_waveform_28);
   fChain->SetBranchAddress("waveform_29", &waveform_29, &b_waveform_29);
   fChain->SetBranchAddress("waveform_30", &waveform_30, &b_waveform_30);
   fChain->SetBranchAddress("waveform_31", &waveform_31, &b_waveform_31);
   fChain->SetBranchAddress("secondshower_U_charge", &secondshower_U_charge, &b_secondshower_U_charge);
   fChain->SetBranchAddress("secondshower_U_nhit", &secondshower_U_nhit, &b_secondshower_U_nhit);
   fChain->SetBranchAddress("secondshower_U_vtxdist", &secondshower_U_vtxdist, &b_secondshower_U_vtxdist);
   fChain->SetBranchAddress("secondshower_U_eigenratio", &secondshower_U_eigenratio, &b_secondshower_U_eigenratio);
   fChain->SetBranchAddress("secondshower_U_dot", &secondshower_U_dot, &b_secondshower_U_dot);
   fChain->SetBranchAddress("secondshower_U_dir", &secondshower_U_dir, &b_secondshower_U_dir);
   fChain->SetBranchAddress("secondshower_V_charge", &secondshower_V_charge, &b_secondshower_V_charge);
   fChain->SetBranchAddress("secondshower_V_nhit", &secondshower_V_nhit, &b_secondshower_V_nhit);
   fChain->SetBranchAddress("secondshower_V_vtxdist", &secondshower_V_vtxdist, &b_secondshower_V_vtxdist);
   fChain->SetBranchAddress("secondshower_V_eigenratio", &secondshower_V_eigenratio, &b_secondshower_V_eigenratio);
   fChain->SetBranchAddress("secondshower_V_dot", &secondshower_V_dot, &b_secondshower_V_dot);
   fChain->SetBranchAddress("secondshower_V_dir", &secondshower_V_dir, &b_secondshower_V_dir);
   fChain->SetBranchAddress("secondshower_Y_charge", &secondshower_Y_charge, &b_secondshower_Y_charge);
   fChain->SetBranchAddress("secondshower_Y_nhit", &secondshower_Y_nhit, &b_secondshower_Y_nhit);
   fChain->SetBranchAddress("secondshower_Y_vtxdist", &secondshower_Y_vtxdist, &b_secondshower_Y_vtxdist);
   fChain->SetBranchAddress("secondshower_Y_eigenratio", &secondshower_Y_eigenratio, &b_secondshower_Y_eigenratio);
   fChain->SetBranchAddress("secondshower_Y_dot", &secondshower_Y_dot, &b_secondshower_Y_dot);
   fChain->SetBranchAddress("secondshower_Y_dir", &secondshower_Y_dir, &b_secondshower_Y_dir);
   fChain->SetBranchAddress("simphoton_number_v", &simphoton_number_v, &b_simphoton_number_v);
   fChain->SetBranchAddress("simphoton_tmin_v", &simphoton_tmin_v, &b_simphoton_tmin_v);
   fChain->SetBranchAddress("simphoton_tmax_v", &simphoton_tmax_v, &b_simphoton_tmax_v);
   fChain->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v, &b_trk_bragg_p_v);
   fChain->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v, &b_trk_bragg_mu_v);
   fChain->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v, &b_trk_bragg_mip_v);
   fChain->SetBranchAddress("trk_pida_v", &trk_pida_v, &b_trk_pida_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_pid_chipi_v", &trk_pid_chipi_v, &b_trk_pid_chipi_v);
   fChain->SetBranchAddress("trk_pid_chika_v", &trk_pid_chika_v, &b_trk_pid_chika_v);
   fChain->SetBranchAddress("trk_pid_chimu_v", &trk_pid_chimu_v, &b_trk_pid_chimu_v);
   fChain->SetBranchAddress("trk_bragg_p_u_v", &trk_bragg_p_u_v, &b_trk_bragg_p_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_u_v", &trk_bragg_mu_u_v, &b_trk_bragg_mu_u_v);
   fChain->SetBranchAddress("trk_bragg_mip_u_v", &trk_bragg_mip_u_v, &b_trk_bragg_mip_u_v);
   fChain->SetBranchAddress("trk_pida_u_v", &trk_pida_u_v, &b_trk_pida_u_v);
   fChain->SetBranchAddress("trk_pid_chipr_u_v", &trk_pid_chipr_u_v, &b_trk_pid_chipr_u_v);
   fChain->SetBranchAddress("trk_pid_chipi_u_v", &trk_pid_chipi_u_v, &b_trk_pid_chipi_u_v);
   fChain->SetBranchAddress("trk_pid_chika_u_v", &trk_pid_chika_u_v, &b_trk_pid_chika_u_v);
   fChain->SetBranchAddress("trk_pid_chimu_u_v", &trk_pid_chimu_u_v, &b_trk_pid_chimu_u_v);
   fChain->SetBranchAddress("trk_bragg_p_v_v", &trk_bragg_p_v_v, &b_trk_bragg_p_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_v_v", &trk_bragg_mu_v_v, &b_trk_bragg_mu_v_v);
   fChain->SetBranchAddress("trk_bragg_mip_v_v", &trk_bragg_mip_v_v, &b_trk_bragg_mip_v_v);
   fChain->SetBranchAddress("trk_pida_v_v", &trk_pida_v_v, &b_trk_pida_v_v);
   fChain->SetBranchAddress("trk_pid_chipr_v_v", &trk_pid_chipr_v_v, &b_trk_pid_chipr_v_v);
   fChain->SetBranchAddress("trk_pid_chipi_v_v", &trk_pid_chipi_v_v, &b_trk_pid_chipi_v_v);
   fChain->SetBranchAddress("trk_pid_chika_v_v", &trk_pid_chika_v_v, &b_trk_pid_chika_v_v);
   fChain->SetBranchAddress("trk_pid_chimu_v_v", &trk_pid_chimu_v_v, &b_trk_pid_chimu_v_v);
   fChain->SetBranchAddress("trk_pfp_id_v", &trk_pfp_id_v, &b_trk_pfp_id_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_start_x_v", &trk_start_x_v, &b_trk_start_x_v);
   fChain->SetBranchAddress("trk_start_y_v", &trk_start_y_v, &b_trk_start_y_v);
   fChain->SetBranchAddress("trk_start_z_v", &trk_start_z_v, &b_trk_start_z_v);
   fChain->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v, &b_trk_sce_start_x_v);
   fChain->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v, &b_trk_sce_start_y_v);
   fChain->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v, &b_trk_sce_start_z_v);
   fChain->SetBranchAddress("trk_end_x_v", &trk_end_x_v, &b_trk_end_x_v);
   fChain->SetBranchAddress("trk_end_y_v", &trk_end_y_v, &b_trk_end_y_v);
   fChain->SetBranchAddress("trk_end_z_v", &trk_end_z_v, &b_trk_end_z_v);
   fChain->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v, &b_trk_sce_end_x_v);
   fChain->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v, &b_trk_sce_end_y_v);
   fChain->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v, &b_trk_sce_end_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_theta_v", &trk_theta_v, &b_trk_theta_v);
   fChain->SetBranchAddress("trk_phi_v", &trk_phi_v, &b_trk_phi_v);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v, &b_trk_range_muon_mom_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_energy_muon_v", &trk_energy_muon_v, &b_trk_energy_muon_v);
   fChain->SetBranchAddress("trk_calo_energy_u_v", &trk_calo_energy_u_v, &b_trk_calo_energy_u_v);
   fChain->SetBranchAddress("trk_calo_energy_v_v", &trk_calo_energy_v_v, &b_trk_calo_energy_v_v);
   fChain->SetBranchAddress("trk_calo_energy_y_v", &trk_calo_energy_y_v, &b_trk_calo_energy_y_v);
   fChain->SetBranchAddress("trk_llr_pid_u_v", &trk_llr_pid_u_v, &b_trk_llr_pid_u_v);
   fChain->SetBranchAddress("trk_llr_pid_v_v", &trk_llr_pid_v_v, &b_trk_llr_pid_v_v);
   fChain->SetBranchAddress("trk_llr_pid_y_v", &trk_llr_pid_y_v, &b_trk_llr_pid_y_v);
   fChain->SetBranchAddress("trk_llr_pid_v", &trk_llr_pid_v, &b_trk_llr_pid_v);
   fChain->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v, &b_trk_llr_pid_score_v);
   Notify();
}

Bool_t mu_sim_mely_class::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mu_sim_mely_class::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mu_sim_mely_class::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mu_sim_mely_class_cxx
