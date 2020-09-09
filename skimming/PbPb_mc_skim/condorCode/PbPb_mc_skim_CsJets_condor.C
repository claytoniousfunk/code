#include "/uscms/home/cmbennet/work/myProcessesEdit/hiforest/plugin/eventMap_hiForest.h"
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <dirent.h>  
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h>


double getPtRel(double MuonPt, double MuonEta, double MuonPhi, double JetPt, double JetEta, double JetPhi){

double Muon_Px = MuonPt*TMath::Cos(MuonPhi);
double Muon_Py = MuonPt*TMath::Sin(MuonPhi);
double Muon_Pz = MuonPt*TMath::SinH(MuonEta);

double Jet_Px = JetPt*TMath::Cos(JetPhi);
double Jet_Py = JetPt*TMath::Sin(JetPhi);
double Jet_Pz = JetPt*TMath::SinH(JetEta);


float lj_x = Jet_Px;
float lj_y = Jet_Py;
float lj_z = Jet_Pz;

// absolute values squared
float lj2 = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;

//float lep2 = lep.px()*lep.px()+lep.py()*lep.py()+lep.pz()*lep.pz();
float lep2 = Muon_Px*Muon_Px + Muon_Py*Muon_Py+Muon_Pz*Muon_Pz;

// projection vec(mu) to lepjet axis
float lepXlj = Muon_Px*lj_x+ Muon_Py*lj_y + Muon_Pz*lj_z;

// absolute value squared and normalized
float pLrel2 = lepXlj*lepXlj/lj2;
float pTrel2 = lep2-pLrel2;

return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}

void PbPb_mc_skim_CsJets_condor(TString input, TString output){
// Define histograms to be filled with data
const int NPhiBins = 300;
const double phiMin = -TMath::Pi();
const double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_cent0to10 = new TH1D("h_jetphi_cent0to10","reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_cent10to30 = new TH1D("h_jetphi_cent10to30","reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_cent30to50 = new TH1D("h_jetphi_cent30to50","reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_cent50to90 = new TH1D("h_jetphi_cent50to90","reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_g = new TH1D("h_jetphi_g","g reco jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_g_cent0to10 = new TH1D("h_jetphi_g_cent0to10","g reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_g_cent10to30 = new TH1D("h_jetphi_g_cent10to30","g reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_g_cent30to50 = new TH1D("h_jetphi_g_cent30to50","g reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_g_cent50to90 = new TH1D("h_jetphi_g_cent50to90","g reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_uds = new TH1D("h_jetphi_uds","uds reco jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_uds_cent0to10 = new TH1D("h_jetphi_uds_cent0to10","uds reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_uds_cent10to30 = new TH1D("h_jetphi_uds_cent10to30","uds reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_uds_cent30to50 = new TH1D("h_jetphi_uds_cent30to50","uds reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_uds_cent50to90 = new TH1D("h_jetphi_uds_cent50to90","uds reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_b = new TH1D("h_jetphi_b","b reco jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_b_cent0to10 = new TH1D("h_jetphi_b_cent0to10","b reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_b_cent10to30 = new TH1D("h_jetphi_b_cent10to30","b reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_b_cent30to50 = new TH1D("h_jetphi_b_cent30to50","b reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_b_cent50to90 = new TH1D("h_jetphi_b_cent50to90","b reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_c = new TH1D("h_jetphi_c","c reco jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_c_cent0to10 = new TH1D("h_jetphi_c_cent0to10","c reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_c_cent10to30 = new TH1D("h_jetphi_c_cent10to30","c reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_c_cent30to50 = new TH1D("h_jetphi_c_cent30to50","c reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_c_cent50to90 = new TH1D("h_jetphi_c_cent50to90","c reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_ghost = new TH1D("h_jetphi_ghost","ghost reco jet #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_ghost_cent0to10 = new TH1D("h_jetphi_ghost_cent0to10","ghost reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_ghost_cent10to30 = new TH1D("h_jetphi_ghost_cent10to30","ghost reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_ghost_cent30to50 = new TH1D("h_jetphi_ghost_cent30to50","ghost reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_jetphi_ghost_cent50to90 = new TH1D("h_jetphi_ghost_cent50to90","ghost reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
TH1D *h_mujetphi = new TH1D("h_mujetphi","reco muon jet #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_cent0to10 = new TH1D("h_mujetphi_cent0to10","reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_cent10to30 = new TH1D("h_mujetphi_cent10to30","reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_cent30to50 = new TH1D("h_mujetphi_cent30to50","reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_cent50to90 = new TH1D("h_mujetphi_cent50to90","reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_g = new TH1D("h_mujetphi_g","g reco muon jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_g_cent0to10 = new TH1D("h_mujetphi_g_cent0to10","g reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_g_cent10to30 = new TH1D("h_mujetphi_g_cent10to30","g reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_g_cent30to50 = new TH1D("h_mujetphi_g_cent30to50","g reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_g_cent50to90 = new TH1D("h_mujetphi_g_cent50to90","g reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_uds = new TH1D("h_mujetphi_uds","uds reco muon jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_uds_cent0to10 = new TH1D("h_mujetphi_uds_cent0to10","uds reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_uds_cent10to30 = new TH1D("h_mujetphi_uds_cent10to30","uds reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_uds_cent30to50 = new TH1D("h_mujetphi_uds_cent30to50","uds reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_uds_cent50to90 = new TH1D("h_mujetphi_uds_cent50to90","uds reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_b = new TH1D("h_mujetphi_b","b reco muon jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_b_cent0to10 = new TH1D("h_mujetphi_b_cent0to10","b reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_b_cent10to30 = new TH1D("h_mujetphi_b_cent10to30","b reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_b_cent30to50 = new TH1D("h_mujetphi_b_cent30to50","b reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_b_cent50to90 = new TH1D("h_mujetphi_b_cent50to90","b reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_c = new TH1D("h_mujetphi_c","c reco muon jet #phi", NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_c_cent0to10 = new TH1D("h_mujetphi_c_cent0to10","c reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_c_cent10to30 = new TH1D("h_mujetphi_c_cent10to30","c reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_c_cent30to50 = new TH1D("h_mujetphi_c_cent30to50","c reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_c_cent50to90 = new TH1D("h_mujetphi_c_cent50to90","c reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_ghost = new TH1D("h_mujetphi_ghost","ghost reco jet #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_ghost_cent0to10 = new TH1D("h_mujetphi_ghost_cent0to10","ghost reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_ghost_cent10to30 = new TH1D("h_mujetphi_ghost_cent10to30","ghost reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_ghost_cent30to50 = new TH1D("h_mujetphi_ghost_cent30to50","ghost reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_mujetphi_ghost_cent50to90 = new TH1D("h_mujetphi_ghost_cent50to90","ghost reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
TH1D *h_genjetphi = new TH1D("h_genjetphi","gen jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_trkPhi = new TH1D("h_trkPhi","track #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_trkPhi_cent0to10 = new TH1D("h_trkPhi_cent0to10","track #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_trkPhi_cent10to30 = new TH1D("h_trkPhi_cent10to30","track #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_trkPhi_cent30to50 = new TH1D("h_trkPhi_cent30to50","track #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_trkPhi_cent50to90 = new TH1D("h_trkPhi_cent50to90","track #phi, 50-90 %",NPhiBins,phiMin,phiMax);
const int NEtaBins = 260;
const double etaMin = -1.3;
const double etaMax = 1.3;
const int NTrkEtaBins = 480;
const double trkEtaMin = -2.4;
const double trkEtaMax = 2.4;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_cent0to10 = new TH1D("h_jeteta_cent0to10","reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_cent10to30 = new TH1D("h_jeteta_cent10to30","reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_cent30to50 = new TH1D("h_jeteta_cent30to50","reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_cent50to90 = new TH1D("h_jeteta_cent50to90","reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_g = new TH1D("h_jeteta_g","g reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_g_cent0to10 = new TH1D("h_jeteta_g_cent0to10","g reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_g_cent10to30 = new TH1D("h_jeteta_g_cent10to30","g reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_g_cent30to50 = new TH1D("h_jeteta_g_cent30to50","g reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_g_cent50to90 = new TH1D("h_jeteta_g_cent50to90","g reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_uds = new TH1D("h_jeteta_uds","uds reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_uds_cent0to10 = new TH1D("h_jeteta_uds_cent0to10","uds reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_uds_cent10to30 = new TH1D("h_jeteta_uds_cent10to30","uds reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_uds_cent30to50 = new TH1D("h_jeteta_uds_cent30to50","uds reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_uds_cent50to90 = new TH1D("h_jeteta_uds_cent50to90","uds reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_b = new TH1D("h_jeteta_b","b reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_b_cent0to10 = new TH1D("h_jeteta_b_cent0to10","b reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_b_cent10to30 = new TH1D("h_jeteta_b_cent10to30","b reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_b_cent30to50 = new TH1D("h_jeteta_b_cent30to50","b reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_b_cent50to90 = new TH1D("h_jeteta_b_cent50to90","b reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_c = new TH1D("h_jeteta_c","c reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_c_cent0to10 = new TH1D("h_jeteta_c_cent0to10","c reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_c_cent10to30 = new TH1D("h_jeteta_c_cent10to30","c reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_c_cent30to50 = new TH1D("h_jeteta_c_cent30to50","c reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_c_cent50to90 = new TH1D("h_jeteta_c_cent50to90","c reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_ghost = new TH1D("h_jeteta_ghost","ghost reco jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_ghost_cent0to10 = new TH1D("h_jeteta_ghost_cent0to10","g reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_ghost_cent10to30 = new TH1D("h_jeteta_ghost_cent10to30","g reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_ghost_cent30to50 = new TH1D("h_jeteta_ghost_cent30to50","g reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_jeteta_ghost_cent50to90 = new TH1D("h_jeteta_ghost_cent50to90","g reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
TH1D *h_mujeteta = new TH1D("h_mujeteta","reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_cent0to10 = new TH1D("h_mujeteta_cent0to10","reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_cent10to30 = new TH1D("h_mujeteta_cent10to30","reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_cent30to50 = new TH1D("h_mujeteta_cent30to50","reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_cent50to90 = new TH1D("h_mujeteta_cent50to90","reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_g = new TH1D("h_mujeteta_g","g reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_g_cent0to10 = new TH1D("h_mujeteta_g_cent0to10","g reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_g_cent10to30 = new TH1D("h_mujeteta_g_cent10to30","g reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_g_cent30to50 = new TH1D("h_mujeteta_g_cent30to50","g reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_g_cent50to90 = new TH1D("h_mujeteta_g_cent50to90","g reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_uds = new TH1D("h_mujeteta_uds","uds reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_uds_cent0to10 = new TH1D("h_mujeteta_uds_cent0to10","uds reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_uds_cent10to30 = new TH1D("h_mujeteta_uds_cent10to30","uds reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_uds_cent30to50 = new TH1D("h_mujeteta_uds_cent30to50","uds reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_uds_cent50to90 = new TH1D("h_mujeteta_uds_cent50to90","uds reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_b = new TH1D("h_mujeteta_b","b reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_b_cent0to10 = new TH1D("h_mujeteta_b_cent0to10","b reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_b_cent10to30 = new TH1D("h_mujeteta_b_cent10to30","b reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_b_cent30to50 = new TH1D("h_mujeteta_b_cent30to50","b reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_b_cent50to90 = new TH1D("h_mujeteta_b_cent50to90","b reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_c = new TH1D("h_mujeteta_c","c reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_c_cent0to10 = new TH1D("h_mujeteta_c_cent0to10","c reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_c_cent10to30 = new TH1D("h_mujeteta_c_cent10to30","c reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_c_cent30to50 = new TH1D("h_mujeteta_c_cent30to50","c reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_c_cent50to90 = new TH1D("h_mujeteta_c_cent50to90","c reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_ghost = new TH1D("h_mujeteta_ghost","ghost reco muon jet #eta",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_ghost_cent0to10 = new TH1D("h_mujeteta_ghost_cent0to10","ghost reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_ghost_cent10to30 = new TH1D("h_mujeteta_ghost_cent10to30","ghost reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_ghost_cent30to50 = new TH1D("h_mujeteta_ghost_cent30to50","ghost reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
                TH1D *h_mujeteta_ghost_cent50to90 = new TH1D("h_mujeteta_ghost_cent50to90","ghost reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
TH1D *h_genjeteta = new TH1D("h_genjeteta","gen jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_trkEta = new TH1D("h_trkEta","track #eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_trkEta_cent0to10 = new TH1D("h_trkEta_cent0to10","track #eta, 0-10 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_trkEta_cent10to30 = new TH1D("h_trkEta_cent10to30","track #eta, 10-30 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_trkEta_cent30to50 = new TH1D("h_trkEta_cent30to50","track #eta, 30-50 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_trkEta_cent50to90 = new TH1D("h_trkEta_cent50to90","track #eta, 50-90 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
const int NPtBins = 450;
const double ptMin = 50.0;
const double ptMax = 500.0;
TH1D *h_jetpt_raw = new TH1D("h_jetpt_raw","raw jet pt",NPtBins,ptMin,ptMax);
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_cent0to10 = new TH1D("h_jetpt_cent0to10","reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_cent10to30 = new TH1D("h_jetpt_cent10to30","reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_cent30to50 = new TH1D("h_jetpt_cent30to50","reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_cent50to90 = new TH1D("h_jetpt_cent50to90","reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_g = new TH1D("h_jetpt_g","g reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_g_cent0to10 = new TH1D("h_jetpt_g_cent0to10","g reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_g_cent10to30 = new TH1D("h_jetpt_g_cent10to30","g reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_g_cent30to50 = new TH1D("h_jetpt_g_cent30to50","g reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_g_cent50to90 = new TH1D("h_jetpt_g_cent50to90","g reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_uds = new TH1D("h_jetpt_uds","uds reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_uds_cent0to10 = new TH1D("h_jetpt_uds_cent0to10","uds reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_uds_cent10to30 = new TH1D("h_jetpt_uds_cent10to30","uds reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_uds_cent30to50 = new TH1D("h_jetpt_uds_cent30to50","uds reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_uds_cent50to90 = new TH1D("h_jetpt_uds_cent50to90","uds reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_b = new TH1D("h_jetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_b_cent0to10 = new TH1D("h_jetpt_b_cent0to10","b reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_b_cent10to30 = new TH1D("h_jetpt_b_cent10to30","b reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_b_cent30to50 = new TH1D("h_jetpt_b_cent30to50","b reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_b_cent50to90 = new TH1D("h_jetpt_b_cent50to90","b reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_c = new TH1D("h_jetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_c_cent0to10 = new TH1D("h_jetpt_c_cent0to10","c reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_c_cent10to30 = new TH1D("h_jetpt_c_cent10to30","c reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_c_cent30to50 = new TH1D("h_jetpt_c_cent30to50","c reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_c_cent50to90 = new TH1D("h_jetpt_c_cent50to90","c reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_ghost = new TH1D("h_jetpt_ghost","ghost reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_ghost_cent0to10 = new TH1D("h_jetpt_ghost_cent0to10","ghost reco jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_ghost_cent10to30 = new TH1D("h_jetpt_ghost_cent10to30","ghost reco jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_ghost_cent30to50 = new TH1D("h_jetpt_ghost_cent30to50","ghost reco jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_jetpt_ghost_cent50to90 = new TH1D("h_jetpt_ghost_cent50to90","ghost reco jet pt, 50-90 %",NPtBins,ptMin,ptMax);
TH1D *h_mujetpt = new TH1D("h_mujetpt","reco muon jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_cent0to10 = new TH1D("h_mujetpt_cent0to10","reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_cent10to30 = new TH1D("h_mujetpt_cent10to30","reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_cent30to50 = new TH1D("h_mujetpt_cent30to50","reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_cent50to90 = new TH1D("h_mujetpt_cent50to90","reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_g = new TH1D("h_mujetpt_g","g reco muon jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_g_cent0to10 = new TH1D("h_mujetpt_g_cent0to10","g reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_g_cent10to30 = new TH1D("h_mujetpt_g_cent10to30","g reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_g_cent30to50 = new TH1D("h_mujetpt_g_cent30to50","g reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_g_cent50to90 = new TH1D("h_mujetpt_g_cent50to90","g reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_uds = new TH1D("h_mujetpt_uds","uds reco muon jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_uds_cent0to10 = new TH1D("h_mujetpt_uds_cent0to10","uds reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_uds_cent10to30 = new TH1D("h_mujetpt_uds_cent10to30","uds reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_uds_cent30to50 = new TH1D("h_mujetpt_uds_cent30to50","uds reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_uds_cent50to90 = new TH1D("h_mujetpt_uds_cent50to90","uds reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_b = new TH1D("h_mujetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_b_cent0to10 = new TH1D("h_mujetpt_b_cent0to10","b reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_b_cent10to30 = new TH1D("h_mujetpt_b_cent10to30","b reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_b_cent30to50 = new TH1D("h_mujetpt_b_cent30to50","b reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_b_cent50to90 = new TH1D("h_mujetpt_b_cent50to90","b reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_c = new TH1D("h_mujetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_c_cent0to10 = new TH1D("h_mujetpt_c_cent0to10","c reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_c_cent10to30 = new TH1D("h_mujetpt_c_cent10to30","c reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_c_cent30to50 = new TH1D("h_mujetpt_c_cent30to50","c reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_c_cent50to90 = new TH1D("h_mujetpt_c_cent50to90","c reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_ghost = new TH1D("h_mujetpt_ghost","ghost reco muon jet pt",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_ghost_cent0to10 = new TH1D("h_mujetpt_ghost_cent0to10","ghost reco muon jet pt, 0-10 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_ghost_cent10to30 = new TH1D("h_mujetpt_ghost_cent10to30","ghost reco muon jet pt, 10-30 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_ghost_cent30to50 = new TH1D("h_mujetpt_ghost_cent30to50","ghost reco muon jet pt, 30-50 %",NPtBins,ptMin,ptMax);
                TH1D *h_mujetpt_ghost_cent50to90 = new TH1D("h_mujetpt_ghost_cent50to90","ghost reco muon jet pt, 50-90 %",NPtBins,ptMin,ptMax);
TH1D *h_genjetpt = new TH1D("h_genjetpt","gen jet pt",NPtBins,ptMin,ptMax);
const int NTrkPtBins = 500;
const double trkPtMin = 0.0;
const double trkPtMax = 500.0;
TH1D *h_trkPt = new TH1D("h_trkPt","track p_{T}",NTrkPtBins,trkPtMin,trkPtMax);
                TH1D *h_trkPt_cent0to10 = new TH1D("h_trkPt_cent0to10","track pt, 0-10 %",NTrkPtBins,trkPtMin,trkPtMax);
                TH1D *h_trkPt_cent10to30 = new TH1D("h_trkPt_cent10to30","track pt, 10-30 %",NTrkPtBins,trkPtMin,trkPtMax);
                TH1D *h_trkPt_cent30to50 = new TH1D("h_trkPt_cent30to50","track pt, 30-50 %",NTrkPtBins,trkPtMin,trkPtMax);
                TH1D *h_trkPt_cent50to90 = new TH1D("h_trkPt_cent50to90","track pt, 50-90 %",NTrkPtBins,trkPtMin,trkPtMax);
// muon data
const int NMuPtBins = 500;
const double muPtMin = 0.0;
const double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muPt_cent0to10 = new TH1D("h_muPt_cent0to10","muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muPt_cent10to30 = new TH1D("h_muPt_cent10to30","muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muPt_cent30to50 = new TH1D("h_muPt_cent30to50","muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muPt_cent50to90 = new TH1D("h_muPt_cent50to90","muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muPt_noCut = new TH1D("h_muPt_noCut","muon pt - no cut",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muInJetPt = new TH1D("h_muInJetPt","in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_cent0to10 = new TH1D("h_muInJetPt_cent0to10","in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_cent10to30 = new TH1D("h_muInJetPt_cent10to30","in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_cent30to50 = new TH1D("h_muInJetPt_cent30to50","in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_cent50to90 = new TH1D("h_muInJetPt_cent50to90","in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muInJetPt_g = new TH1D("h_muInJetPt_g","g in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_g_cent0to10 = new TH1D("h_muInJetPt_g_cent0to10","g in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_g_cent10to30 = new TH1D("h_muInJetPt_g_cent10to30","g in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_g_cent30to50 = new TH1D("h_muInJetPt_g_cent30to50","g in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_g_cent50to90 = new TH1D("h_muInJetPt_g_cent50to90","g in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muInJetPt_uds = new TH1D("h_muInJetPt_uds","uds in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_uds_cent0to10 = new TH1D("h_muInJetPt_uds_cent0to10","uds in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_uds_cent10to30 = new TH1D("h_muInJetPt_uds_cent10to30","uds in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_uds_cent30to50 = new TH1D("h_muInJetPt_uds_cent30to50","uds in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_uds_cent50to90 = new TH1D("h_muInJetPt_uds_cent50to90","uds in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muInJetPt_b = new TH1D("h_muInJetPt_b","b in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_b_cent0to10 = new TH1D("h_muInJetPt_b_cent0to10","b in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_b_cent10to30 = new TH1D("h_muInJetPt_b_cent10to30","b in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_b_cent30to50 = new TH1D("h_muInJetPt_b_cent30to50","b in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_b_cent50to90 = new TH1D("h_muInJetPt_b_cent50to90","b in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muInJetPt_c = new TH1D("h_muInJetPt_c","c in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_c_cent0to10 = new TH1D("h_muInJetPt_c_cent0to10","c in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_c_cent10to30 = new TH1D("h_muInJetPt_c_cent10to30","c in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_c_cent30to50 = new TH1D("h_muInJetPt_c_cent30to50","c in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_c_cent50to90 = new TH1D("h_muInJetPt_c_cent50to90","c in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muInJetPt_ghost = new TH1D("h_muInJetPt_ghost","ghost in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_ghost_cent0to10 = new TH1D("h_muInJetPt_ghost_cent0to10","ghost in-jet muon pt, 0-10 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_ghost_cent10to30 = new TH1D("h_muInJetPt_ghost_cent10to30","ghost in-jet muon pt, 10-30 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_ghost_cent30to50 = new TH1D("h_muInJetPt_ghost_cent30to50","ghost in-jet muon pt, 30-50 %",NMuPtBins,muPtMin,muPtMax);
                TH1D *h_muInJetPt_ghost_cent50to90 = new TH1D("h_muInJetPt_ghost_cent50to90","ghost in-jet muon pt, 50-90 %",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_muEta_cent0to10 = new TH1D("h_muEta_cent0to10","muon eta, 0-10 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_muEta_cent10to30 = new TH1D("h_muEta_cent10to30","muon eta, 10-30 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_muEta_cent30to50 = new TH1D("h_muEta_cent30to50","muon eta, 30-50 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
                TH1D *h_muEta_cent50to90 = new TH1D("h_muEta_cent50to90","muon eta, 50-90 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
                TH1D *h_muPhi_cent0to10 = new TH1D("h_muPhi_cent0to10","muon #phi, 0-10 %",NPhiBins,phiMin,phiMax);
                TH1D *h_muPhi_cent10to30 = new TH1D("h_muPhi_cent10to30","muon #phi, 10-30 %",NPhiBins,phiMin,phiMax);
                TH1D *h_muPhi_cent30to50 = new TH1D("h_muPhi_cent30to50","muon #phi, 30-50 %",NPhiBins,phiMin,phiMax);
                TH1D *h_muPhi_cent50to90 = new TH1D("h_muPhi_cent50to90","muon #phi, 50-90 %",NPhiBins,phiMin,phiMax);
const int NRelPtBins = 50;
const double relPtMin = 0.0;
const double relPtMax = 5.0;
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_cent0to10 = new TH1D("h_muRelPt_cent0to10","muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_cent10to30 = new TH1D("h_muRelPt_cent10to30","muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_cent30to50 = new TH1D("h_muRelPt_cent30to50","muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_cent50to90 = new TH1D("h_muRelPt_cent50to90","muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
	TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","g muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_g_cent0to10 = new TH1D("h_muRelPt_g_cent0to10","g muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_g_cent10to30 = new TH1D("h_muRelPt_g_cent10to30","g muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_g_cent30to50 = new TH1D("h_muRelPt_g_cent30to50","g muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_g_cent50to90 = new TH1D("h_muRelPt_g_cent50to90","g muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
	TH1D *h_muRelPt_uds = new TH1D("h_muRelPt_uds","uds muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_uds_cent0to10 = new TH1D("h_muRelPt_uds_cent0to10","uds muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_uds_cent10to30 = new TH1D("h_muRelPt_uds_cent10to30","uds muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_uds_cent30to50 = new TH1D("h_muRelPt_uds_cent30to50","uds muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_uds_cent50to90 = new TH1D("h_muRelPt_uds_cent50to90","uds muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
	TH1D *h_muRelPt_b = new TH1D("h_muRelPt_b","b muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_b_cent0to10 = new TH1D("h_muRelPt_b_cent0to10","b muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_b_cent10to30 = new TH1D("h_muRelPt_b_cent10to30","b muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_b_cent30to50 = new TH1D("h_muRelPt_b_cent30to50","b muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_b_cent50to90 = new TH1D("h_muRelPt_b_cent50to90","b muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
	TH1D *h_muRelPt_c = new TH1D("h_muRelPt_c","c muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_c_cent0to10 = new TH1D("h_muRelPt_c_cent0to10","g muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_c_cent10to30 = new TH1D("h_muRelPt_c_cent10to30","g muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_c_cent30to50 = new TH1D("h_muRelPt_c_cent30to50","g muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_c_cent50to90 = new TH1D("h_muRelPt_c_cent50to90","g muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
	TH1D *h_muRelPt_ghost = new TH1D("h_muRelPt_ghost","ghost muon rel pt",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_ghost_cent0to10 = new TH1D("h_muRelPt_ghost_cent0to10","ghost muon rel pt, 0-10 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_ghost_cent10to30 = new TH1D("h_muRelPt_ghost_cent10to30","ghost muon rel pt, 10-30 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_ghost_cent30to50 = new TH1D("h_muRelPt_ghost_cent30to50","ghost muon rel pt, 30-50 %",NRelPtBins,relPtMin,relPtMax);
                TH1D *h_muRelPt_ghost_cent50to90 = new TH1D("h_muRelPt_ghost_cent50to90","ghost muon rel pt, 50-90 %",NRelPtBins,relPtMin,relPtMax);
TH1D *h_partonFlavor = new TH1D("h_partonFlavor","parton flavor ID",26,-5.0,21.0);
TH1D *h_hadronFlavor = new TH1D("h_hadronFlavor","hadron flavor ID",26,-5.0,21.0);
        
// event info
TH1D *h_pthat = new TH1D("h_pthat","pthat",NPtBins,ptMin,ptMax);
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
TH1D *h_hiBin = new TH1D("h_hiBin","hiBin (2*centrality)",200,0,200);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",100,0,1);


// Sumw2 commands
h_jetphi->Sumw2();
        h_jetphi_cent0to10->Sumw2();
        h_jetphi_cent10to30->Sumw2();
        h_jetphi_cent30to50->Sumw2();
        h_jetphi_cent50to90->Sumw2();
	h_jetphi_g->Sumw2();
        h_jetphi_g_cent0to10->Sumw2();
        h_jetphi_g_cent10to30->Sumw2();
        h_jetphi_g_cent30to50->Sumw2();
        h_jetphi_g_cent50to90->Sumw2();
	h_jetphi_uds->Sumw2();
        h_jetphi_uds_cent0to10->Sumw2();
        h_jetphi_uds_cent10to30->Sumw2();
        h_jetphi_uds_cent30to50->Sumw2();
        h_jetphi_uds_cent50to90->Sumw2();
	h_jetphi_b->Sumw2();
        h_jetphi_b_cent0to10->Sumw2();
        h_jetphi_b_cent10to30->Sumw2();
        h_jetphi_b_cent30to50->Sumw2();
        h_jetphi_b_cent50to90->Sumw2();
	h_jetphi_c->Sumw2();
        h_jetphi_c_cent0to10->Sumw2();
        h_jetphi_c_cent10to30->Sumw2();
        h_jetphi_c_cent30to50->Sumw2();
        h_jetphi_c_cent50to90->Sumw2();
	h_jetphi_ghost->Sumw2();
        h_jetphi_ghost_cent0to10->Sumw2();
        h_jetphi_ghost_cent10to30->Sumw2();
        h_jetphi_ghost_cent30to50->Sumw2();
        h_jetphi_ghost_cent50to90->Sumw2();
h_mujetphi->Sumw2();
        h_mujetphi_cent0to10->Sumw2();
        h_mujetphi_cent10to30->Sumw2();
        h_mujetphi_cent30to50->Sumw2();
        h_mujetphi_cent50to90->Sumw2();
	h_mujetphi_g->Sumw2();
        h_mujetphi_g_cent0to10->Sumw2();
        h_mujetphi_g_cent10to30->Sumw2();
        h_mujetphi_g_cent30to50->Sumw2();
        h_mujetphi_g_cent50to90->Sumw2();
	h_mujetphi_uds->Sumw2();
        h_mujetphi_uds_cent0to10->Sumw2();
        h_mujetphi_uds_cent10to30->Sumw2();
        h_mujetphi_uds_cent30to50->Sumw2();
        h_mujetphi_uds_cent50to90->Sumw2();
	h_mujetphi_b->Sumw2();
        h_mujetphi_b_cent0to10->Sumw2();
        h_mujetphi_b_cent10to30->Sumw2();
        h_mujetphi_b_cent30to50->Sumw2();
        h_mujetphi_b_cent50to90->Sumw2();
	h_mujetphi_c->Sumw2();
        h_mujetphi_c_cent0to10->Sumw2();
        h_mujetphi_c_cent10to30->Sumw2();
        h_mujetphi_c_cent30to50->Sumw2();
        h_mujetphi_c_cent50to90->Sumw2();
	h_mujetphi_ghost->Sumw2();
        h_mujetphi_ghost_cent0to10->Sumw2();
        h_mujetphi_ghost_cent10to30->Sumw2();
        h_mujetphi_ghost_cent30to50->Sumw2();
        h_mujetphi_ghost_cent50to90->Sumw2();
h_genjetphi->Sumw2();
h_trkPhi->Sumw2();
        h_trkPhi_cent0to10->Sumw2();
        h_trkPhi_cent10to30->Sumw2();
        h_trkPhi_cent30to50->Sumw2();
        h_trkPhi_cent50to90->Sumw2();
h_jeteta->Sumw2();
        h_jeteta_cent0to10->Sumw2();
        h_jeteta_cent10to30->Sumw2();
        h_jeteta_cent30to50->Sumw2();
        h_jeteta_cent50to90->Sumw2();
	h_jeteta_g->Sumw2();
        h_jeteta_g_cent0to10->Sumw2();
        h_jeteta_g_cent10to30->Sumw2();
        h_jeteta_g_cent30to50->Sumw2();
        h_jeteta_g_cent50to90->Sumw2();
	h_jeteta_uds->Sumw2();
        h_jeteta_uds_cent0to10->Sumw2();
        h_jeteta_uds_cent10to30->Sumw2();
        h_jeteta_uds_cent30to50->Sumw2();
        h_jeteta_uds_cent50to90->Sumw2();
	h_jeteta_b->Sumw2();
        h_jeteta_b_cent0to10->Sumw2();
        h_jeteta_b_cent10to30->Sumw2();
        h_jeteta_b_cent30to50->Sumw2();
        h_jeteta_b_cent50to90->Sumw2();
	h_jeteta_c->Sumw2();
        h_jeteta_c_cent0to10->Sumw2();
        h_jeteta_c_cent10to30->Sumw2();
        h_jeteta_c_cent30to50->Sumw2();
        h_jeteta_c_cent50to90->Sumw2();
	h_jeteta_ghost->Sumw2();
        h_jeteta_ghost_cent0to10->Sumw2();
        h_jeteta_ghost_cent10to30->Sumw2();
        h_jeteta_ghost_cent30to50->Sumw2();
        h_jeteta_ghost_cent50to90->Sumw2();
h_mujeteta->Sumw2();
        h_mujeteta_cent0to10->Sumw2();
        h_mujeteta_cent10to30->Sumw2();
        h_mujeteta_cent30to50->Sumw2();
        h_mujeteta_cent50to90->Sumw2();
	h_mujeteta_g->Sumw2();
        h_mujeteta_g_cent0to10->Sumw2();
        h_mujeteta_g_cent10to30->Sumw2();
        h_mujeteta_g_cent30to50->Sumw2();
        h_mujeteta_g_cent50to90->Sumw2();
	h_mujeteta_uds->Sumw2();
        h_mujeteta_uds_cent0to10->Sumw2();
        h_mujeteta_uds_cent10to30->Sumw2();
        h_mujeteta_uds_cent30to50->Sumw2();
        h_mujeteta_uds_cent50to90->Sumw2();
	h_mujeteta_b->Sumw2();
        h_mujeteta_b_cent0to10->Sumw2();
        h_mujeteta_b_cent10to30->Sumw2();
        h_mujeteta_b_cent30to50->Sumw2();
        h_mujeteta_b_cent50to90->Sumw2();
	h_mujeteta_c->Sumw2();
        h_mujeteta_c_cent0to10->Sumw2();
        h_mujeteta_c_cent10to30->Sumw2();
        h_mujeteta_c_cent30to50->Sumw2();
        h_mujeteta_c_cent50to90->Sumw2();
	h_mujeteta_ghost->Sumw2();
        h_mujeteta_ghost_cent0to10->Sumw2();
        h_mujeteta_ghost_cent10to30->Sumw2();
        h_mujeteta_ghost_cent30to50->Sumw2();
        h_mujeteta_ghost_cent50to90->Sumw2();
h_genjeteta->Sumw2();
h_trkEta->Sumw2();
        h_trkEta_cent0to10->Sumw2();
        h_trkEta_cent10to30->Sumw2();
        h_trkEta_cent30to50->Sumw2();
        h_trkEta_cent50to90->Sumw2();
h_jetpt_raw->Sumw2();
h_jetpt->Sumw2();
        h_jetpt_cent0to10->Sumw2();
        h_jetpt_cent10to30->Sumw2();
        h_jetpt_cent30to50->Sumw2();
        h_jetpt_cent50to90->Sumw2();
	h_jetpt_g->Sumw2();
        h_jetpt_g_cent0to10->Sumw2();
        h_jetpt_g_cent10to30->Sumw2();
        h_jetpt_g_cent30to50->Sumw2();
        h_jetpt_g_cent50to90->Sumw2();
	h_jetpt_uds->Sumw2();
        h_jetpt_uds_cent0to10->Sumw2();
        h_jetpt_uds_cent10to30->Sumw2();
        h_jetpt_uds_cent30to50->Sumw2();
        h_jetpt_uds_cent50to90->Sumw2();
	h_jetpt_b->Sumw2();
        h_jetpt_b_cent0to10->Sumw2();
        h_jetpt_b_cent10to30->Sumw2();
        h_jetpt_b_cent30to50->Sumw2();
        h_jetpt_b_cent50to90->Sumw2();
	h_jetpt_c->Sumw2();
        h_jetpt_c_cent0to10->Sumw2();
        h_jetpt_c_cent10to30->Sumw2();
        h_jetpt_c_cent30to50->Sumw2();
        h_jetpt_c_cent50to90->Sumw2();
	h_jetpt_ghost->Sumw2();
        h_jetpt_ghost_cent0to10->Sumw2();
        h_jetpt_ghost_cent10to30->Sumw2();
        h_jetpt_ghost_cent30to50->Sumw2();
        h_jetpt_ghost_cent50to90->Sumw2();
h_mujetpt->Sumw2();
        h_mujetpt_cent0to10->Sumw2();
        h_mujetpt_cent10to30->Sumw2();
        h_mujetpt_cent30to50->Sumw2();
        h_mujetpt_cent50to90->Sumw2();
	h_mujetpt_g->Sumw2();
        h_mujetpt_g_cent0to10->Sumw2();
        h_mujetpt_g_cent10to30->Sumw2();
        h_mujetpt_g_cent30to50->Sumw2();
        h_mujetpt_g_cent50to90->Sumw2();
	h_mujetpt_uds->Sumw2();
        h_mujetpt_uds_cent0to10->Sumw2();
        h_mujetpt_uds_cent10to30->Sumw2();
        h_mujetpt_uds_cent30to50->Sumw2();
        h_mujetpt_uds_cent50to90->Sumw2();
	h_mujetpt_b->Sumw2();
        h_mujetpt_b_cent0to10->Sumw2();
        h_mujetpt_b_cent10to30->Sumw2();
        h_mujetpt_b_cent30to50->Sumw2();
        h_mujetpt_b_cent50to90->Sumw2();
	h_mujetpt_c->Sumw2();
        h_mujetpt_c_cent0to10->Sumw2();
        h_mujetpt_c_cent10to30->Sumw2();
        h_mujetpt_c_cent30to50->Sumw2();
        h_mujetpt_c_cent50to90->Sumw2();
	h_mujetpt_ghost->Sumw2();
        h_mujetpt_ghost_cent0to10->Sumw2();
        h_mujetpt_ghost_cent10to30->Sumw2();
        h_mujetpt_ghost_cent30to50->Sumw2();
        h_mujetpt_ghost_cent50to90->Sumw2();
h_genjetpt->Sumw2();
h_trkPt->Sumw2();
        h_trkPt_cent0to10->Sumw2();
        h_trkPt_cent10to30->Sumw2();
        h_trkPt_cent30to50->Sumw2();
        h_trkPt_cent50to90->Sumw2();
h_muPt->Sumw2();
        h_muPt_cent0to10->Sumw2();
        h_muPt_cent10to30->Sumw2();
        h_muPt_cent30to50->Sumw2();
        h_muPt_cent50to90->Sumw2();
h_muPt_noCut->Sumw2();
h_muInJetPt->Sumw2();
        h_muInJetPt_cent0to10->Sumw2();
        h_muInJetPt_cent10to30->Sumw2();
        h_muInJetPt_cent30to50->Sumw2();
        h_muInJetPt_cent50to90->Sumw2();
    h_muInJetPt_g->Sumw2();
        h_muInJetPt_g_cent0to10->Sumw2();
        h_muInJetPt_g_cent10to30->Sumw2();
        h_muInJetPt_g_cent30to50->Sumw2();
        h_muInJetPt_g_cent50to90->Sumw2();
        h_muInJetPt_uds->Sumw2();
        h_muInJetPt_uds_cent0to10->Sumw2();
        h_muInJetPt_uds_cent10to30->Sumw2();
        h_muInJetPt_uds_cent30to50->Sumw2();
        h_muInJetPt_uds_cent50to90->Sumw2();
	h_muInJetPt_b->Sumw2();
        h_muInJetPt_b_cent0to10->Sumw2();
        h_muInJetPt_b_cent10to30->Sumw2();
        h_muInJetPt_b_cent30to50->Sumw2();
        h_muInJetPt_b_cent50to90->Sumw2();
	h_muInJetPt_c->Sumw2();
        h_muInJetPt_c_cent0to10->Sumw2();
        h_muInJetPt_c_cent10to30->Sumw2();
        h_muInJetPt_c_cent30to50->Sumw2();
        h_muInJetPt_c_cent50to90->Sumw2();
	h_muInJetPt_ghost->Sumw2();
        h_muInJetPt_ghost_cent0to10->Sumw2();
        h_muInJetPt_ghost_cent10to30->Sumw2();
        h_muInJetPt_ghost_cent30to50->Sumw2();
        h_muInJetPt_ghost_cent50to90->Sumw2();
h_muEta->Sumw2();
        h_muEta_cent0to10->Sumw2();
        h_muEta_cent10to30->Sumw2();
        h_muEta_cent30to50->Sumw2();
        h_muEta_cent50to90->Sumw2();
h_muPhi->Sumw2();
        h_muPhi_cent0to10->Sumw2();
        h_muPhi_cent10to30->Sumw2();
        h_muPhi_cent30to50->Sumw2();
        h_muPhi_cent50to90->Sumw2();
h_muRelPt->Sumw2();
        h_muRelPt_cent0to10->Sumw2();
        h_muRelPt_cent10to30->Sumw2();
        h_muRelPt_cent30to50->Sumw2();
        h_muRelPt_cent50to90->Sumw2();
	h_muRelPt_g->Sumw2();
        h_muRelPt_g_cent0to10->Sumw2();
        h_muRelPt_g_cent10to30->Sumw2();
        h_muRelPt_g_cent30to50->Sumw2();
        h_muRelPt_g_cent50to90->Sumw2();
	h_muRelPt_uds->Sumw2();
        h_muRelPt_uds_cent0to10->Sumw2();
        h_muRelPt_uds_cent10to30->Sumw2();
        h_muRelPt_uds_cent30to50->Sumw2();
        h_muRelPt_uds_cent50to90->Sumw2();
	h_muRelPt_b->Sumw2();
        h_muRelPt_b_cent0to10->Sumw2();
        h_muRelPt_b_cent10to30->Sumw2();
        h_muRelPt_b_cent30to50->Sumw2();
        h_muRelPt_b_cent50to90->Sumw2();
	h_muRelPt_c->Sumw2();
        h_muRelPt_c_cent0to10->Sumw2();
        h_muRelPt_c_cent10to30->Sumw2();
        h_muRelPt_c_cent30to50->Sumw2();
        h_muRelPt_c_cent50to90->Sumw2();
	h_muRelPt_ghost->Sumw2();
        h_muRelPt_ghost_cent0to10->Sumw2();
        h_muRelPt_ghost_cent10to30->Sumw2();
        h_muRelPt_ghost_cent30to50->Sumw2();
        h_muRelPt_ghost_cent50to90->Sumw2();
h_pthat->Sumw2();
h_vz->Sumw2();
h_hiBin->Sumw2();
h_deltaR->Sumw2();
h_partonFlavor->Sumw2();
h_hadronFlavor->Sumw2();

// count number of files to be skimmed over

/*
int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("root://cmsxrootd.fnal.gov//store/group/phys_heavyions/jviinika/PbPb2018MC_Dijet_pThat15_CP5_HydjetDrumMB_5p02TeV_HINPbPbAutumn18-mva98_103X_first2k_addedEventPlane/0000");

   	if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    	{ 
        	printf("Could not open current directory" ); 
        	return; 
    	}
		
    	if(dr){ 
		 
		while((de=readdir(dr)) != NULL) {
		if(strcmp(de->d_name,".")!=0 && strcmp(de->d_name,"..")!=0){
		  	NFiles++;
		}
    	    	}
 	closedir(dr);
	}

	cout <<"number of files ="<< NFiles << endl;
*/


TF1 *vz_fit_fxn = new TF1("vz_fit_fxn","[0] + [1]*x + [2]*x*x + [3]*x*x*x",-15.0,15.0);
vz_fit_fxn->SetParameters(1.02937e+00,-5.25572e-03,-8.21959e-04,2.60358e-05);

TF1 *pt_fit_fxn = new TF1("pt_fit_fxn","[0]*exp([1]*x) + [2]*exp([3]*x)",50.0,500.0);
pt_fit_fxn->SetParameters(2.95395e+00,-2.60850e-02,5.71522e-01,-2.90376e-03);


        TFile *f = TFile::Open(input);

	auto em = new eventMap(f);
        em->isMC = 1;
        em->init();
        em->loadTrack();
        em->loadJet("akCs4PFJetAnalyzer");
        em->loadMuon("ggHiNtuplizerGED");
		Long64_t NEvents = em->evtTree->GetEntries();
		
	// event loop
	int evi_frac = 0;
	for(int evi = 0; evi < NEvents; evi++){
	//for(int evi = 1; evi < 5; evi++){	
		em->getEvent(evi);
		//event cuts
                if(em->vz>15.0){continue;}
                if(em->hiBin>180){continue;}
                h_vz->Fill(em->vz);
                h_hiBin->Fill(em->hiBin);
		// vz reweight
		//double w_vz = 1.0/(vz_fit_fxn->Eval(em->vz)); // calculated vz weight
                double w_vz = 1;
		
                double hiBin = em->hiBin;
                double vz = em->vz;
		
		if(vz>15.0){continue;}
		if(em->pthat<30.0){continue;}
		

		h_vz->Fill(em->vz,w_vz);
		h_hiBin->Fill(em->hiBin);
		h_pthat->Fill(em->pthat);
		//cout << "vz = " << em->vz << endl;
		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
                evi_frac = 100*evi/NEvents;
		cout << "Number of jets = " << em->njet << endl;
        
                const double trkPtCut = 1.0;
		// track loop
                
		for(int trki=0; trki < em->ntrk; trki++){
                        

			double trkPti = em->trkpt[trki];
			double trkEtai = em->trketa[trki];
			double trkPhii = em->trkphi[trki];

			

			if(trkPti < trkPtCut || TMath::Abs(trkEtai)>trkEtaMax) {continue;}
			h_trkPt->Fill(trkPti);
			h_trkEta->Fill(trkEtai);
			h_trkPhi->Fill(trkPhii);

                        if(em->hiBin>0 && em->hiBin<20){
                                h_trkPt_cent0to10->Fill(trkPti);
                                h_trkEta_cent0to10->Fill(trkEtai);
                                h_trkPhi_cent0to10->Fill(trkPhii);
                        }
                        if(em->hiBin>20 && em->hiBin<60){
                                h_trkPt_cent10to30->Fill(trkPti);
                                h_trkEta_cent10to30->Fill(trkEtai);
                                h_trkPhi_cent10to30->Fill(trkPhii);
                        }
                        if(em->hiBin>60 && em->hiBin<100){
                                h_trkPt_cent30to50->Fill(trkPti);
                                h_trkEta_cent30to50->Fill(trkEtai);
                                h_trkPhi_cent30to50->Fill(trkPhii);
                        }
                        if(em->hiBin>100 && em->hiBin<180){
                                h_trkPt_cent50to90->Fill(trkPti);
                                h_trkEta_cent50to90->Fill(trkEtai);
                                h_trkPhi_cent50to90->Fill(trkPhii);
                        }

		}
                
		// seperate muon loop
		
		double muPtCut = 10.0;

		for(int muj=0; muj < em->nMu; muj++){

			double muPtj = em->muPt->at(muj);
			double muEtaj = em->muEta->at(muj);
			double muPhij = em->muPhi->at(muj);

			// top level cuts
			if(em->muIsTracker->at(muj)==0 || TMath::Abs(muEtaj)>2.4 || em->muChi2NDF->at(muj)==-99 || em->muChi2NDF->at(muj)>10
					||TMath::Abs(em->muInnerD0->at(muj))>0.2 || TMath::Abs(em->muInnerDz->at(muj))>0.5 || em->muMuonHits->at(muj)<= 0
					|| em->muStations->at(muj)<= 1 || em->muTrkLayers->at(muj)<=5 || em->muPixelHits->at(muj)<=0 || em->muPt->at(muj) < muPtCut){continue;}

			h_muPt->Fill(muPtj);
                        h_muEta->Fill(muEtaj);
                        h_muPhi->Fill(muPhij);
                                        
                        if(em->hiBin>0 && em->hiBin<20){
                                h_muPt_cent0to10->Fill(muPtj);
                                h_muEta_cent0to10->Fill(muEtaj);
                                h_muPhi_cent0to10->Fill(muPhij);
                        }
                        if(em->hiBin>20 && em->hiBin<60){
                                h_muPt_cent10to30->Fill(muPtj);
                                h_muEta_cent10to30->Fill(muEtaj);
                                h_muPhi_cent10to30->Fill(muPhij);
                        }
                        if(em->hiBin>60 && em->hiBin<100){
                                h_muPt_cent30to50->Fill(muPtj);
                                h_muEta_cent30to50->Fill(muEtaj);
                                h_muPhi_cent30to50->Fill(muPhij);
                        }
                        if(em->hiBin>100 && em->hiBin<180){
                                h_muPt_cent50to90->Fill(muPtj);
                                h_muEta_cent50to90->Fill(muEtaj);
                                h_muPhi_cent50to90->Fill(muPhij);
                        }
			
		}
		//reco jet loop
		for(int jetj = 0; jetj < em->njet; jetj++){

			double jetPtj = em->jetpt[jetj];
			double jetEtaj = em->jeteta[jetj];
			double jetPhij = em->jetphi[jetj];

                        if(em->pthat<0.35*jetPtj){continue;}

			double w_pt = 1.0/(pt_fit_fxn->Eval(em->jetpt[jetj]));
			// double w_tot = w_vz*w_pt*em->weight;	
			double w_tot = em->weight;

			if(jetPtj<ptMin || jetPtj>ptMax || fabs(jetEtaj)>etaMax || jetPhij==-999){continue;}
			
			
			h_jetpt_raw->Fill(jetPtj,w_vz*em->weight);						
			h_jetpt->Fill(jetPtj,w_tot);
			h_jeteta->Fill(jetEtaj,w_tot);
			h_jetphi->Fill(jetPhij,w_tot);

                        if(em->hiBin>0 && em->hiBin<20){
                                h_jetpt_cent0to10->Fill(jetPtj,w_tot);
                                h_jeteta_cent0to10->Fill(jetEtaj,w_tot);
                                h_jetphi_cent0to10->Fill(jetPhij,w_tot);
                        }
                        if(em->hiBin>20 && em->hiBin<60){
                                h_jetpt_cent10to30->Fill(jetPtj,w_tot);
                                h_jeteta_cent10to30->Fill(jetEtaj,w_tot);
                                h_jetphi_cent10to30->Fill(jetPhij,w_tot);
                        }
                        if(em->hiBin>60 && em->hiBin<100){
                                h_jetpt_cent30to50->Fill(jetPtj,w_tot);
                                h_jeteta_cent30to50->Fill(jetEtaj,w_tot);
                                h_jetphi_cent30to50->Fill(jetPhij,w_tot);
                        }
                        if(em->hiBin>100 && em->hiBin<180){
                                h_jetpt_cent50to90->Fill(jetPtj,w_tot);
                                h_jeteta_cent50to90->Fill(jetEtaj,w_tot);
                                h_jetphi_cent50to90->Fill(jetPhij,w_tot);
                        }
			
			h_partonFlavor->Fill(em->matchedPartonFlavor[jetj]);
                        h_hadronFlavor->Fill(em->matchedHadronFlavor[jetj]);
			if(em->matchedHadronFlavor[jetj]==4){
                                h_jetpt_c->Fill(jetPtj,w_tot);
                                h_jeteta_c->Fill(jetEtaj,w_tot);
                                h_jetphi_c->Fill(jetPhij,w_tot);
                                if(em->hiBin>0 && em->hiBin<20){
                                        h_jetpt_c_cent0to10->Fill(jetPtj,w_tot);
                                        h_jeteta_c_cent0to10->Fill(jetEtaj,w_tot);
                                        h_jetphi_c_cent0to10->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>20 && em->hiBin<60){
                                        h_jetpt_c_cent10to30->Fill(jetPtj,w_tot);
                                        h_jeteta_c_cent10to30->Fill(jetEtaj,w_tot);
                                        h_jetphi_c_cent10to30->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>60 && em->hiBin<100){
                                        h_jetpt_c_cent30to50->Fill(jetPtj,w_tot);
                                        h_jeteta_c_cent30to50->Fill(jetEtaj,w_tot);
                                        h_jetphi_c_cent30to50->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>100 && em->hiBin<180){
                                        h_jetpt_c_cent50to90->Fill(jetPtj,w_tot);
                                        h_jeteta_c_cent50to90->Fill(jetEtaj,w_tot);
                                        h_jetphi_c_cent50to90->Fill(jetPhij,w_tot);
                                }
			}
			else if(em->matchedHadronFlavor[jetj]==5){
                                h_jetpt_b->Fill(jetPtj,w_tot);
                                h_jeteta_b->Fill(jetEtaj,w_tot);
                                h_jetphi_b->Fill(jetPhij,w_tot);
                                if(em->hiBin>0 && em->hiBin<20){
                                        h_jetpt_b_cent0to10->Fill(jetPtj,w_tot);
                                        h_jeteta_b_cent0to10->Fill(jetEtaj,w_tot);
                                        h_jetphi_b_cent0to10->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>20 && em->hiBin<60){
                                        h_jetpt_b_cent10to30->Fill(jetPtj,w_tot);
                                        h_jeteta_b_cent10to30->Fill(jetEtaj,w_tot);
                                        h_jetphi_b_cent10to30->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>60 && em->hiBin<100){
                                        h_jetpt_b_cent30to50->Fill(jetPtj,w_tot);
                                        h_jeteta_b_cent30to50->Fill(jetEtaj,w_tot);
                                        h_jetphi_b_cent30to50->Fill(jetPhij,w_tot);
                                }
                                if(em->hiBin>100 && em->hiBin<180){
                                        h_jetpt_b_cent50to90->Fill(jetPtj,w_tot);
                                        h_jeteta_b_cent50to90->Fill(jetEtaj,w_tot);
                                        h_jetphi_b_cent50to90->Fill(jetPhij,w_tot);
                                }
			}
			else{
				if(em->matchedPartonFlavor[jetj]==1 || em->matchedPartonFlavor[jetj]==-1 || em->matchedPartonFlavor[jetj]==2 || em->matchedPartonFlavor[jetj]==-2 
                                || em->matchedPartonFlavor[jetj]==3 || em->matchedPartonFlavor[jetj]==-3){
					h_jetpt_uds->Fill(jetPtj,w_tot);
					h_jeteta_uds->Fill(jetEtaj,w_tot);
					h_jetphi_uds->Fill(jetPhij,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_jetpt_uds_cent0to10->Fill(jetPtj,w_tot);
                                                h_jeteta_uds_cent0to10->Fill(jetEtaj,w_tot);
                                                h_jetphi_uds_cent0to10->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_jetpt_uds_cent10to30->Fill(jetPtj,w_tot);
                                                h_jeteta_uds_cent10to30->Fill(jetEtaj,w_tot);
                                                h_jetphi_uds_cent10to30->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_jetpt_uds_cent30to50->Fill(jetPtj,w_tot);
                                                h_jeteta_uds_cent30to50->Fill(jetEtaj,w_tot);
                                                h_jetphi_uds_cent30to50->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_jetpt_uds_cent50to90->Fill(jetPtj,w_tot);
                                                h_jeteta_uds_cent50to90->Fill(jetEtaj,w_tot);
                                                h_jetphi_uds_cent50to90->Fill(jetPhij,w_tot);
                                        }
				}
				if(em->matchedPartonFlavor[jetj]==4 || em->matchedPartonFlavor[jetj]==-4){
					h_jetpt_c->Fill(jetPtj,w_tot);
					h_jeteta_c->Fill(jetEtaj,w_tot);
					h_jetphi_c->Fill(jetPhij,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_jetpt_c_cent0to10->Fill(jetPtj,w_tot);
                                                h_jeteta_c_cent0to10->Fill(jetEtaj,w_tot);
                                                h_jetphi_c_cent0to10->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_jetpt_c_cent10to30->Fill(jetPtj,w_tot);
                                                h_jeteta_c_cent10to30->Fill(jetEtaj,w_tot);
                                                h_jetphi_c_cent10to30->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_jetpt_c_cent30to50->Fill(jetPtj,w_tot);
                                                h_jeteta_c_cent30to50->Fill(jetEtaj,w_tot);
                                                h_jetphi_c_cent30to50->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_jetpt_c_cent50to90->Fill(jetPtj,w_tot);
                                                h_jeteta_c_cent50to90->Fill(jetEtaj,w_tot);
                                                h_jetphi_c_cent50to90->Fill(jetPhij,w_tot);
                                        }
				}
				if(em->matchedPartonFlavor[jetj]==5 || em->matchedPartonFlavor[jetj]==-5){
					h_jetpt_b->Fill(jetPtj,w_tot);
					h_jeteta_b->Fill(jetEtaj,w_tot);
					h_jetphi_b->Fill(jetPhij,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_jetpt_b_cent0to10->Fill(jetPtj,w_tot);
                                                h_jeteta_b_cent0to10->Fill(jetEtaj,w_tot);
                                                h_jetphi_b_cent0to10->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_jetpt_b_cent10to30->Fill(jetPtj,w_tot);
                                                h_jeteta_b_cent10to30->Fill(jetEtaj,w_tot);
                                                h_jetphi_b_cent10to30->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_jetpt_b_cent30to50->Fill(jetPtj,w_tot);
                                                h_jeteta_b_cent30to50->Fill(jetEtaj,w_tot);
                                                h_jetphi_b_cent30to50->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_jetpt_b_cent50to90->Fill(jetPtj,w_tot);
                                                h_jeteta_b_cent50to90->Fill(jetEtaj,w_tot);
                                                h_jetphi_b_cent50to90->Fill(jetPhij,w_tot);
                                        }
				}
				if(em->matchedPartonFlavor[jetj]==21){
					h_jetpt_g->Fill(jetPtj,w_tot);
					h_jeteta_g->Fill(jetEtaj,w_tot);
					h_jetphi_g->Fill(jetPhij,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_jetpt_g_cent0to10->Fill(jetPtj,w_tot);
                                                h_jeteta_g_cent0to10->Fill(jetEtaj,w_tot);
                                                h_jetphi_g_cent0to10->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_jetpt_g_cent10to30->Fill(jetPtj,w_tot);
                                                h_jeteta_g_cent10to30->Fill(jetEtaj,w_tot);
                                                h_jetphi_g_cent10to30->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_jetpt_g_cent30to50->Fill(jetPtj,w_tot);
                                                h_jeteta_g_cent30to50->Fill(jetEtaj,w_tot);
                                                h_jetphi_g_cent30to50->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_jetpt_g_cent50to90->Fill(jetPtj,w_tot);
                                                h_jeteta_g_cent50to90->Fill(jetEtaj,w_tot);
                                                h_jetphi_g_cent50to90->Fill(jetPhij,w_tot);
                                        }
				}
				if(em->matchedPartonFlavor[jetj]==0){
					h_jetpt_ghost->Fill(jetPtj,w_tot);
					h_jeteta_ghost->Fill(jetEtaj,w_tot);
					h_jetphi_ghost->Fill(jetPhij,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_jetpt_ghost_cent0to10->Fill(jetPtj,w_tot);
                                                h_jeteta_ghost_cent0to10->Fill(jetEtaj,w_tot);
                                                h_jetphi_ghost_cent0to10->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_jetpt_ghost_cent10to30->Fill(jetPtj,w_tot);
                                                h_jeteta_ghost_cent10to30->Fill(jetEtaj,w_tot);
                                                h_jetphi_ghost_cent10to30->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_jetpt_ghost_cent30to50->Fill(jetPtj,w_tot);
                                                h_jeteta_ghost_cent30to50->Fill(jetEtaj,w_tot);
                                                h_jetphi_ghost_cent30to50->Fill(jetPhij,w_tot);
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_jetpt_ghost_cent50to90->Fill(jetPtj,w_tot);
                                                h_jeteta_ghost_cent50to90->Fill(jetEtaj,w_tot);
                                                h_jetphi_ghost_cent50to90->Fill(jetPhij,w_tot);
                                        }
				}
			}
			// muon loop
			if(em->nMu==0){continue;}
			double deltaRmin=100000.0;
			double muRelPtMin=-999.0;
			double muJetPtMin = -999.0;
			double muJetPhiMin = -999.0;
			double muJetEtaMin = -999.0;
			double muPtMin = -999.0;
			for(int mui = 0; mui< em->nMu; mui++){
				
				// muon cuts
				if(em->muIsTracker->at(mui)==0 || TMath::Abs(em->muEta->at(mui))>2.4 || em->muChi2NDF->at(mui)==-99 || em->muChi2NDF->at(mui)>10
					||TMath::Abs(em->muInnerD0->at(mui))>0.2 || TMath::Abs(em->muInnerDz->at(mui))>0.5 || em->muMuonHits->at(mui)<= 0
					|| em->muStations->at(mui)<= 1 || em->muTrkLayers->at(mui)<=5 || em->muPixelHits->at(mui)<=0 || em->muPt->at(mui) < muPtCut){continue;}
				//cout << "muon pt =  " << em->muPt->at(mui-1) << endl;
				double muPti = em->muPt->at(mui);
				double muEtai = em->muEta->at(mui);
				double muPhii = em->muPhi->at(mui);

                                double deltaEtaij = muEtai-jetEtaj;
				double deltaPhiij = acos(cos(muPhii-jetPhij));
				double deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				double muRelPt = getPtRel(muPti,muEtai,muPhii,jetPtj,jetEtaj,jetPhij);
				
				if(deltaRij<deltaRmin){
					deltaRmin = deltaRij;
					muRelPtMin = muRelPt;
					muJetPtMin = jetPtj;
                                        muJetEtaMin = jetEtaj;
					muJetPhiMin = jetPhij;
					muPtMin = muPti;
				}
			} // end muon loop
                        if(deltaRmin==100000.0){continue;}	
                        h_deltaR->Fill(deltaRmin,w_tot);
                        h_muPt_noCut->Fill(muPtMin,w_tot);
                        if(muPtMin==-999.0 || muRelPtMin==999.0 || muJetPtMin==999.0 || muJetEtaMin==-999.0 || muJetPhiMin==-999.0){continue;}
                        if(deltaRmin<0.4){
                                h_muRelPt->Fill(muRelPtMin,w_tot);
                                h_muInJetPt->Fill(muPtMin,w_tot);
                                if(em->hiBin>0 && em->hiBin<20){
                                        h_muRelPt_cent0to10->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_cent0to10->Fill(muPtMin,w_tot);
                                }
                                if(em->hiBin>20 && em->hiBin<60){
                                        h_muRelPt_cent10to30->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_cent10to30->Fill(muPtMin,w_tot);
                                }
                                if(em->hiBin>60 && em->hiBin<100){
                                        h_muRelPt_cent30to50->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_cent30to50->Fill(muPtMin,w_tot);
                                }
                                if(em->hiBin>100 && em->hiBin<180){
                                        h_muRelPt_cent50to90->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_cent50to90->Fill(muPtMin,w_tot);
                                }
                                // flavor parse
				if(em->matchedHadronFlavor[jetj]==4){
					h_muRelPt_c->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_c->Fill(muPtMin,w_tot);
				}
				else if(em->matchedHadronFlavor[jetj]==5){
					h_muRelPt_b->Fill(muRelPtMin,w_tot);
                                        h_muInJetPt_b->Fill(muPtMin,w_tot);
				}
				else{
					if(em->matchedPartonFlavor[jetj]==1 || em->matchedPartonFlavor[jetj]==-1 || em->matchedPartonFlavor[jetj]==2 || 
                                        em->matchedPartonFlavor[jetj]==-2 || em->matchedPartonFlavor[jetj]==3 || em->matchedPartonFlavor[jetj]==-3){
                                                h_muRelPt_uds->Fill(muRelPtMin,w_tot);
                                                h_muInJetPt_uds->Fill(muPtMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_muRelPt_uds_cent0to10->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_uds_cent0to10->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_muRelPt_uds_cent10to30->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_uds_cent10to30->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_muRelPt_uds_cent30to50->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_uds_cent30to50->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_muRelPt_uds_cent50to90->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_uds_cent50to90->Fill(muPtMin,w_tot);
                                                }
                                        }
					if(em->matchedPartonFlavor[jetj]==4 || em->matchedPartonFlavor[jetj]==-4){
                                                h_muRelPt_c->Fill(muRelPtMin,w_tot);
                                                h_muInJetPt_c->Fill(muPtMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_muRelPt_c_cent0to10->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_c_cent0to10->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_muRelPt_c_cent10to30->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_c_cent10to30->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_muRelPt_c_cent30to50->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_c_cent30to50->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_muRelPt_c_cent50to90->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_c_cent50to90->Fill(muPtMin,w_tot);
                                                }
                                        }
					if(em->matchedPartonFlavor[jetj]==5 || em->matchedPartonFlavor[jetj]==-5){
                                                h_muRelPt_b->Fill(muRelPtMin,w_tot);
                                                h_muInJetPt_b->Fill(muPtMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_muRelPt_b_cent0to10->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_b_cent0to10->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_muRelPt_b_cent10to30->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_b_cent10to30->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_muRelPt_b_cent30to50->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_b_cent30to50->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_muRelPt_b_cent50to90->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_b_cent50to90->Fill(muPtMin,w_tot);
                                                }
                                        }
					if(em->matchedPartonFlavor[jetj]==21){
                                                h_muRelPt_g->Fill(muRelPtMin,w_tot);
                                                h_muInJetPt_g->Fill(muPtMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_muRelPt_g_cent0to10->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_g_cent0to10->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_muRelPt_g_cent10to30->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_g_cent10to30->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_muRelPt_g_cent30to50->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_g_cent30to50->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_muRelPt_g_cent50to90->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_g_cent50to90->Fill(muPtMin,w_tot);
                                                }
                                        }
					if(em->matchedPartonFlavor[jetj]==0){
                                                h_muRelPt_ghost->Fill(muRelPtMin,w_tot);
                                                h_muInJetPt_ghost->Fill(muPtMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_muRelPt_ghost_cent0to10->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_ghost_cent0to10->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_muRelPt_ghost_cent10to30->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_ghost_cent10to30->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_muRelPt_ghost_cent30to50->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_ghost_cent30to50->Fill(muPtMin,w_tot);
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_muRelPt_ghost_cent50to90->Fill(muRelPtMin,w_tot);
                                                        h_muInJetPt_ghost_cent50to90->Fill(muPtMin,w_tot);
                                                }
                                        }
				}
				
				h_mujetpt->Fill(muJetPtMin,w_tot);
				h_mujeteta->Fill(muJetEtaMin,w_tot);
				h_mujetphi->Fill(muJetPhiMin,w_tot);
                                if(em->hiBin>0 && em->hiBin<20){
                                        h_mujetpt_cent0to10->Fill(muJetPtMin,w_tot);
                                        h_mujeteta_cent0to10->Fill(muJetEtaMin,w_tot);
                                        h_mujetphi_cent0to10->Fill(muJetPhiMin,w_tot);  
                                }
                                if(em->hiBin>20 && em->hiBin<60){
                                        h_mujetpt_cent10to30->Fill(muJetPtMin,w_tot);
                                        h_mujeteta_cent10to30->Fill(muJetEtaMin,w_tot);
                                        h_mujetphi_cent10to30->Fill(muJetPhiMin,w_tot); 
                                }
                                if(em->hiBin>60 && em->hiBin<100){
                                        h_mujetpt_cent30to50->Fill(muJetPtMin,w_tot);
                                        h_mujeteta_cent30to50->Fill(muJetEtaMin,w_tot);
                                        h_mujetphi_cent30to50->Fill(muJetPhiMin,w_tot); 
                                }
                                if(em->hiBin>100 && em->hiBin<180){
                                        h_mujetpt_cent50to90->Fill(muJetPtMin,w_tot);
                                        h_mujeteta_cent50to90->Fill(muJetEtaMin,w_tot);
                                        h_mujetphi_cent50to90->Fill(muJetPhiMin,w_tot); 
                                }
				if(em->matchedHadronFlavor[jetj]==4){
					h_mujetpt_c->Fill(muJetPtMin,w_tot);
					h_mujeteta_c->Fill(muJetEtaMin,w_tot);
					h_mujetphi_c->Fill(muJetPhiMin,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_mujetpt_c_cent0to10->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_c_cent0to10->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_c_cent0to10->Fill(muJetPhiMin,w_tot);  
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_mujetpt_c_cent10to30->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_c_cent10to30->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_c_cent10to30->Fill(muJetPhiMin,w_tot); 
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_mujetpt_c_cent30to50->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_c_cent30to50->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_c_cent30to50->Fill(muJetPhiMin,w_tot); 
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_mujetpt_c_cent50to90->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_c_cent50to90->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_c_cent50to90->Fill(muJetPhiMin,w_tot); 
                                        }
				}
                                else if(em->matchedHadronFlavor[jetj]==5){	
                                        h_mujetpt_b->Fill(muJetPtMin,w_tot);
                                        h_mujeteta_b->Fill(muJetEtaMin,w_tot);
                                        h_mujetphi_b->Fill(muJetPhiMin,w_tot);
                                        if(em->hiBin>0 && em->hiBin<20){
                                                h_mujetpt_b_cent0to10->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_b_cent0to10->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_b_cent0to10->Fill(muJetPhiMin,w_tot);  
                                        }
                                        if(em->hiBin>20 && em->hiBin<60){
                                                h_mujetpt_b_cent10to30->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_b_cent10to30->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_b_cent10to30->Fill(muJetPhiMin,w_tot); 
                                        }
                                        if(em->hiBin>60 && em->hiBin<100){
                                                h_mujetpt_b_cent30to50->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_b_cent30to50->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_b_cent30to50->Fill(muJetPhiMin,w_tot); 
                                        }
                                        if(em->hiBin>100 && em->hiBin<180){
                                                h_mujetpt_b_cent50to90->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_b_cent50to90->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_b_cent50to90->Fill(muJetPhiMin,w_tot); 
                                        }
                                }
                                else{
                                        if(em->matchedPartonFlavor[jetj]==1 || em->matchedPartonFlavor[jetj]==-1 || em->matchedPartonFlavor[jetj]==2 || em->matchedPartonFlavor[jetj]==-2 
                                        || em->matchedPartonFlavor[jetj]==3 || em->matchedPartonFlavor[jetj]==-3){
                                                h_mujetpt_uds->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_uds->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_uds->Fill(muJetPhiMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_mujetpt_uds_cent0to10->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_uds_cent0to10->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_uds_cent0to10->Fill(muJetPhiMin,w_tot);  
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_mujetpt_uds_cent10to30->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_uds_cent10to30->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_uds_cent10to30->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_mujetpt_uds_cent30to50->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_uds_cent30to50->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_uds_cent30to50->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_mujetpt_uds_cent50to90->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_uds_cent50to90->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_uds_cent50to90->Fill(muJetPhiMin,w_tot); 
                                                }
                                        }
                                        if(em->matchedPartonFlavor[jetj]==4 || em->matchedPartonFlavor[jetj]==-4){
                                                h_mujetpt_c->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_c->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_c->Fill(muJetPhiMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_mujetpt_c_cent0to10->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_c_cent0to10->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_c_cent0to10->Fill(muJetPhiMin,w_tot);  
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_mujetpt_c_cent10to30->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_c_cent10to30->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_c_cent10to30->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_mujetpt_c_cent30to50->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_c_cent30to50->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_c_cent30to50->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_mujetpt_c_cent50to90->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_c_cent50to90->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_c_cent50to90->Fill(muJetPhiMin,w_tot); 
                                                }
                                        }
                                        if(em->matchedPartonFlavor[jetj]==5 || em->matchedPartonFlavor[jetj]==-5){
                                                h_mujetpt_b->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_b->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_b->Fill(muJetPhiMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_mujetpt_b_cent0to10->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_b_cent0to10->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_b_cent0to10->Fill(muJetPhiMin,w_tot);  
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_mujetpt_b_cent10to30->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_b_cent10to30->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_b_cent10to30->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_mujetpt_b_cent30to50->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_b_cent30to50->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_b_cent30to50->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_mujetpt_b_cent50to90->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_b_cent50to90->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_b_cent50to90->Fill(muJetPhiMin,w_tot); 
                                                }
                                        }
                                        if(em->matchedPartonFlavor[jetj]==21){
                                                h_mujetpt_g->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_g->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_g->Fill(muJetPhiMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_mujetpt_g_cent0to10->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_g_cent0to10->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_g_cent0to10->Fill(muJetPhiMin,w_tot);  
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_mujetpt_g_cent10to30->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_g_cent10to30->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_g_cent10to30->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_mujetpt_g_cent30to50->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_g_cent30to50->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_g_cent30to50->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_mujetpt_g_cent50to90->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_g_cent50to90->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_g_cent50to90->Fill(muJetPhiMin,w_tot); 
                                                }
                                        }
                                        if(em->matchedPartonFlavor[jetj]==0){
                                                h_mujetpt_ghost->Fill(muJetPtMin,w_tot);
                                                h_mujeteta_ghost->Fill(muJetEtaMin,w_tot);
                                                h_mujetphi_ghost->Fill(muJetPhiMin,w_tot);
                                                if(em->hiBin>0 && em->hiBin<20){
                                                        h_mujetpt_ghost_cent0to10->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_ghost_cent0to10->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_ghost_cent0to10->Fill(muJetPhiMin,w_tot);  
                                                }
                                                if(em->hiBin>20 && em->hiBin<60){
                                                        h_mujetpt_ghost_cent10to30->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_ghost_cent10to30->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_ghost_cent10to30->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>60 && em->hiBin<100){
                                                        h_mujetpt_ghost_cent30to50->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_ghost_cent30to50->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_ghost_cent30to50->Fill(muJetPhiMin,w_tot); 
                                                }
                                                if(em->hiBin>100 && em->hiBin<180){
                                                        h_mujetpt_ghost_cent50to90->Fill(muJetPtMin,w_tot);
                                                        h_mujeteta_ghost_cent50to90->Fill(muJetEtaMin,w_tot);
                                                        h_mujetphi_ghost_cent50to90->Fill(muJetPhiMin,w_tot); 
                                                }
                                        }
                                }
                        }

		} // end reco jet loop
		// gen jet loop
		for(int jeti = 0; jeti < em->ngj; jeti++){
			//cout << "reco jet pt_" << jeti << " = " << em->jetpt[jeti] <<"  ////  " << "gen jet pt_"<<jeti<<" = "<<em->genjetpt[jeti]<<endl;
			
			if(fabs(em->genjeteta[jeti])>etaMax){continue;}
			if(em->genjetpt[jeti]<ptMin || em->genjetpt[jeti]>ptMax){continue;}
			
			//float ref_multiplier = 0.35;
			//if(em->pthat < ref_multiplier*em->ref_jetpt[jeti]){continue;}
			
			h_genjetpt->Fill(em->genjetpt[jeti],em->weight);
			h_genjeteta->Fill(em->genjeteta[jeti],em->weight);
			h_genjetphi->Fill(em->genjetphi[jeti],em->weight);


		} // end gen jet loop

	} // end event loop

auto wf = TFile::Open(output,"recreate");

h_jetphi->Write();
        h_jetphi_cent0to10->Write();
        h_jetphi_cent10to30->Write();
        h_jetphi_cent30to50->Write();
        h_jetphi_cent50to90->Write();
	h_jetphi_g->Write();
        h_jetphi_g_cent0to10->Write();
        h_jetphi_g_cent10to30->Write();
        h_jetphi_g_cent30to50->Write();
        h_jetphi_g_cent50to90->Write();
	h_jetphi_uds->Write();
        h_jetphi_uds_cent0to10->Write();
        h_jetphi_uds_cent10to30->Write();
        h_jetphi_uds_cent30to50->Write();
        h_jetphi_uds_cent50to90->Write();
	h_jetphi_b->Write();
        h_jetphi_b_cent0to10->Write();
        h_jetphi_b_cent10to30->Write();
        h_jetphi_b_cent30to50->Write();
        h_jetphi_b_cent50to90->Write();
	h_jetphi_c->Write();
        h_jetphi_c_cent0to10->Write();
        h_jetphi_c_cent10to30->Write();
        h_jetphi_c_cent30to50->Write();
        h_jetphi_c_cent50to90->Write();
	h_jetphi_ghost->Write();
        h_jetphi_ghost_cent0to10->Write();
        h_jetphi_ghost_cent10to30->Write();
        h_jetphi_ghost_cent30to50->Write();
        h_jetphi_ghost_cent50to90->Write();
h_mujetphi->Write();
        h_mujetphi_cent0to10->Write();
        h_mujetphi_cent10to30->Write();
        h_mujetphi_cent30to50->Write();
        h_mujetphi_cent50to90->Write();
	h_mujetphi_g->Write();
        h_mujetphi_g_cent0to10->Write();
        h_mujetphi_g_cent10to30->Write();
        h_mujetphi_g_cent30to50->Write();
        h_mujetphi_g_cent50to90->Write();
	h_mujetphi_uds->Write();
        h_mujetphi_uds_cent0to10->Write();
        h_mujetphi_uds_cent10to30->Write();
        h_mujetphi_uds_cent30to50->Write();
        h_mujetphi_uds_cent50to90->Write();
	h_mujetphi_b->Write();
        h_mujetphi_b_cent0to10->Write();
        h_mujetphi_b_cent10to30->Write();
        h_mujetphi_b_cent30to50->Write();
        h_mujetphi_b_cent50to90->Write();
	h_mujetphi_c->Write();
        h_mujetphi_c_cent0to10->Write();
        h_mujetphi_c_cent10to30->Write();
        h_mujetphi_c_cent30to50->Write();
        h_mujetphi_c_cent50to90->Write();
	h_mujetphi_ghost->Write();
        h_mujetphi_ghost_cent0to10->Write();
        h_mujetphi_ghost_cent10to30->Write();
        h_mujetphi_ghost_cent30to50->Write();
        h_mujetphi_ghost_cent50to90->Write();
h_genjetphi->Write();
h_trkPhi->Write();
        h_trkPhi_cent0to10->Write();
        h_trkPhi_cent10to30->Write();
        h_trkPhi_cent30to50->Write();
        h_trkPhi_cent50to90->Write();
h_jeteta->Write();
        h_jeteta_cent0to10->Write();
        h_jeteta_cent10to30->Write();
        h_jeteta_cent30to50->Write();
        h_jeteta_cent50to90->Write();
	h_jeteta_g->Write();
        h_jeteta_g_cent0to10->Write();
        h_jeteta_g_cent10to30->Write();
        h_jeteta_g_cent30to50->Write();
        h_jeteta_g_cent50to90->Write();
	h_jeteta_uds->Write();
        h_jeteta_uds_cent0to10->Write();
        h_jeteta_uds_cent10to30->Write();
        h_jeteta_uds_cent30to50->Write();
        h_jeteta_uds_cent50to90->Write();
	h_jeteta_b->Write();
        h_jeteta_b_cent0to10->Write();
        h_jeteta_b_cent10to30->Write();
        h_jeteta_b_cent30to50->Write();
        h_jeteta_b_cent50to90->Write();
	h_jeteta_c->Write();
        h_jeteta_c_cent0to10->Write();
        h_jeteta_c_cent10to30->Write();
        h_jeteta_c_cent30to50->Write();
        h_jeteta_c_cent50to90->Write();
	h_jeteta_ghost->Write();
        h_jeteta_ghost_cent0to10->Write();
        h_jeteta_ghost_cent10to30->Write();
        h_jeteta_ghost_cent30to50->Write();
        h_jeteta_ghost_cent50to90->Sumw2();
h_mujeteta->Sumw2();
        h_mujeteta_cent0to10->Write();
        h_mujeteta_cent10to30->Write();
        h_mujeteta_cent30to50->Write();
        h_mujeteta_cent50to90->Write();
	h_mujeteta_g->Write();
        h_mujeteta_g_cent0to10->Write();
        h_mujeteta_g_cent10to30->Write();
        h_mujeteta_g_cent30to50->Write();
        h_mujeteta_g_cent50to90->Write();
	h_mujeteta_uds->Write();
        h_mujeteta_uds_cent0to10->Write();
        h_mujeteta_uds_cent10to30->Write();
        h_mujeteta_uds_cent30to50->Write();
        h_mujeteta_uds_cent50to90->Write();
	h_mujeteta_b->Write();
        h_mujeteta_b_cent0to10->Write();
        h_mujeteta_b_cent10to30->Write();
        h_mujeteta_b_cent30to50->Write();
        h_mujeteta_b_cent50to90->Write();
	h_mujeteta_c->Write();
        h_mujeteta_c_cent0to10->Write();
        h_mujeteta_c_cent10to30->Write();
        h_mujeteta_c_cent30to50->Write();
        h_mujeteta_c_cent50to90->Write();
	h_mujeteta_ghost->Write();
        h_mujeteta_ghost_cent0to10->Write();
        h_mujeteta_ghost_cent10to30->Write();
        h_mujeteta_ghost_cent30to50->Write();
        h_mujeteta_ghost_cent50to90->Write();
h_genjeteta->Write();
h_trkEta->Write();
        h_trkEta_cent0to10->Write();
        h_trkEta_cent10to30->Write();
        h_trkEta_cent30to50->Write();
        h_trkEta_cent50to90->Write();
h_jetpt_raw->Write();
h_jetpt->Write();
        h_jetpt_cent0to10->Write();
        h_jetpt_cent10to30->Write();
        h_jetpt_cent30to50->Write();
        h_jetpt_cent50to90->Write();
	h_jetpt_g->Write();
        h_jetpt_g_cent0to10->Write();
        h_jetpt_g_cent10to30->Write();
        h_jetpt_g_cent30to50->Write();
        h_jetpt_g_cent50to90->Write();
	h_jetpt_uds->Write();
        h_jetpt_uds_cent0to10->Write();
        h_jetpt_uds_cent10to30->Write();
        h_jetpt_uds_cent30to50->Write();
        h_jetpt_uds_cent50to90->Write();
	h_jetpt_b->Write();
        h_jetpt_b_cent0to10->Write();
        h_jetpt_b_cent10to30->Write();
        h_jetpt_b_cent30to50->Write();
        h_jetpt_b_cent50to90->Write();
	h_jetpt_c->Write();
        h_jetpt_c_cent0to10->Write();
        h_jetpt_c_cent10to30->Write();
        h_jetpt_c_cent30to50->Write();
        h_jetpt_c_cent50to90->Write();
	h_jetpt_ghost->Write();
        h_jetpt_ghost_cent0to10->Write();
        h_jetpt_ghost_cent10to30->Write();
        h_jetpt_ghost_cent30to50->Write();
        h_jetpt_ghost_cent50to90->Write();
h_mujetpt->Write();
        h_mujetpt_cent0to10->Write();
        h_mujetpt_cent10to30->Write();
        h_mujetpt_cent30to50->Write();
        h_mujetpt_cent50to90->Write();
	h_mujetpt_g->Write();
        h_mujetpt_g_cent0to10->Write();
        h_mujetpt_g_cent10to30->Write();
        h_mujetpt_g_cent30to50->Write();
        h_mujetpt_g_cent50to90->Write();
	h_mujetpt_uds->Write();
        h_mujetpt_uds_cent0to10->Write();
        h_mujetpt_uds_cent10to30->Write();
        h_mujetpt_uds_cent30to50->Write();
        h_mujetpt_uds_cent50to90->Write();
	h_mujetpt_b->Write();
        h_mujetpt_b_cent0to10->Write();
        h_mujetpt_b_cent10to30->Write();
        h_mujetpt_b_cent30to50->Write();
        h_mujetpt_b_cent50to90->Write();
	h_mujetpt_c->Write();
        h_mujetpt_c_cent0to10->Write();
        h_mujetpt_c_cent10to30->Write();
        h_mujetpt_c_cent30to50->Write();
        h_mujetpt_c_cent50to90->Write();
	h_mujetpt_ghost->Write();
        h_mujetpt_ghost_cent0to10->Write();
        h_mujetpt_ghost_cent10to30->Write();
        h_mujetpt_ghost_cent30to50->Write();
        h_mujetpt_ghost_cent50to90->Write();
h_genjetpt->Write();
h_trkPt->Write();
        h_trkPt_cent0to10->Write();
        h_trkPt_cent10to30->Write();
        h_trkPt_cent30to50->Write();
        h_trkPt_cent50to90->Write();
h_muPt->Write();
        h_muPt_cent0to10->Write();
        h_muPt_cent10to30->Write();
        h_muPt_cent30to50->Write();
        h_muPt_cent50to90->Write();
h_muPt_noCut->Write();
h_muInJetPt->Write();
        h_muInJetPt_cent0to10->Write();
        h_muInJetPt_cent10to30->Write();
        h_muInJetPt_cent30to50->Write();
        h_muInJetPt_cent50to90->Write();
    h_muInJetPt_g->Write();
        h_muInJetPt_g_cent0to10->Write();
        h_muInJetPt_g_cent10to30->Write();
        h_muInJetPt_g_cent30to50->Write();
        h_muInJetPt_g_cent50to90->Write();
	h_muInJetPt_uds->Write();
        h_muInJetPt_uds_cent0to10->Write();
        h_muInJetPt_uds_cent10to30->Write();
        h_muInJetPt_uds_cent30to50->Write();
        h_muInJetPt_uds_cent50to90->Write();
	h_muInJetPt_b->Write();
        h_muInJetPt_b_cent0to10->Write();
        h_muInJetPt_b_cent10to30->Write();
        h_muInJetPt_b_cent30to50->Write();
        h_muInJetPt_b_cent50to90->Write();
	h_muInJetPt_c->Write();
        h_muInJetPt_c_cent0to10->Write();
        h_muInJetPt_c_cent10to30->Write();
        h_muInJetPt_c_cent30to50->Write();
        h_muInJetPt_c_cent50to90->Write();
	h_muInJetPt_ghost->Write();
        h_muInJetPt_ghost_cent0to10->Write();
        h_muInJetPt_ghost_cent10to30->Write();
        h_muInJetPt_ghost_cent30to50->Write();
        h_muInJetPt_ghost_cent50to90->Write();
h_muEta->Write();
        h_muEta_cent0to10->Write();
        h_muEta_cent10to30->Write();
        h_muEta_cent30to50->Write();
        h_muEta_cent50to90->Write();
h_muPhi->Write();
        h_muPhi_cent0to10->Write();
        h_muPhi_cent10to30->Write();
        h_muPhi_cent30to50->Write();
        h_muPhi_cent50to90->Write();
h_muRelPt->Write();
        h_muRelPt_cent0to10->Write();
        h_muRelPt_cent10to30->Write();
        h_muRelPt_cent30to50->Write();
        h_muRelPt_cent50to90->Write();
	h_muRelPt_g->Write();
        h_muRelPt_g_cent0to10->Write();
        h_muRelPt_g_cent10to30->Write();
        h_muRelPt_g_cent30to50->Write();
        h_muRelPt_g_cent50to90->Write();
	h_muRelPt_uds->Write();
        h_muRelPt_uds_cent0to10->Write();
        h_muRelPt_uds_cent10to30->Write();
        h_muRelPt_uds_cent30to50->Write();
        h_muRelPt_uds_cent50to90->Write();
	h_muRelPt_b->Write();
        h_muRelPt_b_cent0to10->Write();
        h_muRelPt_b_cent10to30->Write();
        h_muRelPt_b_cent30to50->Write();
        h_muRelPt_b_cent50to90->Write();
	h_muRelPt_c->Write();
        h_muRelPt_c_cent0to10->Write();
        h_muRelPt_c_cent10to30->Write();
        h_muRelPt_c_cent30to50->Write();
        h_muRelPt_c_cent50to90->Write();
	h_muRelPt_ghost->Write();
        h_muRelPt_ghost_cent0to10->Write();
        h_muRelPt_ghost_cent10to30->Write();
        h_muRelPt_ghost_cent30to50->Write();
        h_muRelPt_ghost_cent50to90->Write();
h_pthat->Write();
h_vz->Write();
h_hiBin->Write();
h_deltaR->Write();
h_hadronFlavor->Write();
h_partonFlavor->Write();

wf->Close();

return;
} // end program

