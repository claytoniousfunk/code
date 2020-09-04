#include "/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/myProcesses/hiforest/plugin/eventMap_hiForest.h"
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

void PbPb_data_skim(TString input, TString output){
// Define histograms to be filled with data
const int NPhiBins = 300;
const double phiMin = -TMath::Pi();
const double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
    TH1D *h_jetphi_cent0to10 = new TH1D("h_jetphi_cent0to10","reco jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
    TH1D *h_jetphi_cent10to30 = new TH1D("h_jetphi_cent10to30","reco jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
    TH1D *h_jetphi_cent30to50 = new TH1D("h_jetphi_cent30to50","reco jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
    TH1D *h_jetphi_cent50to90 = new TH1D("h_jetphi_cent50to90","reco jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
TH1D *h_mujetphi = new TH1D("h_mujetphi","reco muon jet #phi",NPhiBins,phiMin,phiMax);
    TH1D *h_mujetphi_cent0to10 = new TH1D("h_mujetphi_cent0to10","reco muon jet #phi, 0-10 %",NPhiBins,phiMin,phiMax);
    TH1D *h_mujetphi_cent10to30 = new TH1D("h_mujetphi_cent10to30","reco muon jet #phi, 10-30 %",NPhiBins,phiMin,phiMax);
    TH1D *h_mujetphi_cent30to50 = new TH1D("h_mujetphi_cent30to50","reco muon jet #phi, 30-50 %",NPhiBins,phiMin,phiMax);
    TH1D *h_mujetphi_cent50to90 = new TH1D("h_mujetphi_cent50to90","reco muon jet #phi, 50-90 %",NPhiBins,phiMin,phiMax);
TH1D *h_trkPhi = new TH1D("h_trkPhi","track #phi",NPhiBins,phiMin,phiMax);
    TH1D *h_trkPhi_cent0to10 = new TH1D("h_trkPhi_cent0to10","track #phi, 0-10 %",NPhiBins,phiMin,phiMax);
    TH1D *h_trkPhi_cent10to30 = new TH1D("h_trkPhi_cent10to30","track #phi, 10-30 %",NPhiBins,phiMin,phiMax);
    TH1D *h_trkPhi_cent30to50 = new TH1D("h_trkPhi_cent30to50","track #phi, 30-50 %",NPhiBins,phiMin,phiMax);
    TH1D *h_trkPhi_cent50to90 = new TH1D("h_trkPhi_cent50to90","track #phi, 50-90 %",NPhiBins,phiMin,phiMax);
const int NEtaBins = 300;
const double etaMin = -1.5;
const double etaMax = 1.5;
const int NTrkEtaBins = 450;
const double trkEtaMin = -2.4;
const double trkEtaMax = 2.4;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
    TH1D *h_jeteta_cent0to10 = new TH1D("h_jeteta_cent0to10","reco jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
    TH1D *h_jeteta_cent10to30 = new TH1D("h_jeteta_cent10to30","reco jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
    TH1D *h_jeteta_cent30to50 = new TH1D("h_jeteta_cent30to50","reco jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
    TH1D *h_jeteta_cent50to90 = new TH1D("h_jeteta_cent50to90","reco jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
TH1D *h_mujeteta = new TH1D("h_mujeteta","reco muon jet #eta",NEtaBins,etaMin,etaMax);
    TH1D *h_mujeteta_cent0to10 = new TH1D("h_mujeteta_cent0to10","reco muon jet #eta, 0-10 %",NEtaBins,etaMin,etaMax);
    TH1D *h_mujeteta_cent10to30 = new TH1D("h_mujeteta_cent10to30","reco muon jet #eta, 10-30 %",NEtaBins,etaMin,etaMax);
    TH1D *h_mujeteta_cent30to50 = new TH1D("h_mujeteta_cent30to50","reco muon jet #eta, 30-50 %",NEtaBins,etaMin,etaMax);
    TH1D *h_mujeteta_cent50to90 = new TH1D("h_mujeteta_cent50to90","reco muon jet #eta, 50-90 %",NEtaBins,etaMin,etaMax);
TH1D *h_trkEta = new TH1D("h_trkEta","track #eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
    TH1D *h_trkEta_cent0to10 = new TH1D("h_trkEta_cent0to10","track #eta, 0-10 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
    TH1D *h_trkEta_cent10to30 = new TH1D("h_trkEta_cent10to30","track #eta, 10-30 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
    TH1D *h_trkEta_cent30to50 = new TH1D("h_trkEta_cent30to50","track #eta, 30-50 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
    TH1D *h_trkEta_cent50to90 = new TH1D("h_trkEta_cent50to90","track #eta, 50-90 %",NTrkEtaBins,trkEtaMin,trkEtaMax);
const int NPtBins = 450;
const double ptMin = 50.0;
const double ptMax = 500.0;
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet p_{T}",NPtBins,ptMin,ptMax);
    TH1D *h_jetpt_cent0to10 = new TH1D("h_jetpt_cent0to10","reco jet p_{T}, 0-10 %",NPtBins,ptMin,ptMax);
    TH1D *h_jetpt_cent10to30 = new TH1D("h_jetpt_cent10to30","reco jet p_{T}, 10-30 %",NPtBins,ptMin,ptMax);
    TH1D *h_jetpt_cent30to50 = new TH1D("h_jetpt_cent30to50","reco jet p_{T}, 30-50 %",NPtBins,ptMin,ptMax);
    TH1D *h_jetpt_cent50to90 = new TH1D("h_jetpt_cent50to90","reco jet p_{T}, 50-90 %",NPtBins,ptMin,ptMax);
TH1D *h_mujetpt = new TH1D("h_mujetpt","reco muon jet p_{T}",NPtBins,ptMin,ptMax);
    TH1D *h_mujetpt_cent0to10 = new TH1D("h_mujetpt_cent0to10","reco muon jet p_{T}, 0-10 %",NPtBins,ptMin,ptMax);
    TH1D *h_mujetpt_cent10to30 = new TH1D("h_mujetpt_cent10to30","reco muon jet p_{T}, 10-30 %",NPtBins,ptMin,ptMax);
    TH1D *h_mujetpt_cent30to50 = new TH1D("h_mujetpt_cent30to50","reco muon jet p_{T}, 30-50 %",NPtBins,ptMin,ptMax);
    TH1D *h_mujetpt_cent50to90 = new TH1D("h_mujetpt_cent50to90","reco muon jet p_{T}, 50-90 %",NPtBins,ptMin,ptMax);
const int NTrkPtBins = 500;
const double trkPtMin = 0.0;
const double trkPtMax = 500.0;
TH1D *h_trkPt = new TH1D("h_trkPt","track p_{T}",NTrkPtBins,trkPtMin,trkPtMax);
    TH1D *h_trkPt_cent0to10 = new TH1D("h_trkPt_cent0to10","track p_{T}, 0-10 %",NTrkPtBins,trkPtMin,trkPtMax);
    TH1D *h_trkPt_cent10to30 = new TH1D("h_trkPt_cent10to30","track p_{T}, 10-30 %",NTrkPtBins,trkPtMin,trkPtMax);
    TH1D *h_trkPt_cent30to50 = new TH1D("h_trkPt_cent30to50","track p_{T}, 30-50 %",NTrkPtBins,trkPtMin,trkPtMax);
    TH1D *h_trkPt_cent50to90 = new TH1D("h_trkPt_cent50to90","track p_{T}, 50-90 %",NTrkPtBins,trkPtMin,trkPtMax);
// muon data
int NMuPtBins = 500;
double muPtMin = 0.0;
double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon p_{T}",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muPt_cent0to10 = new TH1D("h_muPt_cent0to10","muon p_{T}, 0-10 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muPt_cent10to30 = new TH1D("h_muPt_cent10to30","muon p_{T}, 10-30 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muPt_cent30to50 = new TH1D("h_muPt_cent30to50","muon p_{T}, 30-50 %",NMuPtBins,muPtMin,muPtMax);
    TH1D *h_muPt_cent50to90 = new TH1D("h_muPt_cent50to90","muon p_{T}, 50-90 %",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon #eta",NEtaBins,etaMin,etaMax);
    TH1D *h_muEta_cent0to10 = new TH1D("h_muEta_cent0to10","muon #eta, 0-10 %",NEtaBins,etaMin,etaMax);
    TH1D *h_muEta_cent10to30 = new TH1D("h_muEta_cent10to30","muon #eta, 10-30 %",NEtaBins,etaMin,etaMax);
    TH1D *h_muEta_cent30to50 = new TH1D("h_muEta_cent30to50","muon #eta, 30-50 %",NEtaBins,etaMin,etaMax);
    TH1D *h_muEta_cent50to90 = new TH1D("h_muEta_cent50to90","muon #eta, 50-90 %",NEtaBins,etaMin,etaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
    TH1D *h_muPhi_cent0to10 = new TH1D("h_muPhi_cent0to10","muon #phi, 0-10 %",NPhiBins,phiMin,phiMax);
    TH1D *h_muPhi_cent10to30 = new TH1D("h_muPhi_cent10to30","muon #phi, 10-30 %",NPhiBins,phiMin,phiMax);
    TH1D *h_muPhi_cent30to50 = new TH1D("h_muPhi_cent30to50","muon #phi, 30-50 %",NPhiBins,phiMin,phiMax);
    TH1D *h_muPhi_cent50to90 = new TH1D("h_muPhi_cent50to90","muon #phi, 50-90 %",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
    TH1D *h_muRelPt_cent0to10 = new TH1D("h_muRelPt_cent0to10","muon rel pt, 0-10 %",50,0.0,5.0);
    TH1D *h_muRelPt_cent10to30 = new TH1D("h_muRelPt_cent10to30","muon rel pt, 10-30 %",50,0.0,5.0);
    TH1D *h_muRelPt_cent30to50 = new TH1D("h_muRelPt_cent30to50","muon rel pt 30-50 %",50,0.0,5.0);
    TH1D *h_muRelPt_cent50to90 = new TH1D("h_muRelPt_cent50to90","muon rel pt, 50-90 %",50,0.0,5.0);
// event info
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
h_mujetphi->Sumw2();
    h_mujetphi_cent0to10->Sumw2();
    h_mujetphi_cent10to30->Sumw2();
    h_mujetphi_cent30to50->Sumw2();
    h_mujetphi_cent50to90->Sumw2();
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
h_mujeteta->Sumw2();
    h_mujeteta_cent0to10->Sumw2();
    h_mujeteta_cent10to30->Sumw2();
    h_mujeteta_cent30to50->Sumw2();
    h_mujeteta_cent50to90->Sumw2();
h_trkEta->Sumw2();
    h_trkEta_cent0to10->Sumw2();
    h_trkEta_cent10to30->Sumw2();
    h_trkEta_cent30to50->Sumw2();
    h_trkEta_cent50to90->Sumw2();
h_jetpt->Sumw2();
    h_jetpt_cent0to10->Sumw2();
    h_jetpt_cent10to30->Sumw2();
    h_jetpt_cent30to50->Sumw2();
    h_jetpt_cent50to90->Sumw2();
h_mujetpt->Sumw2();
    h_mujetpt_cent0to10->Sumw2();
    h_mujetpt_cent10to30->Sumw2();
    h_mujetpt_cent30to50->Sumw2();
    h_mujetpt_cent50to90->Sumw2();
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
h_vz->Sumw2();
h_hiBin->Sumw2();
h_deltaR->Sumw2();

TFile *f = TFile::Open(input);
	auto em = new eventMap(f);
        em->isMC = 0;
        em->init();
        em->loadTrack();
        em->loadJet("akPu4PFJetAnalyzer");
        em->loadMuon("ggHiNtuplizer");
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

		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
           	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
        }
        evi_frac = 100*evi/NEvents;
	    
        // track loop
		for(int trki=0; trki < em->ntrk; trki++){

			double trkPti = em->trkpt[trki];
			double trkEtai = em->trketa[trki];
			double trkPhii = em->trkphi[trki];

			const double trkPtCut = 1.0;

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

        double muPtCut = 10.0;

		// seperate muon loop
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
                    h_muPt_cent0to10->Fill(muPti);
                    h_muEta_cent0to10->Fill(muEtai);
                    h_muPhi_cent0to10->Fill(muPhii);
                }
                if(em->hiBin>20 && em->hiBin<60){
                    h_muPt_cent10to30->Fill(muPti);
                    h_muEta_cent10to30->Fill(muEtai);
                    h_muPhi_cent10to30->Fill(muPhii);
                }
                if(em->hiBin>60 && em->hiBin<100){
                    h_muPt_cent30to50->Fill(muPti);
                    h_muEta_cent30to50->Fill(muEtai);
                    h_muPhi_cent30to50->Fill(muPhii);
                }
                if(em->hiBin>100 && em->hiBin<180){
                    h_muPt_cent50to90->Fill(muPti);
                    h_muEta_cent50to90->Fill(muEtai);
                    h_muPhi_cent50to90->Fill(muPhii);
                }
			
		}
  
		//reco jet loop
		for(int jetj = 0; jetj < em->njet; jetj++){

            double jetPtj = em->jetpt[jetj];
			double jetEtaj = em->jeteta[jetj];
			double jetPhij = em->jetphi[jetj];
			
	
            if(jetPtj<ptMin || jetPtj>ptMax || fabs(jetEtaj)>etaMax || jetPhij==-999){continue;}
				
			h_jetpt->Fill(jetPtj);
			h_jeteta->Fill(jetEtaj);
			h_jetphi->Fill(jetPhij);

            if(em->hiBin>0 && em->hiBin<20){
                h_jetpt_cent0to10->Fill(jetPtj);
                h_jeteta_cent0to10->Fill(jetEtaj);
                h_jetphi_cent0to10->Fill(jetPhij);
            }
            if(em->hiBin>20 && em->hiBin<60){
                h_jetpt_cent10to30->Fill(jetPtj);
                h_jeteta_cent10to30->Fill(jetEtaj);
                h_jetphi_cent10to30->Fill(jetPhij);
            }
            if(em->hiBin>60 && em->hiBin<100){
                h_jetpt_cent30to50->Fill(jetPtj);
                h_jeteta_cent30to50->Fill(jetEtaj);
                h_jetphi_cent30to50->Fill(jetPhij);
            }
            if(em->hiBin>100 && em->hiBin<180){
                h_jetpt_cent50to90->Fill(jetPtj);
                h_jeteta_cent50to90->Fill(jetEtaj);
                h_jetphi_cent50to90->Fill(jetPhij);
            }
	
			// muon loop
			if(em->nMu==0){continue;}
			double deltaRmin=100000.0;
			double muRelPtMin=-999.0;
			double muJetPtMin = -999.0;
			double muJetPhiMin = -999.0;
			double muJetEtaMin = -999.0;
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
				}
			} // end muon loop
				if(deltaRmin==100000.0){continue;}	
				if(muRelPtMin==-999.0 || muJetPtMin==-999.0 || muJetEtaMin==-999.0 || muJetPhiMin==-999.0){continue;}

				h_deltaR->Fill(deltaRmin);
				
				if(deltaRmin<0.4){
					h_muRelPt->Fill(muRelPtMin);
					h_mujetpt->Fill(muJetPtMin);
					h_mujeteta->Fill(muJetEtaMin);
					h_mujetphi->Fill(muJetPhiMin);
                    
                    if(em->hiBin>0 && em->hiBin<20){
                        h_muRelPt_cent0to10->Fill(muRelPtMin);
                        h_mujetpt_cent0to10->Fill(muJetPtMin);
                        h_mujeteta_cent0to10->Fill(muJetEtaMin);
                        h_mujetphi_cent0to10->Fill(muJetPhiMin);
                    }
                    if(em->hiBin>20 && em->hiBin<60){
                        h_muRelPt_cent10to30->Fill(muRelPtMin);
                        h_mujetpt_cent10to30->Fill(muJetPtMin);
                        h_mujeteta_cent10to30->Fill(muJetEtaMin);
                        h_mujetphi_cent10to30->Fill(muJetPhiMin);
                    }
                    if(em->hiBin>60 && em->hiBin<100){
                        h_muRelPt_cent30to50->Fill(muRelPtMin);
                        h_mujetpt_cent30to50->Fill(muJetPtMin);
                        h_mujeteta_cent30to50->Fill(muJetEtaMin);
                        h_mujetphi_cent30to50->Fill(muJetPhiMin);
                    }
                    if(em->hiBin>100 && em->hiBin<180){
                        h_muRelPt_cent50to90->Fill(muRelPtMin);
                        h_mujetpt_cent50to90->Fill(muJetPtMin);
                        h_mujeteta_cent50to90->Fill(muJetEtaMin);
                        h_mujetphi_cent50to90->Fill(muJetPhiMin);
                    }
				}

		} // end reco jet loop

	} // end event loop



auto wf = TFile::Open(output,"recreate");

h_jetphi->Write();
    h_jetphi_cent0to10->Write();
    h_jetphi_cent10to30->Write();
    h_jetphi_cent30to50->Write();
    h_jetphi_cent50to90->Write();
h_mujetphi->Write();
    h_mujetphi_cent0to10->Write();
    h_mujetphi_cent10to30->Write();
    h_mujetphi_cent30to50->Write();
    h_mujetphi_cent50to90->Write();
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
h_mujeteta->Write();
    h_mujeteta_cent0to10->Write();
    h_mujeteta_cent10to30->Write();
    h_mujeteta_cent30to50->Write();
    h_mujeteta_cent50to90->Write();
h_trkEta->Write();
    h_trkEta_cent0to10->Write();
    h_trkEta_cent10to30->Write();
    h_trkEta_cent30to50->Write();
    h_trkEta_cent50to90->Write();
h_jetpt->Write();
    h_jetpt_cent0to10->Write();
    h_jetpt_cent10to30->Write();
    h_jetpt_cent30to50->Write();
    h_jetpt_cent50to90->Write();
h_mujetpt->Write();
    h_mujetpt_cent0to10->Write();
    h_mujetpt_cent10to30->Write();
    h_mujetpt_cent30to50->Write();
    h_mujetpt_cent50to90->Write();
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
h_vz->Write();
h_hiBin->Write();
h_deltaR->Write();

wf->Close();

return;
} // end program

