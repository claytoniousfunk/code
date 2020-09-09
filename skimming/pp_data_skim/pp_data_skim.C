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

void pp_data_skim(TString input, TString output){
// Define histograms to be filled with data
const int NPhiBins = 300;
double phiMin = -TMath::Pi();
double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_mujetphi = new TH1D("h_mujetphi","reco muon jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_trkPhi = new TH1D("h_trkPhi","track #phi",NPhiBins,phiMin,phiMax);
const int NEtaBins = 260;
const double etaMin = -1.3;
const double etaMax = 1.3;
const int NTrkEtaBins = 480;
const double trkEtaMin = -2.4;
const double trkEtaMax = 2.4;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_mujeteta = new TH1D("h_mujeteta","reco muon jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_trkEta = new TH1D("h_trkEta","track #eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
const int NPtBins = 450;
const double ptMin = 50.0;
const double ptMax = 500.0;
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet p_{T}",NPtBins,ptMin,ptMax);
TH1D *h_mujetpt = new TH1D("h_mujetpt","reco muon jet p_{T}",NPtBins,ptMin,ptMax);
const int NTrkPtBins = 500;
const double trkPtMin = 0.0;
const double trkPtMax = 500.0;
TH1D *h_trkPt = new TH1D("h_trkPt","track p_{T}",NTrkPtBins,trkPtMin,trkPtMax);
// muon data
const int NMuPtBins = 500;
const double muPtMin = 0.0;
const double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon p_{T}",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muInJetPt = new TH1D("h_muInJetPt","in-jet muon p_{T}",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon #eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
// event info
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
TH1D *h_hiBin = new TH1D("h_hiBin","hiBin (2*centrality)",200,0,200);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",100,0,1);


// Sumw2 commands
h_jetphi->Sumw2();
h_mujetphi->Sumw2();
h_trkPhi->Sumw2();
h_jeteta->Sumw2();
h_mujeteta->Sumw2();
h_trkEta->Sumw2();
h_jetpt->Sumw2();
h_mujetpt->Sumw2();
h_trkPt->Sumw2();
h_muPt->Sumw2();
h_muInJetPt->Sumw2();
h_muEta->Sumw2();
h_muPhi->Sumw2();
h_muRelPt->Sumw2();
h_vz->Sumw2();
h_hiBin->Sumw2();
h_deltaR->Sumw2();

TFile *f = TFile::Open(input);
	auto em = new eventMap(f);
        em->isMC = 0;
        em->init();
	em->loadTrack();
        em->loadJet("ak4PFJetAnalyzer");
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

		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
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
	
			// muon loop
			if(em->nMu==0){continue;}
			double deltaRmin=100000.0;
			double muRelPtMin=-999.0;
			double muJetPtMin = -999.0;
			double muJetPhiMin = -999.0;
			double muJetEtaMin = -999.0;
			double muPtMin = -999.0;
			// 
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
				if(muPtMin==-999.0 || muRelPtMin==-999.0 || muJetPtMin==-999.0 || muJetEtaMin==-999.0 || muJetPhiMin==-999.0){continue;}

				h_deltaR->Fill(deltaRmin);
				
				if(deltaRmin<0.4){
					h_muRelPt->Fill(muRelPtMin);
					h_mujetpt->Fill(muJetPtMin);
					h_mujeteta->Fill(muJetEtaMin);
					h_mujetphi->Fill(muJetPhiMin);
					h_muInJetPt->Fill(muPtMin);
				}

		} // end reco jet loop

	} // end event loop



auto wf = TFile::Open(output,"recreate");

h_jetphi->Write();
h_mujetphi->Write();
h_trkPhi->Write();
h_jeteta->Write();
h_mujeteta->Write();
h_trkEta->Write();
h_jetpt->Write();
h_mujetpt->Write();
h_trkPt->Write();
h_muPt->Write();
h_muInJetPt->Write();
h_muEta->Write();
h_muPhi->Write();
h_muRelPt->Write();
h_vz->Write();
h_hiBin->Write();
h_deltaR->Write();

wf->Close();

return;
} // end program

