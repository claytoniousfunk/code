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

void data_skim(TString input, TString output){
// Define histograms to be filled with data
int NPhiBins = 300;
double phiMin = -TMath::Pi();
double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
int NEtaBins = 300;
double etaMin = -1.5;
double etaMax = 1.5;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
int NPtBins = 450;
double ptMin = 50.0;
double ptMax = 500.0;
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
// muon data
int NMuPtBins = 500;
double muPtMin = 0.0;
double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muPt_noRcut = new TH1D("h_muPt_noRcut","muon pt, no R cut",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NEtaBins,etaMin,etaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
// event info
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",100,0,1);


// Sumw2 commands
h_jetphi->Sumw2();
h_jeteta->Sumw2();
h_genjeteta->Sumw2();
h_jetpt->Sumw2();
h_muPt->Sumw2();
h_muPt_noRcut->Sumw2();
h_muEta->Sumw2();
h_muPhi->Sumw2();
h_muRelPt->Sumw2();
h_vz->Sumw2();
h_deltaR->Sumw2();

TFile *f = TFile::Open(input);
	auto em = new eventMap(f);
        em->isMC = 0;
        em->init();
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
		h_vz->Fill(em->vz);

		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/NEvents;
		//jet cuts
    double jetEtaCut = 1.5;
	double jetPtCut = 50.0;
		//reco jet loop
		for(int jeti = 0; jeti < em->njet; jeti++){
			
			if(fabs(em->jeteta[jeti])>jetEtaCut){continue;}
			if(em->jetpt[jeti]<jetPtCut){continue;}
				
			h_jetpt->Fill(em->jetpt[jeti]);
			h_jeteta->Fill(em->jeteta[jeti]);
			h_jetphi->Fill(em->jetphi[jeti]);
	
			// muon loop
			if(em->nMu==0){continue;}
			double deltaRmin=100000.0;
			double muRelPtMin=-999.0;
			double muPtMin = -999.0;
			double muPhiMin = -999.0;
			double muEtaMin = -999.0;
			for(int mui = 0; mui< em->nMu; mui++){
				double muPtCut = 10.0;
				// muon cuts
				if(em->muIsTracker->at(mui)==0 || TMath::Abs(em->muEta->at(mui))>2.4 || em->muChi2NDF->at(mui)==-99 || em->muChi2NDF->at(mui)>10
					||TMath::Abs(em->muInnerD0->at(mui))>0.2 || TMath::Abs(em->muInnerDz->at(mui))>0.5 || em->muMuonHits->at(mui)<= 0
					|| em->muStations->at(mui)<= 1 || em->muTrkLayers->at(mui)<=5 || em->muPixelHits->at(mui)<=0 || em->muPt->at(mui) < muPtCut){continue;}
				//cout << "muon pt =  " << em->muPt->at(mui-1) << endl;
				double deltaEtaij = em->muEta->at(mui)-em->jeteta[jeti];
				double deltaPhiij = acos(cos(em->muPhi->at(mui)-em->jetphi[jeti]));
				double deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				double muRelPt = getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]);
				
				if(deltaRij<deltaRmin){
					deltaRmin = deltaRij;
					muRelPtMin = muRelPt;
					muPtMin = em->muPt->at(mui);
					muEtaMin = em->muEta->at(mui);
					muPhiMin = em->muPhi->at(mui);
				}
			} // end muon loop
				if(deltaRmin==100000.0){continue;}	
				if(muRelPtMin==999.0 || muPtMin==999.0 || muEtaMin==-999.0 || muPhiMin==-999.0){continue;}

				h_deltaR->Fill(deltaRmin);
				h_muPt_noRcut->Fill(muPtMin);
				
				if(deltaRmin<0.4){
					h_muRelPt->Fill(muRelPtMin);
					h_muPt->Fill(muPtMin);
					h_muEta->Fill(muEtaMin);
					h_muPhi->Fill(muPhiMin);
				}
				


		} // end reco jet loop

	} // end event loop



auto wf = TFile::Open(output,"recreate");

h_jetpt->Write();
h_jeteta->Write();
h_jetphi->Write();
h_muPt->Write();
h_muPt_noRcut->Write();
h_muEta->Write();
h_muPhi->Write();
h_muRelPt->Write();
h_vz->Write();
h_pthat->Write();
h_deltaR->Write();

wf->Close();

return;
} // end program

