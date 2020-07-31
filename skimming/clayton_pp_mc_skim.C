#include "/home/clayton/Analysis/code/myProcessesEdit/hiforest/plugin/eventMap_hiForest.h"
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

// Define histograms to be filled with data
int NPhiBins = 300;
double phiMin = -TMath::Pi();
double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_genjetphi = new TH1D("h_genjetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
int NEtaBins = 300;
double etaMin = -1.5;
double etaMax = 1.5;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_genjeteta = new TH1D("h_genjeteta","gen jet #eta",NEtaBins,etaMin,etaMax);
int NPtBins = 450;
double ptMin = 50.0;
double ptMax = 500.0;
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
TH1D *h_genjetpt = new TH1D("h_genjetpt","gen jet pt",NPtBins,ptMin,ptMax);
// muon data
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NPtBins,ptMin,ptMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NEtaBins,etaMin,etaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
// event info
TH1D *h_pthat = new TH1D("h_pthat","pthat",NPtBins,ptMin,ptMax);
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);


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

void clayton_pp_mc_skim(){

// Sumw2 commands
h_jetphi->Sumw2();
h_genjetphi->Sumw2();

h_jeteta->Sumw2();
h_genjeteta->Sumw2();

h_jetpt->Sumw2();
h_genjetpt->Sumw2();

h_muPt->Sumw2();
h_muEta->Sumw2();
h_muPhi->Sumw2();
h_muRelPt->Sumw2();

h_pthat->Sumw2();
h_vz->Sumw2();

// count number of files to be skimmed over
int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/data/fnal");

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

// file loop

for(int file = 1; file < NFiles+1; file++){
	cout << "Processing file " << file << "/"<<NFiles<< endl;
	auto f = TFile::Open(Form("/home/clayton/Analysis/data/fnal/HiForestAOD_%d.root",file));

	auto em = new eventMap(f);
        em->isMC = 1;
        em->init();
        em->loadJet("ak4PFJetAnalyzer");
        em->loadMuon("ggHiNtuplizerGED");
		Long64_t NEvents = em->evtTree->GetEntries();
		//cout << "Number of events = " << NEvents << endl;

	// event loop
	int evi_frac = 0;
	for(int evi = 0; evi < NEvents; evi++){
	//for(int evi = 1; evi < 20; evi++){	
		em->getEvent(evi);
		h_vz->Fill(em->vz);
		h_pthat->Fill(em->pthat);
		//cout << "vz = " << em->vz << endl;
		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/NEvents;
		//cout << "Number of jets = " << em->njet << endl;
		for(int jeti = 0; jeti < em->njet; jeti++){
			//cout << "reco jet pt_" << jeti << " = " << em->jetpt[jeti] <<"  ////  " << "gen jet pt_"<<jeti<<" = "<<em->genjetpt[jeti]<<endl;
			h_jetpt->Fill(em->jetpt[jeti],em->weight);
			h_jeteta->Fill(em->jeteta[jeti],em->weight);
			h_jetphi->Fill(em->jetphi[jeti],em->weight);
			h_genjetpt->Fill(em->genjetpt[jeti],em->weight);
			h_genjeteta->Fill(em->genjeteta[jeti],em->weight);
			h_genjetphi->Fill(em->genjetphi[jeti],em->weight);

			for(int mui = 0; mui< em->nMu; mui++){
				//cout << "muon pt =  " << em->muPt->at(mui-1) << endl;
				h_muPt->Fill(em->muPt->at(mui),em->weight);
				h_muEta->Fill(em->muEta->at(mui),em->weight);
				h_muPhi->Fill(em->muPhi->at(mui),em->weight);
				h_muRelPt->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]),em->weight);
			}


		}

	} // end event loop



} // end file loop

auto wf = TFile::Open("/home/clayton/Analysis/code/skimming/rootFiles/clayton_pp_mc_skim.root","recreate");

h_jetpt->Write();
h_jeteta->Write();
h_jetphi->Write();
h_genjetpt->Write();
h_genjeteta->Write();
h_genjetphi->Write();
h_muPt->Write();
h_muEta->Write();
h_muPhi->Write();
h_muRelPt->Write();
h_vz->Write();
h_pthat->Write();

wf->Close();

return;
} // end program

