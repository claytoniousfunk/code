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
	TH1D *h_jetphi_g = new TH1D("h_jetphi_g","g reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_uds = new TH1D("h_jetphi_uds","uds reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_b = new TH1D("h_jetphi_b","b reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_c = new TH1D("h_jetphi_c","c reco jet #phi", NPhiBins,phiMin,phiMax);
TH1D *h_genjetphi = new TH1D("h_genjetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_g = new TH1D("h_genjetphi_g","g gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_uds = new TH1D("h_genjetphi_uds","uds gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_b = new TH1D("h_genjetphi_b","b gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_c = new TH1D("h_genjetphi_c","c gen jet #phi", NPhiBins,phiMin,phiMax);
int NEtaBins = 300;
double etaMin = -1.5;
double etaMax = 1.5;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_g = new TH1D("h_jeteta_g","g reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_uds = new TH1D("h_jeteta_uds","uds reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_b = new TH1D("h_jeteta_b","b reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_c = new TH1D("h_jeteta_c","c reco jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_genjeteta = new TH1D("h_genjeteta","gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_g = new TH1D("h_genjeteta_g","g gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_uds = new TH1D("h_genjeteta_uds","uds gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_b = new TH1D("h_genjeteta_b","b gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_c = new TH1D("h_genjeteta_c","c gen jet #eta",NEtaBins,etaMin,etaMax);
int NPtBins = 450;
double ptMin = 50.0;
double ptMax = 500.0;
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_g = new TH1D("h_jetpt_g","g reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_uds = new TH1D("h_jetpt_uds","uds reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_b = new TH1D("h_jetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_c = new TH1D("h_jetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
TH1D *h_genjetpt = new TH1D("h_genjetpt","gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_g = new TH1D("h_genjetpt_g","g gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_uds = new TH1D("h_genjetpt_uds","uds gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_b = new TH1D("h_genjetpt_b","b gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_c = new TH1D("h_genjetpt_c","c gen jet pt",NPtBins,ptMin,ptMax);
// muon data
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NPtBins,ptMin,ptMax);
	TH1D *h_muPt_g = new TH1D("h_muPt_g","g muon pt",NPtBins,ptMin,ptMax);
	TH1D *h_muPt_uds = new TH1D("h_muPt_uds","uds muon pt",NPtBins,ptMin,ptMax);
	TH1D *h_muPt_b = new TH1D("h_muPt_b","b muon pt",NPtBins,ptMin,ptMax);
	TH1D *h_muPt_c = new TH1D("h_muPt_c","c muon pt",NPtBins,ptMin,ptMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_g = new TH1D("h_muEta_g","g muon eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_uds = new TH1D("h_muEta_uds","uds muon eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_b = new TH1D("h_muEta_b","b muon eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_c = new TH1D("h_muEta_c","c muon eta",NEtaBins,etaMin,etaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_g = new TH1D("h_muPhi_g","g muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_uds = new TH1D("h_muPhi_uds","uds muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_b = new TH1D("h_muPhi_b","b muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_c = new TH1D("h_muPhi_c","c muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","g muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_uds = new TH1D("h_muRelPt_uds","uds muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_b = new TH1D("h_muRelPt_b","b muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_c = new TH1D("h_muRelPt_c","c muon rel pt",50,0.0,5.0);
// event info
TH1D *h_pthat = new TH1D("h_pthat","pthat",NPtBins,ptMin,ptMax);
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",10,0,1);


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
	h_jetphi_g->Sumw2();
	h_jetphi_uds->Sumw2();
	h_jetphi_b->Sumw2();
	h_jetphi_c->Sumw2();
h_genjetphi->Sumw2();
	h_genjetphi_g->Sumw2();
	h_genjetphi_uds->Sumw2();
	h_genjetphi_b->Sumw2();
	h_genjetphi_c->Sumw2();

h_jeteta->Sumw2();
	h_jeteta_g->Sumw2();
	h_jeteta_uds->Sumw2();
	h_jeteta_b->Sumw2();
	h_jeteta_c->Sumw2();
h_genjeteta->Sumw2();
	h_genjeteta_g->Sumw2();
	h_genjeteta_uds->Sumw2();
	h_genjeteta_b->Sumw2();
	h_genjeteta_c->Sumw2();

h_jetpt->Sumw2();
	h_jetpt_g->Sumw2();
	h_jetpt_uds->Sumw2();
	h_jetpt_b->Sumw2();
	h_jetpt_c->Sumw2();
h_genjetpt->Sumw2();
	h_genjetpt_g->Sumw2();
	h_genjetpt_uds->Sumw2();
	h_genjetpt_b->Sumw2();
	h_genjetpt_c->Sumw2();

h_muPt->Sumw2();
	h_muPt_g->Sumw2();
	h_muPt_uds->Sumw2();
	h_muPt_b->Sumw2();
	h_muPt_c->Sumw2();
h_muEta->Sumw2();
	h_muEta_g->Sumw2();
	h_muEta_uds->Sumw2();
	h_muEta_b->Sumw2();
	h_muEta_c->Sumw2();
h_muPhi->Sumw2();
	h_muPhi_g->Sumw2();
	h_muPhi_uds->Sumw2();
	h_muPhi_b->Sumw2();
	h_muPhi_c->Sumw2();
h_muRelPt->Sumw2();
	h_muRelPt_g->Sumw2();
	h_muRelPt_uds->Sumw2();
	h_muRelPt_b->Sumw2();
	h_muRelPt_c->Sumw2();
h_pthat->Sumw2();
h_vz->Sumw2();

h_deltaR->Sumw2();

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
		//event cuts
		if(fabs(em->vz)>15.0){continue;}

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
			//jet cuts
			double jetEtaCut = 1.5;
			double jetPtCut = 50.0;
			if(fabs(em->jeteta[jeti])>jetEtaCut){continue;}
			if(em->jetpt[jeti]<jetPtCut){continue;}
			float ref_multiplier = 0.35;
			if(em->pthat < ref_multiplier*em->ref_jetpt[jeti]){continue;}
			h_jetpt->Fill(em->jetpt[jeti],em->weight);
			h_jeteta->Fill(em->jeteta[jeti],em->weight);
			h_jetphi->Fill(em->jetphi[jeti],em->weight);
			h_genjetpt->Fill(em->genjetpt[jeti],em->weight);
			h_genjeteta->Fill(em->genjeteta[jeti],em->weight);
			h_genjetphi->Fill(em->genjetphi[jeti],em->weight);
			if(em->partonFlavor[jeti]==4){
					h_jetpt_c->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_c->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_c->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_c->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_c->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_c->Fill(em->genjetphi[jeti],em->weight);
			}
			else if(em->partonFlavor[jeti]==5){
					h_jetpt_b->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_b->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_b->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_b->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_b->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_b->Fill(em->genjetphi[jeti],em->weight);
			}
			else{
				if(em->hadronFlavor[jeti]==1 || em->hadronFlavor[jeti]==-1){
					h_jetpt_uds->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_uds->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_uds->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_uds->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_uds->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_uds->Fill(em->genjetphi[jeti],em->weight);
				}
				if(em->hadronFlavor[jeti]==2 || em->hadronFlavor[jeti]==-2){
					h_jetpt_uds->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_uds->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_uds->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_uds->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_uds->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_uds->Fill(em->genjetphi[jeti],em->weight);
				}
				if(em->hadronFlavor[jeti]==3 || em->hadronFlavor[jeti]==-3){
					h_jetpt_uds->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_uds->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_uds->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_uds->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_uds->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_uds->Fill(em->genjetphi[jeti],em->weight);
				}
				if(em->hadronFlavor[jeti]==4 || em->hadronFlavor[jeti]==-4){
					h_jetpt_c->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_c->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_c->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_c->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_c->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_c->Fill(em->genjetphi[jeti],em->weight);
				}
				if(em->hadronFlavor[jeti]==5 || em->hadronFlavor[jeti]==-5){
					h_jetpt_b->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_b->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_b->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_b->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_b->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_b->Fill(em->genjetphi[jeti],em->weight);
				}
				if(em->hadronFlavor[jeti]==21){
					h_jetpt_g->Fill(em->jetpt[jeti],em->weight);
					h_jeteta_g->Fill(em->jeteta[jeti],em->weight);
					h_jetphi_g->Fill(em->jetphi[jeti],em->weight);
					h_genjetpt_g->Fill(em->genjetpt[jeti],em->weight);
					h_genjeteta_g->Fill(em->genjeteta[jeti],em->weight);
					h_genjetphi_g->Fill(em->genjetphi[jeti],em->weight);
				}
			}
			for(int mui = 0; mui< em->nMu; mui++){
				double muPtCut = 5.0;
				// muon cuts
				if(em->muIsTracker->at(mui)==0 || TMath::Abs(em->muEta->at(mui))>2.4 || em->muChi2NDF->at(mui)==-99 || em->muChi2NDF->at(mui)>10
					||TMath::Abs(em->muInnerD0->at(mui))>0.2 || TMath::Abs(em->muInnerDz->at(mui))>0.5 || em->muMuonHits->at(mui)<= 0
					|| em->muStations->at(mui)<= 1 || em->muTrkLayers->at(mui)<=5 || em->muPixelHits->at(mui)<=0 || em->muPt->at(mui) < muPtCut){continue;}
				//cout << "muon pt =  " << em->muPt->at(mui-1) << endl;
				double deltaEtaij = em->muEta->at(mui)-em->jeteta[jeti];
				double deltaPhiij = acos(cos(em->muPhi->at(mui)-em->jetphi[jeti]));
				double deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				h_deltaR->Fill(deltaRij);
				h_muPt->Fill(em->muPt->at(mui));
				h_muEta->Fill(em->muEta->at(mui));
				h_muPhi->Fill(em->muPhi->at(mui));
				h_muRelPt->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
				if(em->partonFlavor[jeti]==4){
					h_muRelPt_c->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));	
					h_muPt_c->Fill(em->muPt->at(mui));
					h_muEta_c->Fill(em->muEta->at(mui));
					h_muPhi_c->Fill(em->muPhi->at(mui));
				}
				else if(em->partonFlavor[jeti]==5){
					h_muRelPt_b->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));	
					h_muPt_b->Fill(em->muPt->at(mui));
					h_muEta_b->Fill(em->muEta->at(mui));
					h_muPhi_b->Fill(em->muPhi->at(mui));
				}
				else{
					if(em->hadronFlavor[jeti]==1 || em->hadronFlavor[jeti]==-1){
						h_muRelPt_uds->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_uds->Fill(em->muPt->at(mui));
						h_muEta_uds->Fill(em->muEta->at(mui));
						h_muPhi_uds->Fill(em->muPhi->at(mui));
					}
					if(em->hadronFlavor[jeti]==2 || em->hadronFlavor[jeti]==-2){
						h_muRelPt_uds->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_uds->Fill(em->muPt->at(mui));
						h_muEta_uds->Fill(em->muEta->at(mui));
						h_muPhi_uds->Fill(em->muPhi->at(mui));
					}
					if(em->hadronFlavor[jeti]==3 || em->hadronFlavor[jeti]==-3){
						h_muRelPt_uds->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_uds->Fill(em->muPt->at(mui));
						h_muEta_uds->Fill(em->muEta->at(mui));
						h_muPhi_uds->Fill(em->muPhi->at(mui));
					}
					if(em->hadronFlavor[jeti]==4 || em->hadronFlavor[jeti]==-4){
						h_muRelPt_c->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_c->Fill(em->muPt->at(mui));
						h_muEta_c->Fill(em->muEta->at(mui));
						h_muPhi_c->Fill(em->muPhi->at(mui));
					}
					if(em->hadronFlavor[jeti]==5 || em->hadronFlavor[jeti]==-5){
						h_muRelPt_b->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_b->Fill(em->muPt->at(mui));
						h_muEta_b->Fill(em->muEta->at(mui));
						h_muPhi_b->Fill(em->muPhi->at(mui));
					}
					if(em->hadronFlavor[jeti]==21){
						h_muRelPt_g->Fill(getPtRel(em->muPt->at(mui),em->muEta->at(mui),em->muPhi->at(mui),em->jetpt[jeti],em->jeteta[jeti],em->jetphi[jeti]));
						h_muPt_g->Fill(em->muPt->at(mui));
						h_muEta_g->Fill(em->muEta->at(mui));
						h_muPhi_g->Fill(em->muPhi->at(mui));
					}
				}

			}


		}

	} // end event loop



} // end file loop

auto wf = TFile::Open("/home/clayton/Analysis/code/skimming/rootFiles/clayton_pp_mc_skim.root","recreate");

h_jetpt->Write();
	h_jetpt_g->Write();
	h_jetpt_uds->Write();
	h_jetpt_b->Write();
	h_jetpt_c->Write();
h_jeteta->Write();
	h_jeteta_g->Write();
	h_jeteta_uds->Write();
	h_jeteta_b->Write();
	h_jeteta_c->Write();
h_jetphi->Write();
	h_jetphi_g->Write();
	h_jetphi_uds->Write();
	h_jetphi_b->Write();
	h_jetphi_c->Write();
h_genjetpt->Write();
	h_genjetpt_g->Write();
	h_genjetpt_uds->Write();
	h_genjetpt_b->Write();
	h_genjetpt_c->Write();
h_genjeteta->Write();
	h_genjeteta_g->Write();
	h_genjeteta_uds->Write();
	h_genjeteta_b->Write();
	h_genjeteta_c->Write();
h_genjetphi->Write();
	h_genjetphi_g->Write();
	h_genjetphi_uds->Write();
	h_genjetphi_b->Write();
	h_genjetphi_c->Write();
h_muPt->Write();
	h_muPt_g->Write();
	h_muPt_uds->Write();
	h_muPt_b->Write();
	h_muPt_c->Write();
h_muEta->Write();
	h_muEta_g->Write();
	h_muEta_uds->Write();
	h_muEta_b->Write();
	h_muEta_c->Write();
h_muPhi->Write();
	h_muPhi_g->Write();
	h_muPhi_uds->Write();
	h_muPhi_b->Write();
	h_muPhi_c->Write();
h_muRelPt->Write();
	h_muRelPt_g->Write();
	h_muRelPt_uds->Write();
	h_muRelPt_b->Write();
	h_muRelPt_c->Write();
h_vz->Write();
h_pthat->Write();
h_deltaR->Write();

wf->Close();

return;
} // end program

