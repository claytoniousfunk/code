#include "/uscms/home/cmbennet/work/myProcesses/hiforest/plugin/eventMap_hiForest.h"
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

void pp_mc_skim(){
// Define histograms to be filled with data
const int NPhiBins = 300;
const double phiMin = -TMath::Pi();
const double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_g = new TH1D("h_jetphi_g","g reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_uds = new TH1D("h_jetphi_uds","uds reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_b = new TH1D("h_jetphi_b","b reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_c = new TH1D("h_jetphi_c","c reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_ghost = new TH1D("h_jetphi_ghost","ghost reco jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_mujetphi = new TH1D("h_mujetphi","reco muon jet #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_g = new TH1D("h_mujetphi_g","g reco muon jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_uds = new TH1D("h_mujetphi_uds","uds reco muon jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_b = new TH1D("h_mujetphi_b","b reco muon jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_c = new TH1D("h_mujetphi_c","c reco muon jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_mujetphi_ghost = new TH1D("h_mujetphi_ghost","ghost reco jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_genjetphi = new TH1D("h_genjetphi","gen jet #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_g = new TH1D("h_genjetphi_g","g gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_uds = new TH1D("h_genjetphi_uds","uds gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_b = new TH1D("h_genjetphi_b","b gen jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_genjetphi_c = new TH1D("h_genjetphi_c","c gen jet #phi", NPhiBins,phiMin,phiMax);
TH1D *h_trkPhi = new TH1D("h_trkPhi","track #phi",NPhiBins,phiMin,phiMax);
const int NEtaBins = 260;
const double etaMin = -1.3;
const double etaMax = 1.3;
const int NTrkEtaBins = 480;
const double trkEtaMin = -2.4;
const double trkEtaMax = 2.4;
TH1D *h_jeteta = new TH1D("h_jeteta","reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_g = new TH1D("h_jeteta_g","g reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_uds = new TH1D("h_jeteta_uds","uds reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_b = new TH1D("h_jeteta_b","b reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_c = new TH1D("h_jeteta_c","c reco jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_jeteta_ghost = new TH1D("h_jeteta_ghost","ghost reco jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_mujeteta = new TH1D("h_mujeteta","reco muon jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_g = new TH1D("h_mujeteta_g","g reco muon jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_uds = new TH1D("h_mujeteta_uds","uds reco muon jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_b = new TH1D("h_mujeteta_b","b reco muon jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_c = new TH1D("h_mujeteta_c","c reco muon jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_mujeteta_ghost = new TH1D("h_mujeteta_ghost","ghost reco muon jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_genjeteta = new TH1D("h_genjeteta","gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_g = new TH1D("h_genjeteta_g","g gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_uds = new TH1D("h_genjeteta_uds","uds gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_b = new TH1D("h_genjeteta_b","b gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_c = new TH1D("h_genjeteta_c","c gen jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_trkEta = new TH1D("h_trkEta","track #eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
const int NPtBins = 450;
const double ptMin = 50.0;
const double ptMax = 500.0;
TH1D *h_jetpt_raw = new TH1D("h_jetpt_raw","raw jet pt",NPtBins,ptMin,ptMax);
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_g = new TH1D("h_jetpt_g","g reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_uds = new TH1D("h_jetpt_uds","uds reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_b = new TH1D("h_jetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_c = new TH1D("h_jetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_ghost = new TH1D("h_jetpt_ghost","ghost reco jet pt",NPtBins,ptMin,ptMax);
TH1D *h_mujetpt = new TH1D("h_mujetpt","reco muon jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_g = new TH1D("h_mujetpt_g","g reco muon jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_uds = new TH1D("h_mujetpt_uds","uds reco muon jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_b = new TH1D("h_mujetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_c = new TH1D("h_mujetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_mujetpt_ghost = new TH1D("h_mujetpt_ghost","ghost reco muon jet pt",NPtBins,ptMin,ptMax);
TH1D *h_genjetpt = new TH1D("h_genjetpt","gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_g = new TH1D("h_genjetpt_g","g gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_uds = new TH1D("h_genjetpt_uds","uds gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_b = new TH1D("h_genjetpt_b","b gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_c = new TH1D("h_genjetpt_c","c gen jet pt",NPtBins,ptMin,ptMax);
const int NTrkPtBins = 500;
const double trkPtMin = 0.0;
const double trkPtMax = 500.0;
TH1D *h_trkPt = new TH1D("h_trkPt","track p_{T}",NTrkPtBins,trkPtMin,trkPtMax);
// muon data
const int NMuPtBins = 500;
const double muPtMin = 0.0;
const double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muPt_noCut = new TH1D("h_muPt_noCut","muon pt - no cut",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muInJetPt = new TH1D("h_muInJetPt","in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muInJetPt_g = new TH1D("h_muPt_g","g in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muInJetPt_uds = new TH1D("h_muPt_uds","uds in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muInJetPt_b = new TH1D("h_muPt_b","b in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muInJetPt_c = new TH1D("h_muPt_c","c in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muInJetPt_ghost = new TH1D("h_muPt_ghost","ghost in-jet muon pt",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NTrkEtaBins,trkEtaMin,trkEtaMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","g muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_uds = new TH1D("h_muRelPt_uds","uds muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_b = new TH1D("h_muRelPt_b","b muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_c = new TH1D("h_muRelPt_c","c muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_ghost = new TH1D("h_muRelPt_ghost","ghost muon rel pt",50,0.0,5.0);
// event info
TH1D *h_pthat = new TH1D("h_pthat","pthat",NPtBins,ptMin,ptMax);
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
TH1D *h_hiBin = new TH1D("h_hiBin","hiBin (2*centrality)",200,0,200);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",100,0,1);


// Sumw2 commands
h_jetphi->Sumw2();
	h_jetphi_g->Sumw2();
	h_jetphi_uds->Sumw2();
	h_jetphi_b->Sumw2();
	h_jetphi_c->Sumw2();
	h_jetphi_ghost->Sumw2();
h_mujetphi->Sumw2();
	h_mujetphi_g->Sumw2();
	h_mujetphi_uds->Sumw2();
	h_mujetphi_b->Sumw2();
	h_mujetphi_c->Sumw2();
	h_mujetphi_ghost->Sumw2();
h_genjetphi->Sumw2();
h_trkPhi->Sumw2();
h_jeteta->Sumw2();
	h_jeteta_g->Sumw2();
	h_jeteta_uds->Sumw2();
	h_jeteta_b->Sumw2();
	h_jeteta_c->Sumw2();
	h_jeteta_ghost->Sumw2();
h_mujeteta->Sumw2();
	h_mujeteta_g->Sumw2();
	h_mujeteta_uds->Sumw2();
	h_mujeteta_b->Sumw2();
	h_mujeteta_c->Sumw2();
	h_mujeteta_ghost->Sumw2();
h_genjeteta->Sumw2();
h_trkEta->Sumw2();
h_jetpt_raw->Sumw2();
h_jetpt->Sumw2();
	h_jetpt_g->Sumw2();
	h_jetpt_uds->Sumw2();
	h_jetpt_b->Sumw2();
	h_jetpt_c->Sumw2();
	h_jetpt_ghost->Sumw2();
h_mujetpt->Sumw2();
	h_mujetpt_g->Sumw2();
	h_mujetpt_uds->Sumw2();
	h_mujetpt_b->Sumw2();
	h_mujetpt_c->Sumw2();
	h_mujetpt_ghost->Sumw2();
h_genjetpt->Sumw2();
h_trkPt->Sumw2();
h_muPt->Sumw2();
h_muPt_noCut->Sumw2();
h_muInJetPt->Sumw2();
	h_muInJetPt_g->Sumw2();
	h_muInJetPt_uds->Sumw2();
	h_muInJetPt_b->Sumw2();
	h_muInJetPt_c->Sumw2();
	h_muInJetPt_ghost->Sumw2();
h_muEta->Sumw2();
h_muPhi->Sumw2();
h_muRelPt->Sumw2();
	h_muRelPt_g->Sumw2();
	h_muRelPt_uds->Sumw2();
	h_muRelPt_b->Sumw2();
	h_muRelPt_c->Sumw2();
	h_muRelPt_ghost->Sumw2();
h_pthat->Sumw2();
h_vz->Sumw2();
h_hiBin->Sumw2();
h_deltaR->Sumw2();

// count number of files to be skimmed over
int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/eos/uscms/store/user/cmbennet/QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/crab_pp2017MC_5TeV_2020/200723_185629/0000");

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
TF1 *vz_fit_fxn = new TF1("vz_fit_fxn","[0] + [1]*x + [2]*x*x + [3]*x*x*x",-15.0,15.0);
vz_fit_fxn->SetParameters(1.02937e+00,-5.25572e-03,-8.21959e-04,2.60358e-05);

TF1 *pt_fit_fxn = new TF1("pt_fit_fxn","[0]*exp([1]*x) + [2]*exp([3]*x)",50.0,500.0);
pt_fit_fxn->SetParameters(2.95395e+00,-2.60850e-02,5.71522e-01,-2.90376e-03);

for(int file = 1; file < NFiles+1; file++){
//for(int file = 82; file < 90; file++){
	cout << "Processing file " << file << "/"<<NFiles<< endl;
        TFile *f = TFile::Open(Form("root://cmsxrootd.fnal.gov//store/user/cmbennet/QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/crab_pp2017MC_5TeV_2020/200723_185629/0000/HiForestAOD_%d.root",file));

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
	//for(int evi = 1; evi < 5; evi++){	
		em->getEvent(evi);
		//event cuts
		// vz reweight
		double w_vz = 1.0/(vz_fit_fxn->Eval(em->vz)); // calculated vz weight
		
		
		if(em->vz>15.0){continue;}
		if(em->pthat<30.0){continue;}
		

		h_vz->Fill(em->vz,w_vz);
		h_hiBin->Fill(em->hiBin);
		h_pthat->Fill(em->pthat);
		//cout << "vz = " << em->vz << endl;
		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/NEvents;
		//cout << "Number of jets = " << em->njet << endl;
        
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
			
			
		}
		//reco jet loop
		for(int jetj = 0; jetj < em->njet; jetj++){

			double jetPtj = em->jetpt[jetj];
			double jetEtaj = em->jeteta[jetj];
			double jetPhij = em->jetphi[jetj];
			//cout << "reco jet pt_" << jeti << " = " << em->jetpt[jeti] <<"  ////  " << "gen jet pt_"<<jeti<<" = "<<em->genjetpt[jeti]<<endl;
			double w_pt = 1.0/(pt_fit_fxn->Eval(em->jetpt[jetj]));
			// double w_tot = w_vz*w_pt*em->weight;	
			double w_tot = em->weight;

			if(jetPtj<ptMin || jetPtj>ptMax || fabs(jetEtaj)>etaMax || jetPhij==-999){continue;}
			
			
			h_jetpt_raw->Fill(jetPtj,w_vz*em->weight);						
			h_jetpt->Fill(jetPtj,w_tot);
			h_jeteta->Fill(jetEtaj,w_tot);
			h_jetphi->Fill(jetPhij,w_tot);
			
			
			if(em->hadronFlavor[jetj]==4){
					h_jetpt_c->Fill(jetPtj,w_tot);
					h_jeteta_c->Fill(jetEtaj,w_tot);
					h_jetphi_c->Fill(jetPhij,w_tot);
			}
			else if(em->hadronFlavor[jetj]==5){
					h_jetpt_b->Fill(jetPtj,w_tot);
					h_jeteta_b->Fill(jetEtaj,w_tot);
					h_jetphi_b->Fill(jetPhij,w_tot);
			}
			else{
				if(em->partonFlavor[jetj]==1 || em->partonFlavor[jetj]==-1){
					h_jetpt_uds->Fill(jetPtj,w_tot);
					h_jeteta_uds->Fill(jetEtaj,w_tot);
					h_jetphi_uds->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==2 || em->partonFlavor[jetj]==-2){
					h_jetpt_uds->Fill(jetPtj,w_tot);
					h_jeteta_uds->Fill(jetEtaj,w_tot);
					h_jetphi_uds->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==3 || em->partonFlavor[jetj]==-3){
					h_jetpt_uds->Fill(jetPtj,w_tot);
					h_jeteta_uds->Fill(jetEtaj,w_tot);
					h_jetphi_uds->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==4 || em->partonFlavor[jetj]==-4){
					h_jetpt_c->Fill(jetPtj,w_tot);
					h_jeteta_c->Fill(jetEtaj,w_tot);
					h_jetphi_c->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==5 || em->partonFlavor[jetj]==-5){
					h_jetpt_b->Fill(jetPtj,w_tot);
					h_jeteta_b->Fill(jetEtaj,w_tot);
					h_jetphi_b->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==21){
					h_jetpt_g->Fill(jetPtj,w_tot);
					h_jeteta_g->Fill(jetEtaj,w_tot);
					h_jetphi_g->Fill(jetPhij,w_tot);
				}
				if(em->partonFlavor[jetj]==0){
					h_jetpt_ghost->Fill(jetPtj,w_tot);
					h_jeteta_ghost->Fill(jetEtaj,w_tot);
					h_jetphi_ghost->Fill(jetPhij,w_tot);
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
					if(em->hadronFlavor[jetj]==4){
						h_muRelPt_c->Fill(muRelPtMin,w_tot);
						h_muInJetPt_c->Fill(muPtMin,w_tot);
					}
					else if(em->hadronFlavor[jetj]==5){
						h_muRelPt_b->Fill(muRelPtMin,w_tot);
						h_muInJetPt_b->Fill(muPtMin,w_tot);
					}
					else{
						if(em->partonFlavor[jetj]==1 || em->partonFlavor[jetj]==-1){h_muRelPt_uds->Fill(muRelPtMin,w_tot);h_muInJetPt_uds->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==2 || em->partonFlavor[jetj]==-2){h_muRelPt_uds->Fill(muRelPtMin,w_tot);h_muInJetPt_uds->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==3 || em->partonFlavor[jetj]==-3){h_muRelPt_uds->Fill(muRelPtMin,w_tot);h_muInJetPt_uds->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==4 || em->partonFlavor[jetj]==-4){h_muRelPt_c->Fill(muRelPtMin,w_tot);h_muInJetPt_c->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==5 || em->partonFlavor[jetj]==-5){h_muRelPt_b->Fill(muRelPtMin,w_tot);h_muInJetPt_b->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==21){h_muRelPt_g->Fill(muRelPtMin,w_tot);h_muInJetPt_g->Fill(muPtMin,w_tot);}
						if(em->partonFlavor[jetj]==0){h_muRelPt_ghost->Fill(muRelPtMin,w_tot);h_muInJetPt_ghost->Fill(muPtMin,w_tot);}
					}
				}
					h_mujetpt->Fill(muJetPtMin,w_tot);
					h_mujeteta->Fill(muJetEtaMin,w_tot);
					h_mujetphi->Fill(muJetPhiMin,w_tot);
					if(em->hadronFlavor[jetj]==4){
						h_mujetpt_c->Fill(muJetPtMin,w_tot);
						h_mujeteta_c->Fill(muJetEtaMin,w_tot);
						h_mujetphi_c->Fill(muJetPhiMin,w_tot);
					}
					else if(em->hadronFlavor[jetj]==5){	
						h_mujetpt_b->Fill(muJetPtMin,w_tot);
						h_mujeteta_b->Fill(muJetEtaMin,w_tot);
						h_mujetphi_b->Fill(muJetPhiMin,w_tot);
					}
					else{
						if(em->partonFlavor[jetj]==1 || em->partonFlavor[jetj]==-1){
							h_mujetpt_uds->Fill(muJetPtMin,w_tot);
							h_mujeteta_uds->Fill(muJetEtaMin,w_tot);
							h_mujetphi_uds->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==2 || em->partonFlavor[jetj]==-2){
							h_mujetpt_uds->Fill(muJetPtMin,w_tot);
							h_mujeteta_uds->Fill(muJetEtaMin,w_tot);
							h_mujetphi_uds->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==3 || em->partonFlavor[jetj]==-3){
							h_mujetpt_uds->Fill(muJetPtMin,w_tot);
							h_mujeteta_uds->Fill(muJetEtaMin,w_tot);
							h_mujetphi_uds->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==4 || em->partonFlavor[jetj]==-4){
							h_mujetpt_c->Fill(muJetPtMin,w_tot);
							h_mujeteta_c->Fill(muJetEtaMin,w_tot);
							h_mujetphi_c->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==5 || em->partonFlavor[jetj]==-5){
							h_mujetpt_b->Fill(muJetPtMin,w_tot);
							h_mujeteta_b->Fill(muJetEtaMin,w_tot);
							h_mujetphi_b->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==21){
							h_mujetpt_g->Fill(muJetPtMin,w_tot);
							h_mujeteta_g->Fill(muJetEtaMin,w_tot);
							h_mujetphi_g->Fill(muJetPhiMin,w_tot);
						}
						if(em->partonFlavor[jetj]==0){
							h_mujetpt_ghost->Fill(muJetPtMin,w_tot);
							h_mujeteta_ghost->Fill(muJetEtaMin,w_tot);
							h_mujetphi_ghost->Fill(muJetPhiMin,w_tot);
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

delete f;
} // end file loop

//auto wf = TFile::Open("/uscms/home/cmbennet/work/pythia_skim_newWeights_pthatcut_30_muptcut_10_16Aug2020.root","recreate");
auto wf = TFile::Open("/uscms/home/cmbennet/work/pp_mc_skim_pthatWeight_pthatcut30_muptcut10_8Sep20.root","recreate");

h_jetpt_raw->Write();
h_jetpt->Write();
	h_jetpt_g->Write();
	h_jetpt_uds->Write();
	h_jetpt_b->Write();
	h_jetpt_c->Write();
	h_jetpt_ghost->Write();
h_mujetpt->Write();
	h_mujetpt_g->Write();
	h_mujetpt_uds->Write();
	h_mujetpt_b->Write();
	h_mujetpt_c->Write();
	h_mujetpt_ghost->Write();
h_jeteta->Write();
	h_jeteta_g->Write();
	h_jeteta_uds->Write();
	h_jeteta_b->Write();
	h_jeteta_c->Write();
	h_jeteta_ghost->Write();
h_mujeteta->Write();
	h_mujeteta_g->Write();
	h_mujeteta_uds->Write();
	h_mujeteta_b->Write();
	h_mujeteta_c->Write();
	h_mujeteta_ghost->Write();
h_jetphi->Write();
	h_jetphi_g->Write();
	h_jetphi_uds->Write();
	h_jetphi_b->Write();
	h_jetphi_c->Write();
	h_jetphi_ghost->Write();
h_mujetphi->Write();
	h_mujetphi_g->Write();
	h_mujetphi_uds->Write();
	h_mujetphi_b->Write();
	h_mujetphi_c->Write();
	h_mujetphi_ghost->Write();
h_genjetpt->Write();
h_genjeteta->Write();
h_genjetphi->Write();
h_muPt->Write();
h_muPt_noCut->Write();
h_muInJetPt->Write();
h_muEta->Write();
h_muPhi->Write();
h_muRelPt->Write();
	h_muRelPt_g->Write();
	h_muRelPt_uds->Write();
	h_muRelPt_b->Write();
	h_muRelPt_c->Write();
	h_muRelPt_ghost->Write();
h_vz->Write();
h_hiBin->Write();
h_pthat->Write();
h_deltaR->Write();
h_trkPt->Write();
h_trkPhi->Write();
h_trkEta->Write();

wf->Close();

return;
} // end program

