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

void pythia_skim(){
// Define histograms to be filled with data
int NPhiBins = 300;
double phiMin = -TMath::Pi();
double phiMax = TMath::Pi();
TH1D *h_jetphi = new TH1D("h_jetphi","reco jet #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_g = new TH1D("h_jetphi_g","g reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_uds = new TH1D("h_jetphi_uds","uds reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_b = new TH1D("h_jetphi_b","b reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_c = new TH1D("h_jetphi_c","c reco jet #phi", NPhiBins,phiMin,phiMax);
	TH1D *h_jetphi_ghost = new TH1D("h_jetphi_ghost","ghost reco jet #phi",NPhiBins,phiMin,phiMax);
TH1D *h_genjetphi = new TH1D("h_genjetphi","gen jet #phi",NPhiBins,phiMin,phiMax);
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
	TH1D *h_jeteta_ghost = new TH1D("h_jeteta_ghost","ghost reco jet #eta",NEtaBins,etaMin,etaMax);
TH1D *h_genjeteta = new TH1D("h_genjeteta","gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_g = new TH1D("h_genjeteta_g","g gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_uds = new TH1D("h_genjeteta_uds","uds gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_b = new TH1D("h_genjeteta_b","b gen jet #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_genjeteta_c = new TH1D("h_genjeteta_c","c gen jet #eta",NEtaBins,etaMin,etaMax);
int NPtBins = 450;
double ptMin = 50.0;
double ptMax = 500.0;
TH1D *h_jetpt_raw = new TH1D("h_jetpt_raw","raw jet pt",NPtBins,ptMin,ptMax);
TH1D *h_jetpt = new TH1D("h_jetpt","reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_g = new TH1D("h_jetpt_g","g reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_uds = new TH1D("h_jetpt_uds","uds reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_b = new TH1D("h_jetpt_b","b reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_c = new TH1D("h_jetpt_c","c reco jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_jetpt_ghost = new TH1D("h_jetpt_ghost","ghost reco jet pt",NPtBins,ptMin,ptMax);
TH1D *h_genjetpt = new TH1D("h_genjetpt","gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_g = new TH1D("h_genjetpt_g","g gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_uds = new TH1D("h_genjetpt_uds","uds gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_b = new TH1D("h_genjetpt_b","b gen jet pt",NPtBins,ptMin,ptMax);
	TH1D *h_genjetpt_c = new TH1D("h_genjetpt_c","c gen jet pt",NPtBins,ptMin,ptMax);
// muon data
int NMuPtBins = 500;
double muPtMin = 0.0;
double muPtMax = 500.0;
TH1D *h_muPt = new TH1D("h_muPt","muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_g = new TH1D("h_muPt_g","g muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_uds = new TH1D("h_muPt_uds","uds muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_b = new TH1D("h_muPt_b","b muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_c = new TH1D("h_muPt_c","c muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_ghost = new TH1D("h_muPt_ghost","ghost muon pt",NMuPtBins,muPtMin,muPtMax);
	TH1D *h_muPt_noCut = new TH1D("h_muPt_noCut","muon pt - no cut",NMuPtBins,muPtMin,muPtMax);
TH1D *h_muEta = new TH1D("h_muEta","muon eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_g = new TH1D("h_muEta_g","g muon #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_uds = new TH1D("h_muEta_uds","uds muon #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_b = new TH1D("h_muEta_b","b muon #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_c = new TH1D("h_muEta_c","c muon #eta",NEtaBins,etaMin,etaMax);
	TH1D *h_muEta_ghost = new TH1D("h_muEta_ghost","ghost muon #eta",NPtBins,ptMin,ptMax);
TH1D *h_muPhi = new TH1D("h_muPhi","muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_g = new TH1D("h_muPhi_g","g muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_uds = new TH1D("h_muPhi_uds","uds muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_b = new TH1D("h_muPhi_b","b muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_c = new TH1D("h_muPhi_c","c muon #phi",NPhiBins,phiMin,phiMax);
	TH1D *h_muPhi_ghost = new TH1D("h_muPhi_ghost","ghost muon #phi",NPhiBins,phiMin,phiMax);
TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","g muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_uds = new TH1D("h_muRelPt_uds","uds muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_b = new TH1D("h_muRelPt_b","b muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_c = new TH1D("h_muRelPt_c","c muon rel pt",50,0.0,5.0);
	TH1D *h_muRelPt_ghost = new TH1D("h_muRelPt_ghost","ghost muon rel pt",50,0.0,5.0);
// event info
TH1D *h_pthat = new TH1D("h_pthat","pthat",NPtBins,ptMin,ptMax);
TH1D *h_vz = new TH1D("h_vz","vz",60,-15.0,15.0);
//calculated values
TH1D *h_deltaR = new TH1D("h_deltaR","#Delta r",100,0,1);


// Sumw2 commands
h_jetphi->Sumw2();
	h_jetphi_g->Sumw2();
	h_jetphi_uds->Sumw2();
	h_jetphi_b->Sumw2();
	h_jetphi_c->Sumw2();
	h_jetphi_ghost->Sumw2();
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
	h_jeteta_ghost->Sumw2();
h_genjeteta->Sumw2();
	h_genjeteta_g->Sumw2();
	h_genjeteta_uds->Sumw2();
	h_genjeteta_b->Sumw2();
	h_genjeteta_c->Sumw2();
h_jetpt_raw->Sumw2();
h_jetpt->Sumw2();
	h_jetpt_g->Sumw2();
	h_jetpt_uds->Sumw2();
	h_jetpt_b->Sumw2();
	h_jetpt_c->Sumw2();
	h_jetpt_ghost->Sumw2();
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
	h_muPt_ghost->Sumw2();
	h_muPt_noCut->Sumw2();
h_muEta->Sumw2();
	h_muEta_g->Sumw2();
	h_muEta_uds->Sumw2();
	h_muEta_b->Sumw2();
	h_muEta_c->Sumw2();
	h_muEta_ghost->Sumw2();
h_muPhi->Sumw2();
	h_muPhi_g->Sumw2();
	h_muPhi_uds->Sumw2();
	h_muPhi_b->Sumw2();
	h_muPhi_c->Sumw2();
	h_muPhi_ghost->Sumw2();
h_muRelPt->Sumw2();
	h_muRelPt_g->Sumw2();
	h_muRelPt_uds->Sumw2();
	h_muRelPt_b->Sumw2();
	h_muRelPt_c->Sumw2();
	h_muRelPt_ghost->Sumw2();
h_pthat->Sumw2();
h_vz->Sumw2();

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
		//if(em->pthat<50.0){continue;}

		h_vz->Fill(em->vz,w_vz);
		h_pthat->Fill(em->pthat);
		//cout << "vz = " << em->vz << endl;
		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
                	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/NEvents;
		//cout << "Number of jets = " << em->njet << endl;
        //jet cuts
        double jetEtaCut = 1.5;
	double jetPtCut = 50.0;
		//reco jet loop
		for(int jeti = 0; jeti < em->njet; jeti++){
			//cout << "reco jet pt_" << jeti << " = " << em->jetpt[jeti] <<"  ////  " << "gen jet pt_"<<jeti<<" = "<<em->genjetpt[jeti]<<endl;
			double w_pt = 1.0/(pt_fit_fxn->Eval(em->jetpt[jeti]));
			double w_tot = w_vz*w_pt*em->weight;	

			if(fabs(em->jeteta[jeti])>jetEtaCut){continue;}
			if(em->jetpt[jeti]<jetPtCut){continue;}
			
			h_jetpt_raw->Fill(em->jetpt[jeti],w_vz*em->weight);						
			h_jetpt->Fill(em->jetpt[jeti],w_tot);
			h_jeteta->Fill(em->jeteta[jeti],w_tot);
			h_jetphi->Fill(em->jetphi[jeti],w_tot);
			
			
			if(em->hadronFlavor[jeti]==4){
					h_jetpt_c->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_c->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_c->Fill(em->jetphi[jeti],w_tot);
			}
			else if(em->hadronFlavor[jeti]==5){
					h_jetpt_b->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_b->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_b->Fill(em->jetphi[jeti],w_tot);
			}
			else{
				if(em->partonFlavor[jeti]==1 || em->partonFlavor[jeti]==-1){
					h_jetpt_uds->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_uds->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_uds->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==2 || em->partonFlavor[jeti]==-2){
					h_jetpt_uds->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_uds->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_uds->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==3 || em->partonFlavor[jeti]==-3){
					h_jetpt_uds->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_uds->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_uds->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==4 || em->partonFlavor[jeti]==-4){
					h_jetpt_c->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_c->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_c->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==5 || em->partonFlavor[jeti]==-5){
					h_jetpt_b->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_b->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_b->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==21){
					h_jetpt_g->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_g->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_g->Fill(em->jetphi[jeti],w_tot);
				}
				if(em->partonFlavor[jeti]==0){
					h_jetpt_ghost->Fill(em->jetpt[jeti],w_tot);
					h_jeteta_ghost->Fill(em->jeteta[jeti],w_tot);
					h_jetphi_ghost->Fill(em->jetphi[jeti],w_tot);
				}
			}
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
				h_deltaR->Fill(deltaRmin,w_tot);
				h_muPt_noCut->Fill(muPtMin,w_tot);
				if(muRelPtMin==999.0 || muPtMin==999.0 || muEtaMin==-999.0 || muPhiMin==-999.0){continue;}
				if(deltaRmin<0.4){
					h_muRelPt->Fill(muRelPtMin,w_tot);
					if(em->hadronFlavor[jeti]==4){
						h_muRelPt_c->Fill(muRelPtMin,w_tot);
					}
					else if(em->hadronFlavor[jeti]==5){
						h_muRelPt_b->Fill(muRelPtMin,w_tot);
					}
					else{
						if(em->partonFlavor[jeti]==1 || em->partonFlavor[jeti]==-1){h_muRelPt_uds->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==2 || em->partonFlavor[jeti]==-2){h_muRelPt_uds->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==3 || em->partonFlavor[jeti]==-3){h_muRelPt_uds->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==4 || em->partonFlavor[jeti]==-4){h_muRelPt_c->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==5 || em->partonFlavor[jeti]==-5){h_muRelPt_b->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==21){h_muRelPt_g->Fill(muRelPtMin,w_tot);}
						if(em->partonFlavor[jeti]==0){h_muRelPt_ghost->Fill(muRelPtMin,w_tot);}
					}
				}
					h_muPt->Fill(muPtMin,w_tot);
					h_muEta->Fill(muEtaMin,w_tot);
					h_muPhi->Fill(muPhiMin,w_tot);
					if(em->hadronFlavor[jeti]==4){
						h_muPt_c->Fill(muPtMin,w_tot);
						h_muEta_c->Fill(muEtaMin,w_tot);
						h_muPhi_c->Fill(muPhiMin,w_tot);
					}
					else if(em->hadronFlavor[jeti]==5){	
						h_muPt_b->Fill(muPtMin,w_tot);
						h_muEta_b->Fill(muEtaMin,w_tot);
						h_muPhi_b->Fill(muPhiMin,w_tot);
					}
					else{
						if(em->partonFlavor[jeti]==1 || em->partonFlavor[jeti]==-1){
							h_muPt_uds->Fill(muPtMin,w_tot);
							h_muEta_uds->Fill(muEtaMin,w_tot);
							h_muPhi_uds->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==2 || em->partonFlavor[jeti]==-2){
							h_muPt_uds->Fill(muPtMin,w_tot);
							h_muEta_uds->Fill(muEtaMin,w_tot);
							h_muPhi_uds->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==3 || em->partonFlavor[jeti]==-3){
							h_muPt_uds->Fill(muPtMin,w_tot);
							h_muEta_uds->Fill(muEtaMin,w_tot);
							h_muPhi_uds->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==4 || em->partonFlavor[jeti]==-4){
							h_muPt_c->Fill(muPtMin,w_tot);
							h_muEta_c->Fill(muEtaMin,w_tot);
							h_muPhi_c->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==5 || em->partonFlavor[jeti]==-5){
							h_muPt_b->Fill(muPtMin,w_tot);
							h_muEta_b->Fill(muEtaMin,w_tot);
							h_muPhi_b->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==21){
							h_muPt_g->Fill(muPtMin,w_tot);
							h_muEta_g->Fill(muEtaMin,w_tot);
							h_muPhi_g->Fill(muPhiMin,w_tot);
						}
						if(em->partonFlavor[jeti]==0){
							h_muPt_ghost->Fill(muPtMin,w_tot);
							h_muEta_ghost->Fill(muEtaMin,w_tot);
							h_muPhi_ghost->Fill(muPhiMin,w_tot);
						}
					}
				


		} // end reco jet loop
		// gen jet loop
		for(int jetj = 0; jetj < em->ngj; jetj++){
			//cout << "reco jet pt_" << jeti << " = " << em->jetpt[jeti] <<"  ////  " << "gen jet pt_"<<jeti<<" = "<<em->genjetpt[jeti]<<endl;
			
			if(fabs(em->genjeteta[jetj])>jetEtaCut){continue;}
			if(em->genjetpt[jetj]<jetPtCut){continue;}
			
			//float ref_multiplier = 0.35;
			//if(em->pthat < ref_multiplier*em->ref_jetpt[jeti]){continue;}
			
			h_genjetpt->Fill(em->genjetpt[jetj],em->weight);
			h_genjeteta->Fill(em->genjeteta[jetj],em->weight);
			h_genjetphi->Fill(em->genjetphi[jetj],em->weight);


		} // end gen jet loop

	} // end event loop

delete f;
} // end file loop

auto wf = TFile::Open("/uscms/home/cmbennet/work/pythia_skim_newWeights_pthatcut_30_muptcut_10_16Aug2020.root","recreate");

h_jetpt_raw->Write();
h_jetpt->Write();
	h_jetpt_g->Write();
	h_jetpt_uds->Write();
	h_jetpt_b->Write();
	h_jetpt_c->Write();
	h_jetpt_ghost->Write();
h_jeteta->Write();
	h_jeteta_g->Write();
	h_jeteta_uds->Write();
	h_jeteta_b->Write();
	h_jeteta_c->Write();
	h_jeteta_ghost->Write();
h_jetphi->Write();
	h_jetphi_g->Write();
	h_jetphi_uds->Write();
	h_jetphi_b->Write();
	h_jetphi_c->Write();
	h_jetphi_ghost->Write();
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
	h_muPt_ghost->Write();
	h_muPt_noCut->Write();
h_muEta->Write();
	h_muEta_g->Write();
	h_muEta_uds->Write();
	h_muEta_b->Write();
	h_muEta_c->Write();
	h_muEta_ghost->Write();
h_muPhi->Write();
	h_muPhi_g->Write();
	h_muPhi_uds->Write();
	h_muPhi_b->Write();
	h_muPhi_c->Write();
	h_muPhi_ghost->Write();
h_muRelPt->Write();
	h_muRelPt_g->Write();
	h_muRelPt_uds->Write();
	h_muRelPt_b->Write();
	h_muRelPt_c->Write();
	h_muRelPt_ghost->Write();
h_vz->Write();
h_pthat->Write();
h_deltaR->Write();

wf->Close();

return;
} // end program

