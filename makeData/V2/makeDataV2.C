
// "etaPtHisto"
// V4: using variable binning
// V5: muon jet filtering, using HydJet_x.root files
// V6: using new weigthing method, pthat.  No longer reweighting according to bin size (delegate to plotting macros).  Including muRelPt calculation
// V7: desinged to run /skims5/ files. Added histogram of untagged jets.
// V8: changed variable binning
// V9: implementing use of hadronFlavor into particle identification.  changed pT max to 500.
// V10 : getting rid of untagged variable

/// "makeData" 
// V1 : adopting from etaPtHisto series.  Got rid of any "_t" variables.  Cleaned up some unused variables.  Most importantly, added vz & hiBin cut
// V2: to be used on data files, no flavor tagging


#include "/home/clayton/Analysis/code/myProcesses/hiforest/plugin/eventMap_hiForest.h"
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




double getPtRel(double MuonPt, double MuonEta, double MuonPhi, double JetPt, double JetEta, double JetPhi)
{

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



void makeDataV2(double muPtCut = 10.0){

	//////////////////////////////////////////////////////////////////////// jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *jetpt=0, *jeteta=0, *jetphi=0, *partonFlavor=0, *hadronFlavor=0; 

	//////////////////////////////////////////////////////////////////////// muon jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *muPt=0, *muEta=0, *muPhi=0, *muChi2NDF=0, *muInnerD0=0, *muInnerDz=0;
	int nMu, hiBin;
	vector <int> *muIsGlobal=0, *muIsTracker=0, *muMuonHits=0, *muStations=0, *muTrkLayers=0, *muPixelHits=0;

	//////////////////////////////////////////////////////////////////////// event variables ////////////////////////////////////////////////////////////////////////

	float vz;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////// DEFINE BINNING ///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////// DEFINE HISTOGRAMS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int NEtaBins = 300;
	double etaMin = -1.5;
	double etaMax = 1.5;
	int NPtBins = 450;
	double ptMin = 50.0;
	double ptMax = 500.0;
	int NPhiBins = 300;
	double phiMin = -TMath::Pi();
	double phiMax = TMath::Pi();
	
		////////////////////////////////////////// /////////////////////////////  Event variables ///////////////////////////////////////////////////////////////////////////////////

		TH1D *h_vz = new TH1D("h_vz","data vz",60,-15.0,15.0);
		TH1D *h_hiBin = new TH1D("h_hiBin","data hiBin",200,0,200);

		////////////////////////////////////////// /////////////////////////////  JETS  ///////////////////////////////////////////////////////////////////////////////////

		TH1D *h_jetpt = new TH1D("h_jetpt","data jet pt",NPtBins,ptMin,ptMax);
			// centrality regions for PbPb
			TH1D *h_jetpt_cent0to10 = new TH1D("h_jetpt_cent0to10","centrality 0-10 %",NPtBins,ptMin,ptMax);
			TH1D *h_jetpt_cent10to30 = new TH1D("h_jetpt_cent10to30","centrality 10-30 %",NPtBins,ptMin,ptMax);
			TH1D *h_jetpt_cent30to50 = new TH1D("h_jetpt_cent30to50","centrality 30-50 %",NPtBins,ptMin,ptMax);
			TH1D *h_jetpt_cent50to90 = new TH1D("h_jetpt_cent50to90","centrality 50-90 %",NPtBins,ptMin,ptMax);
		TH1D *h_jeteta = new TH1D("h_jeteta","data jet #eta",NEtaBins,etaMin,etaMax);
			TH1D *h_jeteta_cent0to10 = new TH1D("h_jeteta_cent0to10","centrality 0-10 %",NEtaBins,etaMin,etaMax);
			TH1D *h_jeteta_cent10to30 = new TH1D("h_jeteta_cent10to30","centrality 10-30 %",NEtaBins,etaMin,etaMax);
			TH1D *h_jeteta_cent30to50 = new TH1D("h_jeteta_cent30to50","centrality 30-50 %",NEtaBins,etaMin,etaMax);
			TH1D *h_jeteta_cent50to90 = new TH1D("h_jeteta_cent50to90","centrality 50-90 %",NEtaBins,etaMin,etaMax);
		TH1D* h_jetphi = new TH1D("h_jetphi","data jet #phi",NPhiBins,phiMin,phiMax);
			TH1D *h_jetphi_cent0to10 = new TH1D("h_jetphi_cent0to10","centrality 0-10 %",NPhiBins,phiMin,phiMax);
			TH1D *h_jetphi_cent10to30 = new TH1D("h_jetphi_cent10to30","centrality 10-30 %",NPhiBins,phiMin,phiMax);
			TH1D *h_jetphi_cent30to50 = new TH1D("h_jetphi_cent30to50","centrality 30-50 %",NPhiBins,phiMin,phiMax);
			TH1D *h_jetphi_cent50to90 = new TH1D("h_jetphi_cent50to90","centrality 50-90 %",NPhiBins,phiMin,phiMax);
		
		TH1D *h_mujetpt = new TH1D("h_mujetpt","data muon-jet pt",NPtBins,ptMin,ptMax);
			// centrality regions for PbPb
			TH1D *h_mujetpt_cent0to10 = new TH1D("h_mujetpt_cent0to10","centrality 0-10 %",NPtBins,ptMin,ptMax);
			TH1D *h_mujetpt_cent10to30 = new TH1D("h_mujetpt_cent10to30","centrality 10-30 %",NPtBins,ptMin,ptMax);
			TH1D *h_mujetpt_cent30to50 = new TH1D("h_mujetpt_cent30to50","centrality 30-50 %",NPtBins,ptMin,ptMax);
			TH1D *h_mujetpt_cent50to90 = new TH1D("h_mujetpt_cent50to90","centrality 50-90 %",NPtBins,ptMin,ptMax);
		TH1D *h_mujeteta = new TH1D("h_mujeteta","data jet #eta",NEtaBins,etaMin,etaMax);
			TH1D *h_mujeteta_cent0to10 = new TH1D("h_mujeteta_cent0to10","centrality 0-10 %",NEtaBins,etaMin,etaMax);
			TH1D *h_mujeteta_cent10to30 = new TH1D("h_mujeteta_cent10to30","centrality 10-30 %",NEtaBins,etaMin,etaMax);
			TH1D *h_mujeteta_cent30to50 = new TH1D("h_mujeteta_cent30to50","centrality 30-50 %",NEtaBins,etaMin,etaMax);
			TH1D *h_mujeteta_cent50to90 = new TH1D("h_mujeteta_cent50to90","centrality 50-90 %",NEtaBins,etaMin,etaMax);
		TH1D* h_mujetphi = new TH1D("h_mujetphi","data jet #phi",NPhiBins,phiMin,phiMax);
			TH1D *h_mujetphi_cent0to10 = new TH1D("h_mujetphi_cent0to10","centrality 0-10 %",NPhiBins,phiMin,phiMax);
			TH1D *h_mujetphi_cent10to30 = new TH1D("h_mujetphi_cent10to30","centrality 10-30 %",NPhiBins,phiMin,phiMax);
			TH1D *h_mujetphi_cent30to50 = new TH1D("h_mujetphi_cent30to50","centrality 30-50 %",NPhiBins,phiMin,phiMax);
			TH1D *h_mujetphi_cent50to90 = new TH1D("h_mujetphi_cent50to90","centrality 50-90 %",NPhiBins,phiMin,phiMax);
		
		
		
		/////////////////////////////////////////////////////////////////   MUON-TAGGED JETS  //////////////////////////////////////////////////////////////////////////////
	int NmuPtBins = 5000;
	double muPtLow = 0.0;
	double muPtHigh = 500.0;
		TH1D *h_muPt = new TH1D("h_muPt","data muon pt",NmuPtBins,muPtLow,muPtHigh);
			TH1D *h_muPt_cent0to10 = new TH1D("h_muPt_cent0to10","centrality 0-10 %",NmuPtBins,muPtLow,muPtHigh);
			TH1D *h_muPt_cent10to30 = new TH1D("h_muPt_cent10to30","centrality 10-30 %",NmuPtBins,muPtLow,muPtHigh);
			TH1D *h_muPt_cent30to50 = new TH1D("h_muPt_cent30to50","centrality 30-50 %",NmuPtBins,muPtLow,muPtHigh);
			TH1D *h_muPt_cent50to90 = new TH1D("h_muPt_cent50to90","centrality 50-90 %",NmuPtBins,muPtLow,muPtHigh);
		TH1D *h_muEta = new TH1D("h_muEta","data muon #eta",NEtaBins,etaMin,etaMax);
			TH1D *h_muEta_cent0to10 = new TH1D("h_muEta_cent0to10","centrality 0-10 %",NEtaBins,etaMin,etaMax);
			TH1D *h_muEta_cent10to30 = new TH1D("h_muEta_cent10to30","centrality 10-30 %",NEtaBins,etaMin,etaMax);
			TH1D *h_muEta_cent30to50 = new TH1D("h_muEta_cent30to50","centrality 30-50 %",NEtaBins,etaMin,etaMax);
			TH1D *h_muEta_cent50to90 = new TH1D("h_muEta_cent50to90","centrality 50-90 %",NEtaBins,etaMin,etaMax);
		TH1D *h_muPhi = new TH1D("h_muPhi","data muon #phi",NPhiBins,phiMin,phiMax);
			TH1D *h_muPhi_cent0to10 = new TH1D("h_muPhi_cent0to10","centrality 0-10 %",NPhiBins,phiMin,phiMax);
			TH1D *h_muPhi_cent10to30 = new TH1D("h_muPhi_cent10to30","centrality 10-30 %",NPhiBins,phiMin,phiMax);
			TH1D *h_muPhi_cent30to50 = new TH1D("h_muPhi_cent30to50","centrality 30-50 %",NPhiBins,phiMin,phiMax);
			TH1D *h_muPhi_cent50to90 = new TH1D("h_muPhi_cent50to90","centrality 50-90 %",NPhiBins,phiMin,phiMax);

		///////////////////////////////////////////////////////////////  CALCULATED VARIABLES  //////////////////////////////////////////////////////////////////////////////

		TH1D *deltaR = new TH1D("deltaR","#Delta r",200,0.0,2.0);
		TH1D *h_muRelPt = new TH1D("h_muRelPt","muon rel pt",50,0.0,5.0); // all muon-jets
			TH1D *h_muRelPt_cent0to10 = new TH1D("h_muRelPt_cent0to10","muon rel pt, centrality 0-10%",50,0.0,5.0);
			TH1D *h_muRelPt_cent10to30 = new TH1D("h_muRelPt_cent10to30","muon rel pt, centrality 10-30%",50,0.0,5.0);
			TH1D *h_muRelPt_cent30to50 = new TH1D("h_muRelPt_cent30to50","muon rel pt, centrality 30-50%",50,0.0,5.0);
			TH1D *h_muRelPt_cent50to90 = new TH1D("h_muRelPt_cent50to90","muon rel pt, centrality 50-90%",50,0.0,5.0);
		

		///////////////////////////////////////////////////////////////  CALCULATED 2D HISTOGRAMS  ////////////////////////////////////////////////////////////////////////////
		
		
		// Sumw2 commands

		h_vz->Sumw2();
		h_hiBin->Sumw2();
		h_jetpt->Sumw2();
			h_jetpt_cent0to10->Sumw2();
			h_jetpt_cent10to30->Sumw2();
			h_jetpt_cent30to50->Sumw2();
			h_jetpt_cent50to90->Sumw2();
		h_jeteta->Sumw2();
			h_jeteta_cent0to10->Sumw2();
			h_jeteta_cent10to30->Sumw2();
			h_jeteta_cent30to50->Sumw2();
			h_jeteta_cent50to90->Sumw2();
		h_jetphi->Sumw2();
			h_jetphi_cent0to10->Sumw2();
			h_jetphi_cent10to30->Sumw2();
			h_jetphi_cent30to50->Sumw2();
			h_jetphi_cent50to90->Sumw2();
		h_mujetpt->Sumw2();
			h_mujetpt_cent0to10->Sumw2();
			h_mujetpt_cent10to30->Sumw2();
			h_mujetpt_cent30to50->Sumw2();
			h_mujetpt_cent50to90->Sumw2();
		h_mujeteta->Sumw2();
			h_mujeteta_cent0to10->Sumw2();
			h_mujeteta_cent10to30->Sumw2();
			h_mujeteta_cent30to50->Sumw2();
			h_mujeteta_cent50to90->Sumw2();
		h_mujetphi->Sumw2();
			h_mujetphi_cent0to10->Sumw2();
			h_mujetphi_cent10to30->Sumw2();
			h_mujetphi_cent30to50->Sumw2();
			h_mujetphi_cent50to90->Sumw2();
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
		deltaR->Sumw2();
		h_muRelPt->Sumw2();
			h_muRelPt_cent0to10->Sumw2();
			h_muRelPt_cent10to30->Sumw2();
			h_muRelPt_cent30to50->Sumw2();
			h_muRelPt_cent50to90->Sumw2();
		
	
	////////////////////////////////////////////////////////////////////////  CUT Values  //////////////////////////////////////////////////////////////////////////////

	
	
	/////////////////////////////////////////////////////////////////////  LOAD DATA  /////////////////////////////////////////////////////////////////////////////
	/*
	int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/skims5/0000");

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

/////////////////////////////////////////////////////////////////////  Begin file loop /////////////////////////////////////////////////////////////////////////////

//for(int file = 1; file < NFiles+1; file++){
//for(int file = 1; file < 20; file++){
	
	//cout << "Processing file " << file << "/"<<NFiles<< endl;
	TFile *f = TFile::Open("/home/clayton/Analysis/data/PbPb_5TeV_SingleMuPD_MuJetTrigger_SuperSlimSkim_21Aug2019.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/data/PPData_5TeV_SingleMuPD_MuJetTrigger_SuperSlimSkim_27Aug2019.root");

    TTree *inp_tree = (TTree*)f->Get("mixing_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();



	////////////////////////////////////////////////////////////////////// RECO JETS ////////////////////////////////////////////////////////////////////
	
	
	inp_tree->SetBranchAddress("pf_corrpt",&jetpt); 
	inp_tree->SetBranchAddress("pf_jteta",&jeteta);
	inp_tree->SetBranchAddress("pf_jtphi",&jetphi);
	
	
	////////////////////////////////////////////////////////////////////// MUON JETS ////////////////////////////////////////////////////////////////////

	inp_tree->SetBranchAddress("mu_pt",&muPt);
	inp_tree->SetBranchAddress("mu_eta",&muEta);
	inp_tree->SetBranchAddress("mu_phi",&muPhi);
	inp_tree->SetBranchAddress("mu_chi2ndf",&muChi2NDF);
	inp_tree->SetBranchAddress("mu_innerD0",&muInnerD0);
	inp_tree->SetBranchAddress("mu_innerDz",&muInnerDz);
	inp_tree->SetBranchAddress("mu_isGlobal",&muIsGlobal);
	inp_tree->SetBranchAddress("mu_isTracker",&muIsTracker);
	inp_tree->SetBranchAddress("mu_muonHits",&muMuonHits);
	inp_tree->SetBranchAddress("mu_stations",&muStations);
	inp_tree->SetBranchAddress("mu_trkLayers",&muTrkLayers);
	inp_tree->SetBranchAddress("mu_pixelHits",&muPixelHits);

	////////////////////////////////////////////////////////////////// EVENT VARIABLES ////////////////////////////////////////////////////////////////////

	inp_tree->SetBranchAddress("vz",&vz);
	inp_tree->SetBranchAddress("hiBin",&hiBin);

	//////////////////////////////////////////////////////////////////  loop variables  ///////////////////////////////////////////////////////////////////////	

	int evi = 0;
	int evi_frac = 0;

	////////////////////////////////////////////////////////////////////  Event loop  //////////////////////////////////////////////////////////////////////////

	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;
        
        
		
        //////////////////////////////////////////////  EVENT CUTS ////////////////////////////////////////////////

		if(fabs(vz)>15.0){continue;}
		if(hiBin > 180){continue;}
        
		h_vz->Fill(vz);
		h_hiBin->Fill(hiBin);
        
        
         
		//////////////////////////////////////////////////////////////////  Jet loop ////////////////////////////////////////////////////////////////////////////
		
		for(int jetj=0; jetj < (int)jeteta->size(); jetj++){
			
			double jetPtj = jetpt->at(jetj);
			double jetEtaj = jeteta->at(jetj);
			double jetPhij = jetphi->at(jetj);

					
			/////////////////////////////////////  JET CUTS /////////////////////////////////////////
			if(fabs(jetEtaj)>etaMax || jetPtj < ptMin || jetPtj>ptMax || jetPhij==-999){continue;}

			h_jetpt->Fill(jetPtj);
				if(hiBin>0 && hiBin<20){h_jetpt_cent0to10->Fill(jetPtj);}
				if(hiBin>20 && hiBin<60){h_jetpt_cent10to30->Fill(jetPtj);}
				if(hiBin>60 && hiBin<100){h_jetpt_cent30to50->Fill(jetPtj);}
				if(hiBin>100 && hiBin<180){h_jetpt_cent50to90->Fill(jetPtj);}
			h_jeteta->Fill(jetEtaj);
				if(hiBin>0 && hiBin<20){h_jeteta_cent0to10->Fill(jetEtaj);}
				if(hiBin>20 && hiBin<60){h_jeteta_cent10to30->Fill(jetEtaj);}
				if(hiBin>60 && hiBin<100){h_jeteta_cent30to50->Fill(jetEtaj);}
				if(hiBin>100 && hiBin<180){h_jeteta_cent50to90->Fill(jetEtaj);}
			h_jetphi->Fill(jetPhij);
				if(hiBin>0 && hiBin<20){h_jetphi_cent0to10->Fill(jetPhij);}
				if(hiBin>20 && hiBin<60){h_jetphi_cent10to30->Fill(jetPhij);}
				if(hiBin>60 && hiBin<100){h_jetphi_cent30to50->Fill(jetPhij);}
				if(hiBin>100 && hiBin<180){h_jetphi_cent50to90->Fill(jetPhij);}
			
			//////////////////////////////////////////////////////////////////  Muon loop ////////////////////////////////////////////////////////////////////////////
			nMu = muPt->size();
			if(nMu==0){continue;}
			double deltaRmin=1000000.0;
			double muRelPtMin = 1000000.0;
			double muJetPtMin = 1000000.0;
			double muJetPhiMin = 1000000.0;
			double muJetEtaMin = 1000000.0;
			for(int mui=0; mui<(int) muPt->size();mui++){
				if(muIsTracker->at(mui)==0 || TMath::Abs(muEta->at(mui))>2.4 || muChi2NDF->at(mui)==-99 || muChi2NDF->at(mui)>10
					||TMath::Abs(muInnerD0->at(mui))>0.2 || TMath::Abs(muInnerDz->at(mui))>0.5 || muMuonHits->at(mui)<= 0
					|| muStations->at(mui)<= 1 || muTrkLayers->at(mui)<=5 || muPixelHits->at(mui)<=0 || muPt->at(mui) < muPtCut){continue;}
			
				double muPti = muPt->at(mui);
				double muEtai = muEta->at(mui);
				double muPhii = muPhi->at(mui);

				h_muPt->Fill(muPti);
					if(hiBin>0 && hiBin<20){h_muPt_cent0to10->Fill(muPti);}
					if(hiBin>20 && hiBin<60){h_muPt_cent10to30->Fill(muPti);}
					if(hiBin>60 && hiBin<100){h_muPt_cent30to50->Fill(muPti);}
					if(hiBin>100 && hiBin<180){h_muPt_cent50to90->Fill(muPti);}
				h_muEta->Fill(muEtai);
					if(hiBin>0 && hiBin<20){h_muEta_cent0to10->Fill(muEtai);}
					if(hiBin>20 && hiBin<60){h_muEta_cent10to30->Fill(muEtai);}
					if(hiBin>60 && hiBin<100){h_muEta_cent30to50->Fill(muEtai);}
					if(hiBin>100 && hiBin<180){h_muEta_cent50to90->Fill(muEtai);}
				h_muPhi->Fill(muPhii);
					if(hiBin>0 && hiBin<20){h_muPhi_cent0to10->Fill(muPhii);}
					if(hiBin>20 && hiBin<60){h_muPhi_cent10to30->Fill(muPhii);}
					if(hiBin>60 && hiBin<100){h_muPhi_cent30to50->Fill(muPhii);}
					if(hiBin>100 && hiBin<180){h_muPhi_cent50to90->Fill(muPhii);}

			
				double deltaEtaij = muEtai-jetEtaj;
				double deltaPhiij = acos(cos(muPhii-jetPhij));
				double deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				double muRelPt = getPtRel(muPti,muEtai,muPhii,jetPtj,jetEtaj,jetPhij);
				
				if(deltaRij<deltaRmin){
					deltaRmin=deltaRij;
					muRelPtMin=muRelPt;
					muJetPtMin=jetPtj;
					muJetEtaMin=jetEtaj;
					muJetPhiMin=jetPhij;
					
				}

			} 
			//////////////////////////////////////////////////////////////////  End muon loop ////////////////////////////////////////////////////////////////////////////
	
			if(deltaRmin != 1000000.0){deltaR->Fill(deltaRmin);}
			if(muJetPtMin != 1000000.0){
				h_mujetpt->Fill(muJetPtMin);
				if(hiBin>0 && hiBin<20){h_mujetpt_cent0to10->Fill(muJetPtMin);}
				if(hiBin>20 && hiBin<60){h_mujetpt_cent10to30->Fill(muJetPtMin);}
				if(hiBin>60 && hiBin<100){h_mujetpt_cent30to50->Fill(muJetPtMin);}
				if(hiBin>100 && hiBin<180){h_mujetpt_cent50to90->Fill(muJetPtMin);}
			}
			if(muJetEtaMin != 1000000.0){
				h_mujeteta->Fill(muJetEtaMin);
				if(hiBin>0 && hiBin<20){h_mujeteta_cent0to10->Fill(muJetEtaMin);}
				if(hiBin>20 && hiBin<60){h_mujeteta_cent10to30->Fill(muJetEtaMin);}
				if(hiBin>60 && hiBin<100){h_mujeteta_cent30to50->Fill(muJetEtaMin);}
				if(hiBin>100 && hiBin<180){h_mujeteta_cent50to90->Fill(muJetEtaMin);}
			}
			if(muJetPhiMin != 1000000.0){
				h_mujetphi->Fill(muJetPhiMin);
				if(hiBin>0 && hiBin<20){h_mujetphi_cent0to10->Fill(muJetPhiMin);}
				if(hiBin>20 && hiBin<60){h_mujetphi_cent10to30->Fill(muJetPhiMin);}
				if(hiBin>60 && hiBin<100){h_mujetphi_cent30to50->Fill(muJetPhiMin);}
				if(hiBin>100 && hiBin<180){h_mujetphi_cent50to90->Fill(muJetPhiMin);}	
			}
			
			if(deltaRmin<0.4){
				
				h_muRelPt->Fill(muRelPtMin);
				if(hiBin>0 && hiBin<20){h_muRelPt_cent0to10->Fill(muRelPtMin);}
				if(hiBin>20 && hiBin<60){h_muRelPt_cent10to30->Fill(muRelPtMin);}
				if(hiBin>60 && hiBin<100){h_muRelPt_cent30to50->Fill(muRelPtMin);}
				if(hiBin>100 && hiBin<180){h_muRelPt_cent50to90->Fill(muRelPtMin);}
				
			}
			
		} 
		//////////////////////////////////////////////////////////////////  End jet loop ////////////////////////////////////////////////////////////////////////////
		
	} 
	//////////////////////////////////////////////////////////////////  End event loop ////////////////////////////////////////////////////////////////////////////

//}
/////////////////////////////////////////////////////////////////////  End file loop /////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////  file creation history /////////////////////////////////////////////////////

//auto wf = TFile::Open("makeDataV2_PbPb_refjets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV2_pp_refjets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_corrpt_muptcut_10.root","recreate");
auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_PbPb_data_corrpt_muptcut_10.root","recreate");


h_vz->Write();
h_hiBin->Write();
h_jetpt->Write();
	h_jetpt_cent0to10->Write();
	h_jetpt_cent10to30->Write();
	h_jetpt_cent30to50->Write();
	h_jetpt_cent50to90->Write();
h_jeteta->Write();
	h_jeteta_cent0to10->Write();
	h_jeteta_cent10to30->Write();
	h_jeteta_cent30to50->Write();
	h_jeteta_cent50to90->Write();
h_jetphi->Write();
	h_jetphi_cent0to10->Write();
	h_jetphi_cent10to30->Write();
	h_jetphi_cent30to50->Write();
	h_jetphi_cent50to90->Write();
h_mujetpt->Write();
	h_mujetpt_cent0to10->Write();
	h_mujetpt_cent10to30->Write();
	h_mujetpt_cent30to50->Write();
	h_mujetpt_cent50to90->Write();
h_mujeteta->Write();
	h_mujeteta_cent0to10->Write();
	h_mujeteta_cent10to30->Write();
	h_mujeteta_cent30to50->Write();
	h_mujeteta_cent50to90->Write();
h_mujetphi->Write();
	h_mujetphi_cent0to10->Write();
	h_mujetphi_cent10to30->Write();
	h_mujetphi_cent30to50->Write();
	h_mujetphi_cent50to90->Write();
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
deltaR->Write();
h_muRelPt->Write();
	h_muRelPt_cent0to10->Write();
	h_muRelPt_cent10to30->Write();
	h_muRelPt_cent30to50->Write();
	h_muRelPt_cent50to90->Write();


wf->Close();

return;

}


 // end program


