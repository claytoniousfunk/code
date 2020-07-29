
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


#include "myProcesses/hiforest/plugin/eventMap_hiForest.h"
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



void makeDataV2(double muPtCut = 5.0, double jetptcut = 50.){

	//////////////////////////////////////////////////////////////////////// jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_jtphi=0, *pf_partonFlavor=0, *pf_hadronFlavor=0; 

	//////////////////////////////////////////////////////////////////////// muon jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *muPt=0, *muEta=0, *muPhi=0, *muChi2NDF=0, *muInnerD0=0, *muInnerDz=0;
	int nMu, hiBin;
	vector <int> *muIsGlobal=0, *muIsTracker=0, *muMuonHits=0, *muStations=0, *muTrkLayers=0, *muPixelHits=0;

	//////////////////////////////////////////////////////////////////////// event variables ////////////////////////////////////////////////////////////////////////

	float vz;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////// DEFINE BINNING ///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const int numUsedEtaBins = 15;
	float usedEtaEdges[numUsedEtaBins+1] = {-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
	const int numUsedPtBins = 15;
	float usedPtEdges[numUsedPtBins+1] = {50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////// DEFINE HISTOGRAMS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	

		////////////////////////////////////////// /////////////////////////////  JETS  ///////////////////////////////////////////////////////////////////////////////////

		TH2D *h2 = new TH2D("h2","All jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		
		
		/////////////////////////////////////////////////////////////////   MUON-TAGGED JETS  //////////////////////////////////////////////////////////////////////////////

		TH2D *h2_MJ = new TH2D("h2_MJ","All muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		// parse jet by species of particles
		
	

		///////////////////////////////////////////////////////////////  CALCULATED VARIABLES  //////////////////////////////////////////////////////////////////////////////

		TH1D *deltaR = new TH1D("deltaR","#Delta r",10,0,1);
		TH1D *h_muRelPt = new TH1D("h_muRelPt","p_{T}^{#mu}",50,0,5); // all muon-jets
		

		///////////////////////////////////////////////////////////////  CALCULATED 2D HISTOGRAMS  ////////////////////////////////////////////////////////////////////////////
		
		const int numUsedRelPtBins = 50;
		float usedRelPtEdges[numUsedRelPtBins+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
			3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0};

		
		TH2D *h2_relPtJetPt = new TH2D("h2_relPtJetPt","All jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		
		


		h2->Sumw2();
		h2_MJ->Sumw2();
		deltaR -> Sumw2();
		h_muRelPt -> Sumw2();
		h2_relPtJetPt->Sumw2();
		
		
	
	////////////////////////////////////////////////////////////////////////  CUT Values  //////////////////////////////////////////////////////////////////////////////

	const double etamaxcut = 1.5; // default = 1.5
	const double pTmaxcut = 500;
	
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
	auto f = TFile::Open("/home/clayton/Analysis/data/PbPb_5TeV_SingleMuPD_MuJetTrigger_SuperSlimSkim_21Aug2019.root");
	//auto f = TFile::Open("/home/clayton/Analysis/data/PPData_5TeV_SingleMuPD_MuJetTrigger_SuperSlimSkim_27Aug2019.root");

    TTree *inp_tree = (TTree*)f->Get("mixing_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();



	////////////////////////////////////////////////////////////////////// RECO JETS ////////////////////////////////////////////////////////////////////
	
	
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_jtphi",&pf_jtphi);
	
	
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
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	double w=0;
	
	int nMuCount=0;

	////////////////////////////////////////////////////////////////////  Event loop  //////////////////////////////////////////////////////////////////////////

	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;
        
        

        //////////////////////////////////////////////  EVENT CUTS ////////////////////////////////////////////////

        

        
        
         
		//////////////////////////////////////////////////////////////////  Jet loop ////////////////////////////////////////////////////////////////////////////
		
		for(int jetj=0; jetj < (int)pf_jteta->size(); jetj++){
			
			double jetPtj = pf_jtpt->at(jetj);
			double jetEtaj = pf_jteta->at(jetj);
			double jetPhij = pf_jtphi->at(jetj);

			// convert eta to theta
			double jetThetaj = 2*atan(exp(-jetEtaj));
			
			double x = pf_jteta->at(jetj);
			double u = pf_jtpt->at(jetj);
			
			/////////////////////////////////////  JET CUTS /////////////////////////////////////////
			if(fabs(x)>etamaxcut || u < jetptcut || u>pTmaxcut || u==-999){continue;}


			h2->Fill(x,u,w);
			
			//////////////////////////////////////////////////////////////////  Muon loop ////////////////////////////////////////////////////////////////////////////
			nMu = muPt->size();
			if(nMu==0){continue;}
			double deltaRmin=100;
			double muRelPtMin = 0;
			for(int mui=0; mui<(int) muPt->size();mui++){
				if(muIsTracker->at(mui)==0 || TMath::Abs(muEta->at(mui))>2.4 || muChi2NDF->at(mui)==-99 || muChi2NDF->at(mui)>10
					||TMath::Abs(muInnerD0->at(mui))>0.2 || TMath::Abs(muInnerDz->at(mui))>0.5 || muMuonHits->at(mui)<= 0
					|| muStations->at(mui)<= 1 || muTrkLayers->at(mui)<=5 || muPixelHits->at(mui)<=0 || muPt->at(mui) < muPtCut){continue;}
			
				double muPti = muPt->at(mui);
				double muEtai = muEta->at(mui);
				double muPhii = muPhi->at(mui);

			
				double deltaEtaij = muEtai-jetEtaj;
				double deltaPhiij = acos(cos(muPhii-jetPhij));
				double deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				double muRelPt = getPtRel(muPti,muEtai,muPhii,jetPtj,jetEtaj,jetPhij);
				
				if(deltaRij<deltaRmin){
					deltaRmin=deltaRij;
					muRelPtMin=muRelPt;
					
				}

			} 
			//////////////////////////////////////////////////////////////////  End muon loop ////////////////////////////////////////////////////////////////////////////
	
			deltaR->Fill(deltaRmin);
			
			if(deltaRmin<0.4){

				h2_MJ->Fill(x,u,w);
				h_muRelPt->Fill(muRelPtMin);
				h2_relPtJetPt->Fill(muRelPtMin,u,w);

			}
			
		} 
		//////////////////////////////////////////////////////////////////  End jet loop ////////////////////////////////////////////////////////////////////////////
		
	} 
	//////////////////////////////////////////////////////////////////  End event loop ////////////////////////////////////////////////////////////////////////////

//}
/////////////////////////////////////////////////////////////////////  End file loop /////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////  file creation history /////////////////////////////////////////////////////

//auto wf = TFile::Open("makeDataV2_PbPb_refjets_pthat_50_muptcut_5.root","recreate");
auto wf = TFile::Open("makeDataV2_pp_refjets_pthat_50_muptcut_5.root","recreate");










h2->Write();

h2_MJ->Write();

deltaR->Write();
h_muRelPt->Write();
h2_relPtJetPt->Write();

wf->Close();

return;

}


 // end program


