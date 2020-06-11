


/*
This version is made to compile skims made by Dhanush, with only relevant data.  No need for eventMap
V5 :  implementing input variables for pT cuts.  Added more histogram definitions (different combinations)
*/





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
#include <array>

//void jetParticleFraction(TString inf){
void jetPtFractionV5(float ptcut1=100, float ptcut2=600){

	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/skims3/PbPbMC2018_hadronFlavor_Apr18.root");
    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_partonFlavor=0; 
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);

	/// DEFINE BINNING ///
	
	
	const Int_t Nbins = 100;
	Float_t edges[Nbins+1];
	double delta = 10.0; // initial width of bin (pT) // 10 GeV up to 100, 20 GeV up to 250, 250 out = 50 GeV, 400-500 = one bin
	///double jtpt_bindouble1 = 350;
	//double jtpt_bindouble2 = 500;
	double jtpt_check=ptcut1;
	edges[0]=ptcut1;
	int numUsedBins = 0;
	for(int zz=1; jtpt_check<ptcut2; zz++){

		if(jtpt_check<100){
			delta=10.;	
			edges[zz]=edges[zz-1]+delta;
		}
		else if (jtpt_check<250){
			delta=20.;
			edges[zz]=edges[zz-1]+delta;
		}
		else if (jtpt_check<400){
			delta=50;
			edges[zz]=edges[zz-1]+delta;
			if(edges[zz]>=400){
				edges[zz]=400;
				edges[zz+1]=ptcut2;
				numUsedBins=zz+1;
				break;
			}
		}
		
		
		
		jtpt_check=jtpt_check+delta;
		numUsedBins=zz+1;
		
	}
	

	Float_t usedEdges[numUsedBins+1];
	for(int jj=0; jj<numUsedBins+1; jj++){
		usedEdges[jj]=edges[jj];
	}
	
	
		/// DEFINE HISTOGRAMS ///
	
	TH1D *h_jetpt = new TH1D("h_jetpt","Jet pT",numUsedBins,usedEdges);
		// parse jetpt by species of particles
		TH1D *h_jetpt_d = new TH1D("h_jetpt_d","",numUsedBins,usedEdges); // down quark
      	TH1D *h_jetpt_u = new TH1D("h_jetpt_u","",numUsedBins,usedEdges); // up quark
      	TH1D *h_jetpt_s = new TH1D("h_jetpt_s","",numUsedBins,usedEdges); // strange quark
      	TH1D *h_jetpt_c = new TH1D("h_jetpt_c","",numUsedBins,usedEdges); // charm quark
      	TH1D *h_jetpt_b = new TH1D("h_jetpt_b","",numUsedBins,usedEdges); // bottom quark
      	TH1D *h_jetpt_dbar = new TH1D("h_jetpt_dbar","",numUsedBins,usedEdges); // down antiquark
      	TH1D *h_jetpt_ubar = new TH1D("h_jetpt_ubar","",numUsedBins,usedEdges); // up antiquark
      	TH1D *h_jetpt_sbar = new TH1D("h_jetpt_sbar","",numUsedBins,usedEdges); // strange antiquark
      	TH1D *h_jetpt_cbar = new TH1D("h_jetpt_cbar","",numUsedBins,usedEdges); // charm antiquark
      	TH1D *h_jetpt_bbar = new TH1D("h_jetpt_bbar","",numUsedBins,usedEdges); // bottom antiquark
      	TH1D *h_jetpt_g = new TH1D("h_jetpt_g","",numUsedBins,usedEdges); // gluon
      	TH1D *h_jetpt_hq = new TH1D("h_jetpt_hq","",numUsedBins,usedEdges); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
      	TH1D *h_jetpt_lq = new TH1D("h_jetpt_lq","",numUsedBins,usedEdges); // all light quarks (u,d,...)
      	TH1D *h_jetpt_qall = new TH1D("h_jetpt_qall","",numUsedBins,usedEdges); // all quarks
      	TH1D *h_jetpt_dall = new TH1D("h_jetpt_dall","",numUsedBins,usedEdges); // d and dbar combined
      	TH1D *h_jetpt_uall = new TH1D("h_jetpt_uall","",numUsedBins,usedEdges); // u and ubar combined
      	TH1D *h_jetpt_sall = new TH1D("h_jetpt_sall","",numUsedBins,usedEdges); // s and sbar combined
      	TH1D *h_jetpt_call = new TH1D("h_jetpt_call","",numUsedBins,usedEdges); // c and cbar combined
      	TH1D *h_jetpt_ball = new TH1D("h_jetpt_ball","",numUsedBins,usedEdges); // b and bbar combined
      	TH1D *h_jetpt_ud = new TH1D("h_jetpt_ud","",numUsedBins,usedEdges); // u and d combined
      	TH1D *h_jetpt_ubardbar = new TH1D("h_jetpt_ubardbar","",numUsedBins,usedEdges); // u and d combined
		TH1D *h_jetpt_other = new TH1D("h_jetpt_other","",numUsedBins,usedEdges); // others
	
	TH1D *h_flavor_forb = new TH1D("h_flavor_forb","",35,-10,25);

	h_jetpt->Sumw2();
	h_jetpt_d->Sumw2();
	h_jetpt_u->Sumw2();
	h_jetpt_s->Sumw2();
	h_jetpt_c->Sumw2();
	h_jetpt_b->Sumw2();
	h_jetpt_dbar->Sumw2();
	h_jetpt_ubar->Sumw2();
	h_jetpt_sbar->Sumw2();
	h_jetpt_cbar->Sumw2();
	h_jetpt_bbar->Sumw2();
	h_jetpt_hq->Sumw2();
	h_jetpt_lq->Sumw2();
	h_jetpt_qall->Sumw2();
	h_jetpt_dall->Sumw2();
	h_jetpt_uall->Sumw2();
	h_jetpt_sall->Sumw2();
	h_jetpt_call->Sumw2();
	h_jetpt_ball->Sumw2();
	h_jetpt_g->Sumw2();
	h_flavor_forb->Sumw2();
	h_jetpt_ud->Sumw2();
	h_jetpt_ubardbar->Sumw2();
	h_jetpt_other -> Sumw2();

	//// CUTS ///////
	const double etamaxcut = 1.5;
	const double etamincut = 0.;
	const double pTmincut = ptcut1; // default = 120.;
	const double pTmaxcut = ptcut2;
	
	// loop variables //	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	//double x = 0;
	//int y = 0;
	
	// Event loop
	for (evi = 0; evi < n_evts; evi++){
	//for (evi = 0; evi < 10; evi++){	
        	inp_tree->GetEntry(evi);
                if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/n_evts;
                                                                                                   

		// Jet loop
		
		for(int jeti = 0; jeti < (int) pf_jteta->size(); jeti++){
			
			double x = pf_jteta->at(jeti);
			double u = pf_jtpt->at(jeti);
			int y = pf_partonFlavor->at(jeti);
			if(fabs(x)>=etamaxcut || u <= pTmincut){continue;}
			h_flavor_forb -> Fill(y);
			h_jetpt -> Fill(u);
			if(y==1){h_jetpt_d->Fill(u);h_jetpt_dall->Fill(u);h_jetpt_lq->Fill(u);h_jetpt_qall->Fill(u);h_jetpt_ud->Fill(u);}
			if(y==2){h_jetpt_u->Fill(u);h_jetpt_uall->Fill(u);h_jetpt_lq->Fill(u);h_jetpt_qall->Fill(u);h_jetpt_ud->Fill(u);}
			if(y==3){h_jetpt_s->Fill(u); h_jetpt_sall->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==4){h_jetpt_c->Fill(u); h_jetpt_call->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==5){h_jetpt_b->Fill(u); h_jetpt_ball->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==-1){h_jetpt_dbar->Fill(u);h_jetpt_dall->Fill(u);h_jetpt_lq->Fill(u);h_jetpt_qall->Fill(u);h_jetpt_ubardbar->Fill(u);}
			if(y==-2){h_jetpt_ubar->Fill(u);h_jetpt_uall->Fill(u);h_jetpt_lq->Fill(u);h_jetpt_qall->Fill(u);h_jetpt_ubardbar->Fill(u);}
			if(y==-3){h_jetpt_sbar->Fill(u); h_jetpt_sall->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==-4){h_jetpt_cbar->Fill(u); h_jetpt_call->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==-5){h_jetpt_bbar->Fill(u); h_jetpt_ball->Fill(u);h_jetpt_hq->Fill(u);h_jetpt_qall->Fill(u);}
			if(y==21){h_jetpt_g->Fill(u);}
			if(y==0){h_jetpt_other -> Fill(u);}
		}

	}
	
	

    
	// Write data	
auto wf = TFile::Open("jetPtFractionV5.root", "recreate");
h_jetpt->Write();
h_jetpt_d->Write();
h_jetpt_u->Write();
h_jetpt_s->Write();
h_jetpt_c->Write();
h_jetpt_b->Write();
h_jetpt_dbar->Write();
h_jetpt_ubar->Write();
h_jetpt_sbar->Write();
h_jetpt_cbar->Write();
h_jetpt_bbar->Write();
h_jetpt_hq->Write();
h_jetpt_lq->Write();
h_jetpt_qall->Write();
h_jetpt_dall->Write();
h_jetpt_uall->Write();
h_jetpt_sall->Write();
h_jetpt_call->Write();
h_jetpt_ball->Write();
h_jetpt_g->Write();
h_flavor_forb->Write();
h_jetpt_ud->Write();
h_jetpt_ubardbar->Write();
h_jetpt_other -> Write();
wf->Close();





} // end program
