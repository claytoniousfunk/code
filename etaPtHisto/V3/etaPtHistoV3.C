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


void etaPtHistoV3(){

	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/skims3/PbPbMC2018_hadronFlavor_Apr18.root");
    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_partonFlavor=0; 
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);

	
	
	/// DEFINE HISTOGRAMS ///
	double eta_min = -1.5;
	double eta_max = 1.5;
	double pTlow = 120.;
	double pThi = 600.;
	int Netabins = 100;
	int NpTbins = 66;
	TH1D *h_jeteta = new TH1D("h_jeteta","Jet eta",Netabins,eta_min,eta_max);
		// parse jetpt by species of particles
			TH1D *h_jeteta_d = new TH1D("h_jeteta_d","",Netabins,eta_min,eta_max); // down quark
      		TH1D *h_jeteta_u = new TH1D("h_jeteta_u","",Netabins,eta_min,eta_max); // up quark
      		TH1D *h_jeteta_s = new TH1D("h_jeteta_s","",Netabins,eta_min,eta_max); // strange quark
      		TH1D *h_jeteta_c = new TH1D("h_jeteta_c","",Netabins,eta_min,eta_max); // charm quark
      		TH1D *h_jeteta_b = new TH1D("h_jeteta_b","",Netabins,eta_min,eta_max); // bottom quark
      		TH1D *h_jeteta_dbar = new TH1D("h_jeteta_dbar","",Netabins,eta_min,eta_max); // down antiquark
      		TH1D *h_jeteta_ubar = new TH1D("h_jeteta_ubar","",Netabins,eta_min,eta_max); // up antiquark
      		TH1D *h_jeteta_sbar = new TH1D("h_jeteta_sbar","",Netabins,eta_min,eta_max); // strange antiquark
      		TH1D *h_jeteta_cbar = new TH1D("h_jeteta_cbar","",Netabins,eta_min,eta_max); // charm antiquark
      		TH1D *h_jeteta_bbar = new TH1D("h_jeteta_bbar","",Netabins,eta_min,eta_max); // bottom antiquark
      		TH1D *h_jeteta_g = new TH1D("h_jeteta_g","",Netabins,eta_min,eta_max); // gluon
      		TH1D *h_jeteta_hq = new TH1D("h_jeteta_hq","",Netabins,eta_min,eta_max); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
			TH1D *h_jeteta_other = new TH1D("h_jeteta_other","",Netabins,eta_min,eta_max); // others

	TH1D *h_jetpt = new TH1D("h_jetpt","Jet pT",NpTbins,pTlow,pThi);
		// parse jetpt by species of particles
			TH1D *h_jetpt_d = new TH1D("h_jetpt_d","",NpTbins,pTlow,pThi); // down quark
      		TH1D *h_jetpt_u = new TH1D("h_jetpt_u","",NpTbins,pTlow,pThi); // up quark
      		TH1D *h_jetpt_s = new TH1D("h_jetpt_s","",NpTbins,pTlow,pThi); // strange quark
      		TH1D *h_jetpt_c = new TH1D("h_jetpt_c","",NpTbins,pTlow,pThi); // charm quark
      		TH1D *h_jetpt_b = new TH1D("h_jetpt_b","",NpTbins,pTlow,pThi); // bottom quark
      		TH1D *h_jetpt_dbar = new TH1D("h_jetpt_dbar","",NpTbins,pTlow,pThi); // down antiquark
      		TH1D *h_jetpt_ubar = new TH1D("h_jetpt_ubar","",NpTbins,pTlow,pThi); // up antiquark
      		TH1D *h_jetpt_sbar = new TH1D("h_jetpt_sbar","",NpTbins,pTlow,pThi); // strange antiquark
      		TH1D *h_jetpt_cbar = new TH1D("h_jetpt_cbar","",NpTbins,pTlow,pThi); // charm antiquark
      		TH1D *h_jetpt_bbar = new TH1D("h_jetpt_bbar","",NpTbins,pTlow,pThi); // bottom antiquark
      		TH1D *h_jetpt_g = new TH1D("h_jetpt_g","",NpTbins,pTlow,pThi); // gluon
      		TH1D *h_jetpt_hq = new TH1D("h_jetpt_hq","",NpTbins,pTlow,pThi); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
			TH1D *h_jetpt_other = new TH1D("h_jetpt_other","",NpTbins,pTlow,pThi); // others
	
	
	TH1D *h_flavor_forb = new TH1D("h_flavor_forb","",35,-10,25);
	
	TH2F *h2 = new TH2F("h2","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_d = new TH2F("h2_d","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_u = new TH2F("h2_u","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_s = new TH2F("h2_s","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_c = new TH2F("h2_c","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_b = new TH2F("h2_b","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_dbar = new TH2F("h2_dbar","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_ubar = new TH2F("h2_ubar","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_sbar = new TH2F("h2_sbar","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_cbar = new TH2F("h2_cbar","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_bbar = new TH2F("h2_bbar","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_g = new TH2F("h2_g","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_hq = new TH2F("h2_hq","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);
		TH2F *h2_other = new TH2F("h2_other","",NpTbins,pTlow,pThi,Netabins,eta_min,eta_max);

	h_jeteta->Sumw2();
	h_jeteta_d->Sumw2();
	h_jeteta_u->Sumw2();
	h_jeteta_s->Sumw2();
	h_jeteta_c->Sumw2();
	h_jeteta_b->Sumw2();
	h_jeteta_dbar->Sumw2();
	h_jeteta_ubar->Sumw2();
	h_jeteta_sbar->Sumw2();
	h_jeteta_cbar->Sumw2();
	h_jeteta_bbar->Sumw2();
	h_jeteta_hq->Sumw2();
	h_jeteta_g->Sumw2();
	h_jeteta_other -> Sumw2();
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
	h_jetpt_g->Sumw2();
	h_flavor_forb->Sumw2();
	h_jetpt_other -> Sumw2();

	//// CUTS ///////
	const double etamaxcut = 1.5; // default = 1.5
	const double etamincut = 0.;
	const double pTmincut = 120.; // default = 120.
	const double pTmaxcut = 600.;
	const double refpTmincut = 50.;
	const double finalpTcut = 120.;



	// loop variables //	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	
	// Event loop
	for (evi = 0; evi < n_evts; evi++){
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
			if(y==0){h_jeteta_other -> Fill(x); continue;}
			if(fabs(x)>=etamaxcut || u <= pTmincut){continue;}
			h_jeteta -> Fill(x);
			h_jetpt -> Fill(u);
			if(y==1){h_jeteta_d->Fill(x); h_jetpt_d->Fill(u); h2_d->Fill(u,x);}
			if(y==2){h_jeteta_u->Fill(x); h_jetpt_u->Fill(u); h2_u->Fill(u,x);}
			if(y==3){h_jeteta_s->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_s->Fill(u); h_jetpt_hq->Fill(u); h2_s->Fill(u,x); h2_hq->Fill(u,x);}
			if(y==4){h_jeteta_c->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_c->Fill(u); h_jetpt_hq->Fill(u); h2_c->Fill(u,x); h2_hq->Fill(u,x);}
			if(y==5){h_jeteta_b->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_b->Fill(u); h_jetpt_hq->Fill(u); h2_b->Fill(u,x);h2_hq->Fill(u,x);}
			if(y==-1){h_jeteta_dbar->Fill(x);h_jetpt_dbar->Fill(u); h2_dbar->Fill(u,x);}
			if(y==-2){h_jeteta_ubar->Fill(x);h_jetpt_ubar->Fill(u); h2_ubar->Fill(u,x);}
			if(y==-3){h_jeteta_sbar->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_sbar->Fill(u); h_jetpt_hq->Fill(u); h2_sbar->Fill(u,x);h2_hq->Fill(u,x);}
			if(y==-4){h_jeteta_cbar->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_cbar->Fill(u); h_jetpt_hq->Fill(u); h2_cbar->Fill(u,x);h2_hq->Fill(u,x);}
			if(y==-5){h_jeteta_bbar->Fill(x); h_jeteta_hq->Fill(x);h_jetpt_bbar->Fill(u); h_jetpt_hq->Fill(u); h2_bbar->Fill(u,x);h2_hq->Fill(u,x);}
			if(y==21){h_jeteta_g->Fill(x);h_jetpt_g->Fill(u); h2_g->Fill(u,x);}
			if(y==0){h_jeteta_other -> Fill(x);h_jetpt_other -> Fill(u); h2_other->Fill(u,x);}

			h2->Fill(u,x);
			
			
			
			
			
			
			
			
			
			
			
			

			
		}

	}

	

	/*
		
	TCanvas *c1 = new TCanvas("c1", "c1",900,900);
   gStyle->SetOptStat(0);
   
   // Create the three pads
   TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
   center_pad->Draw();

   TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
   right_pad->Draw();

   TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
   top_pad->Draw();

   // Create, fill and project a 2D histogram.
   
   
 
   TH1D * projh2X = h2->ProjectionX();
   TH1D * projh2Y = h2->ProjectionY();

   // Drawing
   center_pad->cd();
   gStyle->SetPalette(1);
   h2->Draw("surf2");

   top_pad->cd();
   projh2X->SetFillColor(kBlue+1);
   projh2X->Draw("bar");

   right_pad->cd();
   projh2Y->SetFillColor(kBlue-2);
   projh2Y->Draw("hbar");
   
   c1->cd();
   TLatex *t = new TLatex();
   t->SetTextFont(42);
   t->SetTextSize(0.02);
   t->DrawLatex(0.6,0.88,"Jet pT and jet eta histogram");
   t->DrawLatex(0.6,0.85,"and its two projections.");
	*/
	
/*
TCanvas *c1 = new TCanvas("c1", "c1",900,900);

h2->Draw("surf2");*/
//delete f;



auto wf = TFile::Open("etaPtHistoV3.root", "recreate");
h2->Write();
h2_d->Write();
h2_u->Write();
h2_s->Write();
h2_c->Write();
h2_b->Write();
h2_dbar->Write();
h2_ubar->Write();
h2_sbar->Write();
h2_cbar->Write();
h2_bbar->Write();
h2_g->Write();
h2_other->Write();
h_jeteta->Write();
h_jeteta_d->Write();
h_jeteta_u->Write();
h_jeteta_s->Write();
h_jeteta_c->Write();
h_jeteta_b->Write();
h_jeteta_dbar->Write();
h_jeteta_ubar->Write();
h_jeteta_sbar->Write();
h_jeteta_cbar->Write();
h_jeteta_bbar->Write();
h_jeteta_hq->Write();
h_jeteta_g->Write();
h_jeteta_other -> Write();
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
h_jetpt_g->Write();
h_flavor_forb->Write();
h_jetpt_other -> Write();
wf->Close();


}


 // end program
