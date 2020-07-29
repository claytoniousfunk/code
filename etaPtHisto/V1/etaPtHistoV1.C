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


void etaPtHisto(){

	
	/// DEFINE HISTOGRAMS ///
	double eta_min = -1.5;
	double eta_max = 1.5;
	double pTlow = 0.;
	double pThi = 600.;
	int Nbins = 40;
	TH1D *h_jeteta = new TH1D("h_jeteta","Jet eta",Nbins,eta_min,eta_max);
		// parse jetpt by species of particles
		TH1D *h_jeteta_d = new TH1D("h_jeteta_d","",Nbins,eta_min,eta_max); // down quark
      		TH1D *h_jeteta_u = new TH1D("h_jeteta_u","",Nbins,eta_min,eta_max); // up quark
      		TH1D *h_jeteta_s = new TH1D("h_jeteta_s","",Nbins,eta_min,eta_max); // strange quark
      		TH1D *h_jeteta_c = new TH1D("h_jeteta_c","",Nbins,eta_min,eta_max); // charm quark
      		TH1D *h_jeteta_b = new TH1D("h_jeteta_b","",Nbins,eta_min,eta_max); // bottom quark
      		TH1D *h_jeteta_dbar = new TH1D("h_jeteta_dbar","",Nbins,eta_min,eta_max); // down antiquark
      		TH1D *h_jeteta_ubar = new TH1D("h_jeteta_ubar","",Nbins,eta_min,eta_max); // up antiquark
      		TH1D *h_jeteta_sbar = new TH1D("h_jeteta_sbar","",Nbins,eta_min,eta_max); // strange antiquark
      		TH1D *h_jeteta_cbar = new TH1D("h_jeteta_cbar","",Nbins,eta_min,eta_max); // charm antiquark
      		TH1D *h_jeteta_bbar = new TH1D("h_jeteta_bbar","",Nbins,eta_min,eta_max); // bottom antiquark
      		TH1D *h_jeteta_g = new TH1D("h_jeteta_g","",Nbins,eta_min,eta_max); // gluon
      		TH1D *h_jeteta_hq = new TH1D("h_jeteta_hq","",Nbins,eta_min,eta_max); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
		TH1D *h_jeteta_other = new TH1D("h_jeteta_other","",Nbins,eta_min,eta_max); // others

	TH1D *h_jetpt = new TH1D("h_jetpt","Jet pT",Nbins,pTlow,pThi);
		// parse jetpt by species of particles
		TH1D *h_jetpt_d = new TH1D("h_jetpt_d","",Nbins,pTlow,pThi); // down quark
      		TH1D *h_jetpt_u = new TH1D("h_jetpt_u","",Nbins,pTlow,pThi); // up quark
      		TH1D *h_jetpt_s = new TH1D("h_jetpt_s","",Nbins,pTlow,pThi); // strange quark
      		TH1D *h_jetpt_c = new TH1D("h_jetpt_c","",Nbins,pTlow,pThi); // charm quark
      		TH1D *h_jetpt_b = new TH1D("h_jetpt_b","",Nbins,pTlow,pThi); // bottom quark
      		TH1D *h_jetpt_dbar = new TH1D("h_jetpt_dbar","",Nbins,pTlow,pThi); // down antiquark
      		TH1D *h_jetpt_ubar = new TH1D("h_jetpt_ubar","",Nbins,pTlow,pThi); // up antiquark
      		TH1D *h_jetpt_sbar = new TH1D("h_jetpt_sbar","",Nbins,pTlow,pThi); // strange antiquark
      		TH1D *h_jetpt_cbar = new TH1D("h_jetpt_cbar","",Nbins,pTlow,pThi); // charm antiquark
      		TH1D *h_jetpt_bbar = new TH1D("h_jetpt_bbar","",Nbins,pTlow,pThi); // bottom antiquark
      		TH1D *h_jetpt_g = new TH1D("h_jetpt_g","",Nbins,pTlow,pThi); // gluon
      		TH1D *h_jetpt_hq = new TH1D("h_jetpt_hq","",Nbins,pTlow,pThi); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
		TH1D *h_jetpt_other = new TH1D("h_jetpt_other","",Nbins,pTlow,pThi); // others
	
	
	TH1D *h_flavor_forb = new TH1D("h_flavor_forb","",35,-10,25);
	
	TH2F *h2 = new TH2F("h2","",Nbins,pTlow,pThi,Nbins,eta_min,eta_max);

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

	/// LOAD DATA ///
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/newSkims"); 

   	if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    	{ 
        	printf("Could not open current directory" ); 
        	return 0; 
    	}
	int Nskims = 27;	
	char* filenames[Nskims];
	int i = 0;
    	if(dr){ 
		 
		while((de=readdir(dr)) != NULL) {
            	//printf("%s\n", de->d_name);
	  	filenames[i] = (char*) malloc(strlen(de->d_name));
          	strcpy(filenames[i],de->d_name);
          	//printf("%s", filenames[i]);
          	//printf("\n");
	  	i++;
    	    }
 	closedir(dr);
	}

	
	
for(int ii = 0; ii < Nskims; ii++){
	if(strcmp(filenames[ii],".")==0 || strcmp(filenames[ii],"..")==0 ){continue;}
	cout << filenames[ii] << endl;
	char dir1[33] = "/home/clayton/Analysis/newSkims/";
	char* fileName = strcat(dir1,filenames[ii]);
	//cout << fileName << endl;
        //auto f =TFile::Open(fileName);
	//auto f = TFile::Open(inf);
       // auto em = new eventMap(f);
       // em->isMC = 1;
       // em->init();
       // em->loadJet("akFlowPuCs4PFJetAnalyzer");
       // em->loadTrack();
       // em -> loadGenParticle();
       //Long64_t n_evts = em->evtTree->GetEntries();

/*
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
        	em -> getEvent(evi);
                if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/n_evts;
                                                                                                   

		// Jet loop
		
		for(int jeti = 0; jeti < em->njet; jeti++){
			
			double x = em->jeteta[jeti];
			double u = em->jetpt[jeti];
			int y = em-> flavor_forb[jeti];
			if(y==0){h_jeteta_other -> Fill(x); continue;}
			if(fabs(x)>=etamaxcut || u <= pTmincut){continue;}
			h_jeteta -> Fill(x);
			h_jetpt -> Fill(u);
			if(y==1){h_jeteta_d->Fill(x);}
			if(y==2){h_jeteta_u->Fill(x);}
			if(y==3){h_jeteta_s->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==4){h_jeteta_c->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==5){h_jeteta_b->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==-1){h_jeteta_dbar->Fill(x);}
			if(y==-2){h_jeteta_ubar->Fill(x);}
			if(y==-3){h_jeteta_sbar->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==-4){h_jeteta_cbar->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==-5){h_jeteta_bbar->Fill(x); h_jeteta_hq->Fill(x);}
			if(y==21){h_jeteta_g->Fill(x);}
			if(y==0){h_jeteta_other -> Fill(x);}
			if(y==1){h_jetpt_d->Fill(x);}
			if(y==2){h_jetpt_u->Fill(x);}
			if(y==3){h_jetpt_s->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==4){h_jetpt_c->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==5){h_jetpt_b->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==-1){h_jetpt_dbar->Fill(x);}
			if(y==-2){h_jetpt_ubar->Fill(x);}
			if(y==-3){h_jetpt_sbar->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==-4){h_jetpt_cbar->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==-5){h_jetpt_bbar->Fill(x); h_jetpt_hq->Fill(x);}
			if(y==21){h_jetpt_g->Fill(x);}
			if(y==0){h_jetpt_other -> Fill(x);}

			h2->Fill(u,x);
		}

	}

	*/
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

//TCanvas *c1 = new TCanvas("c1", "c1",900,900);

//h2->Draw("surf2");
//delete f;
}
} // end program
