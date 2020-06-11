


/*
This version is made to compile skims made by Dhanush, with only relevant data.  No need for eventMap

V5: New eta binning of constant 0.2 (as recommended by Dhanush)
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

//void jetParticleFraction(TString inf){
void jetEtaFractionV5(){

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
	int Nbins = 100;
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
	
	TH1D *h_flavor_forb = new TH1D("h_flavor_forb","",35,-10,25);

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
	h_flavor_forb->Sumw2();
	h_jeteta_other -> Sumw2();

	//// CUTS ///////
	const double etamaxcut = 1.5;
	const double etamincut = 0.;
	const double pTmincut = 120.; // default = 120.;
	const double pTmaxcut = 600.;
	const double refpTmincut = 50.;
	const double finalpTcut = 120.;

	// loop variables //	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	//double x = 0;
	//int y = 0;
	
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
			h_flavor_forb -> Fill(y);
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
		}

	}
	
	
		/*
 	// Define the Canvas
        TCanvas *c = new TCanvas("c", "canvas", 800, 800);
        TPad *pad1 = new TPad("pad1", "pad1", 0., 0.3, 1, 0.9);
        pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        //pad1->SetLogy();
       THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
       TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
       TH1D *h2 = (TH1D*)h_jeteta_d->Clone("h2"); // d quarks
       TH1D *h3 = (TH1D*)h_jeteta_dbar->Clone("h3"); // d antiquarks
       TH1D *h4 = (TH1D*)h_jeteta_u->Clone("h4"); // u quarks
       TH1D *h5 = (TH1D*)h_jeteta_ubar->Clone("h5"); // u antiquarks
       TH1D *h6 = (TH1D*)h_jeteta_hq->Clone("h6"); // heavy quarks
       TH1D *h7 = (TH1D*)h_jeteta_other->Clone("h7"); // unidentified particles

       h1 -> SetFillColor(kOrange);
       h2 -> SetFillColor(kBlue);
       h3 -> SetFillColor(kRed-6);
       h4 -> SetFillColor(kMagenta-9);
       h5 -> SetFillColor(kGreen);
       h6 -> SetFillColor(kPink+1);
	h7->SetFillColor(kBlue-10);
       h1 -> SetStats(0);
       h1 -> GetXaxis() -> SetTitle("Jet eta");
       h1 -> Divide(h_jeteta);
       h2 -> Divide(h_jeteta);
       h3 -> Divide(h_jeteta);
       h4 -> Divide(h_jeteta);
       h5 -> Divide(h_jeteta);
       h6 -> Divide(h_jeteta);
	h7->Divide(h_jeteta);
	
       fraction_stack->Add(h1);
       fraction_stack->Add(h2);
       fraction_stack->Add(h3);
       fraction_stack->Add(h4);
       fraction_stack->Add(h5);
       fraction_stack->Add(h6);
	//fraction_stack->Add(h7);
       fraction_stack->Draw("hist");
       auto legend = new TLegend(0.8,0.8,1.0,1.0);
       legend->AddEntry(h1,"g");
       legend->AddEntry(h2,"d");
       legend->AddEntry(h3,"d-bar");
       legend->AddEntry(h4,"u");
       legend->AddEntry(h5,"u-bar");
       legend->AddEntry(h6,"b,s,c,...");
	//legend->AddEntry(h7,"unidentified");
       legend->Draw();
	*/

	// Write data	
auto wf = TFile::Open("jetEtaFractionV5.root", "recreate");
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
h_flavor_forb->Write();
h_jeteta_other -> Write();
wf->Close();





} // end program
