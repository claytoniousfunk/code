


/*
This version is made to compile skims made by Dhanush, with only relevant data.  No need for eventMap
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
void jetPtFractionV3(){

	/// LOAD DATA ///
        TFile* f =TFile::Open("/home/clayton/Analysis/skims3/PbPbMC2018_hadronFlavor_Apr18.root");
        TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_partonFlavor=0; 
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);

	/// DEFINE HISTOGRAMS ///
	double pTlow = 120.;
	double pThi = 600.;
	int Nbins = 100;
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
			if(fabs(x)>=etamaxcut || u <= pTmincut){continue;}
			h_flavor_forb -> Fill(y);
			h_jetpt -> Fill(u);
			if(y==1){h_jetpt_d->Fill(u);}
			if(y==2){h_jetpt_u->Fill(u);}
			if(y==3){h_jetpt_s->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==4){h_jetpt_c->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==5){h_jetpt_b->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==-1){h_jetpt_dbar->Fill(u);}
			if(y==-2){h_jetpt_ubar->Fill(u);}
			if(y==-3){h_jetpt_sbar->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==-4){h_jetpt_cbar->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==-5){h_jetpt_bbar->Fill(u); h_jetpt_hq->Fill(u);}
			if(y==21){h_jetpt_g->Fill(u);}
			if(y==0){h_jetpt_other -> Fill(u);}
		}

	}
	
	
		
 	// Define the Canvas
        TCanvas *c = new TCanvas("c", "canvas", 800, 800);
        TPad *pad1 = new TPad("pad1", "pad1", 0., 0.3, 1, 0.9);
        pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        //pad1->SetLogy();
       THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
       TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
       TH1D *h2 = (TH1D*)h_jetpt_d->Clone("h2"); // d quarks
       TH1D *h3 = (TH1D*)h_jetpt_dbar->Clone("h3"); // d antiquarks
       TH1D *h4 = (TH1D*)h_jetpt_u->Clone("h4"); // u quarks
       TH1D *h5 = (TH1D*)h_jetpt_ubar->Clone("h5"); // u antiquarks
       TH1D *h6 = (TH1D*)h_jetpt_hq->Clone("h6"); // heavy quarks
       TH1D *h7 = (TH1D*)h_jetpt_other->Clone("h7"); // unidentified particles

       h1 -> SetFillColor(kOrange);
       h2 -> SetFillColor(kBlue);
       h3 -> SetFillColor(kRed-6);
       h4 -> SetFillColor(kMagenta-9);
       h5 -> SetFillColor(kGreen);
       h6 -> SetFillColor(kPink+1);
	h7->SetFillColor(kBlue-10);
       h1 -> SetStats(0);
       h1 -> GetXaxis() -> SetTitle("Jet pt");
       h1 -> Divide(h_jetpt);
       h2 -> Divide(h_jetpt);
       h3 -> Divide(h_jetpt);
       h4 -> Divide(h_jetpt);
       h5 -> Divide(h_jetpt);
       h6 -> Divide(h_jetpt);
	h7->Divide(h_jetpt);
	
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
	

	// Write data	
auto wf = TFile::Open("jetPtFractionV3.root", "recreate");
c->SaveAs("jetPtFractionV3.pdf");
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





} // end program
