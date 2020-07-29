


/*
This version is made to compile skims made by Dhanush, with only relevant data.  No need for eventMap
V5 :  implementing input variables for pT cuts.  Added more histogram definitions (different combinations)
V6:   Sumw2() moved to the end, after re-defining bin contents.  This fixed an error bar issue
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
void jetPtFractionV6(Float_t ptcut1=100, Float_t ptcut2=600){

	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/skims3/PbPbMC2018_hadronFlavor_Apr18.root");
    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_partonFlavor=0; 
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);

	/// DEFINE BINNING ///
	
	
	const Int_t Nbins = 1000;
	Float_t edges[Nbins+1];
	Float_t delta = 10.0; // initial width of bin (pT) // 10 GeV up to 100, 20 GeV up to 250, 250 out = 50 GeV, 400-500 = one bin
	Float_t jtpt_check=ptcut1;
	edges[0]=ptcut1;
	int numUsedBins = 0;
	for(int zz=1; jtpt_check<ptcut2; zz++){

		if(jtpt_check<100){
			delta=10.;	
			edges[zz]=edges[zz-1]+delta;
			numUsedBins=numUsedBins+1;
			cout<<numUsedBins<<endl;
		}
		else if (jtpt_check<250){
			delta=20.;
			edges[zz]=edges[zz-1]+delta;
			numUsedBins=numUsedBins+1;
			cout<<numUsedBins<<endl;
		}
		else if (jtpt_check<400){
			delta=50;
			edges[zz]=edges[zz-1]+delta;
			numUsedBins=numUsedBins+1;
			cout<<numUsedBins<<endl;
			if(edges[zz]>=400){
				edges[zz]=400;
				edges[zz+1]=ptcut2;
				numUsedBins=numUsedBins+1;
				cout<<numUsedBins<<endl;
				break;
			}
		}
		
		
		
		jtpt_check=jtpt_check+delta;
		
		
	}
	

	Float_t usedEdges[numUsedBins+1];
	for(int jj=0; jj<numUsedBins+1; jj++){
		usedEdges[jj]=edges[jj];
		cout << usedEdges[jj] << endl;
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
	const double pTmincut = 100; // default = 120.;
	const double pTmaxcut = 600;
	
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
		
		for(int jeti = 0; jeti < (int) pf_jtpt->size(); jeti++){
			
			double x = pf_jteta->at(jeti);
			double u = pf_jtpt->at(jeti);
			int y = pf_partonFlavor->at(jeti);
			if(fabs(x)>etamaxcut || u < pTmincut || u>pTmaxcut){continue;}
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

// generate dN/dp_T plot

	// INCLUDE BIN ERROR!

    Double_t Nbin = h_jetpt->GetSize();

    for (int i=0;i<Nbin;i++){
      Double_t x = h_jetpt->GetBinWidth(i);
      Double_t y_tot = h_jetpt->GetBinContent(i);
      Double_t y_g = h_jetpt_g->GetBinContent(i);
      Double_t y_d = h_jetpt_d->GetBinContent(i);
      Double_t y_dbar = h_jetpt_dbar->GetBinContent(i);
      Double_t y_dall = h_jetpt_dall->GetBinContent(i);
      Double_t y_u = h_jetpt_u->GetBinContent(i);
      Double_t y_ubar = h_jetpt_ubar->GetBinContent(i);
      Double_t y_uall = h_jetpt_uall->GetBinContent(i);
      Double_t y_s = h_jetpt_s->GetBinContent(i);
      Double_t y_sbar = h_jetpt_sbar->GetBinContent(i);
      Double_t y_sall = h_jetpt_sall->GetBinContent(i);
      Double_t y_c = h_jetpt_c->GetBinContent(i);
      Double_t y_cbar = h_jetpt_cbar->GetBinContent(i);
      Double_t y_call = h_jetpt_call->GetBinContent(i);
      Double_t y_b = h_jetpt_b->GetBinContent(i);
      Double_t y_bbar = h_jetpt_bbar->GetBinContent(i);
      Double_t y_ball = h_jetpt_ball->GetBinContent(i);
      Double_t y_hq = h_jetpt_hq->GetBinContent(i);
      Double_t y_lq = h_jetpt_lq->GetBinContent(i);
      Double_t y_qall = h_jetpt_qall->GetBinContent(i);
      Double_t y_other = h_jetpt_other->GetBinContent(i);
      Double_t y_ud = h_jetpt_ud->GetBinContent(i);
      Double_t y_ubardbar = h_jetpt_ubardbar->GetBinContent(i);


      Double_t r_tot = 0;
      Double_t r_g = 0;
      Double_t r_d = 0;
      Double_t r_dbar = 0;
      Double_t r_dall = 0;
      Double_t r_u = 0;
      Double_t r_ubar = 0;
      Double_t r_uall = 0;
      Double_t r_s = 0;
      Double_t r_sbar = 0;
      Double_t r_sall = 0;
      Double_t r_c = 0;
      Double_t r_cbar = 0;
      Double_t r_call = 0;
      Double_t r_b = 0;
      Double_t r_bbar = 0;
      Double_t r_ball = 0;
      Double_t r_hq = 0;
      Double_t r_lq = 0;
      Double_t r_qall = 0;
      Double_t r_other = 0;
      Double_t r_ud = 0;
      Double_t r_ubardbar = 0;

      if(x!=0){
        r_tot = y_tot/x;
        r_g = y_g/x;
        r_d = y_d/x;
        r_dbar = y_dbar/x;
        r_dall = y_dall/x;
        r_u = y_u/x;
        r_ubar = y_ubar/x;
        r_uall = y_uall/x;
        r_s = y_s/x;
        r_sbar = y_sbar/x;
        r_sall = y_sall/x;
        r_c = y_c/x;
        r_cbar = y_cbar/x;
        r_call = y_call/x;
        r_b = y_b/x;
        r_bbar = y_bbar/x;
        r_ball = y_ball/x;
        r_hq = y_hq/x;
        r_lq = y_lq/x;
        r_qall = y_qall/x;
        r_other = y_other/x;
        r_ud = y_ud/x;
        r_ubardbar = y_ubardbar/x;

        h_jetpt->SetBinContent(i,r_tot);
        h_jetpt_g->SetBinContent(i,r_g);
        h_jetpt_d->SetBinContent(i,r_d);
        h_jetpt_dbar->SetBinContent(i,r_dbar);
        h_jetpt_dall->SetBinContent(i,r_dall);
        h_jetpt_u->SetBinContent(i,r_u);
        h_jetpt_ubar->SetBinContent(i,r_ubar);
        h_jetpt_uall->SetBinContent(i,r_uall);
        h_jetpt_s->SetBinContent(i,r_s);
        h_jetpt_sbar->SetBinContent(i,r_sbar);
        h_jetpt_sall->SetBinContent(i,r_sall);
        h_jetpt_c->SetBinContent(i,r_c);
        h_jetpt_cbar->SetBinContent(i,r_cbar);
        h_jetpt_call->SetBinContent(i,r_call);
        h_jetpt_b->SetBinContent(i,r_b);
        h_jetpt_bbar->SetBinContent(i,r_bbar);
        h_jetpt_ball->SetBinContent(i,r_ball);
        h_jetpt_hq->SetBinContent(i,r_hq);
        h_jetpt_lq->SetBinContent(i,r_lq);
        h_jetpt_qall->SetBinContent(i,r_qall);
        h_jetpt_other->SetBinContent(i,r_other);
        h_jetpt_ud->SetBinContent(i,r_ud);
        h_jetpt_ubardbar->SetBinContent(i,r_ubardbar);
      }

      
    }

	


    
	// Write data	
auto wf = TFile::Open("jetPtFractionV6.root", "recreate");
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
