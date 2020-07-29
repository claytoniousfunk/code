


/*
This version is made to compile skims made by Dhanush, with only relevant data.  No need for eventMap

V5: New eta binning of constant 0.2 (as recommended by Dhanush)

V6: copying over and modifying code from jetPtFractionV6.C
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
void jetEtaFractionV6(float etacut1=-1.5, float etacut2=1.5){

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
	Float_t delta = 0.2; // width of bin
	Float_t etamin=-1.5;
	Float_t etamax=1.5;
	Float_t jteta_check=etamin;
	edges[0]=etamin;
	int numUsedBins = 0;
	for(int zz=1; jteta_check<etamax; zz++){
		
		jteta_check=jteta_check+delta;
		edges[zz]=edges[zz-1]+delta;
		numUsedBins=zz;
		
		
	}
	

	Float_t usedEdges[numUsedBins+1];
	for(int jj=0; jj<numUsedBins+1; jj++){
		usedEdges[jj]=edges[jj];
		cout<<usedEdges[jj]<<endl;
	}

	
	
	
		/// DEFINE HISTOGRAMS ///
	
	TH1D *h_jeteta = new TH1D("h_jeteta","Jet #eta",numUsedBins,usedEdges);
		// parse jetpt by species of particles
		TH1D *h_jeteta_d = new TH1D("h_jeteta_d","",numUsedBins,usedEdges); // down quark
      	TH1D *h_jeteta_u = new TH1D("h_jeteta_u","",numUsedBins,usedEdges); // up quark
      	TH1D *h_jeteta_s = new TH1D("h_jeteta_s","",numUsedBins,usedEdges); // strange quark
      	TH1D *h_jeteta_c = new TH1D("h_jeteta_c","",numUsedBins,usedEdges); // charm quark
      	TH1D *h_jeteta_b = new TH1D("h_jeteta_b","",numUsedBins,usedEdges); // bottom quark
      	TH1D *h_jeteta_dbar = new TH1D("h_jeteta_dbar","",numUsedBins,usedEdges); // down antiquark
      	TH1D *h_jeteta_ubar = new TH1D("h_jeteta_ubar","",numUsedBins,usedEdges); // up antiquark
      	TH1D *h_jeteta_sbar = new TH1D("h_jeteta_sbar","",numUsedBins,usedEdges); // strange antiquark
      	TH1D *h_jeteta_cbar = new TH1D("h_jeteta_cbar","",numUsedBins,usedEdges); // charm antiquark
      	TH1D *h_jeteta_bbar = new TH1D("h_jeteta_bbar","",numUsedBins,usedEdges); // bottom antiquark
      	TH1D *h_jeteta_g = new TH1D("h_jeteta_g","",numUsedBins,usedEdges); // gluon
      	TH1D *h_jeteta_hq = new TH1D("h_jeteta_hq","",numUsedBins,usedEdges); // all heavy quarks (s,c,b,s-bar,c-bar,b-bar)
      	TH1D *h_jeteta_lq = new TH1D("h_jeteta_lq","",numUsedBins,usedEdges); // all light quarks (u,d,...)
      	TH1D *h_jeteta_qall = new TH1D("h_jeteta_qall","",numUsedBins,usedEdges); // all quarks
      	TH1D *h_jeteta_dall = new TH1D("h_jeteta_dall","",numUsedBins,usedEdges); // d and dbar combined
      	TH1D *h_jeteta_uall = new TH1D("h_jeteta_uall","",numUsedBins,usedEdges); // u and ubar combined
      	TH1D *h_jeteta_sall = new TH1D("h_jeteta_sall","",numUsedBins,usedEdges); // s and sbar combined
      	TH1D *h_jeteta_call = new TH1D("h_jeteta_call","",numUsedBins,usedEdges); // c and cbar combined
      	TH1D *h_jeteta_ball = new TH1D("h_jeteta_ball","",numUsedBins,usedEdges); // b and bbar combined
      	TH1D *h_jeteta_ud = new TH1D("h_jeteta_ud","",numUsedBins,usedEdges); // u and d combined
      	TH1D *h_jeteta_ubardbar = new TH1D("h_jeteta_ubardbar","",numUsedBins,usedEdges); // u and d combined
		TH1D *h_jeteta_other = new TH1D("h_jeteta_other","",numUsedBins,usedEdges); // others
	
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
	h_jeteta_lq->Sumw2();
	h_jeteta_qall->Sumw2();
	h_jeteta_dall->Sumw2();
	h_jeteta_uall->Sumw2();
	h_jeteta_sall->Sumw2();
	h_jeteta_call->Sumw2();
	h_jeteta_ball->Sumw2();
	h_jeteta_g->Sumw2();
	h_flavor_forb->Sumw2();
	h_jeteta_ud->Sumw2();
	h_jeteta_ubardbar->Sumw2();
	h_jeteta_other -> Sumw2();

	//// CUTS ///////
	double pTmincut=100;
	double pTmaxcut=600;
	
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
			if(fabs(x)>etacut2 || u < pTmincut || u>pTmaxcut){continue;}
			h_flavor_forb -> Fill(y);
			h_jeteta->Fill(x);
			if(y==1){h_jeteta_d->Fill(x);h_jeteta_dall->Fill(x);h_jeteta_lq->Fill(x);h_jeteta_qall->Fill(x);h_jeteta_ud->Fill(x);}
			if(y==2){h_jeteta_u->Fill(x);h_jeteta_uall->Fill(x);h_jeteta_lq->Fill(x);h_jeteta_qall->Fill(x);h_jeteta_ud->Fill(x);}
			if(y==3){h_jeteta_s->Fill(x); h_jeteta_sall->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==4){h_jeteta_c->Fill(x); h_jeteta_call->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==5){h_jeteta_b->Fill(x); h_jeteta_ball->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==-1){h_jeteta_dbar->Fill(x);h_jeteta_dall->Fill(x);h_jeteta_lq->Fill(x);h_jeteta_qall->Fill(x);h_jeteta_ubardbar->Fill(x);}
			if(y==-2){h_jeteta_ubar->Fill(x);h_jeteta_uall->Fill(x);h_jeteta_lq->Fill(x);h_jeteta_qall->Fill(x);h_jeteta_ubardbar->Fill(x);}
			if(y==-3){h_jeteta_sbar->Fill(x); h_jeteta_sall->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==-4){h_jeteta_cbar->Fill(x); h_jeteta_call->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==-5){h_jeteta_bbar->Fill(x); h_jeteta_ball->Fill(x);h_jeteta_hq->Fill(x);h_jeteta_qall->Fill(x);}
			if(y==21){h_jeteta_g->Fill(x);}
			if(y==0){h_jeteta_other -> Fill(x);}
		}

	}

// generate dN/deta plot

	// INCLUDE BIN ERROR!
	
    Double_t Nbin = h_jeteta->GetSize();

    for (int i=0;i<Nbin;i++){
      Double_t x = h_jeteta->GetBinWidth(i);
      Double_t y_tot = h_jeteta->GetBinContent(i);
      Double_t y_g = h_jeteta_g->GetBinContent(i);
      Double_t y_d = h_jeteta_d->GetBinContent(i);
      Double_t y_dbar = h_jeteta_dbar->GetBinContent(i);
      Double_t y_dall = h_jeteta_dall->GetBinContent(i);
      Double_t y_u = h_jeteta_u->GetBinContent(i);
      Double_t y_ubar = h_jeteta_ubar->GetBinContent(i);
      Double_t y_uall = h_jeteta_uall->GetBinContent(i);
      Double_t y_s = h_jeteta_s->GetBinContent(i);
      Double_t y_sbar = h_jeteta_sbar->GetBinContent(i);
      Double_t y_sall = h_jeteta_sall->GetBinContent(i);
      Double_t y_c = h_jeteta_c->GetBinContent(i);
      Double_t y_cbar = h_jeteta_cbar->GetBinContent(i);
      Double_t y_call = h_jeteta_call->GetBinContent(i);
      Double_t y_b = h_jeteta_b->GetBinContent(i);
      Double_t y_bbar = h_jeteta_bbar->GetBinContent(i);
      Double_t y_ball = h_jeteta_ball->GetBinContent(i);
      Double_t y_hq = h_jeteta_hq->GetBinContent(i);
      Double_t y_lq = h_jeteta_lq->GetBinContent(i);
      Double_t y_qall = h_jeteta_qall->GetBinContent(i);
      Double_t y_other = h_jeteta_other->GetBinContent(i);
      Double_t y_ud = h_jeteta_ud->GetBinContent(i);
      Double_t y_ubardbar = h_jeteta_ubardbar->GetBinContent(i);


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

        h_jeteta->SetBinContent(i,r_tot);
        h_jeteta_g->SetBinContent(i,r_g);
        h_jeteta_d->SetBinContent(i,r_d);
        h_jeteta_dbar->SetBinContent(i,r_dbar);
        h_jeteta_dall->SetBinContent(i,r_dall);
        h_jeteta_u->SetBinContent(i,r_u);
        h_jeteta_ubar->SetBinContent(i,r_ubar);
        h_jeteta_uall->SetBinContent(i,r_uall);
        h_jeteta_s->SetBinContent(i,r_s);
        h_jeteta_sbar->SetBinContent(i,r_sbar);
        h_jeteta_sall->SetBinContent(i,r_sall);
        h_jeteta_c->SetBinContent(i,r_c);
        h_jeteta_cbar->SetBinContent(i,r_cbar);
        h_jeteta_call->SetBinContent(i,r_call);
        h_jeteta_b->SetBinContent(i,r_b);
        h_jeteta_bbar->SetBinContent(i,r_bbar);
        h_jeteta_ball->SetBinContent(i,r_ball);
        h_jeteta_hq->SetBinContent(i,r_hq);
        h_jeteta_lq->SetBinContent(i,r_lq);
        h_jeteta_qall->SetBinContent(i,r_qall);
        h_jeteta_other->SetBinContent(i,r_other);
        h_jeteta_ud->SetBinContent(i,r_ud);
        h_jeteta_ubardbar->SetBinContent(i,r_ubardbar);
      }

      
    }
 
	


    
	// Write data	
auto wf = TFile::Open("jetEtaFractionV6.root", "recreate");
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
h_jeteta_lq->Write();
h_jeteta_qall->Write();
h_jeteta_dall->Write();
h_jeteta_uall->Write();
h_jeteta_sall->Write();
h_jeteta_call->Write();
h_jeteta_ball->Write();
h_jeteta_g->Write();
h_flavor_forb->Write();
h_jeteta_ud->Write();
h_jeteta_ubardbar->Write();
h_jeteta_other -> Write();
wf->Close();





} // end program