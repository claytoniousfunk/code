
// V4: using variable binning


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


void etaPtHistoV4(){

	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/skims3/PbPbMC2018_hadronFlavor_Apr18.root");
    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_partonFlavor=0; 
	inp_tree->SetBranchAddress("pf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("pf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("pf_jtpt",&pf_jtpt);
    
	
	/// DEFINE BINNING ///

	// create eta edges array
	const Int_t Nbins = 1000;
	Float_t etaEdges[Nbins+1];
	Float_t deltaEta = 0.2; // width of bin
	Float_t etamin=-1.5;
	Float_t etamax=1.5;
	Float_t jteta_check=etamin;
	etaEdges[0]=etamin;
	int numUsedEtaBins = 0;
	for(int zz=1; jteta_check<etamax; zz++){
		
		jteta_check=jteta_check+deltaEta;
		etaEdges[zz]=etaEdges[zz-1]+deltaEta;
		numUsedEtaBins=zz;
		
		
	}
	

	Float_t usedEtaEdges[numUsedEtaBins+1];
	for(int jj=0; jj<numUsedEtaBins+1; jj++){
		usedEtaEdges[jj]=etaEdges[jj];
	}
	

	// create pT edges array
	
	Float_t ptEdges[Nbins+1];
	Float_t deltaPt = 10.0; // initial width of bin (pT) // 10 GeV up to 100, 20 GeV up to 250, 250 out = 50 GeV, 400-500 = one bin
	Float_t ptmin=100;
	Float_t ptmax=600;
	Float_t jtpt_check=ptmin;
	ptEdges[0]=ptmin;
	int numUsedPtBins = 0;
	for(int yy=1; jtpt_check<ptmax; yy++){

		if(jtpt_check<100){
			deltaPt=10.;	
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<250){
			deltaPt=20.;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<400){
			deltaPt=50;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
			if(ptEdges[yy]>=400){
				ptEdges[yy]=400;
				ptEdges[yy+1]=ptmax;
				numUsedPtBins=numUsedPtBins+1;
				break;
			}
		}
		
		
		
		jtpt_check=jtpt_check+deltaPt;
		
	}
	

	Float_t usedPtEdges[numUsedPtBins+1];
	for(int ii=0; ii<numUsedPtBins+1; ii++){
		usedPtEdges[ii]=ptEdges[ii];
		cout << usedPtEdges[ii] << endl;
	}


	// DEFINE 2D HISTOGRAMS 
		TH1D *h_jetpt = new TH1D("h_jetpt","All jets",numUsedPtBins,usedPtEdges);
		TH1D *h_jeteta = new TH1D("h_jeteta","All jets",numUsedEtaBins,usedEtaEdges);
		TH2D *h2 = new TH2D("h2","All jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		// parse jet by species of particles
		TH2D *h2_d = new TH2D("h2_d","d jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_u = new TH2D("h2_u","u jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_s = new TH2D("h2_s","s jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_c = new TH2D("h2_c","c jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_b = new TH2D("h2_b","b jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dbar = new TH2D("h2_dbar","#bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubar = new TH2D("h2_ubar","#bar{u} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sbar = new TH2D("h2_sbar","#bar{s} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_cbar = new TH2D("h2_cbar","#bar{c} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_bbar = new TH2D("h2_bbar","#bar{b} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_g = new TH2D("h2_g","gluon jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_hq = new TH2D("h2_hq","heavy quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_lq = new TH2D("h2_lq","light quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_qall = new TH2D("h2_qall","All quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dall = new TH2D("h2_dall","d & #bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_uall = new TH2D("h2_uall","u & #bar{u} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sall = new TH2D("h2_sall","s & #bar{s} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_call = new TH2D("h2_call","c & #bar{c} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ball = new TH2D("h2_ball","b & #bar{b} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ud = new TH2D("h2_ud","u & d jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubardbar = new TH2D("h2_ubardbar","#bar{u} & #bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_other = new TH2D("h2_other","ghost jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);

		h2->Sumw2();
	h2_d->Sumw2();
	h2_u->Sumw2();
	h2_s->Sumw2();
	h2_c->Sumw2();
	h2_b->Sumw2();
	h2_dbar->Sumw2();
	h2_ubar->Sumw2();
	h2_sbar->Sumw2();
	h2_cbar->Sumw2();
	h2_bbar->Sumw2();
	h2_hq->Sumw2();
	h2_lq->Sumw2();
	h2_qall->Sumw2();
	h2_dall->Sumw2();
	h2_uall->Sumw2();
	h2_sall->Sumw2();
	h2_call->Sumw2();
	h2_ball->Sumw2();
	h2_g->Sumw2();
	h2_ud->Sumw2();
	h2_ubardbar->Sumw2();
	h2_other -> Sumw2();



	//// CUTS ///////
	const double etamaxcut = 1.5; // default = 1.5
	const double etamincut = 0.;
	const double pTmincut = 100.; // default = 120.
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

			if(fabs(x)>etamaxcut || u < pTmincut || u>pTmaxcut){continue;}

			h_jetpt->Fill(u);
			h_jeteta->Fill(x);

			int y = pf_partonFlavor->at(jeti);
			double w=1; // weight.  will equal inverse area
			double de = 0.2;
			double dp = 10.0;
			if(u<100){dp=10.;}
			else if(u<260){dp=20.;}
			else if(u<360){dp=50.;}
			else if(u<400){dp=40.;}
			else{dp=200;}
			double a = dp*de;
			w=1/a;

			
			
			if(y==0){h2_other -> Fill(x,u,w);}
			if(y==1){h2_d->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_dall->Fill(x,u,w); h2_ud->Fill(x,u,w);}
			if(y==2){h2_u->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_uall->Fill(x,u,w); h2_ud->Fill(x,u,w);}
			if(y==3){h2_s->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_sall->Fill(x,u,w);}
			if(y==4){h2_c->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_call->Fill(x,u,w);}
			if(y==5){h2_b->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_ball->Fill(x,u,w);}
			if(y==-1){h2_dbar->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_dall->Fill(x,u,w); h2_ubardbar->Fill(x,u,w);}
			if(y==-2){h2_ubar->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_uall->Fill(x,u,w); h2_ubardbar->Fill(x,u,w);}
			if(y==-3){h2_sbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_sall->Fill(x,u,w);}
			if(y==-4){h2_cbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_call->Fill(x,u,w);}
			if(y==-5){h2_bbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_ball->Fill(x,u,w);}
			if(y==21){h2_g->Fill(x,u,w);}


			h2->Fill(x,u,w);
		
		}

	}

	

// renormalize
	/*
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

*/ 
	



auto wf = TFile::Open("etaPtHistoV4.root", "recreate");
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
h2_hq->Write();
h2_lq->Write();
h2_qall->Write();
h2_dall->Write();
h2_uall->Write();
h2_sall->Write();
h2_call->Write();
h2_ball->Write();
h2_g->Write();
h2_ud->Write();
h2_ubardbar->Write();
h2_other->Write();
wf->Close();


}


 // end program
