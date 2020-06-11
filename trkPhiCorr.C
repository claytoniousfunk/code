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
 





void trkPhiCorr(){

        /// LOAD DATA ///
        auto f =TFile::Open("/home/clayton/Analysis/skims/PP_2017MC_5TeV_MC_Skim_1.root");
        auto em = new eventMap(f);
        em->isMC = 1;
        em->init();
        em->loadJet("akFlowPuCs4PFJetAnalyzer");
        em->loadTrack();
	em -> loadGenParticle();
	Long64_t n_evts = em->evtTree->GetEntries();
	
	/// DEFINE HISTOGRAMS ///
	TH1D *h_dphi = new TH1D("h_dphi","Track dPhi",500,-1*TMath::Pi(),2*TMath::Pi());

	TH1D *h_phi = new TH1D("h_phi","Track Phi",500,-3*TMath::Pi()/2,3*TMath::Pi()/2);
	//// CUTS ///////
	const double etamaxcut = 1.5;
	const double etamincut = 0.;
	const double pTmincut = 120.; // 
	const double pTmaxcut = 600.;
	const double refpTmincut = 50.;
	const double finalpTcut = 120.;
	


	
        int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;

        /// BEGIN EVENT LOOP ///
        //for(evi = 0; evi < 100; evi++){
	for (evi = 0; evi < n_evts; evi++){
        	em -> getEvent(evi);
                if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
                	cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                 }
        evi_frac = 100*evi/n_evts;
                                                                                                   

		/// BEGIN JET LOOP ///
		//int n_trk  = sizeof(em->trkphi)/sizeof(em->trkphi[0]); // get # of tracks
	/*	for(int jeti = 0; jeti < em->nGenJet(); jeti++){
			if((100*jeti/em->nGenJet())%5==0 && 100*jeti/em->nGenJet() > jeti_frac){
                        	cout<<"trk frac: "<<jeti_frac<<"%"<<endl;
                 	}
       		jeti_frac = 100*jeti/em->nGenJet();
	
			double phi1 = em -> genjetphi[jeti];
	*/		/// BEGIN TRACK LOOP ///
			//for(int trki = 0; trki < 10; trki++){
			for(int trki = 0; trki < em -> ntrk; trki++){
				double phi1 = em -> trkphi[trki];
				h_phi -> Fill(phi1);
				//for(int trkj = 0; trkj < 10; trkj++){
				for(int trkj = 0; trkj < em -> ntrk; trkj++){
					double phi2 = em -> trkphi[trkj];
					double dphi = phi1 - phi2;
					if(trki == trkj) continue;
					
					
					if((phi1-phi2) > -TMath::Pi()/2 && (phi1-phi2) < 3*TMath::Pi()/2){
						h_dphi -> Fill(dphi);
					}
					if((phi1-phi2) > 3*TMath::Pi()/2){
						dphi = dphi - 2*TMath::Pi();
						h_dphi -> Fill(dphi);
					}
					if((phi1-phi2) > -2*TMath::Pi() && (phi1-phi2) < -TMath::Pi()/2){
						dphi = dphi + 2*TMath::Pi();
						h_dphi -> Fill(dphi);
					}
					
				}
			}
			/// END TRACK LOOP ///
		// }
		/// END JET LOOP ///
	}
	/// END EVENT LOOP ///
auto c = new TCanvas("c","", 600, 500);
h_dphi -> Draw("hist");
c->SaveAs("h_dphi.pdf");
auto d = new TCanvas("d","", 600, 500);
h_phi -> Draw("hist");
d->SaveAs("h_phi.pdf");
auto wf = TFile::Open("dphi.root", "recreate");
h_dphi->Write();
h_phi -> Write();
wf->Close();
}
