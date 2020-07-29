
// "etaPtHisto"
// V4 : using variable binning
// V5 : muon jet filtering, using HydJet_x.root files
// V6 : using new weigthing method, pthat.  No longer reweighting according to bin size (delegate to plotting macros).  Including muRelPt calculation
// V7 : desinged to run /skims5/ files. Added histogram of untagged jets.
// V8 : changed variable binning
// V9 : implementing use of hadronFlavor into particle identification.  changed pT max to 500.
// V10 : getting rid of untagged variable

/// "makeData" 
// V1 : adopting from etaPtHisto series.  Got rid of any "_t" variables.  Cleaned up some unused variables.  Most importantly, added vz & hiBin cut
// V2 : to be used on data files, no flavor tagging
// V3 : added reweighting of vz and hibin and overall weight
// V4 : pthat region study
// V5 : add exclusion of pthat<50, ptjet>80 events
// V6 : get rid of variable binning once and for all! (Need to delegate this to the plotting macros).  Creation of h2_lowPthatHighPtjet, weighting investigation histos
// V7 : centrality regions added.  Added h_weight (histo of pthat weight values)
// V8 : Added pthat regions histograms from V4.
// V9 : updated exclusion logic from V5.  Criteria for reco jets updated to: exclude pthat < 50, ptjet>100 events

#include "/home/clayton/Analysis/code/myProcesses/hiforest/plugin/eventMap_hiForest.h"
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



void makeDataV9(bool isReco = 0, double muPtCut = 5.0){

	//////////////////////////////////////////////////////////////////////// jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_jtphi=0, *pf_partonFlavor=0, *pf_hadronFlavor=0; 

	//////////////////////////////////////////////////////////////////////// muon jet variables ////////////////////////////////////////////////////////////////////////

	vector <float> *muPt=0, *muEta=0, *muPhi=0, *muChi2NDF=0, *muInnerD0=0, *muInnerDz=0;
	int nMu, hiBin;
	vector <int> *muIsGlobal=0, *muIsTracker=0, *muMuonHits=0, *muStations=0, *muTrkLayers=0, *muPixelHits=0;

	//////////////////////////////////////////////////////////////////////// event variables ////////////////////////////////////////////////////////////////////////

	float weight, vz;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////// DEFINE BINNING ///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////// DEFINE HISTOGRAMS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int NEtaBins = 300;
	double etaMin = -1.5;
	double etaMax = 1.5;
	int NPtBins = 450;
	double ptMin = 50.0;
	double ptMax = 500.0;


		///////////////////////////////////////////////////////////////////////  EVENTS  ///////////////////////////////////////////////////////////////////////

		TH1D *h_vz = new TH1D("h_vz","vertex position",60,-15.0,15.0);
		TH1D *h_vz_noWeight = new TH1D("h_vz_noWeight","vertex position",60,-15.0,15.0);
		TH1D *h_hiBin = new TH1D("h_hiBin","centrality",200,0.0,200.0);
		TH1D *h_hiBin_noWeight = new TH1D("h_hiBin_noWeight","centrality",200,0.0,200.0);
		TH2D *h2_lowPthatHighPtjet = new TH2D("h2_lowPthatHighPtjet","low pthat, high ptjet",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH1D *h_weight = new TH1D("h_weight","pthat weight parameter",100,0.0,0.05);
	

		///////////////////////////////////////////////////////////////////////  JETS  /////////////////////////////////////////////////////////////////////////

		TH2D *h2 = new TH2D("h2","All jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		// parse jet by species of particles
		TH2D *h2_d = new TH2D("h2_d","d jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_u = new TH2D("h2_u","u jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_s = new TH2D("h2_s","s jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_c = new TH2D("h2_c","c jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_b = new TH2D("h2_b","b jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_dbar = new TH2D("h2_dbar","#bar{d} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ubar = new TH2D("h2_ubar","#bar{u} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_sbar = new TH2D("h2_sbar","#bar{s} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_cbar = new TH2D("h2_cbar","#bar{c} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_bbar = new TH2D("h2_bbar","#bar{b} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_g = new TH2D("h2_g","gluon jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_hq = new TH2D("h2_hq","heavy quark jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_lq = new TH2D("h2_lq","light quark jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_qall = new TH2D("h2_qall","All quark jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_dall = new TH2D("h2_dall","d & #bar{d} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_uall = new TH2D("h2_uall","u & #bar{u} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_sall = new TH2D("h2_sall","s & #bar{s} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_call = new TH2D("h2_call","c & #bar{c} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ball = new TH2D("h2_ball","b & #bar{b} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ud = new TH2D("h2_ud","u & d jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ubardbar = new TH2D("h2_ubardbar","#bar{u} & #bar{d} jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_other = new TH2D("h2_other","ghost jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		// weighting investigation
		TH2D *h2_noVzOrHiBinWeight = new TH2D("h2_noVzOrHiBinWeight","All jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);

		
		/////////////////////////////////////////////////////////////////   MUON-TAGGED JETS  //////////////////////////////////////////////////////////////////////////////

		TH2D *h2_MJ = new TH2D("h2_MJ","All muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		// parse jet by species of particles
		TH2D *h2_d_MJ = new TH2D("h2_d_MJ","d muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_u_MJ = new TH2D("h2_u_MJ","u muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_s_MJ = new TH2D("h2_s_MJ","s muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_c_MJ = new TH2D("h2_c_MJ","c muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_b_MJ = new TH2D("h2_b_MJ","b muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_dbar_MJ = new TH2D("h2_dbar_MJ","#bar{d} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ubar_MJ = new TH2D("h2_ubar_MJ","#bar{u} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_sbar_MJ = new TH2D("h2_sbar_MJ","#bar{s} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_cbar_MJ = new TH2D("h2_cbar_MJ","#bar{c} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_bbar_MJ = new TH2D("h2_bbar_MJ","#bar{b} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_g_MJ = new TH2D("h2_g_MJ","gluon muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_hq_MJ = new TH2D("h2_hq_MJ","heavy quark muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_lq_MJ = new TH2D("h2_lq_MJ","light quark muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_qall_MJ = new TH2D("h2_qall_MJ","All quark muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_dall_MJ = new TH2D("h2_dall_MJ","d & #bar{d} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_uall_MJ = new TH2D("h2_uall_MJ","u & #bar{u} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_sall_MJ = new TH2D("h2_sall_MJ","s & #bar{s} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_call_MJ = new TH2D("h2_call_MJ","c & #bar{c} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ball_MJ = new TH2D("h2_ball_MJ","b & #bar{b} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ud_MJ = new TH2D("h2_ud_MJ","u & d muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_ubardbar_MJ = new TH2D("h2_ubardbar_MJ","#bar{u} & #bar{d} muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_other_MJ = new TH2D("h2_other_MJ","ghost muon-tagged-jets",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);

		
	

		///////////////////////////////////////////////////////////////  CALCULATED VARIABLES  //////////////////////////////////////////////////////////////////////////////

		const int NRelPtBins = 50;
		double relPtMin = 0.0;
		double relPtMax = 5.0;


		TH1D *deltaR = new TH1D("deltaR","#Delta r",10,0,1);
		TH1D *h_muRelPt = new TH1D("h_muRelPt","p_{T}^{#mu}",NRelPtBins,relPtMin,relPtMax); // all muon-jets
		TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","p_{T}^{#mu,g}",NRelPtBins,relPtMin,relPtMax); // all gluon muon-jets
		TH1D *h_muRelPt_lq = new TH1D("h_muRelPt_lq","p_{T}^{#mu,light}",NRelPtBins,relPtMin,relPtMax); // all light-quark muon-jets
		TH1D *h_muRelPt_sall = new TH1D("h_muRelPt_sall","p_{T}^{#mu,s}",NRelPtBins,relPtMin,relPtMax); // all s muon-jets
		TH1D *h_muRelPt_call = new TH1D("h_muRelPt_call","p_{T}^{#mu,c}",NRelPtBins,relPtMin,relPtMax); // all c muon-jets
		TH1D *h_muRelPt_ball = new TH1D("h_muRelPt_ball","p_{T}^{#mu,b}",NRelPtBins,relPtMin,relPtMax); // all b muon-jets

		// centrality dependence
		TH1D *h_muRelPt_centReg1 = new TH1D("h_muRelPt_centReg1","",NRelPtBins,relPtMin,relPtMax); // 0-10 % centrality
		TH1D *h_muRelPt_g_centReg1 = new TH1D("h_muRelPt_g_centReg1","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_lq_centReg1 = new TH1D("h_muRelPt_lq_centReg1","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_sall_centReg1 = new TH1D("h_muRelPt_sall_centReg1","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_call_centReg1 = new TH1D("h_muRelPt_call_centReg1","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_ball_centReg1 = new TH1D("h_muRelPt_ball_centReg1","",NRelPtBins,relPtMin,relPtMax); 

		TH1D *h_muRelPt_centReg2 = new TH1D("h_muRelPt_centReg2","",NRelPtBins,relPtMin,relPtMax); // 10-30 % centrality
		TH1D *h_muRelPt_g_centReg2 = new TH1D("h_muRelPt_g_centReg2","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_lq_centReg2 = new TH1D("h_muRelPt_lq_centReg2","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_sall_centReg2 = new TH1D("h_muRelPt_sall_centReg2","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_call_centReg2 = new TH1D("h_muRelPt_call_centReg2","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_ball_centReg2 = new TH1D("h_muRelPt_ball_centReg2","",NRelPtBins,relPtMin,relPtMax);

		TH1D *h_muRelPt_centReg3 = new TH1D("h_muRelPt_centReg3","",NRelPtBins,relPtMin,relPtMax); // 30-50 % centrality
		TH1D *h_muRelPt_g_centReg3 = new TH1D("h_muRelPt_g_centReg3","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_lq_centReg3 = new TH1D("h_muRelPt_lq_centReg3","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_sall_centReg3 = new TH1D("h_muRelPt_sall_centReg3","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_call_centReg3 = new TH1D("h_muRelPt_call_centReg3","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_ball_centReg3 = new TH1D("h_muRelPt_ball_centReg3","",NRelPtBins,relPtMin,relPtMax);

		TH1D *h_muRelPt_centReg4 = new TH1D("h_muRelPt_centReg4","",NRelPtBins,relPtMin,relPtMax); // 50-90 % centrality
		TH1D *h_muRelPt_g_centReg4 = new TH1D("h_muRelPt_g_centReg4","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_lq_centReg4 = new TH1D("h_muRelPt_lq_centReg4","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_sall_centReg4 = new TH1D("h_muRelPt_sall_centReg4","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_call_centReg4 = new TH1D("h_muRelPt_call_centReg4","",NRelPtBins,relPtMin,relPtMax); 
		TH1D *h_muRelPt_ball_centReg4 = new TH1D("h_muRelPt_ball_centReg4","",NRelPtBins,relPtMin,relPtMax);


		///////////////////////////////////////////////////////////////  CALCULATED 2D HISTOGRAMS  ////////////////////////////////////////////////////////////////////////////
		


		
		TH2D *h2_relPtJetPt = new TH2D("h2_relPtJetPt","All jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);
		TH2D *h2_relPtJetPt_g = new TH2D("h2_relPtJetPt_g","All jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);
		TH2D *h2_relPtJetPt_lq = new TH2D("h2_relPtJetPt_lq","gluon jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);
		TH2D *h2_relPtJetPt_sall = new TH2D("h2_relPtJetPt_sall","s jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);
		TH2D *h2_relPtJetPt_call = new TH2D("h2_relPtJetPt_call","c jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);
		TH2D *h2_relPtJetPt_ball = new TH2D("h2_relPtJetPt_ball","b jets",NRelPtBins,relPtMin,relPtMax,NPtBins,ptMin,ptMax);

		TH2D *h2_jetPtPthatWeight = new TH2D("h2_jetPtPthatWeight","",NPtBins,ptMin,ptMax,100,0.0,0.05);

		TH2D *h2_pthatReg1 = new TH2D("h2_pthatReg1","",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_pthatReg2 = new TH2D("h2_pthatReg2","",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		TH2D *h2_pthatReg3 = new TH2D("h2_pthatReg3","",NEtaBins,etaMin,etaMax,NPtBins,ptMin,ptMax);
		
		

		h_vz->Sumw2();
		h_vz_noWeight->Sumw2();
		h_hiBin->Sumw2();
		h_hiBin_noWeight->Sumw2();
		h2_lowPthatHighPtjet->Sumw2();
		h2_noVzOrHiBinWeight->Sumw2();
		h_weight->Sumw2();

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

		h2_MJ->Sumw2();
		h2_d_MJ->Sumw2();
		h2_u_MJ->Sumw2();
		h2_s_MJ->Sumw2();
		h2_c_MJ->Sumw2();
		h2_b_MJ->Sumw2();
		h2_dbar_MJ->Sumw2();
		h2_ubar_MJ->Sumw2();
		h2_sbar_MJ->Sumw2();
		h2_cbar_MJ->Sumw2();
		h2_bbar_MJ->Sumw2();
		h2_hq_MJ->Sumw2();
		h2_lq_MJ->Sumw2();
		h2_qall_MJ->Sumw2();
		h2_dall_MJ->Sumw2();
		h2_uall_MJ->Sumw2();
		h2_sall_MJ->Sumw2();
		h2_call_MJ->Sumw2();
		h2_ball_MJ->Sumw2();
		h2_g_MJ->Sumw2();
		h2_ud_MJ->Sumw2();
		h2_ubardbar_MJ->Sumw2();
		h2_other_MJ -> Sumw2();

		deltaR -> Sumw2();
		h_muRelPt -> Sumw2();
		h_muRelPt_g -> Sumw2();
		h_muRelPt_lq -> Sumw2();
		h_muRelPt_sall -> Sumw2();
		h_muRelPt_ball -> Sumw2();
		h_muRelPt_call -> Sumw2();
		
		h2_relPtJetPt->Sumw2();
		h2_relPtJetPt_g->Sumw2();
		h2_relPtJetPt_lq->Sumw2();
		h2_relPtJetPt_sall->Sumw2();
		h2_relPtJetPt_ball->Sumw2();
		h2_relPtJetPt_call->Sumw2();

		h2_jetPtPthatWeight->Sumw2();

		h2_pthatReg1->Sumw2();
		h2_pthatReg2->Sumw2();
		h2_pthatReg3->Sumw2();
		

		
		
	//////////////////////////////////////////////////////////////////// weighting functions  //////////////////////////////////////////////////////////////////////////////
		TF1 *fvz = new TF1("fvz","pol6",-15,15);
		fvz->SetParameters(1.00656, -0.0193651, 0.000976851, -1.043e-05, -9.79808e-06, 9.07733e-08, 1.79165e-08);

		TF1* fCentralityWeightFunction= new TF1("fcent","pol6",0,90);
		fCentralityWeightFunction->SetParameters(4.64945,-0.201337, 0.00435794,-7.00799e-05,8.18299e-07,-5.52604e-09,1.54472e-11);

	////////////////////////////////////////////////////////////////////////  CUT Values  //////////////////////////////////////////////////////////////////////////////

	
	/////////////////////////////////////////////////////////////////////  LOAD DATA  /////////////////////////////////////////////////////////////////////////////

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



/////////////////////////////////////////////////////////////////////  Begin file loop /////////////////////////////////////////////////////////////////////////////

for(int file = 1; file < NFiles+1; file++){
//for(int file = 1; file < 20; file++){
	if(file==147 || file==179 || file==184 || file==203 || file==238 || file==3 || file==314 || file==327 || file==360){continue;} // missing these root files [skims5]
	cout << "Processing file " << file << "/"<<NFiles<< endl;
	auto f = TFile::Open(Form("/home/clayton/Analysis/skims5/0000/HydJet_%d.root",file));

    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();



	////////////////////////////////////////////////////////////////////// RECO JETS ////////////////////////////////////////////////////////////////////
	
	if(isReco){
		inp_tree->SetBranchAddress("flowpf_jtpt",&pf_jtpt);
		inp_tree->SetBranchAddress("flowpf_jteta",&pf_jteta);
		inp_tree->SetBranchAddress("flowpf_jtphi",&pf_jtphi);
		inp_tree->SetBranchAddress("flowpf_partonFlavor",&pf_partonFlavor);
		inp_tree->SetBranchAddress("flowpf_hadronFlavor",&pf_hadronFlavor);
	}
	
	////////////////////////////////////////////////////////////////////// GEN JETS ////////////////////////////////////////////////////////////////////

	//inp_tree->SetBranchAddress("gen_jtpt",&pf_jtpt);
	//inp_tree->SetBranchAddress("gen_jteta",&pf_jteta);
	//inp_tree->SetBranchAddress("gen_jtphi",&pf_jtphi);

	////////////////////////////////////////////////////////////////////// REF JETS ////////////////////////////////////////////////////////////////////
	if(!isReco){
		inp_tree->SetBranchAddress("ref_jtpt",&pf_jtpt);
		inp_tree->SetBranchAddress("ref_jteta",&pf_jteta);
		inp_tree->SetBranchAddress("ref_jtphi",&pf_jtphi);
		inp_tree->SetBranchAddress("flowpf_partonFlavor",&pf_partonFlavor);
		inp_tree->SetBranchAddress("flowpf_hadronFlavor",&pf_hadronFlavor);
	}
	
	////////////////////////////////////////////////////////////////////// MUON JETS ////////////////////////////////////////////////////////////////////

	inp_tree->SetBranchAddress("muPt",&muPt);
	inp_tree->SetBranchAddress("muEta",&muEta);
	inp_tree->SetBranchAddress("muPhi",&muPhi);
	inp_tree->SetBranchAddress("muChi2NDF",&muChi2NDF);
	inp_tree->SetBranchAddress("muInnerD0",&muInnerD0);
	inp_tree->SetBranchAddress("muInnerDz",&muInnerDz);
	inp_tree->SetBranchAddress("nMu",&nMu);
	inp_tree->SetBranchAddress("muIsGlobal",&muIsGlobal);
	inp_tree->SetBranchAddress("muIsTracker",&muIsTracker);
	inp_tree->SetBranchAddress("muMuonHits",&muMuonHits);
	inp_tree->SetBranchAddress("muStations",&muStations);
	inp_tree->SetBranchAddress("muTrkLayers",&muTrkLayers);
	inp_tree->SetBranchAddress("muPixelHits",&muPixelHits);

	////////////////////////////////////////////////////////////////// EVENT VARIABLES ////////////////////////////////////////////////////////////////////

	inp_tree->SetBranchAddress("weight",&weight);
	inp_tree->SetBranchAddress("vz",&vz);
	inp_tree->SetBranchAddress("hiBin",&hiBin);

	//////////////////////////////////////////////////////////////////  loop variables  ///////////////////////////////////////////////////////////////////////	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	double w=0;

	////////////////////////////////////////////////////////////////////  Event loop  //////////////////////////////////////////////////////////////////////////

	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;
        
        //////////////////////////////////////////////  EVENT CUTS ////////////////////////////////////////////////

        // pthat / weight cut
        
        if(weight>0.044){continue;} // if pthat < 30
        //if(weight>0.0044){continue;} // if pthat < 50
        //if(weight>0.000535){continue;} // if pthat < 80
        h_weight->Fill(weight);

        // vertex cut
        if (fabs(vz)>15.0){continue;}

        // centrality cut
        if(hiBin>180){continue;}



   

     ///////////////////////////////////////////////////////////////////   Calculate weight (pthat)   //////////////////////////////////////////////////////////////////////////
        
	    
		
        
        double w_vz = fvz->Eval(vz);
        double w_hiBin = fCentralityWeightFunction->Eval(hiBin/2.0); // weights derived from dijet-triggered jets 
    
        w = weight*w_vz*w_hiBin;
          
        double reco_lowPthatHighPtjetCut = 100.0;
	double gen_lowPthatHighPtjetCut = 80.0;
        double lowPthatHighPtjetCut = 0.0;
	if(isReco){lowPthatHighPtjetCut = reco_lowPthatHighPtjetCut;}
	if(!isReco){lowPthatHighPtjetCut = gen_lowPthatHighPtjetCut;}

		//////////////////////////////////////////////////////////////////  Jet loop ////////////////////////////////////////////////////////////////////////////
		
		for(int jetj=0; jetj < (int)pf_jteta->size(); jetj++){
			
			double jetPtj = pf_jtpt->at(jetj);
			double jetEtaj = pf_jteta->at(jetj);
			double jetPhij = pf_jtphi->at(jetj);
			
			double x = pf_jteta->at(jetj);
			double u = pf_jtpt->at(jetj);
			
			/////////////////////////////////////  JET CUTS /////////////////////////////////////////
			if(fabs(x)>etaMax || u < ptMin || u>ptMax || u==-999){continue;}
			
			h2_jetPtPthatWeight->Fill(u,weight);

			

			if(weight>0.0044){
				h2_pthatReg1->Fill(x,u,w);
				if(u>lowPthatHighPtjetCut){
					h2_lowPthatHighPtjet->Fill(x,u,w);
					continue;
				}
			}
			if(weight<0.0044 && weight>0.000535){
				h2_pthatReg2->Fill(x,u,w);
			}
			if(weight<0.000535){
				h2_pthatReg3->Fill(x,u,w);
			}

			

			// weighting investigation

			h_vz->Fill(vz,w);
			h_vz_noWeight->Fill(vz,weight*w_hiBin);

			h_hiBin->Fill(hiBin,w);
			h_hiBin_noWeight->Fill(hiBin,weight*w_vz);




			int z = pf_hadronFlavor->at(jetj);
			int y = pf_partonFlavor->at(jetj);

			if(z==5){
				h2_ball->Fill(x,u,w);
			}
			else if(z==4){
				h2_call->Fill(x,u,w);
			}
			else{
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
			}
			
			
			

			h2->Fill(x,u,w);
			h2_noVzOrHiBinWeight->Fill(x,u,weight);
			
			//////////////////////////////////////////////////////////////////  Muon loop ////////////////////////////////////////////////////////////////////////////

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

				if(z==5){
					h2_ball_MJ->Fill(x,u,w); h_muRelPt_ball->Fill(muRelPtMin);
				}
				else if(z==4){
					h2_call_MJ->Fill(x,u,w); h_muRelPt_call->Fill(muRelPtMin);
				}
				else{
					if(y==0){h2_other_MJ -> Fill(x,u,w); }
					if(y==1){h2_d_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_dall_MJ->Fill(x,u,w); h2_ud_MJ->Fill(x,u,w); h_muRelPt_lq->Fill(muRelPtMin);}
					if(y==2){h2_u_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_uall_MJ->Fill(x,u,w); h2_ud_MJ->Fill(x,u,w); h_muRelPt_lq->Fill(muRelPtMin);}
					if(y==3){h2_s_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_sall_MJ->Fill(x,u,w); h_muRelPt_sall->Fill(muRelPtMin);}
					if(y==4){h2_c_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_call_MJ->Fill(x,u,w); h_muRelPt_call->Fill(muRelPtMin);}
					if(y==5){h2_b_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_ball_MJ->Fill(x,u,w); h_muRelPt_ball->Fill(muRelPtMin);}
					if(y==-1){h2_dbar_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_dall_MJ->Fill(x,u,w); h2_ubardbar_MJ->Fill(x,u,w); h_muRelPt_lq->Fill(muRelPtMin);}
					if(y==-2){h2_ubar_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_uall_MJ->Fill(x,u,w); h2_ubardbar_MJ->Fill(x,u,w); h_muRelPt_lq->Fill(muRelPtMin);}
					if(y==-3){h2_sbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_sall_MJ->Fill(x,u,w); h_muRelPt_sall->Fill(muRelPtMin);}
					if(y==-4){h2_cbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_call_MJ->Fill(x,u,w); h_muRelPt_call->Fill(muRelPtMin);}
					if(y==-5){h2_bbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_ball_MJ->Fill(x,u,w); h_muRelPt_ball->Fill(muRelPtMin);}
					if(y==21){h2_g_MJ->Fill(x,u,w); h_muRelPt_g->Fill(muRelPtMin);}
				}
					
				

				h2_MJ->Fill(x,u,w);
				h_muRelPt->Fill(muRelPtMin);
				if(hiBin>0 && hiBin<20) {
					h_muRelPt_centReg1->Fill(muRelPtMin);
					if(abs(y)==1 || abs(y)==2){h_muRelPt_lq_centReg1->Fill(muRelPtMin);}
					if(abs(y)==3){h_muRelPt_sall_centReg1->Fill(muRelPtMin);}
					if(z==4 || abs(y)==4){h_muRelPt_call_centReg1->Fill(muRelPtMin);}
					if(z==5 || abs(y)==5){h_muRelPt_ball_centReg1->Fill(muRelPtMin);}
					if(y==21){h_muRelPt_g_centReg1->Fill(muRelPtMin);}
				}
				if(hiBin>20 && hiBin<60){
					h_muRelPt_centReg2->Fill(muRelPtMin);
					if(abs(y)==1 || abs(y)==2){h_muRelPt_lq_centReg2->Fill(muRelPtMin);}
					if(abs(y)==3){h_muRelPt_sall_centReg2->Fill(muRelPtMin);}
					if(z==4 || abs(y)==4){h_muRelPt_call_centReg2->Fill(muRelPtMin);}
					if(z==5 || abs(y)==5){h_muRelPt_ball_centReg2->Fill(muRelPtMin);}
					if(y==21){h_muRelPt_g_centReg2->Fill(muRelPtMin);}
				} 
				if(hiBin>60 && hiBin<100){
					h_muRelPt_centReg3->Fill(muRelPtMin);
					if(abs(y)==1 || abs(y)==2){h_muRelPt_lq_centReg3->Fill(muRelPtMin);}
					if(abs(y)==3){h_muRelPt_sall_centReg3->Fill(muRelPtMin);}
					if(z==4 || abs(y)==4){h_muRelPt_call_centReg3->Fill(muRelPtMin);}
					if(z==5 || abs(y)==5){h_muRelPt_ball_centReg3->Fill(muRelPtMin);}
					if(y==21){h_muRelPt_g_centReg3->Fill(muRelPtMin);}
				} 
				if(hiBin>100 && hiBin<180){
					h_muRelPt_centReg4->Fill(muRelPtMin);
					if(abs(y)==1 || abs(y)==2){h_muRelPt_lq_centReg4->Fill(muRelPtMin);}
					if(abs(y)==3){h_muRelPt_sall_centReg4->Fill(muRelPtMin);}
					if(z==4 || abs(y)==4){h_muRelPt_call_centReg4->Fill(muRelPtMin);}
					if(z==5 || abs(y)==5){h_muRelPt_ball_centReg4->Fill(muRelPtMin);}
					if(y==21){h_muRelPt_g_centReg4->Fill(muRelPtMin);}
				}
				
				h2_relPtJetPt->Fill(muRelPtMin,u,w);
				
				if(abs(y)==1 || abs(y)==2){h2_relPtJetPt_lq->Fill(muRelPtMin,u,w);}
				if(abs(y)==3){h2_relPtJetPt_sall->Fill(muRelPtMin,u,w);}
				if(z==4 || abs(y)==4){h2_relPtJetPt_call->Fill(muRelPtMin,u,w);}
				if(z==5 || abs(y)==5){h2_relPtJetPt_ball->Fill(muRelPtMin,u,w);}
				if(y==21){h2_relPtJetPt_g->Fill(muRelPtMin,u,w);}

			}
			
		} 
		//////////////////////////////////////////////////////////////////  End jet loop ////////////////////////////////////////////////////////////////////////////
		
	} 
	//////////////////////////////////////////////////////////////////  End event loop ////////////////////////////////////////////////////////////////////////////

}
/////////////////////////////////////////////////////////////////////  End file loop /////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////  file creation history /////////////////////////////////////////////////////

//auto wf = TFile::Open("makeDataV1_refjets_pthat_80_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV1_recojets_pthat_80_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV1_refjets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV1_recojets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV1_refjets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV1_recojets_pthat_30_muptcut_5.root","recreate");

// mu pt cut dependence study, 6/29/20

//auto wf = TFile::Open("makeDataV1_refjets_pthat_50_muptcut_10.root","recreate");
//auto wf = TFile::Open("makeDataV1_refjets_pthat_50_muptcut_15.root","recreate");
//auto wf = TFile::Open("makeDataV1_refjets_pthat_50_muptcut_20.root","recreate");

// change to version 3
//auto wf = TFile::Open("makeDataV3_recojets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV3_recojets_pthat_50_muptcut_10.root","recreate");

// mu pt cut dependence re-study, 7/6/20
//auto wf = TFile::Open("makeDataV3_refjets_pthat_50_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV3_refjets_pthat_50_muptcut_10.root","recreate");
//auto wf = TFile::Open("makeDataV3_refjets_pthat_50_muptcut_15.root","recreate");
//auto wf = TFile::Open("makeDataV3_refjets_pthat_50_muptcut_20.root","recreate");

// change to version 5
//auto wf = TFile::Open("makeDataV5_recojets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV5_refjets_pthat_30_muptcut_5.root","recreate");

// change to version 6
//auto wf = TFile::Open("makeDataV6_recojets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("makeDataV6_refjets_pthat_30_muptcut_5.root","recreate");

// change to version 7
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_recojets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_refjets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_recojets_pthat_30_muptcut_5_pthatTestLogic.root","recreate");

// change to version 8
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_refjets_pthat_30_muptcut_5.root","recreate");
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root","recreate");

// change to version 9
//auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_refjets_pthat_30_muptcut_5.root","recreate");
auto wf = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_recojets_pthat_30_muptcut_5.root","recreate");

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
h2_MJ->Write();
h2_d_MJ->Write();
h2_u_MJ->Write();
h2_s_MJ->Write();
h2_c_MJ->Write();
h2_b_MJ->Write();
h2_dbar_MJ->Write();
h2_ubar_MJ->Write();
h2_sbar_MJ->Write();
h2_cbar_MJ->Write();
h2_bbar_MJ->Write();
h2_hq_MJ->Write();
h2_lq_MJ->Write();
h2_qall_MJ->Write();
h2_dall_MJ->Write();
h2_uall_MJ->Write();
h2_sall_MJ->Write();
h2_call_MJ->Write();
h2_ball_MJ->Write();
h2_g_MJ->Write();
h2_ud_MJ->Write();
h2_ubardbar_MJ->Write();
h2_other_MJ->Write();
deltaR->Write();
h_muRelPt->Write();
h_muRelPt_g->Write();
h_muRelPt_lq->Write();
h_muRelPt_sall->Write();
h_muRelPt_call->Write();
h_muRelPt_ball->Write();
h2_relPtJetPt->Write();
h2_relPtJetPt_g->Write();
h2_relPtJetPt_lq->Write();
h2_relPtJetPt_ball->Write();
h2_relPtJetPt_call->Write();
h2_relPtJetPt_sall->Write();
h_vz->Write();
h_vz_noWeight->Write();
h_hiBin->Write();
h_hiBin_noWeight->Write();
h2_lowPthatHighPtjet->Write();
h2_noVzOrHiBinWeight->Write();
// centrality dependence
h_muRelPt_centReg1->Write();
h_muRelPt_g_centReg1->Write();
h_muRelPt_lq_centReg1->Write();
h_muRelPt_sall_centReg1->Write();
h_muRelPt_ball_centReg1->Write();
h_muRelPt_call_centReg1->Write();

h_muRelPt_centReg2->Write();
h_muRelPt_g_centReg2->Write();
h_muRelPt_lq_centReg2->Write();
h_muRelPt_sall_centReg2->Write();
h_muRelPt_ball_centReg2->Write();
h_muRelPt_call_centReg2->Write();

h_muRelPt_centReg3->Write();
h_muRelPt_g_centReg3->Write();
h_muRelPt_lq_centReg3->Write();
h_muRelPt_sall_centReg3->Write();
h_muRelPt_ball_centReg3->Write();
h_muRelPt_call_centReg3->Write();

h_muRelPt_centReg4->Write();
h_muRelPt_g_centReg4->Write();
h_muRelPt_lq_centReg4->Write();
h_muRelPt_sall_centReg4->Write();
h_muRelPt_ball_centReg4->Write();
h_muRelPt_call_centReg4->Write();

h_weight->Write();
h2_jetPtPthatWeight->Write();

h2_pthatReg1->Write();
h2_pthatReg2->Write();
h2_pthatReg3->Write();

wf->Close();

return;

}


 // end program


