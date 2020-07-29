
// V4: using variable binning
// V5: muon jet filtering, using HydJet_x.root files
// V6: using new weigthing method, pthat.  No longer reweighting according to bin size (delegate to plotting macros).  Including muRelPt calculation
// V7: desinged to run /skims5/ files. Added histogram of untagged jets.
// V8: no more variable binning, going to delegate to relevant plotting macros


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




Double_t getPtRel(Double_t MuonPt, Double_t MuonEta, Double_t MuonPhi, Double_t JetPt, Double_t JetEta, Double_t JetPhi)
{

Double_t Muon_Px = MuonPt*TMath::Cos(MuonPhi);
Double_t Muon_Py = MuonPt*TMath::Sin(MuonPhi);
Double_t Muon_Pz = MuonPt*TMath::SinH(MuonEta);

Double_t Jet_Px = JetPt*TMath::Cos(JetPhi);
Double_t Jet_Py = JetPt*TMath::Sin(JetPhi);
Double_t Jet_Pz = JetPt*TMath::SinH(JetEta);


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



void etaPtHistoV8(Double_t muPtCut = 2.5, Double_t jetptcut = 120.){

	//////////////////////////////////////////////////////////////////////// jet variables ////////////////////////////////////////////////////////////////////////
	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_jtphi=0, *pf_partonFlavor=0; 
	//////////////////////////////////////////////////////////////////////// muon jet variables ////////////////////////////////////////////////////////////////////////
	vector <float> *muPt=0, *muEta=0, *muPhi=0, *muChi2NDF=0, *muInnerD0=0, *muInnerDz=0;
	int nMu;
	vector <int> *muIsGlobal=0, *muIsTracker=0, *muMuonHits=0, *muStations=0, *muTrkLayers=0, *muPixelHits=0;
	//////////////////////////////////////////////////////////////////////// event variables ////////////////////////////////////////////////////////////////////////
	Float_t weight;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////// DEFINE BINNING ///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////// create eta edges array /////////////////////////////////////////////////////////////////////////
	const Int_t Nbins = 1000;
	


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////// DEFINE HISTOGRAMS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const Int_t numUsedEtaBins = 15;
	Float_t usedEtaEdges[numUsedEtaBins+1] = {-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
	const Int_t numUsedPtBins = 20;
	Float_t usedPtEdges[numUsedPtBins+1] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360, 380, 400, 432, 500};

		////////////////////////////////////////// /////////////////////////////  JETS  ///////////////////////////////////////////////////////////////////////////////////
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

		
		/////////////////////////////////////////////////////////////////   MUON-TAGGED JETS  //////////////////////////////////////////////////////////////////////////////
		TH2D *h2_MJ = new TH2D("h2_MJ","All muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		// parse jet by species of particles
		TH2D *h2_d_MJ = new TH2D("h2_d_MJ","d muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_u_MJ = new TH2D("h2_u_MJ","u muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_s_MJ = new TH2D("h2_s_MJ","s muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_c_MJ = new TH2D("h2_c_MJ","c muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_b_MJ = new TH2D("h2_b_MJ","b muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dbar_MJ = new TH2D("h2_dbar_MJ","#bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubar_MJ = new TH2D("h2_ubar_MJ","#bar{u} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sbar_MJ = new TH2D("h2_sbar_MJ","#bar{s} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_cbar_MJ = new TH2D("h2_cbar_MJ","#bar{c} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_bbar_MJ = new TH2D("h2_bbar_MJ","#bar{b} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_g_MJ = new TH2D("h2_g_MJ","gluon muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_hq_MJ = new TH2D("h2_hq_MJ","heavy quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_lq_MJ = new TH2D("h2_lq_MJ","light quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_qall_MJ = new TH2D("h2_qall_MJ","All quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dall_MJ = new TH2D("h2_dall_MJ","d & #bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_uall_MJ = new TH2D("h2_uall_MJ","u & #bar{u} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sall_MJ = new TH2D("h2_sall_MJ","s & #bar{s} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_call_MJ = new TH2D("h2_call_MJ","c & #bar{c} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ball_MJ = new TH2D("h2_ball_MJ","b & #bar{b} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ud_MJ = new TH2D("h2_ud_MJ","u & d muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubardbar_MJ = new TH2D("h2_ubardbar_MJ","#bar{u} & #bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_other_MJ = new TH2D("h2_other_MJ","ghost muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);

		/////////////////////////////////////////////////////////////////   UNTAGGED JETS  //////////////////////////////////////////////////////////////////////////////
		TH2D *h2_untagged = new TH2D("h2_untagged","untagged jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_untagged_g = new TH2D("h2_untagged_g","untagged gluon jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_untagged_lq = new TH2D("h2_untagged_lq","untagged light-quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_untagged_sall = new TH2D("h2_untagged_sall","untagged s jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_untagged_call = new TH2D("h2_untagged_call","untagged c jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_untagged_ball = new TH2D("h2_untagged_ball","untagged b jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);


		///////////////////////////////////////////////////////////////  CALCULATED VARIABLES  //////////////////////////////////////////////////////////////////////////////
		TH1D *deltaR = new TH1D("deltaR","#Delta r",10,0,1);
		TH1D *h_muRelPt = new TH1D("h_muRelPt","p_{T}^{#mu}",50,0,5); // all muon-jets
		TH1D *h_muRelPt_g = new TH1D("h_muRelPt_g","p_{T}^{#mu,g}",50,0,5); // all gluon muon-jets
		TH1D *h_muRelPt_lq = new TH1D("h_muRelPt_lq","p_{T}^{#mu,light}",50,0,5); // all light-quark muon-jets
		TH1D *h_muRelPt_sall = new TH1D("h_muRelPt_sall","p_{T}^{#mu,s}",50,0,5); // all s muon-jets
		TH1D *h_muRelPt_call = new TH1D("h_muRelPt_call","p_{T}^{#mu,c}",50,0,5); // all c muon-jets
		TH1D *h_muRelPt_ball = new TH1D("h_muRelPt_ball","p_{T}^{#mu,b}",50,0,5); // all b muon-jets

		///////////////////////////////////////////////////////////////  CALCULATED 2D HISTOGRAMS  ////////////////////////////////////////////////////////////////////////////
		
		const Int_t numUsedRelPtBins = 50;
		Float_t usedRelPtEdges[numUsedRelPtBins+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
			3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0};

		
		TH2D *h2_relPtJetPt = new TH2D("h2_relPtJetPt","All jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_relPtJetPt_g = new TH2D("h2_relPtJetPt_g","All jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_relPtJetPt_lq = new TH2D("h2_relPtJetPt_lq","gluon jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_relPtJetPt_sall = new TH2D("h2_relPtJetPt_sall","s jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_relPtJetPt_call = new TH2D("h2_relPtJetPt_call","c jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_relPtJetPt_ball = new TH2D("h2_relPtJetPt_ball","b jets",numUsedRelPtBins,usedRelPtEdges,numUsedPtBins,usedPtEdges);
		


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

		h2_untagged->Sumw2();
		h2_untagged_g->Sumw2();
		h2_untagged_lq->Sumw2();
		h2_untagged_sall->Sumw2();
		h2_untagged_ball->Sumw2();
		h2_untagged_call->Sumw2();

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
		
	
	////////////////////////////////////////////////////////////////////////  CUTS  /////////////////////////////////////////////////////////////////////////////
	const double etamaxcut = 1.6; // default = 1.5
	const double etamincut = 0.;
	const double pTmincut = jetptcut; // default = 120.
	const double pTmaxcut = 500;
	


	/////////////////////////////////////////////////////////////////////  LOAD DATA  /////////////////////////////////////////////////////////////////////////////
	int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/skims5/0000");
	//DIR *dr = opendir("/home/clayton/Analysis/skims4/0000");  

   	if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    	{ 
        	printf("Could not open current directory" ); 
        	return 0; 
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
	
	//inp_tree->SetBranchAddress("flowpf_jtpt",&pf_jtpt);
	//inp_tree->SetBranchAddress("flowpf_jteta",&pf_jteta);
	//inp_tree->SetBranchAddress("flowpf_jtphi",&pf_jtphi);
	//inp_tree->SetBranchAddress("flowpf_partonFlavor",&pf_partonFlavor);
	
	////////////////////////////////////////////////////////////////////// GEN JETS ////////////////////////////////////////////////////////////////////

	//inp_tree->SetBranchAddress("gen_jtpt",&pf_jtpt);
	//inp_tree->SetBranchAddress("gen_jteta",&pf_jteta);
	//inp_tree->SetBranchAddress("gen_jtphi",&pf_jtphi);

	////////////////////////////////////////////////////////////////////// REF JETS ////////////////////////////////////////////////////////////////////
	
	inp_tree->SetBranchAddress("ref_jtpt",&pf_jtpt);
	inp_tree->SetBranchAddress("ref_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("ref_jtphi",&pf_jtphi);
	inp_tree->SetBranchAddress("flowpf_partonFlavor",&pf_partonFlavor);
	
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


	//////////////////////////////////////////////////////////////////  loop variables  ///////////////////////////////////////////////////////////////////////	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	double w=0;
	
	int nMuCount=0;

	////////////////////////////////////////////////////////////////////  Event loop  //////////////////////////////////////////////////////////////////////////

	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;
        
        w=weight;
        if(w>0.044){continue;} // pthat = 30
        //if(w>0.0044){continue;} // pthat = 50
        
         
		//////////////////////////////////////////////////////////////////  Jet loop ////////////////////////////////////////////////////////////////////////////
		
		for(int jetj=0; jetj < (int)pf_jteta->size(); jetj++){
			
			Double_t jetPtj = pf_jtpt->at(jetj);
			Double_t jetEtaj = pf_jteta->at(jetj);
			Double_t jetPhij = pf_jtphi->at(jetj);

			// convert eta to theta
			Double_t jetThetaj = 2*atan(exp(-jetEtaj));
			
			double x = pf_jteta->at(jetj);
			double u = pf_jtpt->at(jetj);
			
			if(fabs(x)>etamaxcut || u < pTmincut || u>pTmaxcut || u==-999){continue;}
			int y = pf_partonFlavor->at(jetj);
			
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
			
			//////////////////////////////////////////////////////////////////  Muon loop ////////////////////////////////////////////////////////////////////////////
			
			if(nMu==0){continue;}
			Double_t deltaRmin=100;
			Double_t muRelPtMin = 0;
			for(int mui=0; mui<(int) muPt->size();mui++){
				if(muIsTracker->at(mui)==0 || TMath::Abs(muEta->at(mui))>2.4 || muChi2NDF->at(mui)==-99 || muChi2NDF->at(mui)>10
					||TMath::Abs(muInnerD0->at(mui))>0.2 || TMath::Abs(muInnerDz->at(mui))>0.5 || muMuonHits->at(mui)<= 0
					|| muStations->at(mui)<= 1 || muTrkLayers->at(mui)<=5 || muPixelHits->at(mui)<=0 || muPt->at(mui) < muPtCut){continue;}
			
				Double_t muPti = muPt->at(mui);
				Double_t muEtai = muEta->at(mui);
				Double_t muPhii = muPhi->at(mui);

				// convert from eta to theta
				Double_t muThetai = 2*atan(exp(-muEtai));
				Double_t muDotJ = sin(jetThetaj)*sin(muThetai)*cos(jetPhij)*cos(muPhii) + sin(jetThetaj)*sin(muThetai)*sin(jetPhij)*sin(muPhii)+cos(jetThetaj)*cos(muThetai);
		
							
				Double_t deltaEtaij = muEtai-jetEtaj;
				Double_t deltaPhiij = acos(cos(muPhii-jetPhij));
				Double_t deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				Double_t muRelPt = getPtRel(muPti,muEtai,muPhii,jetPtj,jetEtaj,jetPhij);
				
				if(deltaRij<deltaRmin){
					deltaRmin=deltaRij;
					muRelPtMin=muRelPt;
					
				}

			} 
			//////////////////////////////////////////////////////////////////  End muon loop ////////////////////////////////////////////////////////////////////////////
	
			deltaR->Fill(deltaRmin);
			
			if(deltaRmin<0.4){
				
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
				

				h2_MJ->Fill(x,u,w);
				h_muRelPt->Fill(muRelPtMin);

				
				h2_relPtJetPt->Fill(muRelPtMin,u,w);
				
				if(abs(y)==1 || abs(y)==2){h2_relPtJetPt_lq->Fill(muRelPtMin,u,w);}
				if(abs(y)==3){h2_relPtJetPt_sall->Fill(muRelPtMin,u,w);}
				if(abs(y)==4){h2_relPtJetPt_call->Fill(muRelPtMin,u,w);}
				if(abs(y)==5){h2_relPtJetPt_ball->Fill(muRelPtMin,u,w);}
				if(y==21){h2_relPtJetPt_g->Fill(muRelPtMin,u,w);}
				
				
				
			

			}
			else{
				h2_untagged->Fill(x,u,w);
				
				if(abs(y)==1 || abs(y)==2){h2_untagged_lq->Fill(x,u,w);}
				if(abs(y)==3){h2_untagged_sall->Fill(x,u,w);}
				if(abs(y)==4){h2_untagged_call->Fill(x,u,w);}
				if(abs(y)==5){h2_untagged_ball->Fill(x,u,w);}
				if(y==21){h2_untagged_g->Fill(x,u,w);}
				
			}
		} 
		//////////////////////////////////////////////////////////////////  End jet loop ////////////////////////////////////////////////////////////////////////////
		
	} 
	//////////////////////////////////////////////////////////////////  End event loop ////////////////////////////////////////////////////////////////////////////

}
/////////////////////////////////////////////////////////////////////  End file loop /////////////////////////////////////////////////////////////////////////////




//auto wf = TFile::Open("etaPtHistoV8_genjets_pthat_50_jetptcut_120.root", "recreate");
auto wf = TFile::Open("etaPtHistoV8_refjets_pthat_30_jetptcut_120_muptcut_2-5.root", "recreate");


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
h2_untagged->Write();
h2_untagged_g->Write();
h2_untagged_lq->Write();
h2_untagged_sall->Write();
h2_untagged_ball->Write();
h2_untagged_call->Write();
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
wf->Close();



}


 // end program


