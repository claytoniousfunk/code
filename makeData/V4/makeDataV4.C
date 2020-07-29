
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
// V3 : added reweighting of vz and hibin and overall weight.
// V4 : pthat region study  


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



void makeDataV4(bool isReco = 0, double muPtCut = 5.0, double jetptcut = 50.0){

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

	const int numUsedEtaBins = 15;
	float usedEtaEdges[numUsedEtaBins+1] = {-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
	const int numUsedPtBins = 15;
	float usedPtEdges[numUsedPtBins+1] = {50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////// DEFINE HISTOGRAMS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	

		////////////////////////////////////////// /////////////////////////////  JETS  ///////////////////////////////////////////////////////////////////////////////////

		///////////////// pthat regions
		TH2D *h2_pthatreg1 = new TH2D("h2_pthatreg1","",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_pthatreg2 = new TH2D("h2_pthatreg2","",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_pthatreg3 = new TH2D("h2_pthatreg3","",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);

		

		h2_pthatreg1->Sumw2();
		h2_pthatreg2->Sumw2();
		h2_pthatreg3->Sumw2();

		
		
		
	//////////////////////////////////////////////////////////////////// weighting functions  //////////////////////////////////////////////////////////////////////////////
		TF1 *fvz = new TF1("fvz","pol6",-15,15);
		fvz->SetParameters(1.00656, -0.0193651, 0.000976851, -1.043e-05, -9.79808e-06, 9.07733e-08, 1.79165e-08);

		TF1* fCentralityWeightFunction= new TF1("fcent","pol6",0,90);
		fCentralityWeightFunction->SetParameters(4.64945,-0.201337, 0.00435794,-7.00799e-05,8.18299e-07,-5.52604e-09,1.54472e-11);

	////////////////////////////////////////////////////////////////////////  CUT Values  //////////////////////////////////////////////////////////////////////////////

	const double etamaxcut = 1.5; // default = 1.5
	const double pTmaxcut = 500;
	
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
	
	int nMuCount=0;

	////////////////////////////////////////////////////////////////////  Event loop  //////////////////////////////////////////////////////////////////////////

	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;



        
        //////////////////////////////////////////////  EVENT CUTS ////////////////////////////////////////////////

        // pthat / weight cut

        if(weight>0.044){continue;} // pthat = 30
        //if(weight>0.0044){continue;} // pthat = 50
        //if(weight>0.000535){continue;} // pthat = 80

        // vertex cut
        if (fabs(vz)>15){continue;}

        // centrality cut
        if(hiBin>180){continue;}



   

     ///////////////////////////////////////////////////////////////////   Calculate weight (pthat)   //////////////////////////////////////////////////////////////////////////
        
	    
		
        
        double w1 = fvz->Eval(vz);
        double w2 = fCentralityWeightFunction->Eval(hiBin/2.0);
    
        w = weight*w1*w2;
          
        
         
		//////////////////////////////////////////////////////////////////  Jet loop ////////////////////////////////////////////////////////////////////////////
		
		for(int jetj=0; jetj < (int)pf_jteta->size(); jetj++){
			
			double jetPtj = pf_jtpt->at(jetj);
			double jetEtaj = pf_jteta->at(jetj);
			double jetPhij = pf_jtphi->at(jetj);

			// convert eta to theta
			double jetThetaj = 2*atan(exp(-jetEtaj));
			
			double x = pf_jteta->at(jetj);
			double u = pf_jtpt->at(jetj);
			
			/////////////////////////////////////  JET CUTS /////////////////////////////////////////
			if(fabs(x)>etamaxcut || u < jetptcut || u>pTmaxcut || u==-999){continue;}

			if(weight<0.044 && weight>0.0044) h2_pthatreg1->Fill(x,u,w);
			if(weight<0.0044 && weight>0.000535) h2_pthatreg2->Fill(x,u,w);
			if(weight<0.000535) h2_pthatreg3->Fill(x,u,w);


			
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

// pthat region study, 7/6/20
auto wf = TFile::Open("makeDataV4_refjets_pthatregions.root","recreate");








h2_pthatreg1->Write();
h2_pthatreg2->Write();
h2_pthatreg3->Write();
wf->Close();

return;

}


 // end program


