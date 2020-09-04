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


void PbPb_data_stitch(){

    // histograms to be filled with all the data (stitched)
    TH1D *h_jetpt, *h_jetpt_cent0to10, *h_jetpt_cent10to30, *h_jetpt_cent30to50, *h_jetpt_cent50to90;
    TH1D *h_jeteta, *h_jeteta_cent0to10, *h_jeteta_cent10to30, *h_jeteta_cent30to50, *h_jeteta_cent50to90;
    TH1D *h_jetphi, *h_jetphi_cent0to10, *h_jetphi_cent10to30, *h_jetphi_cent30to50, *h_jetphi_cent50to90;
    TH1D *h_mujetpt, *h_mujetpt_cent0to10, *h_mujetpt_cent10to30, *h_mujetpt_cent30to50, *h_mujetpt_cent50to90;
    TH1D *h_mujeteta, *h_mujeteta_cent0to10, *h_mujeteta_cent10to30, *h_mujeteta_cent30to50, *h_mujeteta_cent50to90;
    TH1D *h_mujetphi, *h_mujetphi_cent0to10, *h_mujetphi_cent10to30, *h_mujetphi_cent30to50, *h_mujetphi_cent50to90;
    TH1D *h_muPt, *h_muPt_cent0to10, *h_muPt_cent10to30, *h_muPt_cent30to50, *h_muPt_cent50to90;
    TH1D *h_muEta, *h_muEta_cent0to10, *h_muEta_cent10to30, *h_muEta_cent30to50, *h_muEta_cent50to90;
    TH1D *h_muPhi, *h_muPhi_cent0to10, *h_muPhi_cent10to30, *h_muPhi_cent30to50, *h_muPhi_cent50to90;
    TH1D *h_muRelPt, *h_muRelPt_cent0to10, *h_muRelPt_cent10to30, *h_muRelPt_cent30to50, *h_muRelPt_cent50to90;
    TH1D *h_vz, *h_hiBin, *h_deltaR;
    // temporary histograms for each file in file loop
    TH1D *t_jetpt, *t_jetpt_cent0to10, *t_jetpt_cent10to30, *t_jetpt_cent30to50, *t_jetpt_cent50to90;
    TH1D *t_jeteta, *t_jeteta_cent0to10, *t_jeteta_cent10to30, *t_jeteta_cent30to50, *t_jeteta_cent50to90;
    TH1D *t_jetphi, *t_jetphi_cent0to10, *t_jetphi_cent10to30, *t_jetphi_cent30to50, *t_jetphi_cent50to90;
    TH1D *t_mujetpt, *t_mujetpt_cent0to10, *t_mujetpt_cent10to30, *t_mujetpt_cent30to50, *t_mujetpt_cent50to90;
    TH1D *t_mujeteta, *t_mujeteta_cent0to10, *t_mujeteta_cent10to30, *t_mujeteta_cent30to50, *t_mujeteta_cent50to90;
    TH1D *t_mujetphi, *t_mujetphi_cent0to10, *t_mujetphi_cent10to30, *t_mujetphi_cent30to50, *t_mujetphi_cent50to90;
    TH1D *t_muPt, *t_muPt_cent0to10, *t_muPt_cent10to30, *t_muPt_cent30to50, *t_muPt_cent50to90;
    TH1D *t_muEta, *t_muEta_cent0to10, *t_muEta_cent10to30, *t_muEta_cent30to50, *t_muEta_cent50to90;
    TH1D *t_muPhi, *t_muPhi_cent0to10, *t_muPhi_cent10to30, *t_muPhi_cent30to50, *t_muPhi_cent50to90;
    TH1D *t_muRelPt, *t_muRelPt_cent0to10, *t_muRelPt_cent10to30, *t_muRelPt_cent30to50, *t_muRelPt_cent50to90;
    TH1D *t_vz, *t_hiBin, *t_deltaR;
    
    // count number of files in the given directory
    int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/data/PbPbDataSkim_27Aug20/data/");

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

    // begin file loop
    for(int file=0;file<NFiles;file++){
        cout << "Processing file " << file+1 << "/"<<NFiles<< endl;
        TFile *f = TFile::Open(Form("/home/clayton/Analysis/data/PbPbDataSkim_27Aug20/data/job_output_part%d.root",file));

        f->GetObject("h_jetpt",t_jetpt);
        /*
            f->GetObject("h_jetpt_cent0to10",t_jetpt_cent0to10);
            f->GetObject("h_jetpt_cent10to30",t_jetpt_cent10to30);
            f->GetObject("h_jetpt_cent30to50",t_jetpt_cent30to50);
            f->GetObject("h_jetpt_cent50to90",t_jetpt_cent50to90);
        f->GetObject("h_mujetpt",t_mujetpt);
            f->GetObject("h_mujetpt_cent0to10",t_mujetpt_cent0to10);
            f->GetObject("h_mujetpt_cent10to30",t_mujetpt_cent10to30);
            f->GetObject("h_mujetpt_cent30to50",t_mujetpt_cent30to50);
            f->GetObject("h_mujetpt_cent50to90",t_mujetpt_cent50to90);
        f->GetObject("h_jeteta",t_jeteta);
            f->GetObject("h_jeteta_cent0to10",t_jeteta_cent0to10);
            f->GetObject("h_jeteta_cent10to30",t_jeteta_cent10to30);
            f->GetObject("h_jeteta_cent30to50",t_jeteta_cent30to50);
            f->GetObject("h_jeteta_cent50to90",t_jeteta_cent50to90);
        f->GetObject("h_mujeteta",t_mujeteta);
            f->GetObject("h_mujeteta_cent0to10",t_mujeteta_cent0to10);
            f->GetObject("h_mujeteta_cent10to30",t_mujeteta_cent10to30);
            f->GetObject("h_mujeteta_cent30to50",t_mujeteta_cent30to50);
            f->GetObject("h_mujeteta_cent50to90",t_mujeteta_cent50to90);
        f->GetObject("h_jetphi",t_jetphi);
            f->GetObject("h_jetphi_cent0to10",t_jetphi_cent0to10);
            f->GetObject("h_jetphi_cent10to30",t_jetphi_cent10to30);
            f->GetObject("h_jetphi_cent30to50",t_jetphi_cent30to50);
            f->GetObject("h_jetphi_cent50to90",t_jetphi_cent50to90);
        f->GetObject("h_mujetphi",t_mujetphi);
            f->GetObject("h_mujetphi_cent0to10",t_mujetphi_cent0to10);
            f->GetObject("h_mujetphi_cent10to30",t_mujetphi_cent10to30);
            f->GetObject("h_mujetphi_cent30to50",t_mujetphi_cent30to50);
            f->GetObject("h_mujetphi_cent50to90",t_mujetphi_cent50to90);
        f->GetObject("h_muPt",t_muPt);
            f->GetObject("h_muPt_cent0to10",t_muPt_cent0to10);
            f->GetObject("h_muPt_cent10to30",t_muPt_cent10to30);
            f->GetObject("h_muPt_cent30to50",t_muPt_cent30to50);
            f->GetObject("h_muPt_cent50to90",t_muPt_cent50to90);
        f->GetObject("h_muEta",t_muEta);
            f->GetObject("h_muEta_cent0to10",t_muEta_cent0to10);
            f->GetObject("h_muEta_cent10to30",t_muEta_cent10to30);
            f->GetObject("h_muEta_cent30to50",t_muEta_cent30to50);
            f->GetObject("h_muEta_cent50to90",t_muEta_cent50to90);
        f->GetObject("h_muPhi",t_muPhi);
            f->GetObject("h_muPhi_cent0to10",t_muPhi_cent0to10);
            f->GetObject("h_muPhi_cent10to30",t_muPhi_cent10to30);
            f->GetObject("h_muPhi_cent30to50",t_muPhi_cent30to50);
            f->GetObject("h_muPhi_cent50to90",t_muPhi_cent50to90);
        f->GetObject("h_muRelPt",t_muRelPt);
            f->GetObject("h_muRelPt_cent0to10",t_muRelPt_cent0to10);
            f->GetObject("h_muRelPt_cent10to30",t_muRelPt_cent10to30);
            f->GetObject("h_muRelPt_cent30to50",t_muRelPt_cent30to50);
            f->GetObject("h_muRelPt_cent50to90",t_muRelPt_cent50to90);
        f->GetObject("h_vz",t_vz);
        f->GetObject("h_hiBin",t_hiBin);
        f->GetObject("h_deltaR",t_deltaR);
        */

        if(file==0){
            h_jetpt = (TH1D*) t_jetpt->Clone("h_jetpt");
            /*
                h_jetpt_cent0to10 = (TH1D*) t_jetpt_cent0to10->Clone("h_jetpt_cent0to10");
                h_jetpt_cent10to30 = (TH1D*) t_jetpt_cent10to30->Clone("h_jetpt_cent10to30");
                h_jetpt_cent30to50 = (TH1D*) t_jetpt_cent30to50->Clone("h_jetpt_cent30to50");
                h_jetpt_cent50to90 = (TH1D*) t_jetpt_cent50to90->Clone("h_jetpt_cent50to90");
            h_mujetpt = (TH1D*) t_mujetpt->Clone("h_mujetpt");
                h_mujetpt_cent0to10 = (TH1D*) t_mujetpt_cent0to10->Clone("h_mujetpt_cent0to10");
                h_mujetpt_cent10to30 = (TH1D*) t_mujetpt_cent10to30->Clone("h_mujetpt_cent10to30");
                h_mujetpt_cent30to50 = (TH1D*) t_mujetpt_cent30to50->Clone("h_mujetpt_cent30to50");
                h_mujetpt_cent50to90 = (TH1D*) t_mujetpt_cent50to90->Clone("h_mujetpt_cent50to90");
            h_jeteta = (TH1D*) t_jeteta->Clone("h_jeteta");
                h_jeteta_cent0to10 = (TH1D*) t_jeteta_cent0to10->Clone("h_jeteta_cent0to10");
                h_jeteta_cent10to30 = (TH1D*) t_jeteta_cent10to30->Clone("h_jeteta_cent10to30");
                h_jeteta_cent30to50 = (TH1D*) t_jeteta_cent30to50->Clone("h_jeteta_cent30to50");
                h_jeteta_cent50to90 = (TH1D*) t_jeteta_cent50to90->Clone("h_jeteta_cent50to90");
            h_mujeteta = (TH1D*) t_mujeteta->Clone("h_mujeteta");
                h_mujeteta_cent0to10 = (TH1D*) t_mujeteta_cent0to10->Clone("h_mujeteta_cent0to10");
                h_mujeteta_cent10to30 = (TH1D*) t_mujeteta_cent10to30->Clone("h_mujeteta_cent10to30");
                h_mujeteta_cent30to50 = (TH1D*) t_mujeteta_cent30to50->Clone("h_mujeteta_cent30to50");
                h_mujeteta_cent50to90 = (TH1D*) t_mujeteta_cent50to90->Clone("h_mujeteta_cent50to90");
            h_jetphi = (TH1D*) t_jetphi->Clone("h_jetphi");
                h_jetphi_cent0to10 = (TH1D*) t_jetphi_cent0to10->Clone("h_jetphi_cent0to10");
                h_jetphi_cent10to30 = (TH1D*) t_jetphi_cent10to30->Clone("h_jetphi_cent10to30");
                h_jetphi_cent30to50 = (TH1D*) t_jetphi_cent30to50->Clone("h_jetphi_cent30to50");
                h_jetphi_cent50to90 = (TH1D*) t_jetphi_cent50to90->Clone("h_jetphi_cent50to90");
            h_mujetphi = (TH1D*) t_mujetphi->Clone("h_mujetphi");
                h_mujetphi_cent0to10 = (TH1D*) t_mujetphi_cent0to10->Clone("h_mujetphi_cent0to10");
                h_mujetphi_cent10to30 = (TH1D*) t_mujetphi_cent10to30->Clone("h_mujetphi_cent10to30");
                h_mujetphi_cent30to50 = (TH1D*) t_mujetphi_cent30to50->Clone("h_mujetphi_cent30to50");
                h_mujetphi_cent50to90 = (TH1D*) t_mujetphi_cent50to90->Clone("h_mujetphi_cent50to90");
            h_muPt = (TH1D*) t_muPt->Clone("h_muPt");
                h_muPt_cent0to10 = (TH1D*) t_muPt_cent0to10->Clone("h_muPt_cent0to10");
                h_muPt_cent10to30 = (TH1D*) t_muPt_cent10to30->Clone("h_muPt_cent10to30");
                h_muPt_cent30to50 = (TH1D*) t_muPt_cent30to50->Clone("h_muPt_cent30to50");
                h_muPt_cent50to90 = (TH1D*) t_muPt_cent50to90->Clone("h_muPt_cent50to90");
            h_muEta = (TH1D*) t_muEta->Clone("h_muEta");
                h_muEta_cent0to10 = (TH1D*) t_muEta_cent0to10->Clone("h_muEta_cent0to10");
                h_muEta_cent10to30 = (TH1D*) t_muEta_cent10to30->Clone("h_muEta_cent10to30");
                h_muEta_cent30to50 = (TH1D*) t_muEta_cent30to50->Clone("h_muEta_cent30to50");
                h_muEta_cent50to90 = (TH1D*) t_muEta_cent50to90->Clone("h_muEta_cent50to90");
            h_muPhi = (TH1D*) t_muPhi->Clone("h_muPhi");
                h_muPhi_cent0to10 = (TH1D*) t_muPhi_cent0to10->Clone("h_muPhi_cent0to10");
                h_muPhi_cent10to30 = (TH1D*) t_muPhi_cent10to30->Clone("h_muPhi_cent10to30");
                h_muPhi_cent30to50 = (TH1D*) t_muPhi_cent30to50->Clone("h_muPhi_cent30to50");
                h_muPhi_cent50to90 = (TH1D*) t_muPhi_cent50to90->Clone("h_muPhi_cent50to90");
            h_vz = (TH1D*) t_vz->Clone("h_vz");
            h_hiBin = (TH1D*) t_hiBin->Clone("h_hiBin");
            h_deltaR = (TH1D*) t_deltaR->Clone("h_deltaR");
            */
        }
        else{
            h_jetpt->Add(t_jetpt);
            /*
                h_jetpt_cent0to10->Add(t_jetpt_cent0to10);
                h_jetpt_cent10to30->Add(t_jetpt_cent10to30);
                h_jetpt_cent30to50->Add(t_jetpt_cent30to50);
                h_jetpt_cent50to90->Add(t_jetpt_cent50to90);
            h_mujetpt->Add(t_mujetpt);
                h_mujetpt_cent0to10->Add(t_mujetpt_cent0to10);
                h_mujetpt_cent10to30->Add(t_mujetpt_cent10to30);
                h_mujetpt_cent30to50->Add(t_mujetpt_cent30to50);
                h_mujetpt_cent50to90->Add(t_mujetpt_cent50to90);
            h_jeteta->Add(t_jeteta);
                h_jeteta_cent0to10->Add(t_jeteta_cent0to10);
                h_jeteta_cent10to30->Add(t_jeteta_cent10to30);
                h_jeteta_cent30to50->Add(t_jeteta_cent30to50);
                h_jeteta_cent50to90->Add(t_jeteta_cent50to90);
            h_mujeteta->Add(t_mujeteta);
                h_mujeteta_cent0to10->Add(t_mujeteta_cent0to10);
                h_mujeteta_cent10to30->Add(t_mujeteta_cent10to30);
                h_mujeteta_cent30to50->Add(t_mujeteta_cent30to50);
                h_mujeteta_cent50to90->Add(t_mujeteta_cent50to90);
            h_jetphi->Add(t_jetphi);
                h_jetphi_cent0to10->Add(t_jetphi_cent0to10);
                h_jetphi_cent10to30->Add(t_jetphi_cent10to30);
                h_jetphi_cent30to50->Add(t_jetphi_cent30to50);
                h_jetphi_cent50to90->Add(t_jetphi_cent50to90);
            h_mujetphi->Add(t_mujetphi);
                h_mujetphi_cent0to10->Add(t_mujetphi_cent0to10);
                h_mujetphi_cent10to30->Add(t_mujetphi_cent10to30);
                h_mujetphi_cent30to50->Add(t_mujetphi_cent30to50);
                h_mujetphi_cent50to90->Add(t_mujetphi_cent50to90);
            h_muPt->Add(t_muPt);
                h_muPt_cent0to10->Add(t_muPt_cent0to10);
                h_muPt_cent10to30->Add(t_muPt_cent10to30);
                h_muPt_cent30to50->Add(t_muPt_cent30to50);
                h_muPt_cent50to90->Add(t_muPt_cent50to90);
            h_muEta->Add(t_muEta);
                h_muEta_cent0to10->Add(t_muEta_cent0to10);
                h_muEta_cent10to30->Add(t_muEta_cent10to30);
                h_muEta_cent30to50->Add(t_muEta_cent30to50);
                h_muEta_cent50to90->Add(t_muEta_cent50to90);
            h_muPhi->Add(t_muPhi);
                h_muPhi_cent0to10->Add(t_muPhi_cent0to10);
                h_muPhi_cent10to30->Add(t_muPhi_cent10to30);
                h_muPhi_cent30to50->Add(t_muPhi_cent30to50);
                h_muPhi_cent50to90->Add(t_muPhi_cent50to90);
            h_vz->Add(t_vz);
            h_hiBin->Add(t_hiBin);
            h_deltaR->Add(t_deltaR);
            */
        }

    t_jetpt->Reset("ICESM");
        
    delete f;

    } // end file loop

    auto wf = TFile::Open("/home/clayton/Analysis/code/stitch/PbPb_data_stitch/PbPb_data_stitch.root","recreate");
    //h_jetpt->Write();
    wf->Close();


} // end program