
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

void muRelPtCentralityRegions(){

	TH1D *h_muRelPt_centReg1, *h_muRelPt_centReg2, *h_muRelPt_centReg3, *h_muRelPt_centReg4;
	TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV7_recojets_pthat_30_muptcut_5.root");

	f->GetObject("h_muRelPt_centReg1",h_muRelPt_centReg1);
	f->GetObject("h_muRelPt_centReg2",h_muRelPt_centReg2);
	f->GetObject("h_muRelPt_centReg3",h_muRelPt_centReg3);
	f->GetObject("h_muRelPt_centReg4",h_muRelPt_centReg4);






}