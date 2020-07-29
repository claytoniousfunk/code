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
#include <array>

void RejectAcceptRatio(float etacut1=-1.5, float etacut2=1.5){

	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_recojets_pthat_30_jetptcut_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_recojets_pthat_50_jetptcut_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_recojets_pthat_80_jetptcut_50_muptcut_5.root");

	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_30_jetptcut_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_50_jetptcut_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_80_jetptcut_50_muptcut_5.root");

	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_80_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_80_muptcut_5.root");	
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_30_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_30_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_80_muptcut_5.root");
	
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");
	//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_refjets_pthat_30_muptcut_5.root");
	TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_recojets_pthat_30_muptcut_5.root");

	TH2D *h2_reject, *h2_accept;
	TH1D *h_reject, *h_accept;

	
	f->GetObject("h2",h2_accept);
	f->GetObject("h2_lowPthatHighPtjet",h2_reject);

	TH1D* hx_test=h2_accept->ProjectionX();
	TAxis *xaxis = hx_test->GetXaxis();
	int firstxbin = xaxis->FindBin(etacut1);
	cout << "first x bin = " << firstxbin << endl;
	int lastxbin = xaxis->FindBin(etacut2);
	cout << "last x bin = " << lastxbin << endl;
	
	h_reject = h2_reject->ProjectionY("h_reject",firstxbin,lastxbin);
	h_accept = h2_accept->ProjectionY("h_accept",firstxbin,lastxbin);
	

	
	

	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0);
	pad1->Draw();
	pad1->SetLogy();
	pad1->SetLeftMargin(0.15);
	pad1->cd();
	
	


	
	
	const int NPtAxisBins = 19;
	double pt_axis[NPtAxisBins+1] = {50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
	
	
	TH1D *h_accept_rebin = (TH1D*) h_accept->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *h_reject_rebin = (TH1D*) h_reject->Rebin(NPtAxisBins,"",pt_axis);
	
	TH1D *r_rebin = (TH1D*) h_reject_rebin->Clone("r");
	r_rebin->Divide(h_accept_rebin);

	


	r_rebin->SetMarkerStyle(8);
    r_rebin->SetMarkerColor(kAzure+9);
    r_rebin->SetMarkerSize(1.2);
    r_rebin->SetFillColorAlpha(kAzure+10,0.5);
    r_rebin->SetStats(0);
    r_rebin->Draw("e2");

    r_rebin->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
    r_rebin->GetYaxis()->SetTitle("rejected jets / accepted jets");
    r_rebin->SetTitle("");



      

    




    c1->SaveAs("/home/clayton/Analysis/code/RejectAcceptRatio/figures/RejectAcceptRatio.png");


}
