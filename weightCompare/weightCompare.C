



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

//TFile *f = TFile::Open("makeDataV6_recojets_pthat_30_muptcut_5.root");
TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");

TH1D *h_vz, *h_vz_noWeight, *h_hiBin, *h_hiBin_noWeight;
TH2D *h2, *h2_noVzOrHiBinWeight, *h2_lowPthatHighPtjet;



void weightCompare(){

	f->GetObject("h_vz",h_vz);
	f->GetObject("h_vz_noWeight",h_vz_noWeight);
	f->GetObject("h_hiBin",h_hiBin);
	f->GetObject("h_hiBin_noWeight",h_hiBin_noWeight);
	f->GetObject("h2_lowPthatHighPtjet",h2_lowPthatHighPtjet);
	f->GetObject("h2",h2);
	f->GetObject("h2_noVzOrHiBinWeight",h2_noVzOrHiBinWeight);

	double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
	double vz_axis[31] = {-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,
		2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};

	TH1D *h_vz_rebin = (TH1D*)h_vz->Rebin(30,"",vz_axis);
	TH1D *h_vz_noWeight_rebin = (TH1D*)h_vz_noWeight->Rebin(30,"",vz_axis);
	TH1D *vz_ratio = (TH1D*)h_vz_rebin->Clone("vz_ratio");
	vz_ratio->Divide(h_vz_noWeight_rebin);

	TH1D *hiBin_ratio = (TH1D*)h_hiBin->Clone("hiBin_ratio");
	hiBin_ratio->Divide(h_hiBin_noWeight);

	TH1D *h = h2->ProjectionY();
	TH1D *h_noWeight = h2_noVzOrHiBinWeight->ProjectionY();
	TH1D *h_rebin = (TH1D*)h->Rebin(15,"",pt_axis);
	TH1D *h_noWeight_rebin = (TH1D*)h_noWeight->Rebin(15,"",pt_axis);
	TH1D *pt_ratio = (TH1D*) h_rebin->Clone("pt_ratio");
	pt_ratio->Divide(h_noWeight_rebin);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.3,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	pad1->Draw();
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad1->cd();

	double markerSize = 1.2;

	h_vz_rebin->SetMarkerStyle(33);
	h_vz_rebin->SetMarkerSize(markerSize);
	h_vz_rebin->SetMarkerColor(kAzure+3);
	h_vz_rebin->SetFillColorAlpha(kAzure+3,0.5);
	h_vz_rebin->SetStats(0);
	
	h_vz_rebin->SetTitle("");
	h_vz_rebin->GetYaxis()->SetTitle("Entries");

	h_vz_rebin->Draw("e2");

	h_vz_noWeight_rebin->SetMarkerStyle(20);
	h_vz_noWeight_rebin->SetMarkerSize(markerSize);
	h_vz_noWeight_rebin->SetMarkerColor(kOrange-7);
	h_vz_noWeight_rebin->SetFillColorAlpha(kOrange-7,0.5);
	

	h_vz_noWeight_rebin->Draw("e2 same");

	auto legend = new TLegend(0.8,0.7,0.95,0.85);
	legend->AddEntry(h_vz_rebin,"with weights","fp");
	legend->AddEntry(h_vz_noWeight_rebin,"no weights","fp");
	legend->Draw();

	pad2->cd();
	vz_ratio->SetMarkerStyle(8);
	vz_ratio->SetMarkerSize(1.0);
	vz_ratio->SetMarkerColor(kBlack);
	vz_ratio->Draw("ep");
	vz_ratio->GetXaxis()->SetTitle("Vertex position [cm]");
	vz_ratio->GetYaxis()->SetTitle("weights / no weights");
	vz_ratio->GetYaxis()->SetTitleSize(0.1);
	vz_ratio->GetYaxis()->SetTitleOffset(0.5);
	vz_ratio->GetXaxis()->SetTitleOffset(1.1);
	vz_ratio->GetXaxis()->SetTitleSize(.1);
	vz_ratio->SetTitle("");
	vz_ratio->SetStats(0);

	c1->SaveAs("/home/clayton/Analysis/code/weightCompare/figures/vzWeightCompare.pdf");
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TCanvas *c2 = new TCanvas("c2","c2",900,900);
	TPad *pad3 = new TPad("pad3","pad3",0.,0.3,1.0,1.0);
	TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,0.3);
	pad3->Draw();
	pad4->SetBottomMargin(0.3);
	pad4->Draw();
	pad3->cd();

	h_hiBin->SetMarkerStyle(33);
	h_hiBin->SetMarkerSize(markerSize);
	h_hiBin->SetMarkerColor(kAzure+3);
	h_hiBin->SetFillColorAlpha(kAzure+3,0.5);
	h_hiBin->SetStats(0);
	
	h_hiBin->SetTitle("");
	h_hiBin->GetYaxis()->SetTitle("Entries");

	h_hiBin->Draw("e2");

	h_hiBin_noWeight->SetMarkerStyle(20);
	h_hiBin_noWeight->SetMarkerSize(markerSize);
	h_hiBin_noWeight->SetMarkerColor(kOrange-7);
	h_hiBin_noWeight->SetFillColorAlpha(kOrange-7,0.5);
	

	h_hiBin_noWeight->Draw("e2 same");

	auto legend2 = new TLegend(0.8,0.7,0.95,0.85);
	legend2->AddEntry(h_hiBin,"with weights","fp");
	legend2->AddEntry(h_hiBin_noWeight,"no weights","fp");
	legend2->Draw();

	pad4->cd();
	hiBin_ratio->SetMarkerStyle(8);
	hiBin_ratio->SetMarkerSize(1.0);
	hiBin_ratio->SetMarkerColor(kBlack);
	hiBin_ratio->Draw("ep");
	hiBin_ratio->GetXaxis()->SetTitle("hiBin (centrality*2)");
	hiBin_ratio->GetYaxis()->SetTitle("weights / no weights");
	hiBin_ratio->GetYaxis()->SetTitleSize(0.1);
	hiBin_ratio->GetYaxis()->SetTitleOffset(0.5);
	hiBin_ratio->GetXaxis()->SetTitleOffset(1.1);
	hiBin_ratio->GetXaxis()->SetTitleSize(.1);
	hiBin_ratio->SetTitle("");
	hiBin_ratio->SetStats(0);

	c2->SaveAs("/home/clayton/Analysis/code/weightCompare/figures/hiBinWeightCompare.pdf");
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TCanvas *c3 = new TCanvas("c3","c3",900,900);
	TPad *pad5 = new TPad("pad5","pad5",0.,0.3,1.0,1.0);
	TPad *pad6 = new TPad("pad6","pad6",0.0,0.0,1.0,0.3);
	pad5->Draw();
	pad5->SetLogy();
	pad6->SetBottomMargin(0.3);
	pad6->Draw();
	pad5->cd();

	h_rebin->SetMarkerStyle(33);
	h_rebin->SetMarkerSize(markerSize);
	h_rebin->SetMarkerColor(kAzure+3);
	h_rebin->SetFillColorAlpha(kAzure+3,0.5);
	h_rebin->SetStats(0);
	
	h_rebin->SetTitle("");
	h_rebin->GetYaxis()->SetTitle("Entries");

	h_rebin->Draw("e2");

	h_noWeight_rebin->SetMarkerStyle(20);
	h_noWeight_rebin->SetMarkerSize(markerSize);
	h_noWeight_rebin->SetMarkerColor(kOrange-7);
	h_noWeight_rebin->SetFillColorAlpha(kOrange-7,0.5);
	

	h_noWeight_rebin->Draw("e2 same");

	auto legend3 = new TLegend(0.8,0.7,0.95,0.85);
	legend3->AddEntry(h_rebin,"with weights","fp");
	legend3->AddEntry(h_noWeight_rebin,"no weights","fp");
	legend3->Draw();

	pad6->cd();
	pt_ratio->SetMarkerStyle(8);
	pt_ratio->SetMarkerSize(1.0);
	pt_ratio->SetMarkerColor(kBlack);
	pt_ratio->Draw("ep");
	pt_ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	pt_ratio->GetYaxis()->SetTitle("weights / no weights");
	pt_ratio->GetYaxis()->SetTitleSize(0.1);
	pt_ratio->GetYaxis()->SetTitleOffset(0.5);
	pt_ratio->GetXaxis()->SetTitleOffset(1.1);
	pt_ratio->GetXaxis()->SetTitleSize(.1);
	pt_ratio->SetTitle("");
	pt_ratio->SetStats(0);

	c3->SaveAs("/home/clayton/Analysis/code/weightCompare/pthatWeightCompare.pdf");
	





}
