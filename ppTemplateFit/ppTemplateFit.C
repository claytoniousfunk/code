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

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_80_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_80_muptcut_5.root");
// mu pt cut dependence study
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_10.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_15.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_20.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV2_PbPb_refjets_pthat_50_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV2_pp_refjets_pthat_50_muptcut_5.root");


// rel pt muon study, 7/6/20
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_5.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_10.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_10.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_15.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_15.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_20.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_20.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");

TFile *f = TFile::Open("/home/clayton/Analysis/data/ppskim/pythia_muptcut_10.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/data/ppskim/clayton_pp_mc_skim.root");
TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_muptcut_5.root");

TH1D *h_muRelPt, *h_muRelPt_g, *h_muRelPt_c, *h_muRelPt_b, *h_muRelPt_uds, *h_muRelPt_q, *h_muRelPt_udsgc, *h_muRelPt_udsg, *h_muRelPt_ghost;

double func_temp_0(double *x, double *par){
	double xx = x[0];
	int bin = h_muRelPt_g->FindBin(xx);
	double q = (par[0])*h_muRelPt_q->GetBinContent(bin);
	double g = (1-par[0])*h_muRelPt_g->GetBinContent(bin);
	return q+g;
}

double func_temp_1(double *x, double *par){
	double xx = x[0];
	int bin = h_muRelPt_b->FindBin(xx);
	double b = (par[0])*h_muRelPt_b->GetBinContent(bin);
	double e = (1-par[0])*h_muRelPt_udsgc->GetBinContent(bin);
	return b+e;
}
double func_temp_2(double *x, double *par){
	double xx = x[0];
	int bin = h_muRelPt_b->FindBin(xx);
	double b = (par[0])*h_muRelPt_b->GetBinContent(bin);
	double c = (par[1])*h_muRelPt_c->GetBinContent(bin);
	double e = (1-par[0]-par[1])*h_muRelPt_udsg->GetBinContent(bin);
	return b+c+e;
}

void ppTemplateFit(int plotID = 2, bool isStack = 1){
	
	

	f_data->GetObject("h_muRelPt",h_muRelPt);
	f->GetObject("h_muRelPt_g",h_muRelPt_g);
	f->GetObject("h_muRelPt_uds",h_muRelPt_uds);
	f->GetObject("h_muRelPt_c",h_muRelPt_c);
	f->GetObject("h_muRelPt_b",h_muRelPt_b);
	f->GetObject("h_muRelPt_ghost",h_muRelPt_ghost);
	

	
///////////////////////////////////////////////////////////////  all quarks & gluons ///////////////////////////////////////////////////////////////
	if(plotID==0){

		// build the all-quarks template
		h_muRelPt_q = (TH1D*)h_muRelPt_c->Clone("h_muRelPt_q");
		h_muRelPt_q->Add(h_muRelPt_b);
		h_muRelPt_q->Add(h_muRelPt_uds);
		
		// normalize the templates
		h_muRelPt_g->Scale(1./h_muRelPt_g->Integral());
		h_muRelPt_q->Scale(1./h_muRelPt_q->Integral());
		// normalize the data
		h_muRelPt->Scale(1./h_muRelPt->Integral());

		TH1D *fitRatio = (TH1D*) h_muRelPt->Clone("fitRatio");

		double low_x = 0.0;
		double high_x = 5.0;
		int numPar = 1;

		TCanvas *c1 = new TCanvas("c1","c1",800,800);
		c1->cd();
		TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
		TPad *pad2 = new TPad("pad2","pad2",0.,0.0,1.0,0.3);
		pad1->Draw();
		pad1->cd();
		
		TF1 *func = new TF1("func",func_temp_0,low_x,high_x,numPar);
		h_muRelPt->Fit(func,"M R","N",low_x,high_x);
		gStyle->SetOptFit(1);
		h_muRelPt->SetTitle("");
		
		h_muRelPt->SetStats(1);
		h_muRelPt->SetMarkerStyle(22);
		h_muRelPt->SetMarkerColor(kBlue);
		h_muRelPt->SetFillColorAlpha(kBlue,0.5);

		h_muRelPt_q->SetMarkerStyle(8);
		h_muRelPt_q->SetMarkerColor(kRed);
		h_muRelPt_q->SetFillColorAlpha(kRed,0.5);

		h_muRelPt_g->SetMarkerStyle(33);
		h_muRelPt_g->SetMarkerColor(kGreen);
		h_muRelPt_g->SetFillColorAlpha(kGreen,0.5);

		func->SetLineStyle(10);

		h_muRelPt->Draw("e2");
		h_muRelPt_q->Draw("e2 same");
		h_muRelPt_g->Draw("e2 same");
		func->Draw("same");

		auto legend = new TLegend(0.6,0.3,0.8,0.6);
		legend->AddEntry(h_muRelPt,"inclusive","fp");
		legend->AddEntry(h_muRelPt_q,"quark template","fp");
		legend->AddEntry(h_muRelPt_g,"gluon template","fp");
		legend->AddEntry(func,"fit");
		legend->Draw();

		c1->cd();
		pad2->SetBottomMargin(0.3);
		pad2->Draw();
		pad2->cd();
		
		fitRatio->Divide(func);
		fitRatio->SetStats(0);
		fitRatio->SetTitle("");
		fitRatio->SetMarkerStyle(8);
		fitRatio->SetMarkerColor(kAzure);
		fitRatio->SetFillColorAlpha(kAzure,0.5);
		fitRatio->Draw("e2");
		fitRatio->GetYaxis()->SetTitle("data / fit");
		fitRatio->GetYaxis()->SetTitleSize(0.1);
		fitRatio->GetYaxis()->SetTitleOffset(0.5);
		fitRatio->GetXaxis()->SetTitle("muon rel p_{T} [GeV/c]");
		fitRatio->GetXaxis()->SetTitleOffset(1.1);
		fitRatio->GetXaxis()->SetTitleSize(.1);
		fitRatio->SetMinimum(0.8);
		fitRatio->SetMaximum(1.2);

		c1->SaveAs("/home/clayton/Analysis/code/ppTemplateFit/figures/templateFit_gVsElse.pdf");
	}

	///////////////////////////////////////////////////////////////  b quarks vs everything else ///////////////////////////////////////////////////////////////
	if(plotID==1){
		//build the "everything else" template
		h_muRelPt_udsgc = (TH1D*)h_muRelPt_g->Clone("h_muRelPt_udsgc");
		h_muRelPt_udsgc->Add(h_muRelPt_uds);
		h_muRelPt_udsgc->Add(h_muRelPt_c);
		// normalize the templates
		h_muRelPt_udsgc->Scale(1./h_muRelPt_udsgc->Integral());
		h_muRelPt_b->Scale(1./h_muRelPt_b->Integral());
		// normalize the data
		h_muRelPt->Scale(1./h_muRelPt->Integral());

		TH1D *fitRatio = (TH1D*) h_muRelPt->Clone("fitRatio");

		double low_x = 0.0;
		double high_x = 5.0;
		int numPar = 1;

		TCanvas *c1 = new TCanvas("c1","c1",800,800);
		c1->cd();
		TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
		TPad *pad2 = new TPad("pad2","pad2",0.,0.0,1.0,0.3);
		pad1->Draw();
		pad1->cd();

		TF1 *func = new TF1("func",func_temp_1,low_x,high_x,numPar);
		h_muRelPt->Fit(func,"M R","N",low_x,high_x);
		gStyle->SetOptFit(1);
		h_muRelPt->SetTitle("");
		
		h_muRelPt->SetStats(1);
		h_muRelPt->SetMarkerStyle(22);
		h_muRelPt->SetMarkerColor(kBlue);
		h_muRelPt->SetFillColorAlpha(kBlue,0.5);

		h_muRelPt_b->SetMarkerStyle(33);
		h_muRelPt_b->SetMarkerColor(kOrange);
		h_muRelPt_b->SetFillColorAlpha(kOrange,0.5);

		h_muRelPt_udsgc->SetMarkerStyle(43);
		h_muRelPt_udsgc->SetMarkerColor(kPink);
		h_muRelPt_udsgc->SetFillColorAlpha(kPink,0.5);

		func->SetLineStyle(10);

		h_muRelPt->Draw("e2");
		h_muRelPt_b->Draw("e2 same");
		h_muRelPt_udsgc->Draw("e2 same");
		func->Draw("same");

		auto legend = new TLegend(0.6,0.3,0.8,0.6);
		legend->AddEntry(h_muRelPt,"inclusive","fp");
		legend->AddEntry(h_muRelPt_b,"b template","fp");
		legend->AddEntry(h_muRelPt_udsgc,"udscg template","fp");
		legend->AddEntry(func,"fit");
		legend->Draw();

		c1->cd();
		pad2->SetBottomMargin(0.3);
		pad2->Draw();
		pad2->cd();
		
		fitRatio->Divide(func);
		fitRatio->SetStats(0);
		fitRatio->SetTitle("");
		fitRatio->SetMarkerStyle(8);
		fitRatio->SetMarkerColor(kAzure);
		fitRatio->SetFillColorAlpha(kAzure,0.5);
		fitRatio->Draw("e2");
		fitRatio->GetYaxis()->SetTitle("data / fit");
		fitRatio->GetYaxis()->SetTitleSize(0.1);
		fitRatio->GetYaxis()->SetTitleOffset(0.5);
		fitRatio->GetXaxis()->SetTitle("muon rel p_{T} [GeV/c]");
		fitRatio->GetXaxis()->SetTitleOffset(1.1);
		fitRatio->GetXaxis()->SetTitleSize(.1);
		fitRatio->SetMinimum(0.8);
		fitRatio->SetMaximum(1.2);

		c1->SaveAs("/home/clayton/Analysis/code/ppTemplateFit/figures/templateFit_bVsElse.pdf");

	}
	///////////////////////////////////////////////////////////////  b quarks vs c quarks vs everything else ///////////////////////////////////////////////////////////////
	if(plotID==2){
		// build h_udsg
		h_muRelPt_udsg = (TH1D*) h_muRelPt_uds->Clone("h_muRelPt_udsg");
		h_muRelPt_udsg->Add(h_muRelPt_g);
		
		// normalize the templates
		h_muRelPt_udsg->Scale(1./h_muRelPt_udsg->Integral());
		h_muRelPt_b->Scale(1./h_muRelPt_b->Integral());
		h_muRelPt_c->Scale(1./h_muRelPt_c->Integral());
		// normalize the data
		h_muRelPt->Scale(1./h_muRelPt->Integral());
		// define the stack
		THStack *h_stack = new THStack("h_stack","");



		TH1D *fitRatio = (TH1D*) h_muRelPt->Clone("fitRatio");

		double low_x = 0.0;
		double high_x = 5.0;
		int numPar = 2;

		TCanvas *c1 = new TCanvas("c1","c1",800,800);
		c1->cd();
		TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
		TPad *pad2 = new TPad("pad2","pad2",0.,0.0,1.0,0.3);
		pad1->Draw();
		pad1->cd();

		TF1 *func = new TF1("func",func_temp_2,low_x,high_x,numPar);
		//func->SetParameter(0,0.3);
		//func->SetParameter(1,0.3);
		func->SetParName(0,"b");
		func->SetParName(1,"c");
		func->SetParLimits(0,0.0,1.0);
		func->SetParLimits(1,0.0,1.0);
		h_muRelPt->Fit(func,"M R","N",low_x,high_x);

		gStyle->SetOptFit(1);
		h_muRelPt->SetTitle("");

		h_muRelPt->SetMarkerStyle(8);
		h_muRelPt->SetMarkerColor(kBlack);
		h_muRelPt->SetFillColorAlpha(kBlack,0.5);

		h_muRelPt_b->SetMarkerStyle(8);
		h_muRelPt_b->SetMarkerColor(kBlue);
		h_muRelPt_b->SetFillColorAlpha(kBlue,0.5);

		h_muRelPt_c->SetMarkerStyle(8);
		h_muRelPt_c->SetMarkerColor(kRed);
		h_muRelPt_c->SetFillColorAlpha(kRed,0.5);

		h_muRelPt_udsg->SetMarkerStyle(8);
		h_muRelPt_udsg->SetMarkerColor(kGreen+2);
		h_muRelPt_udsg->SetFillColorAlpha(kGreen+2,0.5);

		double p0 = func->GetParameter(0);
		double p1 = func->GetParameter(1);


		h_muRelPt->SetStats(1);

		TH1D *h_muRelPt_udsg_scaled = (TH1D*) h_muRelPt_udsg->Clone("h_muRelPt_udsg_scaled");
		TH1D *h_muRelPt_c_scaled = (TH1D*) h_muRelPt_c->Clone("h_muRelPt_c_scaled");
		TH1D *h_muRelPt_b_scaled = (TH1D*) h_muRelPt_b->Clone("h_muRelPt_b_scaled");

		h_muRelPt_udsg_scaled->Scale(1.0-p0-p1);
		h_muRelPt_c_scaled->Scale(p1);
		h_muRelPt_b_scaled->Scale(p0);

		h_stack->Add(h_muRelPt_udsg_scaled);
		h_stack->Add(h_muRelPt_c_scaled);
		h_stack->Add(h_muRelPt_b_scaled);

		func->SetLineStyle(10);

		if(!isStack){
			h_muRelPt->Draw("e2");
			h_muRelPt_b->Draw("e2 same");
			h_muRelPt_c->Draw("e2 same");
			h_muRelPt_udsg->Draw("e2 same");
			func->Draw("same");

			auto legend = new TLegend(0.6,0.25,0.8,0.6);
			legend->AddEntry(h_muRelPt,"inclusive","fp");
			legend->AddEntry(h_muRelPt_b,"b template","fp");
			legend->AddEntry(h_muRelPt_c,"c template","fp");
			legend->AddEntry(h_muRelPt_udsg,"udsg template","fp");
			legend->AddEntry(func,"fit");
			legend->Draw();
		}

		if(isStack){
			
			h_stack->Draw("hist");
			h_muRelPt->Draw("ep same");
			func->Draw("same");
			auto legend = new TLegend(0.6,0.25,0.8,0.6);
			legend->AddEntry(h_muRelPt,"inclusive","p");
			legend->AddEntry(h_muRelPt_b_scaled,"b template","f");
			legend->AddEntry(h_muRelPt_c_scaled,"c template","f");
			legend->AddEntry(h_muRelPt_udsg_scaled,"udsg template","f");
			legend->AddEntry(func,"fit");
			legend->Draw();

		}

		

		c1->cd();
		pad2->SetBottomMargin(0.3);
		pad2->Draw();
		pad2->cd();
		
		fitRatio->Divide(func);
		fitRatio->SetStats(0);
		fitRatio->SetTitle("");
		fitRatio->SetMarkerStyle(8);
		fitRatio->SetMarkerColor(kAzure);
		fitRatio->SetFillColorAlpha(kAzure,0.5);
		fitRatio->Draw("e2");
		fitRatio->GetYaxis()->SetTitle("data / fit");
		fitRatio->GetYaxis()->SetTitleSize(0.1);
		fitRatio->GetYaxis()->SetTitleOffset(0.5);
		fitRatio->GetXaxis()->SetTitle("muon rel p_{T} [GeV/c]");
		fitRatio->GetXaxis()->SetTitleOffset(1.1);
		fitRatio->GetXaxis()->SetTitleSize(.1);
		fitRatio->SetMinimum(0.8);
		fitRatio->SetMaximum(1.2);

		if(isStack){c1->SaveAs("/home/clayton/Analysis/code/ppTemplateFit/figures/ppTemplateFit_bVsCVsElse_stack.pdf");}
		if(!isStack){c1->SaveAs("/home/clayton/Analysis/code/ppTemplateFit/figures/ppTemplateFit_bVsCVsElse.pdf");}

		
		// check fit against count
		
		double n_tot = h_muRelPt_udsg->GetEntries() + h_muRelPt_b->GetEntries() + h_muRelPt_c->GetEntries();
		double n_udsg = h_muRelPt_udsg->GetEntries();
		double n_b = h_muRelPt_b->GetEntries();
		double n_c = h_muRelPt_c->GetEntries();
		double n_ghost = h_muRelPt_ghost->GetEntries();

		cout << "Total entries = " << n_tot << endl;
		cout << "udsg entries = " << n_udsg << endl;
		cout << "b entries = " << n_b << endl;
		cout << "c entries = " << n_c << endl;
		cout << "ghost entries = " << n_ghost << endl;

		double nr_udsg = n_udsg/n_tot;
		double nr_b = n_b/n_tot;
		double nr_c = n_c/n_tot;
		
		cout << "b fraction = " << nr_b << endl;
		cout << "c fraction = " << nr_c << endl;
		cout << "udsg fraction = " << nr_udsg << endl;


		

	}
	


}

