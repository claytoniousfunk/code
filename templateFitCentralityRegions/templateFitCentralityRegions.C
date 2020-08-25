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

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_recojets_pthat_30_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_recojets_pthat_30_muptcut_5.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_refjets_pthat_30_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V7/rootFiles/makeDataV7_refjets_pthat_30_muptcut_5.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");

//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_recojets_pthat_30_muptcut_5.root");
//TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_recojets_pthat_30_muptcut_5.root");

TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_refjets_pthat_30_muptcut_5.root");
TFile *f_data = TFile::Open("/home/clayton/Analysis/code/makeData/V9/rootFiles/makeDataV9_refjets_pthat_30_muptcut_5.root");



TH1D *h_muRelPt_centReg1, *h_muRelPt_centReg2, *h_muRelPt_centReg3, *h_muRelPt_centReg4;
TH1D *h_muRelPt_g_centReg1, *h_muRelPt_lq_centReg1, *h_muRelPt_sall_centReg1, *h_muRelPt_ball_centReg1, *h_muRelPt_call_centReg1, *h_muRelPt_ghost_centReg1;
TH1D *h_muRelPt_g_centReg2, *h_muRelPt_lq_centReg2, *h_muRelPt_sall_centReg2, *h_muRelPt_ball_centReg2, *h_muRelPt_call_centReg2, *h_muRelPt_ghost_centReg2;
TH1D *h_muRelPt_g_centReg3, *h_muRelPt_lq_centReg3, *h_muRelPt_sall_centReg3, *h_muRelPt_ball_centReg3, *h_muRelPt_call_centReg3, *h_muRelPt_ghost_centReg3;
TH1D *h_muRelPt_g_centReg4, *h_muRelPt_lq_centReg4, *h_muRelPt_sall_centReg4, *h_muRelPt_ball_centReg4, *h_muRelPt_call_centReg4, *h_muRelPt_ghost_centReg4;
TH1D *h_muRelPt, *h_muRelPt_g, *h_muRelPt_call, *h_muRelPt_ball, *h_muRelPt_sall, *h_muRelPt_lq, *h_muRelPt_q, *h_muRelPt_ee, *h_muRelPt_ghost;


double func_temp_2(double *x, double *par){
	double xx = x[0];
	int bin = h_muRelPt_ball->FindBin(xx);
	double b = (par[0])*h_muRelPt_ball->GetBinContent(bin);
	double c = (par[1])*h_muRelPt_call->GetBinContent(bin);
	double e = (1-par[0]-par[1])*h_muRelPt_ee->GetBinContent(bin);
	return b+c+e;
}

void templateFitCentralityRegions(bool isStack = 1, int centReg = 1){
	
	

	f_data->GetObject("h_muRelPt_centReg1",h_muRelPt_centReg1);
	f_data->GetObject("h_muRelPt_g_centReg1",h_muRelPt_g_centReg1);
	f_data->GetObject("h_muRelPt_lq_centReg1",h_muRelPt_lq_centReg1);
	f_data->GetObject("h_muRelPt_sall_centReg1",h_muRelPt_sall_centReg1);
	f_data->GetObject("h_muRelPt_ball_centReg1",h_muRelPt_ball_centReg1);
	f_data->GetObject("h_muRelPt_call_centReg1",h_muRelPt_call_centReg1);
	f_data->GetObject("h_muRelPt_ghost_centReg1",h_muRelPt_ghost_centReg1);


	f_data->GetObject("h_muRelPt_centReg2",h_muRelPt_centReg2);
	f_data->GetObject("h_muRelPt_g_centReg2",h_muRelPt_g_centReg2);
	f_data->GetObject("h_muRelPt_lq_centReg2",h_muRelPt_lq_centReg2);
	f_data->GetObject("h_muRelPt_sall_centReg2",h_muRelPt_sall_centReg2);
	f_data->GetObject("h_muRelPt_ball_centReg2",h_muRelPt_ball_centReg2);
	f_data->GetObject("h_muRelPt_call_centReg2",h_muRelPt_call_centReg2);
	f_data->GetObject("h_muRelPt_ghost_centReg2",h_muRelPt_ghost_centReg2);


	f_data->GetObject("h_muRelPt_centReg3",h_muRelPt_centReg3);
	f_data->GetObject("h_muRelPt_g_centReg3",h_muRelPt_g_centReg3);
	f_data->GetObject("h_muRelPt_lq_centReg3",h_muRelPt_lq_centReg3);
	f_data->GetObject("h_muRelPt_sall_centReg3",h_muRelPt_sall_centReg3);
	f_data->GetObject("h_muRelPt_ball_centReg3",h_muRelPt_ball_centReg3);
	f_data->GetObject("h_muRelPt_call_centReg3",h_muRelPt_call_centReg3);
	f_data->GetObject("h_muRelPt_ghost_centReg3",h_muRelPt_ghost_centReg3);


	f_data->GetObject("h_muRelPt_centReg4",h_muRelPt_centReg4);
	f_data->GetObject("h_muRelPt_g_centReg4",h_muRelPt_g_centReg4);
	f_data->GetObject("h_muRelPt_lq_centReg4",h_muRelPt_lq_centReg4);
	f_data->GetObject("h_muRelPt_sall_centReg4",h_muRelPt_sall_centReg4);
	f_data->GetObject("h_muRelPt_ball_centReg4",h_muRelPt_ball_centReg4);
	f_data->GetObject("h_muRelPt_call_centReg4",h_muRelPt_call_centReg4);
	f_data->GetObject("h_muRelPt_ghost_centReg4",h_muRelPt_ghost_centReg4);




	///////////////////////////////////////////////////////////////  b quarks vs c quarks vs everything else ///////////////////////////////////////////////////////////////
	
		// define the histo to be fit to
		if(centReg==1) {
			h_muRelPt = (TH1D*) h_muRelPt_centReg1->Clone("h_muRelPt");
			h_muRelPt_g = (TH1D*) h_muRelPt_g_centReg1->Clone("h_muRelPt_g");
			h_muRelPt_lq = (TH1D*) h_muRelPt_lq_centReg1->Clone("h_muRelPt_lq");
			h_muRelPt_sall = (TH1D*) h_muRelPt_sall_centReg1->Clone("h_muRelPt_sall");
			h_muRelPt_ball = (TH1D*) h_muRelPt_ball_centReg1->Clone("h_muRelPt_ball");
			h_muRelPt_call = (TH1D*) h_muRelPt_call_centReg1->Clone("h_muRelPt_call");
			h_muRelPt_ghost = (TH1D*) h_muRelPt_ghost_centReg1->Clone("h_muRelPt_ghost");
		}
		if(centReg==2){
			h_muRelPt = (TH1D*) h_muRelPt_centReg2->Clone("h_muRelPt");
			h_muRelPt_g = (TH1D*) h_muRelPt_g_centReg2->Clone("h_muRelPt_g");
			h_muRelPt_lq = (TH1D*) h_muRelPt_lq_centReg2->Clone("h_muRelPt_lq");
			h_muRelPt_sall = (TH1D*) h_muRelPt_sall_centReg2->Clone("h_muRelPt_sall");
			h_muRelPt_ball = (TH1D*) h_muRelPt_ball_centReg2->Clone("h_muRelPt_ball");
			h_muRelPt_call = (TH1D*) h_muRelPt_call_centReg2->Clone("h_muRelPt_call");
			h_muRelPt_ghost = (TH1D*) h_muRelPt_ghost_centReg2->Clone("h_muRelPt_ghost");
		}
		if(centReg==3){
			h_muRelPt = (TH1D*) h_muRelPt_centReg3->Clone("h_muRelPt");
			h_muRelPt_g = (TH1D*) h_muRelPt_g_centReg3->Clone("h_muRelPt_g");
			h_muRelPt_lq = (TH1D*) h_muRelPt_lq_centReg3->Clone("h_muRelPt_lq");
			h_muRelPt_sall = (TH1D*) h_muRelPt_sall_centReg3->Clone("h_muRelPt_sall");
			h_muRelPt_ball = (TH1D*) h_muRelPt_ball_centReg3->Clone("h_muRelPt_ball");
			h_muRelPt_call = (TH1D*) h_muRelPt_call_centReg3->Clone("h_muRelPt_call");
			h_muRelPt_ghost = (TH1D*) h_muRelPt_ghost_centReg3->Clone("h_muRelPt_ghost");
		}
		if(centReg==4){
			h_muRelPt = (TH1D*) h_muRelPt_centReg4->Clone("h_muRelPt");
			h_muRelPt_g = (TH1D*) h_muRelPt_g_centReg4->Clone("h_muRelPt_g");
			h_muRelPt_lq = (TH1D*) h_muRelPt_lq_centReg4->Clone("h_muRelPt_lq");
			h_muRelPt_sall = (TH1D*) h_muRelPt_sall_centReg4->Clone("h_muRelPt_sall");
			h_muRelPt_ball = (TH1D*) h_muRelPt_ball_centReg4->Clone("h_muRelPt_ball");
			h_muRelPt_call = (TH1D*) h_muRelPt_call_centReg4->Clone("h_muRelPt_call");
			h_muRelPt_ghost = (TH1D*) h_muRelPt_ghost_centReg4->Clone("h_muRelPt_ghost");
		}



		//build the "everything else" template
		h_muRelPt_ee = (TH1D*)h_muRelPt_g->Clone("h_muRelPt_ee");
		h_muRelPt_ee->Add(h_muRelPt_lq);
		h_muRelPt_ee->Add(h_muRelPt_sall);
		// normalize the templates
		h_muRelPt_ee->Scale(1./h_muRelPt_ee->Integral());
		h_muRelPt_ball->Scale(1./h_muRelPt_ball->Integral());
		h_muRelPt_call->Scale(1./h_muRelPt_call->Integral());
		// subtract the ghosts
		h_muRelPt->Add(h_muRelPt_ghost,-1.0);
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

		h_muRelPt->SetMarkerStyle(22);
		h_muRelPt->SetMarkerColor(kBlue);
		h_muRelPt->SetFillColorAlpha(kBlue,0.5);

		h_muRelPt_ball->SetMarkerStyle(33);
		h_muRelPt_ball->SetMarkerColor(kOrange);
		h_muRelPt_ball->SetFillColorAlpha(kOrange,0.5);

		h_muRelPt_call->SetMarkerStyle(34);
		h_muRelPt_call->SetMarkerColor(kCyan+3);
		h_muRelPt_call->SetFillColorAlpha(kCyan+3,0.5);

		h_muRelPt_ee->SetMarkerStyle(43);
		h_muRelPt_ee->SetMarkerColor(kPink-4);
		h_muRelPt_ee->SetFillColorAlpha(kPink-4,0.5);

		double p0 = func->GetParameter(0); // b parameter
		double p1 = func->GetParameter(1); // c parameter
		double e0 = func->GetParError(0); // b parameter error
		double e1 = func->GetParError(1); // c parameter error


		h_muRelPt->SetStats(1);

		TH1D *h_muRelPt_ee_scaled = (TH1D*) h_muRelPt_ee->Clone("h_muRelPt_ee_scaled");
		TH1D *h_muRelPt_call_scaled = (TH1D*) h_muRelPt_call->Clone("h_muRelPt_call_scaled");
		TH1D *h_muRelPt_ball_scaled = (TH1D*) h_muRelPt_ball->Clone("h_muRelPt_ball_scaled");

		h_muRelPt_ee_scaled->Scale(1.0-p0-p1);
		h_muRelPt_call_scaled->Scale(p1);
		h_muRelPt_ball_scaled->Scale(p0);

		h_stack->Add(h_muRelPt_ee_scaled);
		h_stack->Add(h_muRelPt_call_scaled);
		h_stack->Add(h_muRelPt_ball_scaled);

		func->SetLineStyle(10);

		if(!isStack){
			h_muRelPt->Draw("e2");
			h_muRelPt_ball->Draw("e2 same");
			h_muRelPt_call->Draw("e2 same");
			h_muRelPt_ee->Draw("e2 same");
			func->Draw("same");

			if(centReg==1){h_muRelPt->SetTitle("0-10");}
			if(centReg==2){h_muRelPt->SetTitle("10-30");}
			if(centReg==3){h_muRelPt->SetTitle("30-50");}
			if(centReg==4){h_muRelPt->SetTitle("50-90");}





			auto legend = new TLegend(0.6,0.25,0.8,0.6);
			legend->AddEntry(h_muRelPt,"inclusive","fp");
			legend->AddEntry(h_muRelPt_ball,"b template","fp");
			legend->AddEntry(h_muRelPt_call,"c template","fp");
			legend->AddEntry(h_muRelPt_ee,"udsg template","fp");
			legend->AddEntry(func,"fit");
			legend->Draw();
		}

		if(isStack){
			
			h_stack->Draw("hist");
			h_muRelPt->Draw("ep same");
			func->Draw("same");
			auto legend = new TLegend(0.6,0.25,0.8,0.6);
			legend->AddEntry(h_muRelPt,"inclusive","p");
			legend->AddEntry(h_muRelPt_ball,"b template","f");
			legend->AddEntry(h_muRelPt_call,"c template","f");
			legend->AddEntry(h_muRelPt_ee,"udsg template","f");
			legend->AddEntry(func,"fit");
			legend->Draw();

			if(centReg==1){h_stack->SetTitle("0-10 % centrality");}
			if(centReg==2){h_stack->SetTitle("10-30 % centrality");}
			if(centReg==3){h_stack->SetTitle("30-50 % centrality");}
			if(centReg==4){h_stack->SetTitle("50-90 % centrality");}

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

		if(isStack){c1->SaveAs("/home/clayton/Analysis/code/templateFitCentralityRegions/figures/templateFitCentralityRegions_bVsCVsElse_stack.pdf");}
		if(!isStack){c1->SaveAs("/home/clayton/Analysis/code/templateFitCentralityRegions/figures/templateFitCentralityRegions_bVsCVsElse.pdf");}

		// check fit against count

		double n_tot = h_muRelPt_ee->GetEntries() + h_muRelPt_ball->GetEntries() + h_muRelPt_call->GetEntries();
		double n_udsg = h_muRelPt_ee->GetEntries();
		double n_b = h_muRelPt_ball->GetEntries();
		double n_c = h_muRelPt_call->GetEntries();

		double nr_udsg = n_udsg/n_tot;
		double nr_b = n_b/n_tot;
		double nr_c = n_c/n_tot;

		
		cout << "b fraction = " << nr_b << endl;
		cout << "c fraction = " << nr_c << endl;
		cout << "udsg fraction = " << nr_udsg << endl;


		

		
	
	


}

