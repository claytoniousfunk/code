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

void RecoGenCompare(Float_t etacut1=-1.5, Float_t etacut2=1.5){

	//TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_30_muptcut_5.root");
	//TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_30_muptcut_5.root");

	//TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_50_muptcut_5.root");
	//TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_5.root");

	//TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_80_muptcut_5.root");
	//TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_80_muptcut_5.root");

	//TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeDataV5_recojets_pthat_30_muptcut_5.root");
	//TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeDataV5_refjets_pthat_30_muptcut_5.root");

	//TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeDataV6_recojets_pthat_30_muptcut_5.root");
	//TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeDataV6_refjets_pthat_30_muptcut_5.root");
	
	TFile *f_reco = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");
	TFile *f_gen = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_refjets_pthat_30_muptcut_5.root");


	double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};

	TH2D *h2_reco, *h2_reco_g, *h2_reco_lq, *h2_reco_sall, *h2_reco_ball, *h2_reco_call;
	TH2D *h2_gen, *h2_gen_g, *h2_gen_lq, *h2_gen_sall, *h2_gen_ball, *h2_gen_call;
	//TH1D *h_reco, *h_reco_g, *h_reco_lq, *h_reco_sall, *h_reco_ball, *h_reco_call;
	//TH1D *h_gen, *h_gen_g, *h_gen_lq, *h_gen_sall, *h_gen_ball, *h_gen_call;

	f_reco->GetObject("h2",h2_reco);
	f_reco->GetObject("h2_g",h2_reco_g);
	f_reco->GetObject("h2_lq",h2_reco_lq);
	f_reco->GetObject("h2_sall",h2_reco_sall);
	f_reco->GetObject("h2_ball",h2_reco_ball);
	f_reco->GetObject("h2_call",h2_reco_call);

	TH1D* hx_test=h2_reco->ProjectionX();
	TAxis *xaxis = hx_test->GetXaxis();
	Int_t firstxbin = xaxis->FindBin(etacut1);
	cout << "first x bin = " << firstxbin << endl;
	Int_t lastxbin = xaxis->FindBin(etacut2);
	cout << "last x bin = " << lastxbin << endl;
	
	TH1D *h_reco_raw = h2_reco->ProjectionY("h_reco_raw",firstxbin,lastxbin);
	TH1D *h_reco_g_raw = h2_reco_g->ProjectionY("h_reco_g_raw",firstxbin,lastxbin);
	TH1D *h_reco_lq_raw = h2_reco_lq->ProjectionY("h_reco_lq_raw",firstxbin,lastxbin);
	TH1D *h_reco_sall_raw = h2_reco_sall->ProjectionY("h_reco_sall_raw",firstxbin,lastxbin);
	TH1D *h_reco_ball_raw = h2_reco_ball->ProjectionY("h_reco_ball_raw",firstxbin,lastxbin);
	TH1D *h_reco_call_raw = h2_reco_call->ProjectionY("h_reco_call_raw",firstxbin,lastxbin);

	TH1D *h_reco = (TH1D*) h_reco_raw->Rebin(15,"",pt_axis);
	TH1D *h_reco_g = (TH1D*) h_reco_g_raw->Rebin(15,"",pt_axis);
	TH1D *h_reco_lq = (TH1D*) h_reco_lq_raw->Rebin(15,"",pt_axis);
	TH1D *h_reco_sall = (TH1D*) h_reco_sall_raw->Rebin(15,"",pt_axis);
	TH1D *h_reco_ball = (TH1D*) h_reco_ball_raw->Rebin(15,"",pt_axis);
	TH1D *h_reco_call = (TH1D*) h_reco_call_raw->Rebin(15,"",pt_axis);

	f_gen->GetObject("h2",h2_gen);
	f_gen->GetObject("h2_g",h2_gen_g);
	f_gen->GetObject("h2_lq",h2_gen_lq);
	f_gen->GetObject("h2_sall",h2_gen_sall);
	f_gen->GetObject("h2_ball",h2_gen_ball);
	f_gen->GetObject("h2_call",h2_gen_call);
	
	TH1D *h_gen_raw = h2_gen->ProjectionY("h_gen_raw",firstxbin,lastxbin);
	TH1D *h_gen_g_raw = h2_gen_g->ProjectionY("h_gen_g_raw",firstxbin,lastxbin);
	TH1D *h_gen_lq_raw = h2_gen_lq->ProjectionY("h_gen_lq_raw",firstxbin,lastxbin);
	TH1D *h_gen_sall_raw = h2_gen_sall->ProjectionY("h_gen_sall_raw",firstxbin,lastxbin);
	TH1D *h_gen_ball_raw = h2_gen_ball->ProjectionY("h_gen_ball_raw",firstxbin,lastxbin);
	TH1D *h_gen_call_raw = h2_gen_call->ProjectionY("h_gen_call_raw",firstxbin,lastxbin);

	TH1D *h_gen = (TH1D*) h_gen_raw->Rebin(15,"",pt_axis);
	TH1D *h_gen_g = (TH1D*) h_gen_g_raw->Rebin(15,"",pt_axis);
	TH1D *h_gen_lq = (TH1D*) h_gen_lq_raw->Rebin(15,"",pt_axis);
	TH1D *h_gen_sall = (TH1D*) h_gen_sall_raw->Rebin(15,"",pt_axis);
	TH1D *h_gen_ball = (TH1D*) h_gen_ball_raw->Rebin(15,"",pt_axis);
	TH1D *h_gen_call = (TH1D*) h_gen_call_raw->Rebin(15,"",pt_axis);


	TH1D *h_reco_uds = (TH1D*)h_reco_lq->Clone("h_reco_uds");
	h_reco_uds->Add(h_reco_sall);

	TH1D *h_gen_uds = (TH1D*)h_gen_lq->Clone("h_gen_uds");
	h_gen_uds->Add(h_gen_sall);

	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
	pad1->Draw();

	TH1D *r_raw = (TH1D*)h_reco->Clone("r_raw");
	TH1D *r_g_raw = (TH1D*)h_reco_g->Clone("r_g_raw");
	TH1D *r_uds_raw = (TH1D*)h_reco_uds->Clone("r_uds_raw");
	TH1D *r_ball_raw = (TH1D*)h_reco_ball->Clone("r_ball_raw");
	TH1D *r_call_raw = (TH1D*)h_reco_call->Clone("r_call_raw");

	


	r_raw->Divide(h_gen);
	r_g_raw->Divide(h_gen_g);
	r_uds_raw->Divide(h_gen_uds);
	r_ball_raw->Divide(h_gen_ball);
	r_call_raw->Divide(h_gen_call);
	
	TH1D *r = (TH1D*)r_raw->Clone("r");
	TH1D *r_g = (TH1D*)r_g_raw->Clone("r_g");
	TH1D *r_uds = (TH1D*)r_uds_raw->Clone("r_uds");
	TH1D *r_ball = (TH1D*)r_ball_raw->Clone("r_ball");
	TH1D *r_call = (TH1D*)r_call_raw->Clone("r_call");
	
	/*
	TH1D *r = (TH1D*)r_raw->Rebin(15,"",pt_axis);
	TH1D *r_g = (TH1D*)r_g_raw->Rebin(15,"",pt_axis);
	TH1D *r_uds = (TH1D*)r_uds_raw->Rebin(15,"",pt_axis);
	TH1D *r_ball = (TH1D*)r_ball_raw->Rebin(15,"",pt_axis);
	TH1D *r_call = (TH1D*)r_call_raw->Rebin(15,"",pt_axis);
	*/
	r->SetMarkerStyle(8);
    r->SetMarkerColor(kAzure+9);
    r->SetMarkerSize(1.2);
    r->SetFillColorAlpha(kAzure+10,0.5);
    r->SetStats(0);
    //r->Draw("e2");
    

    r_g->SetMarkerStyle(33);
    r_g->SetMarkerColor(kBlue+4);
    r_g->SetMarkerSize(1.4);
    r_g->SetFillColorAlpha(kBlue-1,0.5);
    r_g->SetStats(0);
    r_g->Draw("e2");

    r_g->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
    r_g->GetYaxis()->SetTitle("reco jets / gen jets");
    r_g->SetTitle("");
    r_g->SetMinimum(0.0);
    r_g->SetMaximum(2.0);

    

    r_uds->SetMarkerStyle(45);
    r_uds->SetMarkerColor(kGreen-1);
    r_uds->SetMarkerSize(1.4);
    r_uds->SetFillColorAlpha(kGreen-5,0.5);
    r_uds->SetStats(0);
    r_uds->Draw("e2 same");

    r_ball->SetMarkerStyle(29);
    r_ball->SetMarkerColor(kRed+1);
    r_ball->SetMarkerSize(1.4);
    r_ball->SetFillColorAlpha(kRed,0.5);
    r_ball->SetStats(0);
    r_ball->Draw("e2 same");

    r_call->SetMarkerStyle(47);
    r_call->SetMarkerColor(kAzure+9);
    r_call->SetMarkerSize(1.4);
    r_call->SetFillColorAlpha(kAzure+10,0.5);
    r_call->SetStats(0);
    r_call->Draw("e2 same");

    auto legend = new TLegend(0.905,0.4,0.995,0.7);


    
    legend->AddEntry(r_call,"c","fp");
    legend->AddEntry(r_ball,"b","fp");
    legend->AddEntry(r_uds,"uds","fp");
    legend->AddEntry(r_g,"g","fp");
    legend->Draw();




    c1->SaveAs("RecoGenCompare.pdf");


}
