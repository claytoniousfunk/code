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

void taggedUntagged(Float_t etacut1=-1.5, Float_t etacut2=1.5){

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
	TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_refjets_pthat_30_muptcut_5.root");

	TH2D *h2_tagged, *h2_tagged_g, *h2_tagged_lq, *h2_tagged_sall, *h2_tagged_ball, *h2_tagged_call;
	TH2D *h2_untagged, *h2_untagged_g, *h2_untagged_lq, *h2_untagged_sall, *h2_untagged_ball, *h2_untagged_call;
	TH1D *h_tagged, *h_tagged_g, *h_tagged_lq, *h_tagged_sall, *h_tagged_ball, *h_tagged_call;
	TH1D *h_untagged, *h_untagged_g, *h_untagged_lq, *h_untagged_sall, *h_untagged_ball, *h_untagged_call;

	f->GetObject("h2_MJ",h2_tagged);
	f->GetObject("h2_g_MJ",h2_tagged_g);
	f->GetObject("h2_lq_MJ",h2_tagged_lq);
	f->GetObject("h2_sall_MJ",h2_tagged_sall);
	f->GetObject("h2_ball_MJ",h2_tagged_ball);
	f->GetObject("h2_call_MJ",h2_tagged_call);

	TH1D* hx_test=h2_tagged->ProjectionX();
	TAxis *xaxis = hx_test->GetXaxis();
	Int_t firstxbin = xaxis->FindBin(etacut1);
	cout << "first x bin = " << firstxbin << endl;
	Int_t lastxbin = xaxis->FindBin(etacut2);
	cout << "last x bin = " << lastxbin << endl;
	
	h_tagged = h2_tagged->ProjectionY("h_tagged",firstxbin,lastxbin);
	h_tagged_g = h2_tagged_g->ProjectionY("h_tagged_g",firstxbin,lastxbin);
	h_tagged_lq = h2_tagged_lq->ProjectionY("h_tagged_lq",firstxbin,lastxbin);
	h_tagged_sall = h2_tagged_sall->ProjectionY("h_tagged_sall",firstxbin,lastxbin);
	h_tagged_ball = h2_tagged_ball->ProjectionY("h_tagged_ball",firstxbin,lastxbin);
	h_tagged_call = h2_tagged_call->ProjectionY("h_tagged_call",firstxbin,lastxbin);

	f->GetObject("h2",h2_untagged);
	f->GetObject("h2_g",h2_untagged_g);
	f->GetObject("h2_lq",h2_untagged_lq);
	f->GetObject("h2_sall",h2_untagged_sall);
	f->GetObject("h2_ball",h2_untagged_ball);
	f->GetObject("h2_call",h2_untagged_call);
	
	h_untagged = h2_untagged->ProjectionY("h_untagged",firstxbin,lastxbin);
	h_untagged_g = h2_untagged_g->ProjectionY("h_untagged_g",firstxbin,lastxbin);
	h_untagged_lq = h2_untagged_lq->ProjectionY("h_untagged_lq",firstxbin,lastxbin);
	h_untagged_sall = h2_untagged_sall->ProjectionY("h_untagged_sall",firstxbin,lastxbin);
	h_untagged_ball = h2_untagged_ball->ProjectionY("h_untagged_ball",firstxbin,lastxbin);
	h_untagged_call = h2_untagged_call->ProjectionY("h_untagged_call",firstxbin,lastxbin);

	TH1D *h_tagged_uds = (TH1D*)h_tagged_lq->Clone("h_tagged_uds");
	h_tagged_uds->Add(h_tagged_sall);

	TH1D *h_untagged_uds = (TH1D*)h_untagged_lq->Clone("h_untagged_uds");
	h_untagged_uds->Add(h_untagged_sall);

	TCanvas *c1 = new TCanvas("c1","c1",1200,500);
	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,0.5,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.5,0.0,1.0,1.0);
	pad1->Draw();
	pad2->Draw();

	pad1->SetLeftMargin(0.15);
	pad2->SetLeftMargin(0.15);

	TH1D *r = (TH1D*)h_tagged->Clone("r");
	TH1D *r_g = (TH1D*)h_tagged_g->Clone("r_g");
	TH1D *r_uds = (TH1D*)h_tagged_uds->Clone("r_uds");
	TH1D *r_ball = (TH1D*)h_tagged_ball->Clone("r_ball");
	TH1D *r_call = (TH1D*)h_tagged_call->Clone("r_call");

	const int NPtAxisBins = 19;
	double pt_axis[NPtAxisBins+1] = {50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
	
	TH1D *h_untagged_rebin = (TH1D*) h_untagged->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *h_untagged_g_rebin = (TH1D*) h_untagged_g->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *h_untagged_uds_rebin = (TH1D*) h_untagged_uds->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *h_untagged_ball_rebin = (TH1D*) h_untagged_ball->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *h_untagged_call_rebin = (TH1D*) h_untagged_call->Rebin(NPtAxisBins,"",pt_axis);

	TH1D *r_rebin = (TH1D*) r->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *r_g_rebin = (TH1D*) r_g->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *r_uds_rebin = (TH1D*) r_uds->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *r_ball_rebin = (TH1D*) r_ball->Rebin(NPtAxisBins,"",pt_axis);
	TH1D *r_call_rebin = (TH1D*) r_call->Rebin(NPtAxisBins,"",pt_axis);


	r_rebin->Divide(h_untagged_rebin);
	r_g_rebin->Divide(h_untagged_g_rebin);
	r_uds_rebin->Divide(h_untagged_uds_rebin);
	r_ball_rebin->Divide(h_untagged_ball_rebin);
	r_call_rebin->Divide(h_untagged_call_rebin);

	


	r_rebin->SetMarkerStyle(8);
    r_rebin->SetMarkerColor(kAzure+9);
    r_rebin->SetMarkerSize(1.2);
    r_rebin->SetFillColorAlpha(kAzure+10,0.5);
    r_rebin->SetStats(0);
    //r->Draw("e2");

    pad1->cd();
    r_g_rebin->SetMarkerStyle(33);
    r_g_rebin->SetMarkerColor(kBlue+4);
    r_g_rebin->SetMarkerSize(1.4);
    r_g_rebin->SetFillColorAlpha(kBlue-1,0.5);
    r_g_rebin->SetStats(0);
    r_g_rebin->Draw("e2");

    r_g_rebin->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r_g_rebin->GetYaxis()->SetTitle("#mu-jets / inclusive jets");
    r_g_rebin->SetTitle("");
    r_g_rebin->SetMinimum(0.0);
    r_g_rebin->SetMaximum(0.02);

    r_uds_rebin->SetMarkerStyle(45);
    r_uds_rebin->SetMarkerColor(kGreen-1);
    r_uds_rebin->SetMarkerSize(1.4);
    r_uds_rebin->SetFillColorAlpha(kGreen-5,0.5);
    r_uds_rebin->SetStats(0);
    r_uds_rebin->Draw("e2 same");

    pad2->cd();
    r_ball_rebin->SetMarkerStyle(29);
    r_ball_rebin->SetMarkerColor(kRed+1);
    r_ball_rebin->SetMarkerSize(1.4);
    r_ball_rebin->SetFillColorAlpha(kRed,0.5);
    r_ball_rebin->SetStats(0);
    r_ball_rebin->Draw("e2");

    r_ball_rebin->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r_ball_rebin->GetYaxis()->SetTitle("#mu-jets / inclusive jets");
    r_ball_rebin->SetTitle("");
    r_ball_rebin->SetMinimum(0.0);
    r_ball_rebin->SetMaximum(0.3);

    r_call_rebin->SetMarkerStyle(47);
    r_call_rebin->SetMarkerColor(kAzure+9);
    r_call_rebin->SetMarkerSize(1.4);
    r_call_rebin->SetFillColorAlpha(kAzure+10,0.5);
    r_call_rebin->SetStats(0);
    r_call_rebin->Draw("e2 same");


    

    auto legend = new TLegend(0.905,0.4,0.995,0.7);

    
    legend->AddEntry(r_ball_rebin,"b","fp");
    legend->AddEntry(r_call_rebin,"c","fp");
    legend->AddEntry(r_g_rebin,"g","fp");
    legend->AddEntry(r_uds_rebin,"uds","fp");

    legend->Draw();




    c1->SaveAs("/home/clayton/Analysis/code/taggedUntagged/figures/taggedUntagged.pdf");


}
