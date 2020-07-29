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

void taggedUntaggedRecoGen(Float_t etacut1=-1.5, Float_t etacut2=1.5){

	//TFile *f1 = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_10.root");
	//TFile *f2 = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_refjets_pthat_30_jetptcut_50_muptcut_10.root");
	TFile *f1 = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_50_jetptcut_50_muptcut_5.root");
	TFile *f2 = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_refjets_pthat_50_jetptcut_50_muptcut_5.root");

	TH2D *h2_tagged, *h2_tagged_g, *h2_tagged_lq, *h2_tagged_sall, *h2_tagged_ball, *h2_tagged_call;
	TH2D *h2_untagged, *h2_untagged_g, *h2_untagged_lq, *h2_untagged_sall, *h2_untagged_ball, *h2_untagged_call;
	TH1D *h_tagged, *h_tagged_g, *h_tagged_lq, *h_tagged_sall, *h_tagged_ball, *h_tagged_call;
	TH1D *h_untagged, *h_untagged_g, *h_untagged_lq, *h_untagged_sall, *h_untagged_ball, *h_untagged_call;

	TH2D *h2_untagged_reco;
	TH1D *h_untagged_reco;

	f1->GetObject("h2_MJ",h2_tagged);
	f1->GetObject("h2_g_MJ",h2_tagged_g);
	f1->GetObject("h2_lq_MJ",h2_tagged_lq);
	f1->GetObject("h2_sall_MJ",h2_tagged_sall);
	f1->GetObject("h2_ball_MJ",h2_tagged_ball);
	f1->GetObject("h2_call_MJ",h2_tagged_call);

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

	f1->GetObject("h2",h2_untagged_reco);
	f2->GetObject("h2",h2_untagged);
	f2->GetObject("h2_g",h2_untagged_g);
	f2->GetObject("h2_lq",h2_untagged_lq);
	f2->GetObject("h2_sall",h2_untagged_sall);
	f2->GetObject("h2_ball",h2_untagged_ball);
	f2->GetObject("h2_call",h2_untagged_call);
	
	h_untagged_reco = h2_untagged_reco->ProjectionY("h_untagged_reco",firstxbin,lastxbin);
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

	TCanvas *c1 = new TCanvas("c1","c1",900,400);
	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
	pad1->Draw();

	TH1D *r = (TH1D*)h_tagged->Clone("r");
	TH1D *r_g = (TH1D*)h_tagged_g->Clone("r_g");
	TH1D *r_uds = (TH1D*)h_tagged_uds->Clone("r_uds");
	TH1D *r_ball = (TH1D*)h_tagged_ball->Clone("r_ball");
	TH1D *r_call = (TH1D*)h_tagged_call->Clone("r_call");
	TH1D *r_untagged_recoToGen = (TH1D*)h_untagged_reco->Clone("r_untagged_recoToGen");

	r->Divide(h_untagged);
	r_g->Divide(h_untagged_g);
	r_uds->Divide(h_untagged_uds);
	r_ball->Divide(h_untagged_ball);
	r_call->Divide(h_untagged_call);
	r_untagged_recoToGen->Divide(h_untagged);
	

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
    //r_g->Draw("e2");

    r_g->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r_g->GetYaxis()->SetTitle("#frac{number of tagged (muon) reco jets}{number of untagged gen jets}");
    r_g->SetTitle("");
    //r_g->SetMinimum(0.0);
    //r_g->SetMaximum(1.0);

    r_uds->SetMarkerStyle(45);
    r_uds->SetMarkerColor(kGreen-1);
    r_uds->SetMarkerSize(1.4);
    r_uds->SetFillColorAlpha(kGreen-5,0.5);
    r_uds->SetStats(0);
    //r_uds->Draw("e2 same");

    r_ball->SetMarkerStyle(29);
    r_ball->SetMarkerColor(kRed+1);
    r_ball->SetMarkerSize(1.4);
    r_ball->SetFillColorAlpha(kRed,0.5);
    r_ball->SetStats(0);
    //r_ball->Draw("e2");

    r_ball->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r_ball->GetYaxis()->SetTitle("#frac{number of tagged (muon) reco jets}{number of untagged gen jets}");
    r_ball->SetTitle("");
    r_ball->SetMinimum(0.0);
    //r_g->SetMaximum(1.0);

    r_call->SetMarkerStyle(47);
    r_call->SetMarkerColor(kAzure+9);
    r_call->SetMarkerSize(1.4);
    r_call->SetFillColorAlpha(kAzure+10,0.5);
    r_call->SetStats(0);
    //r_call->Draw("e2 same");

    r_untagged_recoToGen->SetMarkerStyle(43);
    r_untagged_recoToGen->SetMarkerColor(kAzure-1);
    r_untagged_recoToGen->SetMarkerSize(1.4);
    r_untagged_recoToGen->SetFillColorAlpha(kAzure-2,0.5);
    r_untagged_recoToGen->SetStats(0);
    r_untagged_recoToGen->Draw("e2");

    r_untagged_recoToGen->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r_untagged_recoToGen->GetYaxis()->SetTitle("#frac{number of untagged reco jets}{number of untagged gen jets}");
    r_untagged_recoToGen->SetTitle("");
    r_ball->SetMinimum(0.0);


    

    auto legend = new TLegend(0.905,0.4,0.995,0.7);

    //legend->AddEntry(r_call,"c","fp");
    //legend->AddEntry(r_ball,"b","fp");
    //legend->AddEntry(r_uds,"uds","fp");
    //legend->AddEntry(r_g,"g","fp");
    //legend->Draw();




    c1->SaveAs("taggedUntaggedRecoGen.pdf");


}